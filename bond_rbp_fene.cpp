/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */
/* ----------------------------------------------------------------------
   Contributing author: Enrico Skoruppa (School of Physics and Astronomy, University of Edinburgh, Edinburgh)
------------------------------------------------------------------------- */


#include "bond_rbp_fene.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "neighbor.h"
#include "error.h"
#include "so3.h"
#include "lamath.h"
#include "domain.h"
#include "math_extra.h"
#include "update.h"
#include "atom_vec_ellipsoid.h"
#include "utils.h"

#include <cstring>    
#include <algorithm> 

using namespace LAMMPS_NS;

/* ----------------------------------------------------------------------
   Constructor: initializes internal data structures and sets defaults
------------------------------------------------------------------------- */

BondRBPFene::BondRBPFene(LAMMPS *lmp) : Bond(lmp) {
  params = nullptr;
}

/* ----------------------------------------------------------------------
   Destructor: deallocates dynamically allocated memory
------------------------------------------------------------------------- */

BondRBPFene::~BondRBPFene() {
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(params);
    allocated = 0;
  }
}

/* ----------------------------------------------------------------------
   Compute forces and energy contribution for all bonds of this style
------------------------------------------------------------------------- */

void BondRBPFene::compute(int eflag, int vflag) {

  int id1, id2, bond_type;
  
  double T1[3][3],T2[3][3];
  // double T1tp[3][3];
  double r1[3], r2[3], dr[3];
  double R[3][3];
  double Om[3], w[3];
  double Omd[3], wd[3]; 
  double Jinvtp[3][3];

  double A[3], B[3];
  double torque_1[3], torque_2[3];
  double force_1[3];

  double tmp1[3],tmp2[3];

  double **x = atom->x;                
  double **f = atom->f;             
  double **torque = atom->torque;   
  double *quat1,*quat2;

  // fene constant
  // const double rlogarg_min = 0.2;
  const double rlogarg_min = 0.1;
  
  auto avec = dynamic_cast<AtomVecEllipsoid *>(atom->style_match("ellipsoid"));
  AtomVecEllipsoid::Bonus *bonus = avec->bonus;
  int *ellipsoid = atom->ellipsoid;

  ev_init(eflag, vflag);

  // Access bond list
  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  // @Oliver: Can you confirm that this is safe to use with hybrid bond style? Will bondlist only contain the bonds of this bond style?  
  for (int bid = 0; bid < nbondlist; bid++) {
    id1 = bondlist[bid][0];
    id2 = bondlist[bid][1];
    bond_type = bondlist[bid][2];

    //-------------------------------------------------------------//
    // Compute triads and positions
    //-------------------------------------------------------------//

    // get quaternions
    quat1=bonus[ellipsoid[id1]].quat;
    quat2=bonus[ellipsoid[id2]].quat;

    // transform quat to triads [SO(3)]
    MathExtra::quat_to_mat(quat1, T1);
    MathExtra::quat_to_mat(quat2, T2);

    // compute R
    lamath::mul_AtB(T1,T2,R);

    // get positions
    r1[0] = x[id1][0];
    r1[1] = x[id1][1];
    r1[2] = x[id1][2];
    r2[0] = x[id2][0];
    r2[1] = x[id2][1];
    r2[2] = x[id2][2];

    // compute dr
    lamath::subtract(r2,r1,dr);
    // check for domain mismatch of special neighbor
    if (domain->minimum_image_check(dr[0],dr[1],dr[2])) {
      domain->minimum_image(FLERR,dr[0],dr[1],dr[2]);
    }

    // compute w
    lamath::mul_Atx(T1,dr,w);

    if (params[bond_type].subtract_groundstate) {
      //-------------------------------------------------------------//
      // Compute force wrench for se(3) (X) convention
      //-------------------------------------------------------------//

      // compute Omega
      so3::rotmat2euler(R,Om);
      // compute Omega_Delta
      lamath::subtract(Om,params[bond_type].srot,Omd);
      // compute w_Delta
      lamath::subtract(w,params[bond_type].svec,wd);
    
      // compute partial E / partial Omega_Delta
      lamath::mul(params[bond_type].Mr,Omd,A);
      lamath::mul(params[bond_type].Mtr_tr,wd,tmp1);
      lamath::add_to(A,tmp1);

      // compute partial E / partial Omega_Delta
      lamath::mul(params[bond_type].Mtr_bl,Omd,B);
      lamath::mul(params[bond_type].Mt,wd,tmp1);
      lamath::add_to(B,tmp1);

      // compute transposed inverse left Jacobian
      so3::leftJacobianInverseTransposed(Om,Jinvtp);

      // compute force
      lamath::mul(T1,B,force_1);
      
      // compute torque
      lamath::mul(Jinvtp,A,tmp1);
      lamath::mul(T1,tmp1,torque_2);

      MathExtra::cross3(w,B,tmp1);
      lamath::mul(T1,tmp1,torque_1);
      lamath::add_to(torque_1,torque_2);

    }
    else{
      //-------------------------------------------------------------//
      // Compute force wrench for SE(3) (Y) convention
      //-------------------------------------------------------------//
      
      double Dmat[3][3];
      // lamath::mul(params[bond_type].Smat_tp,R,Dmat);
      lamath::mul_AtB(params[bond_type].Smat,R,Dmat);

      // compute Phi_delta (reuse Omd for this)
      so3::rotmat2euler(Dmat,Omd);

      // compute d (reuse wd for this)
      lamath::subtract(w,params[bond_type].svec,tmp1);
      lamath::mul_Atx(params[bond_type].Smat,tmp1,wd);
      // lamath::mul(params[bond_type].Smat_tp,tmp1,wd);

      // compute partial E / partial Omega_Delta
      lamath::mul(params[bond_type].Mr,Omd,A);
      lamath::mul(params[bond_type].Mtr_tr,wd,tmp1);
      lamath::add_to(A,tmp1);

      // compute partial E / partial Omega_Delta
      lamath::mul(params[bond_type].Mtr_bl,Omd,B);
      lamath::mul(params[bond_type].Mt,wd,tmp1);
      lamath::add_to(B,tmp1);

      // compute transposed inverse left Jacobian
      so3::leftJacobianInverseTransposed(Omd,Jinvtp);

      // compute torque 2
      lamath::mul(Jinvtp,A,tmp1);
      lamath::mul(params[bond_type].Smat,tmp1,tmp2);
      lamath::mul(T1,tmp2,torque_2);

      // compute force 1
      lamath::mul(params[bond_type].Smat,B,tmp1);
      lamath::mul(T1,tmp1,force_1);

      // compute torque 1
      MathExtra::cross3(w,tmp1,tmp2);
      lamath::mul(T1,tmp2,torque_1);
      lamath::add_to(torque_1,torque_2);
    }


    // ///////////////////////////////////////////////////////////////
    // ///////////////////////////////////////////////////////////////
    // // DEBUG CODE 
    // double be = 0;
    // double Yv[6]; 
    // Yv[0] = Omd[0]; 
    // Yv[1] = Omd[1];
    // Yv[2] = Omd[2];
    // Yv[3] = wd[0];   
    // Yv[4] = wd[1];
    // Yv[5] = wd[2];
    // // Calculate energy: 0.5 * Y_vec^T * Mmat * Y_vec
    // for (int i_mat = 0; i_mat < 6; ++i_mat) { 
    //   for (int j_mat = 0; j_mat < 6; ++j_mat) {
    //     be += Yv[i_mat] * params[bond_type].Mmat[i_mat][j_mat] * Yv[j_mat];
    //   }
    // }
    // be *= 0.5;
    // if (be > 24) {
    //   error->warning(FLERR, "High elastic energy: {} {} {} {}", update->ntimestep, atom->tag[id1],
    //     atom->tag[id2], be);
    //     fprintf(screen," High elastic energy: %.3f\n",be);
    //     fprintf(screen," %.4f %.4f %.4f\n",Omd[0],Omd[1],Omd[2]);
    //     fprintf(screen," %.4f %.4f %.4f\n",wd[0],wd[1],wd[2]);
    // }
    // ///////////////////////////////////////////////////////////////
    // ///////////////////////////////////////////////////////////////


    //-------------------------------------------------------------//
    // Compute force due to FENE potential
    //-------------------------------------------------------------//
    
    // FENE force (lab frame) and energy for this bond
    double Flab_fene[3] = {0.0, 0.0, 0.0};
    double E_fene = 0.0;
    
    if (params[bond_type].fene_active) {
      
      // r = ||wd|| in T1 frame
      double r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
      double r  = sqrt(r2);
      
      // FENE only acts for r >= Rc
      if (r >= params[bond_type].Rc && r > 0) {
        
        // shorthand
        double Rc      = params[bond_type].Rc;      // plays role of r0
        double Rspan   = params[bond_type].Rspan;   // Delta
        double Rspan2  = params[bond_type].Rspan2;  // Delta^2
        double Kf      = params[bond_type].K;       // K in your definition
        double keff    = Kf * Rspan2;               // k_eff = K * (R0 - Rc)^2
        
        double rdel     = r - Rc;
        double rdelsq   = rdel * rdel;
        
        // rlogarg = 1 - (r - Rc)^2 / (R0 - Rc)^2
        double rlogarg = 1.0 - rdelsq / Rspan2;
        
        // Overstretch protection logic, adapted from your snippet
        if (rlogarg < rlogarg_min) {
          error->warning(FLERR, "RBP FENE bond too long: {} {} {} {}",
            update->ntimestep, atom->tag[id1],
            atom->tag[id2], r);
          if (rlogarg <= -3.0) {
          // if (rlogarg <= -5.0) {
            fprintf(screen, "\nRc = %.3f, R0 = %.3f rlogarg = %.3f\n",params[bond_type].Rc,params[bond_type].R0,rlogarg);
            fprintf(screen, "keff = %.3f\n",keff);
            fprintf(screen, "Kf = %.3f\n",Kf);
            fprintf(screen, "r = %.3f\n",r);
            fprintf(screen, "w  = %.3f %.3f %.3f\n",w[0],w[1],w[2]);
            fprintf(screen, "wd = %.3f %.3f %.3f\n",wd[0],wd[1],wd[2]);
            fprintf(screen, "ws = %.3f %.3f %.3f\n",params[bond_type].svec[0],params[bond_type].svec[1],params[bond_type].svec[2]);
            for (int ii=0;ii<6;ii++) {
              for (int jj=0;jj<6;jj++) {
                fprintf(screen, "%.3f ",params[bond_type].Mmat[ii][jj]);
              }
              fprintf(screen, "\n");
            }
            if (!lamath::is_positive_definite(params[bond_type].Mmat)) {
              error->all(FLERR, "Stiffness matrix M is not positive definite");
            }
            error->one(FLERR, "Bad RBP FENE bond: r = {}, w = {} {} {}",r,w[0],w[1],w[2]);
          }
          rlogarg = rlogarg_min;
          
          // if overstretched F(r)=F(r_max)=F_max, E(r)=E(r_max)+F_max*(r-r_max)
          rdel = Rspan * sqrt(1.0 - rlogarg);
          if (eflag) {
            E_fene = -0.5 * keff * log(rlogarg)
            + keff * sqrt(1.0 - rlogarg) / (rlogarg * Rspan)
            * (r - Rc - Rspan * sqrt(1.0 - rlogarg));
            // fprintf(screen, "E_fene_in = %.3f\n",E_fene);
          }
        }
        else if (eflag) {
          E_fene = -0.5 * keff * log(rlogarg);
        }
        
        double fbond = keff * rdel / (rlogarg * Rspan2 * r);
        force_1[0] += fbond * dr[0];
        force_1[1] += fbond * dr[1];
        force_1[2] += fbond * dr[2];
      }
    }

    //-------------------------------------------------------------//
    // Apply forces and torques to each atom
    //-------------------------------------------------------------//

    if (newton_bond || id1 < nlocal) {

      f[id1][0] += force_1[0];
      f[id1][1] += force_1[1];
      f[id1][2] += force_1[2];

      torque[id1][0] += torque_1[0];
      torque[id1][1] += torque_1[1];
      torque[id1][2] += torque_1[2];
      
    }

    if (newton_bond || id2 < nlocal) {

      f[id2][0] -= force_1[0];
      f[id2][1] -= force_1[1];
      f[id2][2] -= force_1[2];

      torque[id2][0] -= torque_2[0];
      torque[id2][1] -= torque_2[1];
      torque[id2][2] -= torque_2[2];
    }

    if (evflag) {

      double bond_energy = 0.0;
      const double Yv[6] = {Omd[0], Omd[1], Omd[2], wd[0], wd[1], wd[2]};
      const auto& M = params[bond_type].Mmat;
      for (int i = 0; i < 6; ++i) {
          bond_energy += Yv[i] * M[i][i] * Yv[i];
          for (int j = i + 1; j < 6; ++j) {
              bond_energy += 2.0 * Yv[i] * M[i][j] * Yv[j];
          }
      }
      bond_energy *= 0.5;

      if (params[bond_type].fene_active) {
        bond_energy += E_fene;
      }

      ev_tally_xyz(id1, id2, nlocal, newton_bond, bond_energy,
                  force_1[0], force_1[1], force_1[2],
                  -dr[0], -dr[1], -dr[2]);  // -(r2-r1) = r1-r2
    }
  }
}


/* ----------------------------------------------------------------------
   Compute single
------------------------------------------------------------------------- */

double BondRBPFene::single(int i, double rsq, int itype, int jtype, double &fforce) {
  fforce = 0.0;
  return 0.0;
}


void BondRBPFene::coeff(int narg, char **arg)
{

  // check conditions for properly formated arg !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // may not be necessary

  if (!allocated) allocate();

  int ilo, ihi;
  utils::bounds(FLERR, arg[0], 1, atom->nbondtypes, ilo, ihi, error);

  // load from database file
  if (narg == 4 && strcmp(arg[1], "dbfile") == 0) {

    RBPDatabase db(lmp, error);
    db.read(arg[2]);
    
    // Validate that bond style in database matches this style
    if (db.metadata().bond_style != RBPBOND_FENE_STYLE) {
      error->all(FLERR, 
        "Bond style mismatch: database specifies '" + db.metadata().bond_style + 
        "' but using bond_style rbp/fene");
    }

    int dbid = utils::inumeric(FLERR, arg[3], false, lmp);
    for (int bond_type=ilo;bond_type<=ihi;bond_type++) {
      assign_coeffs(bond_type,db.bond(dbid++).coeffs,db.metadata().subtract_groundstate);
    }
    return;
  }
  
  // assign from passed args
  std::vector<double> coeffs;
  coeffs.reserve(narg - 1);
  for (int i = 1; i < narg; i++) {
    double val = utils::numeric(FLERR, arg[i], false, lmp);  
    coeffs.push_back(val);
  }
  for (int bond_type=ilo;bond_type<=ihi;bond_type++) {
    assign_coeffs(bond_type,coeffs,RBPFENE_BOND_DEFAULT_SUBTRACT_GROUNDSTATE);
  }
}


/* ----------------------------------------------------------------------
   Initialization: verify atom properties and prepare style
------------------------------------------------------------------------- */

void BondRBPFene::init_style() {

  // Ensure ellipsoidal atoms are being used
  if (!atom->ellipsoid_flag)
    error->all(FLERR, "Bond style rbp requires atom style with ellipsoids");

  if (domain->dimension != 3)
    error->all(FLERR, "Bond style rbp requires a 3D simulation");

  if (force->special_lj[1] != 0.0) {
    if (comm->me == 0) error->warning(FLERR, "Use special bonds = 0,x,x with bond style rbp");
  }
}


/* ----------------------------------------------------------------------
   Return equilibrium distance for bond type i
------------------------------------------------------------------------- */

double BondRBPFene::equilibrium_distance(int bond_type) 
{
  return params[bond_type].equidist;
}

/* ----------------------------------------------------------------------
   Write bond style parameters to restart file
------------------------------------------------------------------------- */

void BondRBPFene::write_restart(FILE *fp) {

  fwrite(&params[1], sizeof(RBPParams), atom->nbondtypes, fp);

}

/* ----------------------------------------------------------------------
   Read bond style parameters from restart file
------------------------------------------------------------------------- */

void BondRBPFene::read_restart(FILE *fp) {

  allocate();

  if (comm->me == 0) {
    utils::sfread(FLERR,
                  &params[1],
                  sizeof(RBPParams),
                  atom->nbondtypes,
                  fp,
                  nullptr,
                  error);
  }

  // broadcast the whole params block to all ranks
  MPI_Bcast(&params[1],
            atom->nbondtypes * static_cast<int>(sizeof(RBPParams)),
            MPI_BYTE,
            0,
            world);

  // mark all types as having coefficients
  for (int i = 1; i <= atom->nbondtypes; i++) setflag[i] = 1;

}

/* ----------------------------------------------------------------------
   Write bond style data (coefficients) to LAMMPS data file
------------------------------------------------------------------------- */

void BondRBPFene::write_data(FILE *fp) {

  for (int i = 1; i <= atom->nbondtypes; i++) {
    if (!setflag[i]) continue;

    // type index
    fprintf(fp, "%d", i);

    // 6 ground-state components (Ystatic)
    for (int k = 0; k < 6; k++)
      fprintf(fp, " %g", params[i].Ystatic[k]);

    // 21 upper-triangular stiffness components (Mmat)
    for (int r = 0; r < 6; r++)
      for (int c = r; c < 6; c++)
        fprintf(fp, " %g", params[i].Mmat[r][c]);

    fprintf(fp, "\n");
  }

}

/* ----------------------------------------------------------------------
   Allocate memory for internal data structures
------------------------------------------------------------------------- */

void BondRBPFene::allocate() {
  allocated = 1;
  int n = atom->nbondtypes;

  memory->create(params, n + 1, "bond:rbp:params");
  memory->create(setflag, n+1, "bond:rbp:setflag");
  for (int i = 0; i <= n; i++) setflag[i] = 0;
}


void BondRBPFene::assign_coeffs(int bond_type, const std::vector<double> &args, bool subtract_groundstate) {
  // ------------------------------------------------------------
  // Inline numeric definitions (legacy / TWLC / full RBP)
  // Supported forms:
  //   6 + 6  = 12 args (diagonal M)
  //   6 + 12 = 18 args (block-diagonal M)
  //   6 + 21 = 27 args (full upper-triangular M)
  //
  // these are the arguments following the leading bond_type id
  // ------------------------------------------------------------

  constexpr int opt1 = 15;
  constexpr int opt2 = 21;
  constexpr int opt3 = 30;

  int narg = static_cast<int>(args.size());

  if (narg != opt1 && narg != opt2 && narg != opt3) {
    std::string msg =
      "Incorrect number of arguments for bond_coeff (" + std::to_string(narg) + ")\n" +
      "Expected:\n" +
      "  3+6+6  (3 fene, 6 groundstate and 6 diagonal components of M)\n" +
      "  3+6+12 (3 fene, 6 groundstate and block-diagonal M: Upper triangular assignment for each 3x3 block. Assignment by first iterating through the rows)\n" +
      "  3+6+21 (3 fene, 6 groundstate and full symmetric M assignment. Assignment by first iterating through the rows)";
    error->all(FLERR, msg.c_str());
  }

  params[bond_type].subtract_groundstate = subtract_groundstate;

  // -------------------------------------------------------------------
  // fene parameters
  int argid = 0;
  params[bond_type].K     = args[argid++];
  params[bond_type].Rc    = args[argid++];
  params[bond_type].R0    = args[argid++];

  if (params[bond_type].Rc < 0) {
    error->warning(FLERR, "Passed Rc < 0. Set to zero.");
    params[bond_type].Rc = 0;
  }


  params[bond_type].Rspan  = params[bond_type].R0 - params[bond_type].Rc;
  params[bond_type].Rspan2 = params[bond_type].Rspan * params[bond_type].Rspan;
  if (params[bond_type].K > 0 && params[bond_type].Rspan > 0) {
    params[bond_type].fene_active = true;
  }
  else {
    error->warning(FLERR, "Invalid FENE arguments. FENE deactivated.");
  }


  // -------------------------------------------------------------------
  // read groundstate
  for (int i = 0; i < 6; i++)
    params[bond_type].Ystatic[i] = args[argid++];

  set_static_(bond_type);
  set_equidist(bond_type);
  zero_stiffmat_(bond_type);

  // -------------------------------------------------------------------
  // set diagnoal only
  if (narg == opt1) {
    for (int i = 0; i < 6; i++)
      params[bond_type].Mmat[i][i] = args[argid++];
  }

  // -------------------------------------------------------------------
  // set rotation and translation blocks independently (block diagonal)
  if (narg == opt2) {
    params[bond_type].Mmat[0][0] = args[argid++];
    params[bond_type].Mmat[0][1] = args[argid++];
    params[bond_type].Mmat[0][2] = args[argid++];
    params[bond_type].Mmat[1][1] = args[argid++];
    params[bond_type].Mmat[1][2] = args[argid++];
    params[bond_type].Mmat[2][2] = args[argid++];
    params[bond_type].Mmat[3][3] = args[argid++];
    params[bond_type].Mmat[3][4] = args[argid++];
    params[bond_type].Mmat[3][5] = args[argid++];
    params[bond_type].Mmat[4][4] = args[argid++];
    params[bond_type].Mmat[4][5] = args[argid++];
    params[bond_type].Mmat[5][5] = args[argid++];
  }

  
  
  // -------------------------------------------------------------------
  // set everything by assigning upper triangular coefficients
  if (narg == opt3) {
    // int k = 6;
    for (int i = 0; i < 6; i++) {
      for (int j = i; j < 6; j++) {
        params[bond_type].Mmat[i][j] = args[argid++];
      }
    }
  }
    
  set_lower_triangle_(bond_type);
  assign_blocks_(bond_type);
  
  if (!lamath::is_positive_definite(params[bond_type].Mmat)) {
    error->all(FLERR, "Stiffness matrix M is not positive definite");
  }

  setflag[bond_type] = 1;
}

void BondRBPFene::set_static_(int bond_type) {
  // set rotation and translation components and static rotation matrix
  for (int i = 0; i < 3; i++) {
    params[bond_type].srot[i] = params[bond_type].Ystatic[i];
    params[bond_type].svec[i] = params[bond_type].Ystatic[i + 3];
  }
  so3::euler2rotmat(params[bond_type].srot, params[bond_type].Smat);
}

void BondRBPFene::set_equidist(int bond_type) {
  params[bond_type].equidist = lamath::norm(params[bond_type].svec);
}

void BondRBPFene::zero_stiffmat_(int bond_type) {
  // initialize stiffness matrix to zero
  for (int i = 0; i < 6; i++)
    for (int j = 0; j < 6; j++)
      params[bond_type].Mmat[i][j] = 0.0;
}

void BondRBPFene::set_lower_triangle_(int bond_type) {
  // symmetrize by setting lower triangular components
  for (int i = 0; i < 6; i++)
    for (int j = i + 1; j < 6; j++)
      params[bond_type].Mmat[j][i] = params[bond_type].Mmat[i][j];
}

void BondRBPFene::assign_blocks_(int bond_type) {
  // assign rotation, translation and coupling blocks
  for (int r = 0; r < 3; r++) {
    for (int c = 0; c < 3; c++) {
      params[bond_type].Mr[r][c]     = params[bond_type].Mmat[r][c];
      params[bond_type].Mtr_tr[r][c] = params[bond_type].Mmat[r][c + 3];
      params[bond_type].Mtr_bl[r][c] = params[bond_type].Mmat[r + 3][c];
      params[bond_type].Mt[r][c]     = params[bond_type].Mmat[r + 3][c + 3];
    }
  }
}
