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
   Contributing author: Enrico Skoruppa (Physics of Life, TU Dresden, Dresden)
------------------------------------------------------------------------- */

#include "bond_rbp.h"
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

BondRBP::BondRBP(LAMMPS *lmp) : Bond(lmp) {
  params = nullptr;
}

/* ----------------------------------------------------------------------
   Destructor: deallocates dynamically allocated memory
------------------------------------------------------------------------- */

BondRBP::~BondRBP() {
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(params);
    allocated = 0;
  }
}

/* ----------------------------------------------------------------------
   Compute forces and energy contribution for all bonds of this style
------------------------------------------------------------------------- */

void BondRBP::compute(int eflag, int vflag) {

  int id1, id2, bond_type;
  
  double T1[3][3],T1tp[3][3],T2[3][3];
  double r1[3], r2[3], dr[3];
  double R[3][3];
  double Om[3], w[3];
  double Omd[3], wd[3]; 
  double Jinvtp[3][3];

  double gamma[3],mu[3];

  double tau1[3],tau2[3];
  double tau1lab[3],tau2lab[3];
  double Flab[3];
  double tau_pure[3];
  double levertorque[3];
  double tmp[3];

  double **x = atom->x;                
  double **f = atom->f;             
  double **torque = atom->torque;   
  double *quat1,*quat2;
  
  auto avec = dynamic_cast<AtomVecEllipsoid *>(atom->style_match("ellipsoid"));
  AtomVecEllipsoid::Bonus *bonus = avec->bonus;
  int *ellipsoid = atom->ellipsoid;

  ev_init(eflag, vflag);

  // Access bond list
  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;
 
  for (int bid = 0; bid < nbondlist; bid++) {
    id1 = bondlist[bid][0];
    id2 = bondlist[bid][1];
    bond_type = bondlist[bid][2];

    // get quaternions
    quat1=bonus[ellipsoid[id1]].quat;
    quat2=bonus[ellipsoid[id2]].quat;

    // get positions
    r1[0] = x[id1][0];
    r1[1] = x[id1][1];
    r1[2] = x[id1][2];
    r2[0] = x[id2][0];
    r2[1] = x[id2][1];
    r2[2] = x[id2][2];

    // transform quat to triads [SO(3)]
    MathExtra::quat_to_mat(quat1, T1);
    MathExtra::quat_to_mat(quat2, T2);

    // compute T1^T
    lamath::transpose(T1,T1tp);

    // compute R
    lamath::mul(T1tp,T2,R);
    
    // compute Omega
    so3::rotmat2euler(R,Om);
    
    // compute w
    dr[0] = r2[0] - r1[0]; 
    dr[1] = r2[1] - r1[1]; 
    dr[2] = r2[2] - r1[2]; 

    // check for domain mismatch of special neighbor
    if (domain->minimum_image_check(dr[0],dr[1],dr[2])) {
      domain->minimum_image(dr[0],dr[1],dr[2]);
    }

    lamath::mul(T1tp,dr,w);

    // compute Omega_Delta
    Omd[0] = Om[0] - params[bond_type].srot[0];
    Omd[1] = Om[1] - params[bond_type].srot[1];
    Omd[2] = Om[2] - params[bond_type].srot[2];

    // compute w_Delta
    wd[0] = w[0] - params[bond_type].svec[0];
    wd[1] = w[1] - params[bond_type].svec[1];
    wd[2] = w[2] - params[bond_type].svec[2];
    
    // compute transposed inverse left Jacobian
    so3::leftJacobianInverseTransposed(Om,Jinvtp);
    
    // double tmpmat[3][3];
    // so3::leftJacobianInverseTransposed(Om,tmpmat);
    // lamath::transpose(tmpmat,Jinvtp);
    
    // compute gamma
    lamath::mul(params[bond_type].Mr,Omd,gamma);
    lamath::mul(params[bond_type].Mtr_tr,wd,tmp);
    gamma[0] += tmp[0];
    gamma[1] += tmp[1];
    gamma[2] += tmp[2];
    
    // compute mu
    lamath::mul(params[bond_type].Mtr_bl,Omd,mu);
    lamath::mul(params[bond_type].Mt,wd,tmp);
    mu[0] += tmp[0];
    mu[1] += tmp[1];
    mu[2] += tmp[2];

    // compute the torque component in frame T1 without the lever arm torque
    lamath::mul(Jinvtp,gamma,tau_pure);

    // compute Force in lab frame acting on T1 (minus for T2)
    lamath::mul(T1,mu,Flab);

    // compute torque on T1 in lab frame
    MathExtra::cross3(w,mu,levertorque);

    tmp[0] = tau_pure[0] + levertorque[0];
    tmp[1] = tau_pure[1] + levertorque[1];
    tmp[2] = tau_pure[2] + levertorque[2];

    lamath::mul(T1,tmp,tau1lab);
    lamath::mul(T1,tau_pure,tau2lab);

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
    //   }
    // ///////////////////////////////////////////////////////////////
    // ///////////////////////////////////////////////////////////////

    ////////////////////////////////////////////
    ////////////////////////////////////////////

    // apply force and torque to each of 2 atoms
    if (newton_bond || id1 < nlocal) {

      f[id1][0] += Flab[0];
      f[id1][1] += Flab[1];
      f[id1][2] += Flab[2];

      torque[id1][0] += tau1lab[0];
      torque[id1][1] += tau1lab[1];
      torque[id1][2] += tau1lab[2];
      
    }

    if (newton_bond || id2 < nlocal) {

      f[id2][0] -= Flab[0];
      f[id2][1] -= Flab[1];
      f[id2][2] -= Flab[2];

      torque[id2][0] -= tau2lab[0];
      torque[id2][1] -= tau2lab[1];
      torque[id2][2] -= tau2lab[2];
    }

    if (evflag) {

      double bond_energy = 0.0;
      double Y_vec[6]; 

      Y_vec[0] = Omd[0]; 
      Y_vec[1] = Omd[1];
      Y_vec[2] = Omd[2];
      Y_vec[3] = wd[0];   
      Y_vec[4] = wd[1];
      Y_vec[5] = wd[2];

      // Calculate energy: 0.5 * Y_vec^T * Mmat * Y_vec
      for (int i_mat = 0; i_mat < 6; ++i_mat) { 
        for (int j_mat = 0; j_mat < 6; ++j_mat) {
          bond_energy += Y_vec[i_mat] * params[bond_type].Mmat[i_mat][j_mat] * Y_vec[j_mat];
        }
      }
      bond_energy *= 0.5;

      // For virial, delx,dely,delz should be vector from j to i (id2 to id1)
      double vir_delx = x[id1][0] - x[id2][0];
      double vir_dely = x[id1][1] - x[id2][1];
      double vir_delz = x[id1][2] - x[id2][2];

      // Call the base class's ev_tally_xyz function
      // Force on atom i (id1) from atom j (id2) is -Flab
      ev_tally_xyz(id1, id2, nlocal, newton_bond, bond_energy, 
                    Flab[0], Flab[1], Flab[2], // Force on id1 from id2
                    vir_delx, vir_dely, vir_delz);
    }
  }
}


/* ----------------------------------------------------------------------
   Compute single
------------------------------------------------------------------------- */

double BondRBP::single(int i, double rsq, int itype, int jtype, double &fforce) {
  fforce = 0.0;
  return 0.0;
}


void BondRBP::coeff(int narg, char **arg)
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
    if (db.metadata().bond_style != RBPBOND_STYLE) {
      error->all(FLERR, 
        "Bond style mismatch: database specifies '" + db.metadata().bond_style + 
        "' but using bond_style rbp");
    }

    int dbid = utils::inumeric(FLERR, arg[3], false, lmp);
    for (int bond_type=ilo;bond_type<=ihi;bond_type++) {
      assign_coeffs(bond_type,db.bond(dbid++).coeffs);
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
    assign_coeffs(bond_type,coeffs);
  }
}


/* ----------------------------------------------------------------------
   Initialization: verify atom properties and prepare style
------------------------------------------------------------------------- */

void BondRBP::init_style() {

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

double BondRBP::equilibrium_distance(int bond_type) 
{
  return params[bond_type].equidist;
}

/* ----------------------------------------------------------------------
   Write bond style parameters to restart file
------------------------------------------------------------------------- */

void BondRBP::write_restart(FILE *fp) {

  fwrite(&params[1], sizeof(RBPParams), atom->nbondtypes, fp);

}

/* ----------------------------------------------------------------------
   Read bond style parameters from restart file
------------------------------------------------------------------------- */

void BondRBP::read_restart(FILE *fp) {

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

void BondRBP::write_data(FILE *fp) {

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

void BondRBP::allocate() {
  allocated = 1;
  int n = atom->nbondtypes;

  memory->create(params, n + 1, "bond:rbp:params");
  memory->create(setflag, n+1, "bond:rbp:setflag");
  for (int i = 0; i <= n; i++) setflag[i] = 0;
}


void BondRBP::assign_coeffs(int bond_type, const std::vector<double> &args) {
  // ------------------------------------------------------------
  // Inline numeric definitions (legacy / TWLC / full RBP)
  // Supported forms:
  //   6 + 6  = 12 args (diagonal M)
  //   6 + 12 = 18 args (block-diagonal M)
  //   6 + 21 = 27 args (full upper-triangular M)
  //
  // these are the arguments following the leading bond_type id
  // ------------------------------------------------------------
  constexpr int opt1 = 12;
  constexpr int opt2 = 18;
  constexpr int opt3 = 27;

  int narg = static_cast<int>(args.size());

  if (narg != opt1 && narg != opt2 && narg != opt3) {
    std::string msg =
      "Incorrect number of arguments for bond_coeff (" + std::to_string(narg) + ")\n" +
      "Expected:\n" +
      "  6+6  (6 groundstate and 6 diagonal components of M)\n" +
      "  6+12 (6 groundstate and block-diagonal M: Upper triangular assignment for each 3x3 block. Assignment by first iterating through the rows)\n" +
      "  6+21 (6 groundstate and full symmetric M assignment. Assignment by first iterating through the rows)";
    error->all(FLERR, msg.c_str());
  }

  // -------------------------------------------------------------------
  // read groundstate
  for (int i = 0; i < 6; i++)
    params[bond_type].Ystatic[i] = args[i];

  set_static_(bond_type);
  set_equidist(bond_type);
  zero_stiffmat_(bond_type);

  // -------------------------------------------------------------------
  // set diagnoal only
  if (narg == opt1) {
    for (int i = 0; i < 6; i++)
      params[bond_type].Mmat[i][i] = args[6+i];
  }

  // -------------------------------------------------------------------
  // set rotation and translation blocks independently (block diagonal)
  if (narg == opt2) {
    int k = 6;
    params[bond_type].Mmat[0][0] = args[k++];
    params[bond_type].Mmat[0][1] = args[k++];
    params[bond_type].Mmat[0][2] = args[k++];
    params[bond_type].Mmat[1][1] = args[k++];
    params[bond_type].Mmat[1][2] = args[k++];
    params[bond_type].Mmat[2][2] = args[k++];
    params[bond_type].Mmat[3][3] = args[k++];
    params[bond_type].Mmat[3][4] = args[k++];
    params[bond_type].Mmat[3][5] = args[k++];
    params[bond_type].Mmat[4][4] = args[k++];
    params[bond_type].Mmat[4][5] = args[k++];
    params[bond_type].Mmat[5][5] = args[k++];
  }

  // -------------------------------------------------------------------
  // set everything by assigning upper triangular coefficients
  if (narg == opt3) {
    int k = 6;
    for (int i = 0; i < 6; i++)
      for (int j = i; j < 6; j++)
        params[bond_type].Mmat[i][j] = args[k++];
  }

  set_lower_triangle_(bond_type);
  assign_blocks_(bond_type);

  if (!lamath::is_positive_definite(params[bond_type].Mmat)) {
    error->all(FLERR, "Stiffness matrix M is not positive definite");
  }
  
  setflag[bond_type] = 1;
}

void BondRBP::set_static_(int bond_type) {
  // set rotation and translation components and static rotation matrix
  for (int i = 0; i < 3; i++) {
    params[bond_type].srot[i] = params[bond_type].Ystatic[i];
    params[bond_type].svec[i] = params[bond_type].Ystatic[i + 3];
  }
  so3::euler2rotmat(params[bond_type].srot, params[bond_type].Smat);
}

void BondRBP::set_equidist(int bond_type) {
  params[bond_type].equidist = lamath::norm(params[bond_type].svec);
}

void BondRBP::zero_stiffmat_(int bond_type) {
  // initialize stiffness matrix to zero
  for (int i = 0; i < 6; i++)
    for (int j = 0; j < 6; j++)
      params[bond_type].Mmat[i][j] = 0.0;
}

void BondRBP::set_lower_triangle_(int bond_type) {
  // symmetrize by setting lower triangular components
  for (int i = 0; i < 6; i++)
    for (int j = i + 1; j < 6; j++)
      params[bond_type].Mmat[j][i] = params[bond_type].Mmat[i][j];
}

void BondRBP::assign_blocks_(int bond_type) {
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
