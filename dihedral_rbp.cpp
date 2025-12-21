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

#include "dihedral_rbp.h"
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


using namespace LAMMPS_NS;

/* ----------------------------------------------------------------------
   Constructor: initializes internal data structures and sets defaults
------------------------------------------------------------------------- */

DihedralRBP::DihedralRBP(LAMMPS *lmp) : Dihedral(lmp) {
  params = nullptr;
}

/* ----------------------------------------------------------------------
   Destructor: deallocates dynamically allocated memory
------------------------------------------------------------------------- */

DihedralRBP::~DihedralRBP() {
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(params);
    allocated = 0;
  }
}

/* ----------------------------------------------------------------------
   Compute forces and energy contribution for all dihedrals of this style
------------------------------------------------------------------------- */


void DihedralRBP::compute(int eflag, int vflag) {


  int id1,id2,id3,id4;
  int dihedral_type;
  
  int **dihedrallist = neighbor->dihedrallist;
  int ndihedrallist  = neighbor->ndihedrallist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;
  
  auto avec = dynamic_cast<AtomVecEllipsoid *>(atom->style_match("ellipsoid"));
  AtomVecEllipsoid::Bonus *bonus = avec->bonus;
  int *ellipsoid = atom->ellipsoid;
  
  double **x = atom->x;                
  double **f = atom->f;             
  double **torque = atom->torque;   

  // Rotation variables
  double *quat1,*quat2,*quat3,*quat4;
  double T1[3][3],T2[3][3],T3[3][3],T4[3][3];
  double T1tp[3][3],T3tp[3][3];

  double R_a[3][3],R_b[3][3];
  double Om_a[3],Om_b[3];
  double Omd_a[3],Omd_b[3];
  double Jinvtp_a[3][3];
  double Jinvtp_b[3][3];

  // Translation variables
  double r1[3],r2[3],r3[3],r4[3];
  double dr_a[3],dr_b[3];
  double w_a[3],w_b[3];
  double wd_a[3],wd_b[3];

  // temp
  double tmp1[3],tmp2[3],tmp3[3];

  // wrench components
  double gamma_a[3];
  double gamma_b[3];
  double mu_a[3];
  double mu_b[3];

  // torques
  double tau1[3],tau2[3],tau3[3],tau4[3];
  // forces
  double f1[3],f2[3],f3[3],f4[3];
  
  ev_init(eflag, vflag);
  
  for (int did = 0; did < ndihedrallist; did++) {
    // unpack dihedral endpoints and type
    id1 = dihedrallist[did][0];
    id2 = dihedrallist[did][1];
    id3 = dihedrallist[did][2];
    id4 = dihedrallist[did][3];
    dihedral_type = dihedrallist[did][4];

    // fprintf(screen,"DihedralRBP %d %d %d\n",id1,id2,id3);
    
    ////////////////////////////////////////
    // Rotational Part
    // get quaternions
    quat1=bonus[ellipsoid[id1]].quat;
    quat2=bonus[ellipsoid[id2]].quat;
    quat3=bonus[ellipsoid[id3]].quat;
    quat4=bonus[ellipsoid[id4]].quat;
    
    // transform quat to triads [SO(3)]
    MathExtra::quat_to_mat(quat1, T1);
    MathExtra::quat_to_mat(quat2, T2);
    MathExtra::quat_to_mat(quat3, T3);
    MathExtra::quat_to_mat(quat4, T4);
    
    // compute transposed
    lamath::transpose(T1,T1tp);
    lamath::transpose(T3,T3tp);
    
    // compute R
    lamath::mul(T1tp,T2,R_a);
    lamath::mul(T3tp,T4,R_b);
    
    // compute Omega_1
    so3::rotmat2euler(R_a,Om_a);
    so3::rotmat2euler(R_b,Om_b);
    
    // compute Omega_Delta
    Omd_a[0] = Om_a[0] - params[dihedral_type].srot1[0];
    Omd_a[1] = Om_a[1] - params[dihedral_type].srot1[1];
    Omd_a[2] = Om_a[2] - params[dihedral_type].srot1[2];
    Omd_b[0] = Om_b[0] - params[dihedral_type].srot2[0];
    Omd_b[1] = Om_b[1] - params[dihedral_type].srot2[1];
    Omd_b[2] = Om_b[2] - params[dihedral_type].srot2[2];
    
    // compute transposed inverse left Jacobian
    so3::leftJacobianInverseTransposed(Om_a,Jinvtp_a);
    so3::leftJacobianInverseTransposed(Om_b,Jinvtp_b);
    
    ////////////////////////////////////////
    // Translational Part
    // get positions
    r1[0] = x[id1][0];
    r1[1] = x[id1][1];
    r1[2] = x[id1][2];
    r2[0] = x[id2][0];
    r2[1] = x[id2][1];
    r2[2] = x[id2][2];
    r3[0] = x[id3][0];
    r3[1] = x[id3][1];
    r3[2] = x[id3][2];
    r4[0] = x[id4][0];
    r4[1] = x[id4][1];
    r4[2] = x[id4][2];

    // compute w
    dr_a[0] = r2[0] - r1[0]; 
    dr_a[1] = r2[1] - r1[1]; 
    dr_a[2] = r2[2] - r1[2]; 
    dr_b[0] = r4[0] - r3[0]; 
    dr_b[1] = r4[1] - r3[1]; 
    dr_b[2] = r4[2] - r3[2]; 


    // check for domain mismatch of special neighbor
    if (domain->minimum_image_check(dr_a[0],dr_a[1],dr_a[2])) {

      fprintf(screen," dr_a too large!\n");
      std::exit(0);
      
      domain->minimum_image(dr_a[0],dr_a[1],dr_a[2]);
    }
    if (domain->minimum_image_check(dr_b[0],dr_b[1],dr_b[2])) {
      fprintf(screen," dr_a too large!\n");
      std::exit(0);
      domain->minimum_image(dr_b[0],dr_b[1],dr_b[2]);
    }
    
    // fprintf(screen," Terminated\n");
    // std::exit(0);

    lamath::mul(T1tp,dr_a,w_a);
    lamath::mul(T3tp,dr_b,w_b);
  
    // compute w_Delta
    wd_a[0] = w_a[0] - params[dihedral_type].svec1[0];
    wd_a[1] = w_a[1] - params[dihedral_type].svec1[1];
    wd_a[2] = w_a[2] - params[dihedral_type].svec1[2];
    wd_b[0] = w_b[0] - params[dihedral_type].svec2[0];
    wd_b[1] = w_b[1] - params[dihedral_type].svec2[1];
    wd_b[2] = w_b[2] - params[dihedral_type].svec2[2];

    //////////////////////////////////////////////////////
    // Compute gamma and mu

    // gamma1
    lamath::mul(params[dihedral_type].Mrr,Omd_b,tmp1);
    lamath::mul(params[dihedral_type].Mrt,wd_b,tmp2);
    lamath::add(tmp1,tmp2,gamma_a);
    // gamma1[0] = tmp1[0] + tmp2[0];
    // gamma1[1] = tmp1[1] + tmp2[1];
    // gamma1[2] = tmp1[2] + tmp2[2];

    // mu1
    lamath::mul(params[dihedral_type].Mtr,Omd_b,tmp1);
    lamath::mul(params[dihedral_type].Mtt,wd_b,tmp2);
    lamath::add(tmp1,tmp2,mu_a);
    // mu1[0] = tmp1[0] + tmp2[0];
    // mu1[1] = tmp1[1] + tmp2[1];
    // mu1[2] = tmp1[2] + tmp2[2];

    // gamma2
    lamath::mul(params[dihedral_type].Mrr_tp,Omd_a,tmp1);
    lamath::mul(params[dihedral_type].Mtr_tp,wd_a,tmp2);
    lamath::add(tmp1,tmp2,gamma_b);
    // gamma2[0] = tmp1[0] + tmp2[0];
    // gamma2[1] = tmp1[1] + tmp2[1];
    // gamma2[2] = tmp1[2] + tmp2[2];

    // mu1
    lamath::mul(params[dihedral_type].Mrt_tp,Omd_a,tmp1);
    lamath::mul(params[dihedral_type].Mtt_tp,wd_a,tmp2);
    lamath::add(tmp1,tmp2,mu_b);
    // mu2[0] = tmp1[0] + tmp2[0];
    // mu2[1] = tmp1[1] + tmp2[1];
    // mu2[2] = tmp1[2] + tmp2[2];

    // Triad 1 ////////
    // tau1 
    MathExtra::cross3(w_a,mu_a,tmp1);
    lamath::mul(Jinvtp_a,gamma_a,tmp2);
    lamath::add(tmp1,tmp2,tmp3);
    lamath::mul(T1,tmp3,tau1);
    // f1
    lamath::mul(T1,mu_a,f1);
    
    // Triad 2 ////////
    // tau2
    lamath::mul(T1,tmp2,tmp3);
    lamath::signflip(tmp3,tau2);
    
    // f2
    lamath::mul(T1,mu_a,tmp3);
    lamath::signflip(tmp3,f2);
    
    // Triad 3 ////////
    // tau3 
    MathExtra::cross3(w_b,mu_b,tmp1);
    lamath::mul(Jinvtp_b,gamma_b,tmp2);
    lamath::add(tmp1,tmp2,tmp3);
    lamath::mul(T3,tmp3,tau3);
    // f3
    lamath::mul(T3,mu_b,f3);

    // Triad 4 ////////
    // tau4
    lamath::mul(T3,tmp2,tmp3);
    lamath::signflip(tmp3,tau4);
    
    // f4
    lamath::mul(T3,mu_b,tmp3);
    lamath::signflip(tmp3,f4);
    
    ////////////////////////////////////////////
    ////////////////////////////////////////////

    // apply force and torque to each of 4 atoms
    if (newton_bond || id1 < nlocal) {

      f[id1][0] += f1[0];
      f[id1][1] += f1[1];
      f[id1][2] += f1[2];

      torque[id1][0] += tau1[0];
      torque[id1][1] += tau1[1];
      torque[id1][2] += tau1[2];
    }

    if (newton_bond || id2 < nlocal) {

      f[id2][0] += f2[0];
      f[id2][1] += f2[1];
      f[id2][2] += f2[2];

      torque[id2][0] += tau2[0];
      torque[id2][1] += tau2[1];
      torque[id2][2] += tau2[2];
    }

    if (newton_bond || id3 < nlocal) {

      f[id3][0] += f3[0];
      f[id3][1] += f3[1];
      f[id3][2] += f3[2];

      torque[id3][0] += tau3[0];
      torque[id3][1] += tau3[1];
      torque[id3][2] += tau3[2];
      
    }

    if (newton_bond || id4 < nlocal) {

      f[id4][0] += f4[0];
      f[id4][1] += f4[1];
      f[id4][2] += f4[2];

      torque[id4][0] += tau4[0];
      torque[id4][1] += tau4[1];
      torque[id4][2] += tau4[2];
      
    }

    if (evflag) {

      double dihedral_energy = 0.0;
      double Yd1[6]; 
      double Yd2[6]; 

      Yd1[0] = Omd_a[0]; 
      Yd1[1] = Omd_a[1];
      Yd1[2] = Omd_a[2];
      Yd1[3] = wd_a[0];   
      Yd1[4] = wd_a[1];
      Yd1[5] = wd_a[2];

      Yd2[0] = Omd_b[0]; 
      Yd2[1] = Omd_b[1];
      Yd2[2] = Omd_b[2];
      Yd2[3] = wd_b[0];   
      Yd2[4] = wd_b[1];
      Yd2[5] = wd_b[2];

      // Calculate energy: 0.5 * Y_vec^T * Mmat * Y_vec
      for (int i_mat = 0; i_mat < 6; ++i_mat) { 
        for (int j_mat = 0; j_mat < 6; ++j_mat) {
          dihedral_energy += Yd1[i_mat] * params[dihedral_type].Mmat[i_mat][j_mat] * Yd2[j_mat];
        }
      }
      dihedral_energy *= 0.5;

      // bond vectors for virial:
      // vb1: 1 -> 2
      // vb2: 3 -> 2
      // vb3: 4 -> 3
      double vb1x = x[id1][0] - x[id2][0];
      double vb1y = x[id1][1] - x[id2][1];
      double vb1z = x[id1][2] - x[id2][2];

      double vb2x = x[id3][0] - x[id2][0];
      double vb2y = x[id3][1] - x[id2][1];
      double vb2z = x[id3][2] - x[id2][2];

      double vb3x = x[id4][0] - x[id3][0];
      double vb3y = x[id4][1] - x[id3][1];
      double vb3z = x[id4][2] - x[id3][2];

      // Tally energy + virial.
      // IMPORTANT: f1, f3, f4 must be the *net* forces on atoms id1, id3, id4
      // from THIS dihedral.
      ev_tally(id1, id2, id3, id4,
               nlocal, newton_bond,
               dihedral_energy,
               f1,   // force on atom id1
               f3,   // force on atom id3
               f4,   // force on atom id4
               vb1x, vb1y, vb1z,
               vb2x, vb2y, vb2z,
               vb3x, vb3y, vb3z);
    }

  }
}


/* ----------------------------------------------------------------------
   Parse coefficients for each dihedral type
------------------------------------------------------------------------- */

void DihedralRBP::coeff(int narg, char **arg) {

  // check conditions for properly formated arg !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // may not be necessary

  if (!allocated) allocate();

  int ilo, ihi;
  utils::bounds(FLERR, arg[0], 1, atom->ndihedraltypes, ilo, ihi, error);

  // load from database file
  if (narg == 4 && strcmp(arg[1], "dbfile") == 0) {

    RBPDatabase db(lmp, error);
    db.read(arg[2]);

    int dbid = utils::inumeric(FLERR, arg[3], false, lmp);
    for (int dihedral_type=ilo;dihedral_type<=ihi;dihedral_type++) {
      assign_coeffs(dihedral_type,db.dihedral(dbid++).coeffs);
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
  for (int dihedral_type=ilo;dihedral_type<=ihi;dihedral_type++) {
    assign_coeffs(dihedral_type,coeffs);
  }
}

/* ----------------------------------------------------------------------
   Initialization: verify atom properties and prepare style
------------------------------------------------------------------------- */

void DihedralRBP::init_style() {

  // Ensure ellipsoidal atoms are being used
  if (!atom->ellipsoid_flag)
    error->all(FLERR, "Dihedral style rbp requires atom style with ellipsoids");

  if (domain->dimension != 3)
    error->all(FLERR, "Dihedral style rbp requires a 3D simulation");

  // for (int i = 1; i <= atom->ndihedraltypes; i++) {
  //   if (!setflag[i])
  //     error->all(FLERR, "Not all dihedral coefficients are set for dihedral style rbp");
  // }
}


/* ----------------------------------------------------------------------
   Write dihedral style parameters to restart file
------------------------------------------------------------------------- */

void DihedralRBP::write_restart(FILE *fp) {
  // - Loop over ndihedraltypes
  // - Write 36 values for Mmat[type][6][6]
  // - Write 12 values for ss[type] (3x3 rotation + 3 translation)
  // - Use fwrite for binary output
}

/* ----------------------------------------------------------------------
   Read dihedral style parameters from restart file
------------------------------------------------------------------------- */

void DihedralRBP::read_restart(FILE *fp) {
  // - Allocate Mmat and ss if not already done
  // - Read 36 values for Mmat[type]
  // - Read 12 values for ss[type]
  // - Optional: validate input and sanity-check matrix
}

/* ----------------------------------------------------------------------
   Write dihedral style data (coefficients) to LAMMPS data file
------------------------------------------------------------------------- */

void DihedralRBP::write_data(FILE *fp) {
  // - Loop over all dihedral types
  // - Write dihedral_coeff line: 21 M values + 6 SE(3) values
  // - Match the format expected by coeff()
  // - Use fprintf for clean output formatting
}

/* ----------------------------------------------------------------------
   Allocate memory for internal data structures  //   // compute gamma
  //   lamath::mul(params[dihedral_type].Mr,Omd,gamma);
  //   lamath::mul(params[dihedral_type].Mtr_tr,wd,tmp);
  //   gamma[0] += tmp[0];
  //   gamma[1] += tmp[1];
  //   gamma[2] += tmp[2];
    
  //   // compute mu
  //   lamath::mul(params[dihedral_type].Mtr_bl,Omd,mu);
  //   lamath::mul(params[dihedral_type].Mt,wd,tmp);
  //   mu[0] += tmp[0];
  //   mu[1] += tmp[1];
  //   mu[2] += tmp[2];
------------------------------------------------------------------------- */

void DihedralRBP::allocate() {
  // - Allocate Mmat[ndihedraltypes+1][6][6]
  // - Allocate ss[ndihedraltypes+1]
  // - Initialize values to zero or identity where appropriate
  // - Use LAMMPS memory->create interface

  allocated = 1;
  int n = atom->ndihedraltypes;

  memory->create(params,  n+1, "dihedral:rbp:params");
  memory->create(setflag, n+1, "dihedral:rbp:setflag");
  for (int i = 0; i <= n; i++) setflag[i] = 0;
}


void DihedralRBP::assign_coeffs(int dihedral_type, const std::vector<double> &args) {


  if (args.size() != 48) {
    error->all(FLERR, "Invalid number of coefficients found for dihedral style rbp. Requires 48 coefficients: X0_1 (6) X0_2 (6) stiffmat (6x6).");
  }

  // assign Ystatic
  for (int i=0;i<6;i++) {
    params[dihedral_type].Ystatic1[i] = args[i];
    params[dihedral_type].Ystatic2[i] = args[i+6];
  }

  // assign partial static
  for (int ii=0;ii<3;ii++) {
    params[dihedral_type].srot1[ii] = params[dihedral_type].Ystatic1[ii];
    params[dihedral_type].svec1[ii] = params[dihedral_type].Ystatic1[ii+3];
    params[dihedral_type].srot2[ii] = params[dihedral_type].Ystatic2[ii];
    params[dihedral_type].svec2[ii] = params[dihedral_type].Ystatic2[ii+3];
  }

  // assigne Smat
  so3::euler2rotmat(params[dihedral_type].srot1,params[dihedral_type].Smat1);
  so3::euler2rotmat(params[dihedral_type].srot2,params[dihedral_type].Smat2);

  int argid = 12;
  for (int ii=0;ii<6;ii++) {
    for (int jj=0;jj<6;jj++) {
      params[dihedral_type].Mmat[ii][jj] = args[argid++];
    }
  }

  for (int r = 0; r < 3; r++) {
    for (int c = 0; c < 3; c++) {
      params[dihedral_type].Mrr[r][c] = params[dihedral_type].Mmat[r][c];
      params[dihedral_type].Mrt[r][c] = params[dihedral_type].Mmat[r][c+3];
      params[dihedral_type].Mtr[r][c] = params[dihedral_type].Mmat[r+3][c];
      params[dihedral_type].Mtt[r][c] = params[dihedral_type].Mmat[r+3][c+3];
    }
  }

  lamath::transpose(params[dihedral_type].Mrr,params[dihedral_type].Mrr_tp);
  lamath::transpose(params[dihedral_type].Mtt,params[dihedral_type].Mtt_tp);
  lamath::transpose(params[dihedral_type].Mtr,params[dihedral_type].Mtr_tp);
  lamath::transpose(params[dihedral_type].Mrt,params[dihedral_type].Mrt_tp);

  // THESE OFF-DIAGONAL COMPONENTS NEED NOT NECESSARILY BE POSITIVE DEFINITE!!!
  // // check if matrix is positive definite
  // if (!lamath::is_positive_definite(params[dihedral_type].Mmat)) {
  //   error->all(FLERR, "Stiffness matrix M is not positive definite.");
  // }

  setflag[dihedral_type] = 1;
}