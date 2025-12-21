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

#include "atom.h"
#include "atom_vec_ellipsoid.h"
#include "bond_rbp_fene.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "math_extra.h"
#include "memory.h"
#include "neighbor.h"
#include "update.h"
#include "so3.h"
#include "lamath.h"


using namespace LAMMPS_NS;
using MathConst::MY_CUBEROOT2;

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

  // FENE vars
  double delx, dely, delz, ebond, fbond;
  double rsq, r0sq, rlogarg, sr2, sr6;
  ebond = 0.0;

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

  // Loop over all bonds
  for (int bid = 0; bid < nbondlist; bid++) {
    id1 = bondlist[bid][0];
    id2 = bondlist[bid][1];
    bond_type = bondlist[bid][2];

    ///////////////////////////////////////////////
    // FENE Force /////////////////////////////////
    delx = x[id1][0] - x[id2][0];
    dely = x[id1][1] - x[id2][1];
    delz = x[id1][2] - x[id2][2];

    // force from log term

    rsq = delx * delx + dely * dely + delz * delz;
    r0sq = params[bond_type].r0 * params[bond_type].r0;
    rlogarg = 1.0 - rsq / r0sq;

    // if r -> r0, then rlogarg < 0.0 which is an error
    // issue a warning and reset rlogarg = epsilon
    // if r > 2*r0 something serious is wrong, abort

    if (rlogarg < 0.1) {
      error->warning(FLERR, "FENE bond too long: {} {} {} {}", update->ntimestep, atom->tag[id1],
                     atom->tag[id2], sqrt(rsq));
      if (rlogarg <= -3.0) error->one(FLERR, "Bad FENE bond");
      rlogarg = 0.1;
    }
    fbond = -params[bond_type].K / rlogarg;

    // force from LJ term

    if (rsq < MY_CUBEROOT2 * params[bond_type].sigma * params[bond_type].sigma) {
      sr2 = params[bond_type].sigma * params[bond_type].sigma / rsq;
      sr6 = sr2 * sr2 * sr2;
      fbond += 48.0 * params[bond_type].epsilon * sr6 * (sr6 - 0.5) / rsq;
    }

    ///////////////////////////////////////////////
    // FBP FORCES AND TORQUES /////////////////////

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
      
      // Fene
      f[id1][0] += delx * fbond;
      f[id1][1] += dely * fbond;
      f[id1][2] += delz * fbond;
    }

    if (newton_bond || id2 < nlocal) {

      f[id2][0] -= Flab[0];
      f[id2][1] -= Flab[1];
      f[id2][2] -= Flab[2];

      torque[id2][0] -= tau2lab[0];
      torque[id2][1] -= tau2lab[1];
      torque[id2][2] -= tau2lab[2];

      // Fene
      f[id2][0] -= delx * fbond;
      f[id2][1] -= dely * fbond;
      f[id2][2] -= delz * fbond;

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

      // Fene Energy
      bond_energy += -0.5 * params[bond_type].K * r0sq * log(rlogarg);

      // For virial, delx,dely,delz should be vector from j to i (id2 to id1)
      double vir_delx = x[id1][0] - x[id2][0];
      double vir_dely = x[id1][1] - x[id2][1];
      double vir_delz = x[id1][2] - x[id2][2];

      // Call the base class's ev_tally_xyz function
      // Force on atom i (id1) from atom j (id2) is -Flab
      ev_tally_xyz(id1, id2, nlocal, newton_bond, bond_energy, 
                    // -Flab[0], -Flab[1], -Flab[2], // Force on id1 from id2
                    Flab[0], Flab[1], Flab[2], // Force on id1 from id2
                    vir_delx, vir_dely, vir_delz);
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

/* ----------------------------------------------------------------------
   Parse coefficients for each bond type
------------------------------------------------------------------------- */

void BondRBPFene::coeff(int narg, char **arg) {

  int opt1,opt2,opt3;
  opt1 = 17;
  opt2 = 23;
  opt3 = 32;

  if (narg != opt1 && narg != opt2 && narg != opt3) {
    error->all(FLERR, "Incorrect number of arguments for bond_coeff (expected 1 + 2 + 6 + 6/12/21 args)");
  }

  int bond_type = utils::inumeric(FLERR, arg[0], false, lmp);
  if (bond_type < 1 || bond_type > atom->nbondtypes )
    error->all(FLERR, "Invalid bond type index in bond_coeff");

  // allocate memory
  if (!allocated) allocate();

  int add = 1;

  // Init Fene parameters ///////////////
  params[bond_type].K       = utils::numeric(FLERR, arg[add+0], false, lmp);
  params[bond_type].r0      = utils::numeric(FLERR, arg[add+1], false, lmp);
  params[bond_type].epsilon = utils::numeric(FLERR, arg[add+2], false, lmp);
  params[bond_type].sigma   = utils::numeric(FLERR, arg[add+3], false, lmp);
  add += 4;

  ///////////////////////////////////////

  // assign Ystatic
  for (int dd=0;dd<6;dd++) {
    params[bond_type].Ystatic[dd] = utils::numeric(FLERR, arg[add+dd], false, lmp);
  }
  add += 6;

  // assign partial static
  for (int ii=0;ii<3;ii++) {
    params[bond_type].srot[ii] = params[bond_type].Ystatic[ii];
    params[bond_type].svec[ii] = params[bond_type].Ystatic[ii+3];
  }

  // assigne Smat
  so3::euler2rotmat(params[bond_type].srot,params[bond_type].Smat);

  // set equilibrium distance
  params[bond_type].equidist = lamath::norm(params[bond_type].svec);

  // initiate M to zero
  for (int ii = 0; ii < 6; ii++) {
    for (int jj = 0; jj < 6; jj++) {
      params[bond_type].Mmat[ii][jj] = 0.0;
    }
  }
  
  if (narg == opt1) {
    for (int cr = 0; cr < 6; cr++) {
      params[bond_type].Mmat[cr][cr] = utils::numeric(FLERR, arg[add+cr], false, lmp);
    }
    add += 6;
  }

  if (narg == opt2) {

    params[bond_type].Mmat[0][0] = utils::numeric(FLERR, arg[add+0], false, lmp);
    params[bond_type].Mmat[0][1] = utils::numeric(FLERR, arg[add+1], false, lmp);
    params[bond_type].Mmat[0][2] = utils::numeric(FLERR, arg[add+2], false, lmp);
    params[bond_type].Mmat[1][1] = utils::numeric(FLERR, arg[add+3], false, lmp);
    params[bond_type].Mmat[1][2] = utils::numeric(FLERR, arg[add+4], false, lmp);
    params[bond_type].Mmat[2][2] = utils::numeric(FLERR, arg[add+5], false, lmp);
    params[bond_type].Mmat[3][3] = utils::numeric(FLERR, arg[add+6], false, lmp);
    params[bond_type].Mmat[3][4] = utils::numeric(FLERR, arg[add+7], false, lmp);
    params[bond_type].Mmat[3][5] = utils::numeric(FLERR, arg[add+8], false, lmp);
    params[bond_type].Mmat[4][4] = utils::numeric(FLERR, arg[add+9], false, lmp);
    params[bond_type].Mmat[4][5] = utils::numeric(FLERR, arg[add+10], false, lmp);
    params[bond_type].Mmat[5][5] = utils::numeric(FLERR, arg[add+11], false, lmp);
    add += 12;
  }

  if (narg == opt3) {
    params[bond_type].Mmat[0][0] = utils::numeric(FLERR, arg[add+0], false, lmp);
    params[bond_type].Mmat[0][1] = utils::numeric(FLERR, arg[add+1], false, lmp);
    params[bond_type].Mmat[0][2] = utils::numeric(FLERR, arg[add+2], false, lmp);
    params[bond_type].Mmat[0][3] = utils::numeric(FLERR, arg[add+3], false, lmp);
    params[bond_type].Mmat[0][4] = utils::numeric(FLERR, arg[add+4], false, lmp);
    params[bond_type].Mmat[0][5] = utils::numeric(FLERR, arg[add+5], false, lmp);
    params[bond_type].Mmat[1][1] = utils::numeric(FLERR, arg[add+6], false, lmp);
    params[bond_type].Mmat[1][2] = utils::numeric(FLERR, arg[add+7], false, lmp);
    params[bond_type].Mmat[1][3] = utils::numeric(FLERR, arg[add+8], false, lmp);
    params[bond_type].Mmat[1][4] = utils::numeric(FLERR, arg[add+9], false, lmp);
    params[bond_type].Mmat[1][5] = utils::numeric(FLERR, arg[add+10], false, lmp);
    params[bond_type].Mmat[2][2] = utils::numeric(FLERR, arg[add+11], false, lmp);
    params[bond_type].Mmat[2][3] = utils::numeric(FLERR, arg[add+12], false, lmp);
    params[bond_type].Mmat[2][4] = utils::numeric(FLERR, arg[add+13], false, lmp);
    params[bond_type].Mmat[2][5] = utils::numeric(FLERR, arg[add+14], false, lmp);
    params[bond_type].Mmat[3][3] = utils::numeric(FLERR, arg[add+15], false, lmp);
    params[bond_type].Mmat[3][4] = utils::numeric(FLERR, arg[add+16], false, lmp);
    params[bond_type].Mmat[3][5] = utils::numeric(FLERR, arg[add+17], false, lmp);
    params[bond_type].Mmat[4][4] = utils::numeric(FLERR, arg[add+18], false, lmp);
    params[bond_type].Mmat[4][5] = utils::numeric(FLERR, arg[add+19], false, lmp);
    params[bond_type].Mmat[5][5] = utils::numeric(FLERR, arg[add+20], false, lmp);
    add += 21;
  }

  // mirror entries
  for (int ii = 0; ii < 6; ii++) {
    for (int jj = ii+1; jj < 6; jj++) {
      params[bond_type].Mmat[jj][ii] = params[bond_type].Mmat[ii][jj];
    }
  }

  for (int r = 0; r < 3; r++) {
    for (int c = 0; c < 3; c++) {
      params[bond_type].Mr[r][c] = params[bond_type].Mmat[r][c];
      params[bond_type].Mtr_tr[r][c] = params[bond_type].Mmat[r][c+3];
      params[bond_type].Mtr_bl[r][c] = params[bond_type].Mmat[r+3][c];
      params[bond_type].Mt[r][c] = params[bond_type].Mmat[r+3][c+3];
    }
  }

  // check if matrix is positive definite
  if (!lamath::is_positive_definite(params[bond_type].Mmat)) {
    error->all(FLERR, "Stiffness matrix M is not positive definite.");
  }
  
  setflag[bond_type] = 1;
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

  for (int i = 1; i <= atom->nbondtypes; i++) {
    if (!setflag[i])
      error->all(FLERR, "Not all bond coefficients are set for bond style rbp");
  }

  if (force->special_lj[1] != 0.0) {
    if (comm->me == 0) error->warning(FLERR, "Use special bonds = 0,x,x with bond style rbpfene");
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
  // - Loop over nbondtypes
  // - Write 36 values for Mmat[type][6][6]
  // - Write 12 values for ss[type] (3x3 rotation + 3 translation)
  // - Use fwrite for binary output
}

/* ----------------------------------------------------------------------
   Read bond style parameters from restart file
------------------------------------------------------------------------- */

void BondRBPFene::read_restart(FILE *fp) {
  // - Allocate Mmat and ss if not already done
  // - Read 36 values for Mmat[type]
  // - Read 12 values for ss[type]
  // - Optional: validate input and sanity-check matrix
}

/* ----------------------------------------------------------------------
   Write bond style data (coefficients) to LAMMPS data file
------------------------------------------------------------------------- */

void BondRBPFene::write_data(FILE *fp) {
  // - Loop over all bond types
  // - Write bond_coeff line: 21 M values + 6 SE(3) values
  // - Match the format expected by coeff()
  // - Use fprintf for clean output formatting
}

/* ----------------------------------------------------------------------
   Allocate memory for internal data structures
------------------------------------------------------------------------- */

void BondRBPFene::allocate() 
{
  allocated = 1;
  const int np1 = atom->nbondtypes + 1;

  memory->create(params, np1, "bond:rbp:params");
  memory->create(setflag, np1, "bond:rbp:setflag");
  for (int i = 1; i < np1; i++) setflag[i] = 0;
}
