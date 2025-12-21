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

#include "angle_rbp.h"
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

AngleRBP::AngleRBP(LAMMPS *lmp) : Angle(lmp) {
  params = nullptr;
}

/* ----------------------------------------------------------------------
   Destructor: deallocates dynamically allocated memory
------------------------------------------------------------------------- */

AngleRBP::~AngleRBP() {
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(params);
    allocated = 0;
  }
}

/* ----------------------------------------------------------------------
   Compute forces and energy contribution for all angles of this style
------------------------------------------------------------------------- */


void AngleRBP::compute(int eflag, int vflag) {


  int id1,id2,id3;
  int angle_type;
  
  int **anglelist = neighbor->anglelist;
  int nanglelist  = neighbor->nanglelist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;
  
  auto avec = dynamic_cast<AtomVecEllipsoid *>(atom->style_match("ellipsoid"));
  AtomVecEllipsoid::Bonus *bonus = avec->bonus;
  int *ellipsoid = atom->ellipsoid;
  
  double **x = atom->x;                
  double **f = atom->f;             
  double **torque = atom->torque;   

  // Rotation variables
  double *quat1,*quat2,*quat3;
  double T1[3][3],T2[3][3],T3[3][3];
  double T1tp[3][3],T2tp[3][3];

  double R1[3][3],R2[3][3];
  double Om1[3],Om2[3];
  double Omd1[3],Omd2[3];
  double Jinvtp1[3][3];
  double Jinvtp2[3][3];

  // Translation variables
  double r1[3],r2[3],r3[3];
  double dr1[3],dr2[3];
  double w1[3],w2[3];
  double wd1[3],wd2[3];

  // temp
  double tmp1[3],tmp2[3],tmp3[3];

  // wrench components
  double gamma1[3];
  double gamma2[3];
  double mu1[3];
  double mu2[3];

  // torques
  double tau1[3],tau2[3],tau3[3],tau4[3];
  // forces
  double f1[3],f2[3],f3[3],f4[3];
  
  ev_init(eflag, vflag);
  
  for (int aid = 0; aid < nanglelist; aid++) {
    // unpack angle endpoints and type
    id1 = anglelist[aid][0];
    id2 = anglelist[aid][1];
    id3 = anglelist[aid][2];
    angle_type = anglelist[aid][3];

    // fprintf(screen,"AngleRBP %d %d %d\n",id1,id2,id3);
    
    ////////////////////////////////////////
    // Rotational Part
    // get quaternions
    quat1=bonus[ellipsoid[id1]].quat;
    quat2=bonus[ellipsoid[id2]].quat;
    quat3=bonus[ellipsoid[id3]].quat;
    
    // transform quat to triads [SO(3)]
    MathExtra::quat_to_mat(quat1, T1);
    MathExtra::quat_to_mat(quat2, T2);
    MathExtra::quat_to_mat(quat3, T3);
    
    // compute transposed
    lamath::transpose(T1,T1tp);
    lamath::transpose(T2,T2tp);
    
    // compute R
    lamath::mul(T1tp,T2,R1);
    lamath::mul(T2tp,T3,R2);
    
    // compute Omega_1
    so3::rotmat2euler(R1,Om1);
    so3::rotmat2euler(R2,Om2);
    
    // compute Omega_Delta
    Omd1[0] = Om1[0] - params[angle_type].srot1[0];
    Omd1[1] = Om1[1] - params[angle_type].srot1[1];
    Omd1[2] = Om1[2] - params[angle_type].srot1[2];
    Omd2[0] = Om2[0] - params[angle_type].srot2[0];
    Omd2[1] = Om2[1] - params[angle_type].srot2[1];
    Omd2[2] = Om2[2] - params[angle_type].srot2[2];
    
    // compute transposed inverse left Jacobian
    so3::leftJacobianInverseTransposed(Om1,Jinvtp1);
    so3::leftJacobianInverseTransposed(Om2,Jinvtp2);
    
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

    // compute w
    dr1[0] = r2[0] - r1[0]; 
    dr1[1] = r2[1] - r1[1]; 
    dr1[2] = r2[2] - r1[2]; 
    dr2[0] = r3[0] - r2[0]; 
    dr2[1] = r3[1] - r2[1]; 
    dr2[2] = r3[2] - r2[2]; 

    // check for domain mismatch of special neighbor
    if (domain->minimum_image_check(dr1[0],dr1[1],dr1[2])) {
      domain->minimum_image(dr1[0],dr1[1],dr1[2]);
    }
    if (domain->minimum_image_check(dr2[0],dr2[1],dr2[2])) {
      domain->minimum_image(dr2[0],dr2[1],dr2[2]);
    }

    lamath::mul(T1tp,dr1,w1);
    lamath::mul(T2tp,dr2,w2);
  
    // compute w_Delta
    wd1[0] = w1[0] - params[angle_type].svec1[0];
    wd1[1] = w1[1] - params[angle_type].svec1[1];
    wd1[2] = w1[2] - params[angle_type].svec1[2];
    wd2[0] = w2[0] - params[angle_type].svec2[0];
    wd2[1] = w2[1] - params[angle_type].svec2[1];
    wd2[2] = w2[2] - params[angle_type].svec2[2];

    //////////////////////////////////////////////////////
    // Compute gamma and mu

    // gamma1
    lamath::mul(params[angle_type].Mrr,Omd2,tmp1);
    lamath::mul(params[angle_type].Mrt,wd2,tmp2);
    lamath::add(tmp1,tmp2,gamma1);
    // gamma1[0] = tmp1[0] + tmp2[0];
    // gamma1[1] = tmp1[1] + tmp2[1];
    // gamma1[2] = tmp1[2] + tmp2[2];

    // mu1
    lamath::mul(params[angle_type].Mtr,Omd2,tmp1);
    lamath::mul(params[angle_type].Mtt,wd2,tmp2);
    lamath::add(tmp1,tmp2,mu1);
    // mu1[0] = tmp1[0] + tmp2[0];
    // mu1[1] = tmp1[1] + tmp2[1];
    // mu1[2] = tmp1[2] + tmp2[2];

    // gamma2
    lamath::mul(params[angle_type].Mrr_tp,Omd1,tmp1);
    lamath::mul(params[angle_type].Mtr_tp,wd1,tmp2);
    lamath::add(tmp1,tmp2,gamma2);
    // gamma2[0] = tmp1[0] + tmp2[0];
    // gamma2[1] = tmp1[1] + tmp2[1];
    // gamma2[2] = tmp1[2] + tmp2[2];

    // mu1
    lamath::mul(params[angle_type].Mrt_tp,Omd1,tmp1);
    lamath::mul(params[angle_type].Mtt_tp,wd1,tmp2);
    lamath::add(tmp1,tmp2,mu2);
    // mu2[0] = tmp1[0] + tmp2[0];
    // mu2[1] = tmp1[1] + tmp2[1];
    // mu2[2] = tmp1[2] + tmp2[2];

    // Triad 1 ////////
    // tau1 
    MathExtra::cross3(w1,mu1,tmp1);
    lamath::mul(Jinvtp1,gamma1,tmp2);
    lamath::add(tmp1,tmp2,tmp3);
    lamath::mul(T1,tmp3,tau1);
    // f1
    lamath::mul(T1,mu1,f1);
    
    // Triad 2 ////////
    // tau2
    lamath::mul(T1,tmp2,tmp3);
    lamath::signflip(tmp3,tau2);
    
    // f2
    lamath::mul(T1,mu1,tmp3);
    lamath::signflip(tmp3,f2);
    
    // Triad 3 ////////
    // tau3 
    MathExtra::cross3(w2,mu2,tmp1);
    lamath::mul(Jinvtp2,gamma2,tmp2);
    lamath::add(tmp1,tmp2,tmp3);
    lamath::mul(T2,tmp3,tau3);
    // f3
    lamath::mul(T2,mu2,f3);

    // Triad 4 ////////
    // tau4
    lamath::mul(T2,tmp2,tmp3);
    lamath::signflip(tmp3,tau4);
    
    // f4
    lamath::mul(T2,mu2,tmp3);
    lamath::signflip(tmp3,f4);
    


    ////////////////////////////////////////////
    ////////////////////////////////////////////

    // apply force and torque to each of 3 atoms
    if (newton_bond || id1 < nlocal) {

      f[id1][0] += f1[0];
      f[id1][1] += f1[1];
      f[id1][2] += f1[2];

      torque[id1][0] += tau1[0];
      torque[id1][1] += tau1[1];
      torque[id1][2] += tau1[2];
      
    }

    if (newton_bond || id2 < nlocal) {

      f[id2][0] += f2[0] + f3[0];
      f[id2][1] += f2[1] + f3[1];
      f[id2][2] += f2[2] + f3[2];

      torque[id2][0] += tau2[0] + tau3[0];
      torque[id2][1] += tau2[1] + tau3[1];
      torque[id2][2] += tau2[2] + tau3[2];
    }

    if (newton_bond || id3 < nlocal) {

      f[id3][0] += f4[0];
      f[id3][1] += f4[1];
      f[id3][2] += f4[2];

      torque[id3][0] += tau4[0];
      torque[id3][1] += tau4[1];
      torque[id3][2] += tau4[2];
      
    }

    if (evflag) {

      double angle_energy = 0.0;
      double Yd1[6]; 
      double Yd2[6]; 

      Yd1[0] = Omd1[0]; 
      Yd1[1] = Omd1[1];
      Yd1[2] = Omd1[2];
      Yd1[3] = wd1[0];   
      Yd1[4] = wd1[1];
      Yd1[5] = wd1[2];

      Yd2[0] = Omd2[0]; 
      Yd2[1] = Omd2[1];
      Yd2[2] = Omd2[2];
      Yd2[3] = wd2[0];   
      Yd2[4] = wd2[1];
      Yd2[5] = wd2[2];

      // Calculate energy: 0.5 * Y_vec^T * Mmat * Y_vec
      for (int i_mat = 0; i_mat < 6; ++i_mat) { 
        for (int j_mat = 0; j_mat < 6; ++j_mat) {
          angle_energy += Yd1[i_mat] * params[angle_type].Mmat[i_mat][j_mat] * Yd2[j_mat];
        }
      }
      angle_energy *= 0.5;

      // bond vectors for virial: from center (id2) to ends (id1, id3)
      double delx1 = -dr1[0];
      double dely1 = -dr1[1];
      double delz1 = -dr1[2];

      double delx2 = dr2[0];
      double dely2 = dr2[1];
      double delz2 = dr2[2];

      // Tally energy + virial (forces on atoms 1 and 3)
      ev_tally(id1, id2, id3,
               nlocal, newton_bond,
               angle_energy,
               f1,    // force on atom 1
               f4,    // force on atom 3
               delx1, dely1, delz1,
               delx2, dely2, delz2);

    }
  }
}

/* ----------------------------------------------------------------------
   Compute single
------------------------------------------------------------------------- */

double AngleRBP::single(int type, int i1, int i2, int i3) {
  return 0.0;
}

/* ----------------------------------------------------------------------
   Parse coefficients for each angle type
------------------------------------------------------------------------- */


void AngleRBP::coeff(int narg, char **arg) {

  // check conditions for properly formated arg !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // may not be necessary

  if (!allocated) allocate();

  int ilo, ihi;
  utils::bounds(FLERR, arg[0], 1, atom->nangletypes, ilo, ihi, error);

  // load from database file
  if (narg == 4 && strcmp(arg[1], "dbfile") == 0) {

    RBPDatabase db(lmp, error);
    db.read(arg[2]);

    int dbid = utils::inumeric(FLERR, arg[3], false, lmp);
    for (int angle_type=ilo;angle_type<=ihi;angle_type++) {
      assign_coeffs(angle_type,db.angle(dbid++).coeffs);
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
  for (int angle_type=ilo;angle_type<=ihi;angle_type++) {
    assign_coeffs(angle_type,coeffs);
  }
}

/* ----------------------------------------------------------------------
   Initialization: verify atom properties and prepare style
------------------------------------------------------------------------- */

void AngleRBP::init_style() {

  // Ensure ellipsoidal atoms are being used
  if (!atom->ellipsoid_flag)
    error->all(FLERR, "Angle style rbp requires atom style with ellipsoids");

  if (domain->dimension != 3)
    error->all(FLERR, "Angle style rbp requires a 3D simulation");
}


/* ----------------------------------------------------------------------
   Return equilibrium angle for angle type i
------------------------------------------------------------------------- */

double AngleRBP::equilibrium_angle(int angle_type) 
{
  return 0.0;
}

/* ----------------------------------------------------------------------
   Write angle style parameters to restart file
------------------------------------------------------------------------- */

void AngleRBP::write_restart(FILE *fp) {
  // - Loop over nangletypes
  // - Write 36 values for Mmat[type][6][6]
  // - Write 12 values for ss[type] (3x3 rotation + 3 translation)
  // - Use fwrite for binary output
}

/* ----------------------------------------------------------------------
   Read angle style parameters from restart file
------------------------------------------------------------------------- */

void AngleRBP::read_restart(FILE *fp) {
  // - Allocate Mmat and ss if not already done
  // - Read 36 values for Mmat[type]
  // - Read 12 values for ss[type]
  // - Optional: validate input and sanity-check matrix
}

/* ----------------------------------------------------------------------
   Write angle style data (coefficients) to LAMMPS data file
------------------------------------------------------------------------- */

void AngleRBP::write_data(FILE *fp) {
  // - Loop over all angle types
  // - Write angle_coeff line: 21 M values + 6 SE(3) values
  // - Match the format expected by coeff()
  // - Use fprintf for clean output formatting
}

/* ----------------------------------------------------------------------
   Allocate memory for internal data structures  //   // compute gamma
  //   lamath::mul(params[angle_type].Mr,Omd,gamma);
  //   lamath::mul(params[angle_type].Mtr_tr,wd,tmp);
  //   gamma[0] += tmp[0];
  //   gamma[1] += tmp[1];
  //   gamma[2] += tmp[2];
    
  //   // compute mu
  //   lamath::mul(params[angle_type].Mtr_bl,Omd,mu);
  //   lamath::mul(params[angle_type].Mt,wd,tmp);
  //   mu[0] += tmp[0];
  //   mu[1] += tmp[1];
  //   mu[2] += tmp[2];
------------------------------------------------------------------------- */

void AngleRBP::allocate() {
  // - Allocate Mmat[nangletypes+1][6][6]
  // - Allocate ss[nangletypes+1]
  // - Initialize values to zero or identity where appropriate
  // - Use LAMMPS memory->create interface

  allocated = 1;
  int n = atom->nangletypes;

  memory->create(params,  n+1, "angle:rbp:params");
  memory->create(setflag, n+1, "angle:rbp:setflag");
  for (int i = 0; i <= n; i++) setflag[i] = 0;
}


void AngleRBP::assign_coeffs(int angle_type, const std::vector<double> &args) {


  if (args.size() != 48) {
    error->all(FLERR, "Invalid number of coefficients found for angle style rbp. Requires 48 coefficients: X0_1 (6) X0_2 (6) stiffmat (6x6).");
  }

  // assign Ystatic
  for (int i=0;i<6;i++) {
    params[angle_type].Ystatic1[i] = args[i];
    params[angle_type].Ystatic2[i] = args[i+6];
  }

  // assign partial static
  for (int ii=0;ii<3;ii++) {
    params[angle_type].srot1[ii] = params[angle_type].Ystatic1[ii];
    params[angle_type].svec1[ii] = params[angle_type].Ystatic1[ii+3];
    params[angle_type].srot2[ii] = params[angle_type].Ystatic2[ii];
    params[angle_type].svec2[ii] = params[angle_type].Ystatic2[ii+3];
  }

  // assigne Smat
  so3::euler2rotmat(params[angle_type].srot1,params[angle_type].Smat1);
  so3::euler2rotmat(params[angle_type].srot2,params[angle_type].Smat2);

  int argid = 12;
  for (int ii=0;ii<6;ii++) {
    for (int jj=0;jj<6;jj++) {
      params[angle_type].Mmat[ii][jj] = args[argid++];
    }
  }

  for (int r = 0; r < 3; r++) {
    for (int c = 0; c < 3; c++) {
      params[angle_type].Mrr[r][c] = params[angle_type].Mmat[r][c];
      params[angle_type].Mrt[r][c] = params[angle_type].Mmat[r][c+3];
      params[angle_type].Mtr[r][c] = params[angle_type].Mmat[r+3][c];
      params[angle_type].Mtt[r][c] = params[angle_type].Mmat[r+3][c+3];
    }
  }

  lamath::transpose(params[angle_type].Mrr,params[angle_type].Mrr_tp);
  lamath::transpose(params[angle_type].Mtt,params[angle_type].Mtt_tp);
  lamath::transpose(params[angle_type].Mtr,params[angle_type].Mtr_tp);
  lamath::transpose(params[angle_type].Mrt,params[angle_type].Mrt_tp);

  // THESE OFF-DIAGONAL COMPONENTS NEED NOT NECESSARILY BE POSITIVE DEFINITE!!!
  // // check if matrix is positive definite
  // if (!lamath::is_positive_definite(params[angle_type].Mmat)) {
  //   error->all(FLERR, "Stiffness matrix M is not positive definite.");
  // }

  setflag[angle_type] = 1;
}