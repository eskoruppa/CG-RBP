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
  
  double A[3], B[3], C[3], D[3];
  double torque_1[3], torque_2[3], torque_3[3], torque_4[3];
  double force_1[3], force_2[3], force_3[3], force_4[3];
  
  // temp
  double tmp1[3],tmp2[3],tmp3[3];
  
  ev_init(eflag, vflag);
  
  for (int aid = 0; aid < nanglelist; aid++) {
    // unpack angle endpoints and type
    id1 = anglelist[aid][0];
    id2 = anglelist[aid][1];
    id3 = anglelist[aid][2];
    angle_type = anglelist[aid][3];

    //-------------------------------------------------------------//
    // Compute triads and positions
    //-------------------------------------------------------------//

    // get quaternions
    quat1=bonus[ellipsoid[id1]].quat;
    quat2=bonus[ellipsoid[id2]].quat;
    quat3=bonus[ellipsoid[id3]].quat;
    
    // transform quat to triads [SO(3)]
    MathExtra::quat_to_mat(quat1, T1);
    MathExtra::quat_to_mat(quat2, T2);
    MathExtra::quat_to_mat(quat3, T3);
    
    // compute R
    lamath::mul_AtB(T1,T2,R_a);
    lamath::mul_AtB(T2,T3,R_b);

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

      // compute dr_a and dr_b
    lamath::subtract(r2,r1,dr_a);
    lamath::subtract(r3,r2,dr_b);
    // check for domain mismatch of special neighbor
    if (domain->minimum_image_check(dr_a[0],dr_a[1],dr_a[2])) {
      domain->minimum_image(FLERR,dr_a[0],dr_a[1],dr_a[2]);
    }
    if (domain->minimum_image_check(dr_b[0],dr_b[1],dr_b[2])) {
      domain->minimum_image(FLERR,dr_b[0],dr_b[1],dr_b[2]);
    }

    // compute w_a and w_b
    lamath::mul_Atx(T1,dr_a,w_a);
    lamath::mul_Atx(T2,dr_b,w_b);


    if (params[angle_type].subtract_groundstate) {
      //-------------------------------------------------------------//
      // Compute force wrench for se(3) (X) convention
      //-------------------------------------------------------------//

      // compute Omega_a and Omega_b
      so3::rotmat2euler(R_a,Om_a);
      so3::rotmat2euler(R_b,Om_b);

      // compute Omega_delta (Omega_d = Omega - Omega_s)
      lamath::subtract(Om_a,params[angle_type].srot1,Omd_a);
      lamath::subtract(Om_b,params[angle_type].srot2,Omd_b);

      // compute w_delta (w_d = w - w_s)
      lamath::subtract(w_a,params[angle_type].svec1,wd_a);
      lamath::subtract(w_b,params[angle_type].svec2,wd_b);

      // compute partial E / partial Omega_Delta,a (A)
      lamath::mul(params[angle_type].Mrr,Omd_b,tmp1);
      lamath::mul(params[angle_type].Mrt,wd_b,tmp2);
      lamath::add(tmp1,tmp2,A);

      // compute partial E / partial w_Delta,a (B)
      lamath::mul(params[angle_type].Mtt,wd_b,tmp1);
      lamath::mul(params[angle_type].Mtr,Omd_b,tmp2);
      lamath::add(tmp1,tmp2,B);

      // compute partial E / partial Omega_Delta,b (C)
      // lamath::mul(params[angle_type].Mrr_tp,Omd_a,tmp1);
      // lamath::mul(params[angle_type].Mtr_tp,wd_a,tmp2);
      lamath::mul_Atx(params[angle_type].Mrr,Omd_a,tmp1);
      lamath::mul_Atx(params[angle_type].Mtr,wd_a,tmp2);
      lamath::add(tmp1,tmp2,C);

      // compute partial E / partial w_Delta,b (D)
      // lamath::mul(params[angle_type].Mrt_tp,Omd_a,tmp1);
      // lamath::mul(params[angle_type].Mtt_tp,wd_a,tmp2);
      lamath::mul_Atx(params[angle_type].Mrt,Omd_a,tmp1);
      lamath::mul_Atx(params[angle_type].Mtt,wd_a,tmp2);
      lamath::add(tmp1,tmp2,D);

      // force 1 and 2
      lamath::mul(T1,B,force_1);
      lamath::signflip(force_1,force_2);
      
      // torques 1 and 2
      so3::leftJacobianInverseTransposed(Om_a,Jinvtp_a);
      lamath::mul(Jinvtp_a,A,tmp1);
      MathExtra::cross3(w_a,B,tmp2);
      lamath::add(tmp1,tmp2,tmp3);
      
      lamath::mul(T1,tmp3,torque_1);
      lamath::mul(T1,tmp1,tmp2);
      lamath::signflip(tmp2,torque_2);
      
      // force 2
      lamath::mul(T2,D,force_3);
      lamath::signflip(force_3,force_4);
      
      // torques 3 and 4
      so3::leftJacobianInverseTransposed(Om_b,Jinvtp_b);
      lamath::mul(Jinvtp_b,C,tmp1);
      MathExtra::cross3(w_b,D,tmp2);
      lamath::add(tmp1,tmp2,tmp3);

      lamath::mul(T2,tmp3,torque_3);
      lamath::mul(T2,tmp1,tmp2);
      lamath::signflip(tmp2,torque_4);

    }
    else {
      //-------------------------------------------------------------//
      // Compute force wrench for SE(3) (Y) convention
      //-------------------------------------------------------------//

      // compute Phi_delta_a (use Jinvtp_a for D_a and Omd_a for Phi_delta_a)
      lamath::mul_AtB(params[angle_type].Smat1,R_a,Jinvtp_a);
      so3::rotmat2euler(Jinvtp_a,Omd_a);
      
      // compute Phi_delta_b (use Jinvtp_b for D_b and Omd_b for Phi_delta_b)
      lamath::mul_AtB(params[angle_type].Smat2,R_b,Jinvtp_b);
      so3::rotmat2euler(Jinvtp_b,Omd_b);

      // compute d_a (reuse wd_a for d_a)
      lamath::subtract(w_a,params[angle_type].svec1,tmp1);
      lamath::mul_Atx(params[angle_type].Smat1,tmp1,wd_a);

      // compute d_b (reuse wd_b for d_b)
      lamath::subtract(w_b,params[angle_type].svec2,tmp1);
      lamath::mul_Atx(params[angle_type].Smat2,tmp1,wd_b);

      // compute partial E / partial Phi_Delta,a (A)
      lamath::mul(params[angle_type].Mrr,Omd_b,tmp1);
      lamath::mul(params[angle_type].Mrt,wd_b,tmp2);
      lamath::add(tmp1,tmp2,A);

      // compute partial E / partial d_a (B)
      lamath::mul(params[angle_type].Mtt,wd_b,tmp1);
      lamath::mul(params[angle_type].Mtr,Omd_b,tmp2);
      lamath::add(tmp1,tmp2,B);

      // compute partial E / partial Phi_Delta,b (C)
      // lamath::mul(params[angle_type].Mrr_tp,Omd_a,tmp1);
      // lamath::mul(params[angle_type].Mtr_tp,wd_a,tmp2);
      lamath::mul_Atx(params[angle_type].Mrr,Omd_a,tmp1);
      lamath::mul_Atx(params[angle_type].Mtr,wd_a,tmp2);
      lamath::add(tmp1,tmp2,C);

      // compute partial E / partial d_b (D)
      // lamath::mul(params[angle_type].Mrt_tp,Omd_a,tmp1);
      // lamath::mul(params[angle_type].Mtt_tp,wd_a,tmp2);
      lamath::mul_Atx(params[angle_type].Mrt,Omd_a,tmp1);
      lamath::mul_Atx(params[angle_type].Mtt,wd_a,tmp2);
      lamath::add(tmp1,tmp2,D);

      // compute transposed inverse left Jacobians
      so3::leftJacobianInverseTransposed(Omd_a,Jinvtp_a);
      so3::leftJacobianInverseTransposed(Omd_b,Jinvtp_b);

      // compute torque 2
      lamath::mul(Jinvtp_a,A,tmp1);
      lamath::mul(params[angle_type].Smat1,tmp1,tmp2);
      lamath::mul(T1,tmp2,tmp3);
      lamath::signflip(tmp3,torque_2);
      
      // compute force 1 and force 2
      lamath::mul(params[angle_type].Smat1,B,tmp1);
      lamath::mul(T1,tmp1,force_1);
      lamath::signflip(force_1,force_2);
      
      // compute torque 1
      MathExtra::cross3(w_a,tmp1,tmp2);
      lamath::mul(T1,tmp2,torque_1);
      lamath::add_to(torque_1,tmp3);

      // compute torque 4
      lamath::mul(Jinvtp_b,C,tmp1);
      lamath::mul(params[angle_type].Smat2,tmp1,tmp2);
      lamath::mul(T2,tmp2,tmp3);
      lamath::signflip(tmp3,torque_4);

      // compute force 3 and force 4
      lamath::mul(params[angle_type].Smat2,D,tmp1);
      lamath::mul(T2,tmp1,force_3);
      lamath::signflip(force_3,force_4);

      // compute torque 3
      MathExtra::cross3(w_b,tmp1,tmp2);
      lamath::mul(T2,tmp2,torque_3);
      lamath::add_to(torque_3,tmp3);
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

      f[id2][0] += force_2[0] + force_3[0];
      f[id2][1] += force_2[1] + force_3[1];
      f[id2][2] += force_2[2] + force_3[2];

      torque[id2][0] += torque_2[0] + torque_3[0];
      torque[id2][1] += torque_2[1] + torque_3[1];
      torque[id2][2] += torque_2[2] + torque_3[2];
    }

    if (newton_bond || id3 < nlocal) {

      f[id3][0] += force_4[0];
      f[id3][1] += force_4[1];
      f[id3][2] += force_4[2];

      torque[id3][0] += torque_4[0];
      torque[id3][1] += torque_4[1];
      torque[id3][2] += torque_4[2];
      
    }

    if (evflag) {

      //-------------------------------------------------------------//
      // ev tally 
      //-------------------------------------------------------------//

      double angle_energy = 0.0;
      const double Yd1[6] = {Omd_a[0], Omd_a[1], Omd_a[2], wd_a[0], wd_a[1], wd_a[2]};
      const double Yd2[6] = {Omd_b[0], Omd_b[1], Omd_b[2], wd_b[0], wd_b[1], wd_b[2]};
      for (int i_mat = 0; i_mat < 6; ++i_mat) {
        for (int j_mat = 0; j_mat < 6; ++j_mat) {
          angle_energy += Yd1[i_mat] * params[angle_type].Mmat[i_mat][j_mat] * Yd2[j_mat];
        }
      }
      angle_energy *= 0.5;

      // Tally energy + virial (forces on atoms 1 and 3)
      ev_tally(id1, id2, id3,
               nlocal, newton_bond,
               angle_energy,
               force_1,    // force on atom 1
               force_4,    // force on atom 3
               -dr_a[0], -dr_a[1], -dr_a[2],
               dr_b[0], dr_b[1], dr_b[2]);
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
    
    // Validate that angle style in database matches this style
    if (db.metadata().angle_style != RBPANGLE_STYLE) {
      error->all(FLERR, 
        "Angle style mismatch: database specifies '" + db.metadata().angle_style + 
        "' but using angle_style rbp");
    }

    int dbid = utils::inumeric(FLERR, arg[3], false, lmp);
    for (int angle_type=ilo;angle_type<=ihi;angle_type++) {
      assign_coeffs(angle_type,db.angle(dbid++).coeffs,db.metadata().subtract_groundstate);
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
    assign_coeffs(angle_type,coeffs,RBP_ANGLE_DEFAULT_SUBTRACT_GROUNDSTATE);
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
  fwrite(&params[1], sizeof(RBPParams), atom->nangletypes, fp);
}

/* ----------------------------------------------------------------------
   Read angle style parameters from restart file
------------------------------------------------------------------------- */

void AngleRBP::read_restart(FILE *fp) {
  allocate();

  if (comm->me == 0) {
    utils::sfread(FLERR,
                  &params[1],
                  sizeof(RBPParams),
                  atom->nangletypes,
                  fp,
                  nullptr,
                  error);
  }

  // broadcast the whole params block to all ranks
  MPI_Bcast(&params[1],
            atom->nangletypes * static_cast<int>(sizeof(RBPParams)),
            MPI_BYTE,
            0,
            world);

  // mark all types as having coefficients
  for (int i = 1; i <= atom->nangletypes; i++) setflag[i] = 1;
}

/* ----------------------------------------------------------------------
   Write angle style data (coefficients) to LAMMPS data file
------------------------------------------------------------------------- */

void AngleRBP::write_data(FILE *fp) {
  for (int i = 1; i <= atom->nangletypes; i++) {
    if (!setflag[i]) continue;

    // type index
    fprintf(fp, "%d", i);

    // 6 ground-state components (Ystatic1)
    for (int k = 0; k < 6; k++)
      fprintf(fp, " %g", params[i].Ystatic1[k]);

    // 6 ground-state components (Ystatic2)
    for (int k = 0; k < 6; k++)
      fprintf(fp, " %g", params[i].Ystatic2[k]);

    // 36 full stiffness matrix components (Mmat)
    for (int r = 0; r < 6; r++)
      for (int c = 0; c < 6; c++)
        fprintf(fp, " %g", params[i].Mmat[r][c]);

    fprintf(fp, "\n");
  }
}


void AngleRBP::allocate() {
  allocated = 1;
  int n = atom->nangletypes;

  memory->create(params,  n+1, "angle:rbp:params");
  memory->create(setflag, n+1, "angle:rbp:setflag");
  for (int i = 0; i <= n; i++) setflag[i] = 0;
}


void AngleRBP::assign_coeffs(int angle_type, const std::vector<double> &args, bool subtract_groundstate) {


  if (args.size() != 48) {
    error->all(FLERR, "Invalid number of coefficients found for angle style rbp. Requires 48 coefficients: X0_1 (6) X0_2 (6) stiffmat (6x6).");
  }

  params[angle_type].subtract_groundstate = subtract_groundstate;

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

  setflag[angle_type] = 1;
}