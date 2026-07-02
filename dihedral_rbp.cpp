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

#ifdef DIHEDRAL_RBP_PRECOMPUTE_ACTIVE
#include "fix_rbp_lrf.h"
#include "modify.h"
#endif


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
  
  #ifndef DIHEDRAL_RBP_PRECOMPUTE_ACTIVE
  auto avec = dynamic_cast<AtomVecEllipsoid *>(atom->style_match("ellipsoid"));
  AtomVecEllipsoid::Bonus *bonus = avec->bonus;
  int *ellipsoid = atom->ellipsoid;
  #endif
  
  double **x = atom->x;                
  double **f = atom->f;             
  double **torque = atom->torque;   

  // Rotation variables
  #ifndef DIHEDRAL_RBP_PRECOMPUTE_ACTIVE
  double *quat1,*quat2,*quat3,*quat4;
  double T1_arr[3][3], T2_arr[3][3], T3_arr[3][3], T4_arr[3][3];
  double R_a[3][3], R_b[3][3];
  double Om_a[3], Om_b[3];
  double Jinvtp_a[3][3];
  double Jinvtp_b[3][3];
  #endif
  
  // Translation variables
  double r1[3],r2[3],r3[3],r4[3];
  double dr_a[3],dr_b[3];
  double w_a[3],w_b[3];
  double wd_a[3],wd_b[3];
  
  double Omd_a[3],Omd_b[3];
  double A[3], B[3], C[3], D[3];
  double torque_1[3], torque_2[3], torque_3[3], torque_4[3];
  double force_1[3], force_2[3], force_3[3], force_4[3];
  
  // temp
  double tmp1[3],tmp2[3],tmp3[3];
  
  ev_init(eflag, vflag);
  
  for (int did = 0; did < ndihedrallist; did++) {
    // unpack dihedral endpoints and type
    id1 = dihedrallist[did][0];
    id2 = dihedrallist[did][1];
    id3 = dihedrallist[did][2];
    id4 = dihedrallist[did][3];
    dihedral_type = dihedrallist[did][4];

    //-------------------------------------------------------------//
    // Compute triads and positions
    //-------------------------------------------------------------//

    #ifdef DIHEDRAL_RBP_PRECOMPUTE_ACTIVE
    // load precomputed triads (zero-copy pointer cast)
    double (*T1)[3] = (double(*)[3]) fix_lrf->triads[id1];
    double (*T2)[3] = (double(*)[3]) fix_lrf->triads[id2];
    double (*T3)[3] = (double(*)[3]) fix_lrf->triads[id3];
    double (*T4)[3] = (double(*)[3]) fix_lrf->triads[id4];
    #else
    // get quaternions
    quat1=bonus[ellipsoid[id1]].quat;
    quat2=bonus[ellipsoid[id2]].quat;
    quat3=bonus[ellipsoid[id3]].quat;
    quat4=bonus[ellipsoid[id4]].quat;
    
    // transform quat to triads [SO(3)]
    MathExtra::quat_to_mat(quat1, T1_arr);
    MathExtra::quat_to_mat(quat2, T2_arr);
    MathExtra::quat_to_mat(quat3, T3_arr);
    MathExtra::quat_to_mat(quat4, T4_arr);
    double (*T1)[3] = T1_arr;
    double (*T2)[3] = T2_arr;
    double (*T3)[3] = T3_arr;
    double (*T4)[3] = T4_arr;

    // compute R
    lamath::mul_AtB(T1,T2,R_a);
    lamath::mul_AtB(T3,T4,R_b);
    #endif

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

    // compute dr_a and dr_b
    lamath::subtract(r2,r1,dr_a);
    lamath::subtract(r4,r3,dr_b);
    // check for domain mismatch of special neighbor
    if (domain->minimum_image_check(dr_a[0],dr_a[1],dr_a[2])) {
      domain->minimum_image(FLERR,dr_a[0],dr_a[1],dr_a[2]);
    }
    if (domain->minimum_image_check(dr_b[0],dr_b[1],dr_b[2])) {
      domain->minimum_image(FLERR,dr_b[0],dr_b[1],dr_b[2]);
    }

    // compute w_a and w_b
    lamath::mul_Atx(T1,dr_a,w_a);
    lamath::mul_Atx(T3,dr_b,w_b);

    if (params[dihedral_type].subtract_groundstate) {
      //-------------------------------------------------------------//
      // Compute force wrench for se(3) (X) convention
      //-------------------------------------------------------------//

      #ifdef DIHEDRAL_RBP_PRECOMPUTE_ACTIVE
      // Omd = Om - srot, where Om is precomputed in fwd_euler
      lamath::subtract(fix_lrf->fwd_euler[id1],params[dihedral_type].srot1,Omd_a);
      lamath::subtract(fix_lrf->fwd_euler[id3],params[dihedral_type].srot2,Omd_b);
      #else
      // compute Omega_a and Omega_b
      so3::rotmat2euler(R_a,Om_a);
      so3::rotmat2euler(R_b,Om_b);
      lamath::subtract(Om_a,params[dihedral_type].srot1,Omd_a);
      lamath::subtract(Om_b,params[dihedral_type].srot2,Omd_b);
      #endif

      // compute w_delta (w_d = w - w_s)
      lamath::subtract(w_a,params[dihedral_type].svec1,wd_a);
      lamath::subtract(w_b,params[dihedral_type].svec2,wd_b);

      // compute partial E / partial Omega_Delta,a (A)
      lamath::mul(params[dihedral_type].Mrr,Omd_b,tmp1);
      lamath::mul(params[dihedral_type].Mrt,wd_b,tmp2);
      lamath::add(tmp1,tmp2,A);

      // compute partial E / partial w_Delta,a (B)
      lamath::mul(params[dihedral_type].Mtt,wd_b,tmp1);
      lamath::mul(params[dihedral_type].Mtr,Omd_b,tmp2);
      lamath::add(tmp1,tmp2,B);

      // compute partial E / partial Omega_Delta,b (C)
      // lamath::mul(params[dihedral_type].Mrr_tp,Omd_a,tmp1);
      // lamath::mul(params[dihedral_type].Mtr_tp,wd_a,tmp2);
      lamath::mul_Atx(params[dihedral_type].Mrr,Omd_a,tmp1);
      lamath::mul_Atx(params[dihedral_type].Mtr,wd_a,tmp2);
      lamath::add(tmp1,tmp2,C);

      // compute partial E / partial w_Delta,b (D)
      // lamath::mul(params[dihedral_type].Mrt_tp,Omd_a,tmp1);
      // lamath::mul(params[dihedral_type].Mtt_tp,wd_a,tmp2);
      lamath::mul_Atx(params[dihedral_type].Mrt,Omd_a,tmp1);
      lamath::mul_Atx(params[dihedral_type].Mtt,wd_a,tmp2);
      lamath::add(tmp1,tmp2,D);

      // compute transposed inverse left Jacobians
      #ifdef DIHEDRAL_RBP_PRECOMPUTE_ACTIVE
      double (*Jinvtp_a)[3] = (double(*)[3]) fix_lrf->fwd_Jinvtp[id1];
      double (*Jinvtp_b)[3] = (double(*)[3]) fix_lrf->fwd_Jinvtp[id3];
      #else
      so3::leftJacobianInverseTransposed(Om_a,Jinvtp_a);
      so3::leftJacobianInverseTransposed(Om_b,Jinvtp_b);
      #endif

      // force 1 and 2
      lamath::mul(T1,B,force_1);
      lamath::signflip(force_1,force_2);

      // torques 1 and 2
      lamath::mul(Jinvtp_a,A,tmp1);
      MathExtra::cross3(w_a,B,tmp2);
      lamath::add(tmp1,tmp2,tmp3);

      lamath::mul(T1,tmp3,torque_1);
      lamath::mul(T1,tmp1,tmp2);
      lamath::signflip(tmp2,torque_2);

      // force 2
      lamath::mul(T3,D,force_3);
      lamath::signflip(force_3,force_4);

      // torques 3 and 4
      lamath::mul(Jinvtp_b,C,tmp1);
      MathExtra::cross3(w_b,D,tmp2);
      lamath::add(tmp1,tmp2,tmp3);

      lamath::mul(T3,tmp3,torque_3);
      lamath::mul(T3,tmp1,tmp2);
      lamath::signflip(tmp2,torque_4);

      // lamath::print_mat3_lammps(screen, "Jinvtp_a", Jinvtp_a);
      // lamath::print_mat3_lammps(screen, "Jinvtp_b", Jinvtp_b);

    }
    else {
      //-------------------------------------------------------------//
      // Compute force wrench for SE(3) (Y) convention
      //-------------------------------------------------------------//

      #ifdef DIHEDRAL_RBP_PRECOMPUTE_ACTIVE
      // Phi_delta precomputed in fwd_euler (already the delta for Y)
      Omd_a[0] = fix_lrf->fwd_euler[id1][0];
      Omd_a[1] = fix_lrf->fwd_euler[id1][1];
      Omd_a[2] = fix_lrf->fwd_euler[id1][2];
      Omd_b[0] = fix_lrf->fwd_euler[id3][0];
      Omd_b[1] = fix_lrf->fwd_euler[id3][1];
      Omd_b[2] = fix_lrf->fwd_euler[id3][2];
      #else
      // compute Phi_delta_a (use Jinvtp_a for D_a and Omd_a for Phi_delta_a)
      lamath::mul_AtB(params[dihedral_type].Smat1,R_a,Jinvtp_a);
      so3::rotmat2euler(Jinvtp_a,Omd_a);
      
      // compute Phi_delta_b (use Jinvtp_b for D_b and Omd_b for Phi_delta_b)
      lamath::mul_AtB(params[dihedral_type].Smat2,R_b,Jinvtp_b);
      so3::rotmat2euler(Jinvtp_b,Omd_b);
      #endif

      // compute d_a (reuse wd_a for d_a)
      lamath::subtract(w_a,params[dihedral_type].svec1,tmp1);
      lamath::mul_Atx(params[dihedral_type].Smat1,tmp1,wd_a);

      // compute d_b (reuse wd_b for d_b)
      lamath::subtract(w_b,params[dihedral_type].svec2,tmp1);
      lamath::mul_Atx(params[dihedral_type].Smat2,tmp1,wd_b);

      // compute partial E / partial Phi_Delta,a (A)
      lamath::mul(params[dihedral_type].Mrr,Omd_b,tmp1);
      lamath::mul(params[dihedral_type].Mrt,wd_b,tmp2);
      lamath::add(tmp1,tmp2,A);

      // compute partial E / partial d_a (B)
      lamath::mul(params[dihedral_type].Mtt,wd_b,tmp1);
      lamath::mul(params[dihedral_type].Mtr,Omd_b,tmp2);
      lamath::add(tmp1,tmp2,B);

      // compute partial E / partial Phi_Delta,b (C)
      // lamath::mul(params[dihedral_type].Mrr_tp,Omd_a,tmp1);
      // lamath::mul(params[dihedral_type].Mtr_tp,wd_a,tmp2);
      lamath::mul_Atx(params[dihedral_type].Mrr,Omd_a,tmp1);
      lamath::mul_Atx(params[dihedral_type].Mtr,wd_a,tmp2);
      lamath::add(tmp1,tmp2,C);

      // compute partial E / partial d_b (D)
      // lamath::mul(params[dihedral_type].Mrt_tp,Omd_a,tmp1);
      // lamath::mul(params[dihedral_type].Mtt_tp,wd_a,tmp2);
      lamath::mul_Atx(params[dihedral_type].Mrt,Omd_a,tmp1);
      lamath::mul_Atx(params[dihedral_type].Mtt,wd_a,tmp2);
      lamath::add(tmp1,tmp2,D);

      // compute torque 2
      #ifdef DIHEDRAL_RBP_PRECOMPUTE_ACTIVE
      double (*Jinvtp_a)[3] = (double(*)[3]) fix_lrf->fwd_Jinvtp[id1];
      #else
      so3::leftJacobianInverseTransposed(Omd_a,Jinvtp_a);
      #endif

      lamath::mul(Jinvtp_a,A,tmp1);
      lamath::mul(params[dihedral_type].Smat1,tmp1,tmp2);
      lamath::mul(T1,tmp2,tmp3);
      lamath::signflip(tmp3,torque_2);
      
      // compute force 1 and force 2
      lamath::mul(params[dihedral_type].Smat1,B,tmp1);
      lamath::mul(T1,tmp1,force_1);
      lamath::signflip(force_1,force_2);
      
      // compute torque 1
      MathExtra::cross3(w_a,tmp1,tmp2);
      lamath::mul(T1,tmp2,torque_1);
      lamath::add_to(torque_1,tmp3);
      
      // compute torque 4
      #ifdef DIHEDRAL_RBP_PRECOMPUTE_ACTIVE
      double (*Jinvtp_b)[3] = (double(*)[3]) fix_lrf->fwd_Jinvtp[id3];
      #else
      so3::leftJacobianInverseTransposed(Omd_b,Jinvtp_b);
      #endif

      lamath::mul(Jinvtp_b,C,tmp1);
      lamath::mul(params[dihedral_type].Smat2,tmp1,tmp2);
      lamath::mul(T3,tmp2,tmp3);
      lamath::signflip(tmp3,torque_4);

      // compute force 3 and force 4
      lamath::mul(params[dihedral_type].Smat2,D,tmp1);
      lamath::mul(T3,tmp1,force_3);
      lamath::signflip(force_3,force_4);

      // compute torque 3
      MathExtra::cross3(w_b,tmp1,tmp2);
      lamath::mul(T3,tmp2,torque_3);
      lamath::add_to(torque_3,tmp3);

      // lamath::print_mat3_lammps(screen, "Jinvtp_a", Jinvtp_a);
      // lamath::print_mat3_lammps(screen, "Jinvtp_b", Jinvtp_b);

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

      f[id2][0] += force_2[0];
      f[id2][1] += force_2[1];
      f[id2][2] += force_2[2];

      torque[id2][0] += torque_2[0];
      torque[id2][1] += torque_2[1];
      torque[id2][2] += torque_2[2];
    }

    if (newton_bond || id3 < nlocal) {

      f[id3][0] += force_3[0];
      f[id3][1] += force_3[1];
      f[id3][2] += force_3[2];

      torque[id3][0] += torque_3[0];
      torque[id3][1] += torque_3[1];
      torque[id3][2] += torque_3[2];
      
    }

    if (newton_bond || id4 < nlocal) {

      f[id4][0] += force_4[0];
      f[id4][1] += force_4[1];
      f[id4][2] += force_4[2];

      torque[id4][0] += torque_4[0];
      torque[id4][1] += torque_4[1];
      torque[id4][2] += torque_4[2];
      
    }

    // // ---- DEBUG: verify ev_tally equivalence (delete this block when done) ----
    // {
    //   double dr32_dbg[3];
    //   dr32_dbg[0] = x[id3][0] - x[id2][0];
    //   dr32_dbg[1] = x[id3][1] - x[id2][1];
    //   dr32_dbg[2] = x[id3][2] - x[id2][2];
    //   if (domain->minimum_image_check(dr32_dbg[0], dr32_dbg[1], dr32_dbg[2]))
    //     domain->minimum_image(FLERR, dr32_dbg[0], dr32_dbg[1], dr32_dbg[2]);

    //   verify_ev_tally(id1, id2, id3, id4,
    //                   nlocal, newton_bond,
    //                   0.0,
    //                   force_1, force_2, force_3, force_4,
    //                   dr_a, dr_b, dr32_dbg);
    // }
    // // ---- END DEBUG ----

    if (evflag) {

      //-------------------------------------------------------------//
      // ev tally 
      //-------------------------------------------------------------//
    
      // Coupling energy E = Yd_a^T M Yd_b (no 1/2: M is the full off-diagonal
      // coupling block M_{i,i+2} of the global 1/2 sum_ij Y_i^T M_ij Y_j, so the
      // M_{i,i+2} and M_{i+2,i}=M^T terms add and cancel the 1/2). This matches
      // the forces, which use the full gradient (A,B)=M Yd_b, (C,D)=M^T Yd_a.
      double dihedral_energy = 0.0;
      const double Yd1[6] = {Omd_a[0], Omd_a[1], Omd_a[2], wd_a[0], wd_a[1], wd_a[2]};
      const double Yd2[6] = {Omd_b[0], Omd_b[1], Omd_b[2], wd_b[0], wd_b[1], wd_b[2]};
      for (int i_mat = 0; i_mat < 6; ++i_mat) {
        for (int j_mat = 0; j_mat < 6; ++j_mat) {
          dihedral_energy += Yd1[i_mat] * params[dihedral_type].Mmat[i_mat][j_mat] * Yd2[j_mat];
        }
      }
      
      #ifdef RBP_DIHEDRAL_USE_CUSTOM_EV_TALLY
      double dr32[3];
      dr32[0] = x[id3][0] - x[id2][0];
      dr32[1] = x[id3][1] - x[id2][1];
      dr32[2] = x[id3][2] - x[id2][2];
      if (domain->minimum_image_check(dr32[0], dr32[1], dr32[2]))
        domain->minimum_image(FLERR, dr32[0], dr32[1], dr32[2]);

      ev_tally_rbp(id1, id2, id3, id4,
                   nlocal, newton_bond,
                   dihedral_energy,
                   force_1, force_3,
                   dr_a, dr_b, dr32);
      #else
      // vb2 = r3 - r2  (needs its own minimum image)
      double vb2x = x[id3][0] - x[id2][0];
      double vb2y = x[id3][1] - x[id2][1];
      double vb2z = x[id3][2] - x[id2][2];
      if (domain->minimum_image_check(vb2x, vb2y, vb2z)) {
          domain->minimum_image(FLERR, vb2x, vb2y, vb2z);
      }

      ev_tally(id1, id2, id3, id4,
              nlocal, newton_bond,
              dihedral_energy,
              force_1, force_3, force_4,
              -dr_a[0], -dr_a[1], -dr_a[2],
              vb2x, vb2y, vb2z,
              dr_b[0], dr_b[1], dr_b[2]);
      #endif
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
    
    // Validate that dihedral style in database matches this style
    if (db.metadata().dihedral_style != RBPDIHEDRAL_STYLE) {
      error->all(FLERR, 
        "Dihedral style mismatch: database specifies '" + db.metadata().dihedral_style + 
        "' but using dihedral_style rbp");
    }

    int dbid = utils::inumeric(FLERR, arg[3], false, lmp);
    for (int dihedral_type=ilo;dihedral_type<=ihi;dihedral_type++) {
      assign_coeffs(dihedral_type,db.dihedral(dbid++).coeffs,db.metadata().subtract_groundstate);
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
    assign_coeffs(dihedral_type,coeffs,RBP_DIHEDRAL_DEFAULT_SUBTRACT_GROUNDSTATE);
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

  #ifdef DIHEDRAL_RBP_PRECOMPUTE_ACTIVE
  // auto-create or find fix rbp/lrf
  auto fixes = modify->get_fix_by_style("^rbp/lrf");
  if (fixes.empty())
    fix_lrf = dynamic_cast<FixRBPLRF*>(modify->add_fix("rbp_lrf all rbp/lrf"));
  else
    fix_lrf = dynamic_cast<FixRBPLRF*>(fixes[0]);

  // Register dihedral sub-junction params
  for (int i = 1; i <= atom->ndihedraltypes; i++) {
    if (setflag[i]) {
      fix_lrf->register_dihedral_junction(i, 0, params[i].subtract_groundstate,
                                          params[i].srot1, params[i].svec1,
                                          "dihedral_rbp");
      fix_lrf->register_dihedral_junction(i, 1, params[i].subtract_groundstate,
                                          params[i].srot2, params[i].svec2,
                                          "dihedral_rbp");
    }
  }
  #endif

  // for (int i = 1; i <= atom->ndihedraltypes; i++) {
  //   if (!setflag[i])
  //     error->all(FLERR, "Not all dihedral coefficients are set for dihedral style rbp");
  // }
}


/* ----------------------------------------------------------------------
   Write dihedral style parameters to restart file
------------------------------------------------------------------------- */

void DihedralRBP::write_restart(FILE *fp) {
  // Write only primary fields (Ystatic1, Ystatic2, Mmat, subtract_groundstate);
  // derived fields (Smat, srot, svec, M-blocks and their transposes) are
  // recomputed on read.
  for (int i = 1; i <= atom->ndihedraltypes; i++) {
    int sg = params[i].subtract_groundstate ? 1 : 0;
    fwrite(&sg, sizeof(int), 1, fp);
    fwrite(params[i].Ystatic1, sizeof(double), 6, fp);
    fwrite(params[i].Ystatic2, sizeof(double), 6, fp);
    fwrite(&params[i].Mmat[0][0], sizeof(double), 36, fp);
  }
}

/* ----------------------------------------------------------------------
   Read dihedral style parameters from restart file
------------------------------------------------------------------------- */

void DihedralRBP::read_restart(FILE *fp) {
  allocate();

  for (int i = 1; i <= atom->ndihedraltypes; i++) {
    int sg = 0;
    if (comm->me == 0) {
      utils::sfread(FLERR, &sg, sizeof(int), 1, fp, nullptr, error);
      utils::sfread(FLERR, params[i].Ystatic1, sizeof(double), 6, fp, nullptr, error);
      utils::sfread(FLERR, params[i].Ystatic2, sizeof(double), 6, fp, nullptr, error);
      utils::sfread(FLERR, &params[i].Mmat[0][0], sizeof(double), 36, fp, nullptr, error);
    }
    MPI_Bcast(&sg, 1, MPI_INT, 0, world);
    MPI_Bcast(params[i].Ystatic1, 6, MPI_DOUBLE, 0, world);
    MPI_Bcast(params[i].Ystatic2, 6, MPI_DOUBLE, 0, world);
    MPI_Bcast(&params[i].Mmat[0][0], 36, MPI_DOUBLE, 0, world);

    params[i].subtract_groundstate = (sg != 0);

    compute_derived_(i);

    setflag[i] = 1;
  }
}

/* ----------------------------------------------------------------------
   Write dihedral style data (coefficients) to LAMMPS data file
------------------------------------------------------------------------- */

void DihedralRBP::write_data(FILE *fp) {
  for (int i = 1; i <= atom->ndihedraltypes; i++) {
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


void DihedralRBP::allocate() {
  allocated = 1;
  int n = atom->ndihedraltypes;

  memory->create(params,  n+1, "dihedral:rbp:params");
  memory->create(setflag, n+1, "dihedral:rbp:setflag");
  for (int i = 0; i <= n; i++) setflag[i] = 0;
}


void DihedralRBP::assign_coeffs(int dihedral_type, const std::vector<double> &args, bool subtract_groundstate) {


  if (args.size() != 48) {
    error->all(FLERR, "Invalid number of coefficients found for dihedral style rbp. Requires 48 coefficients: X0_1 (6) X0_2 (6) stiffmat (6x6).");
  }

  params[dihedral_type].subtract_groundstate = subtract_groundstate;

  // assign Ystatic
  for (int i=0;i<6;i++) {
    params[dihedral_type].Ystatic1[i] = args[i];
    params[dihedral_type].Ystatic2[i] = args[i+6];
  }

  int argid = 12;
  for (int ii=0;ii<6;ii++) {
    for (int jj=0;jj<6;jj++) {
      params[dihedral_type].Mmat[ii][jj] = args[argid++];
    }
  }

  compute_derived_(dihedral_type);

  setflag[dihedral_type] = 1;
}

/* ----------------------------------------------------------------------
   Recompute derived fields (srot, svec, Smat, M-blocks and transposes)
   from the primary fields (Ystatic1, Ystatic2, Mmat).
------------------------------------------------------------------------- */

void DihedralRBP::compute_derived_(int dihedral_type) {
  for (int ii=0;ii<3;ii++) {
    params[dihedral_type].srot1[ii] = params[dihedral_type].Ystatic1[ii];
    params[dihedral_type].svec1[ii] = params[dihedral_type].Ystatic1[ii+3];
    params[dihedral_type].srot2[ii] = params[dihedral_type].Ystatic2[ii];
    params[dihedral_type].svec2[ii] = params[dihedral_type].Ystatic2[ii+3];
  }

  so3::euler2rotmat(params[dihedral_type].srot1,params[dihedral_type].Smat1);
  so3::euler2rotmat(params[dihedral_type].srot2,params[dihedral_type].Smat2);

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
}


void DihedralRBP::ev_tally_rbp(int i1, int i2, int i3, int i4,
                                int nlocal, int newton_bond,
                                double edihedral,
                                const double *force_1, const double *force_3,
                                const double *dr_a, const double *dr_b,
                                const double *dr32)
{
  // Energy tallying (identical to standard)
  if (eflag_either) {
    if (eflag_global) {
      if (newton_bond)
        energy += edihedral;
      else {
        double eq = 0.25 * edihedral;
        if (i1 < nlocal) energy += eq;
        if (i2 < nlocal) energy += eq;
        if (i3 < nlocal) energy += eq;
        if (i4 < nlocal) energy += eq;
      }
    }
    if (eflag_atom) {
      double eq = 0.25 * edihedral;
      if (newton_bond || i1 < nlocal) eatom[i1] += eq;
      if (newton_bond || i2 < nlocal) eatom[i2] += eq;
      if (newton_bond || i3 < nlocal) eatom[i3] += eq;
      if (newton_bond || i4 < nlocal) eatom[i4] += eq;
    }
  }

  if (vflag_either) {
    // Virial: W = (r1-r2) ⊗ f1 + (r3-r4) ⊗ f3
    //           = -dr_a ⊗ force_1 - dr_b ⊗ force_3
    double v[6];
    v[0] = -(dr_a[0]*force_1[0] + dr_b[0]*force_3[0]);
    v[1] = -(dr_a[1]*force_1[1] + dr_b[1]*force_3[1]);
    v[2] = -(dr_a[2]*force_1[2] + dr_b[2]*force_3[2]);
    v[3] = -(dr_a[0]*force_1[1] + dr_b[0]*force_3[1]);
    v[4] = -(dr_a[0]*force_1[2] + dr_b[0]*force_3[2]);
    v[5] = -(dr_a[1]*force_1[2] + dr_b[1]*force_3[2]);

    if (vflag_global) {
      if (newton_bond) {
        for (int i = 0; i < 6; i++) virial[i] += v[i];
      } else {
        double vq[6];
        for (int i = 0; i < 6; i++) vq[i] = 0.25 * v[i];
        if (i1 < nlocal) for (int i = 0; i < 6; i++) virial[i] += vq[i];
        if (i2 < nlocal) for (int i = 0; i < 6; i++) virial[i] += vq[i];
        if (i3 < nlocal) for (int i = 0; i < 6; i++) virial[i] += vq[i];
        if (i4 < nlocal) for (int i = 0; i < 6; i++) virial[i] += vq[i];
      }
    }

    if (vflag_atom) {
      double vq[6];
      for (int i = 0; i < 6; i++) vq[i] = 0.25 * v[i];
      if (newton_bond || i1 < nlocal)
        for (int i = 0; i < 6; i++) vatom[i1][i] += vq[i];
      if (newton_bond || i2 < nlocal)
        for (int i = 0; i < 6; i++) vatom[i2][i] += vq[i];
      if (newton_bond || i3 < nlocal)
        for (int i = 0; i < 6; i++) vatom[i3][i] += vq[i];
      if (newton_bond || i4 < nlocal)
        for (int i = 0; i < 6; i++) vatom[i4][i] += vq[i];
    }
  }

  // Per-atom centroid virial: (r_a - r_centroid) ⊗ f_a
  // Needs dr32 = minimum_image(r3 - r2) to locate all atoms relative to r2
  if (cvflag_atom) {
    // All positions relative to r2 (using minimum-image displacements)
    double dr12[3], dr42[3];
    for (int i = 0; i < 3; i++) {
      dr12[i] = -dr_a[i];                  // r1 - r2
      dr42[i] = dr32[i] + dr_b[i];         // r4 - r2 = (r3-r2) + (r4-r3)
    }

    // Centroid offset from r2: c = ((r1-r2) + 0 + (r3-r2) + (r4-r2)) / 4
    double c[3];
    for (int i = 0; i < 3; i++)
      c[i] = 0.25 * (dr12[i] + dr32[i] + dr42[i]);

    // r_a - r_centroid for each atom
    double a1[3], a2[3], a3[3], a4[3];
    for (int i = 0; i < 3; i++) {
      a1[i] = dr12[i] - c[i];
      a2[i] =         - c[i];
      a3[i] = dr32[i] - c[i];
      a4[i] = dr42[i] - c[i];
    }

    // Forces: f2 = -force_1, f4 = -force_3
    if (newton_bond || i1 < nlocal) {
      cvatom[i1][0] += a1[0]*force_1[0];
      cvatom[i1][1] += a1[1]*force_1[1];
      cvatom[i1][2] += a1[2]*force_1[2];
      cvatom[i1][3] += a1[0]*force_1[1];
      cvatom[i1][4] += a1[0]*force_1[2];
      cvatom[i1][5] += a1[1]*force_1[2];
      cvatom[i1][6] += a1[1]*force_1[0];
      cvatom[i1][7] += a1[2]*force_1[0];
      cvatom[i1][8] += a1[2]*force_1[1];
    }
    if (newton_bond || i2 < nlocal) {
      cvatom[i2][0] += -a2[0]*force_1[0];
      cvatom[i2][1] += -a2[1]*force_1[1];
      cvatom[i2][2] += -a2[2]*force_1[2];
      cvatom[i2][3] += -a2[0]*force_1[1];
      cvatom[i2][4] += -a2[0]*force_1[2];
      cvatom[i2][5] += -a2[1]*force_1[2];
      cvatom[i2][6] += -a2[1]*force_1[0];
      cvatom[i2][7] += -a2[2]*force_1[0];
      cvatom[i2][8] += -a2[2]*force_1[1];
    }
    if (newton_bond || i3 < nlocal) {
      cvatom[i3][0] += a3[0]*force_3[0];
      cvatom[i3][1] += a3[1]*force_3[1];
      cvatom[i3][2] += a3[2]*force_3[2];
      cvatom[i3][3] += a3[0]*force_3[1];
      cvatom[i3][4] += a3[0]*force_3[2];
      cvatom[i3][5] += a3[1]*force_3[2];
      cvatom[i3][6] += a3[1]*force_3[0];
      cvatom[i3][7] += a3[2]*force_3[0];
      cvatom[i3][8] += a3[2]*force_3[1];
    }
    if (newton_bond || i4 < nlocal) {
      cvatom[i4][0] += -a4[0]*force_3[0];
      cvatom[i4][1] += -a4[1]*force_3[1];
      cvatom[i4][2] += -a4[2]*force_3[2];
      cvatom[i4][3] += -a4[0]*force_3[1];
      cvatom[i4][4] += -a4[0]*force_3[2];
      cvatom[i4][5] += -a4[1]*force_3[2];
      cvatom[i4][6] += -a4[1]*force_3[0];
      cvatom[i4][7] += -a4[2]*force_3[0];
      cvatom[i4][8] += -a4[2]*force_3[1];
    }
  }
}


void DihedralRBP::verify_ev_tally(int i1, int i2, int i3, int i4,
                                   int nlocal, int newton_bond,
                                   double edihedral,
                                   const double *force_1, const double *force_2,
                                   const double *force_3, const double *force_4,
                                   const double *dr_a, const double *dr_b,
                                   const double *dr32)
{
  constexpr double tol = 1e-10;
  bool mismatch = false;

  // ===================================================================
  // Build the inputs that the original ev_tally expects
  // ===================================================================

  // vb1 = r1 - r2 = -dr_a
  double vb1[3] = {-dr_a[0], -dr_a[1], -dr_a[2]};
  // vb2 = r3 - r2 = dr32
  double vb2[3] = {dr32[0], dr32[1], dr32[2]};
  // vb3 = r4 - r3 = dr_b
  double vb3[3] = {dr_b[0], dr_b[1], dr_b[2]};

  // ===================================================================
  // 1. VIRIAL: compare the two formulations
  // ===================================================================

  // --- Original LAMMPS dihedral virial ---
  // v[0] = vb1x*f1[0] + vb2x*f3[0] + (vb3x+vb2x)*f4[0]
  // etc.
  double v_orig[6];
  v_orig[0] = vb1[0]*force_1[0] + vb2[0]*force_3[0] + (vb3[0]+vb2[0])*force_4[0];
  v_orig[1] = vb1[1]*force_1[1] + vb2[1]*force_3[1] + (vb3[1]+vb2[1])*force_4[1];
  v_orig[2] = vb1[2]*force_1[2] + vb2[2]*force_3[2] + (vb3[2]+vb2[2])*force_4[2];
  v_orig[3] = vb1[0]*force_1[1] + vb2[0]*force_3[1] + (vb3[0]+vb2[0])*force_4[1];
  v_orig[4] = vb1[0]*force_1[2] + vb2[0]*force_3[2] + (vb3[0]+vb2[0])*force_4[2];
  v_orig[5] = vb1[1]*force_1[2] + vb2[1]*force_3[2] + (vb3[1]+vb2[1])*force_4[2];

  // --- New simplified virial ---
  // W = -(dr_a ⊗ force_1 + dr_b ⊗ force_3)
  double v_new[6];
  v_new[0] = -(dr_a[0]*force_1[0] + dr_b[0]*force_3[0]);
  v_new[1] = -(dr_a[1]*force_1[1] + dr_b[1]*force_3[1]);
  v_new[2] = -(dr_a[2]*force_1[2] + dr_b[2]*force_3[2]);
  v_new[3] = -(dr_a[0]*force_1[1] + dr_b[0]*force_3[1]);
  v_new[4] = -(dr_a[0]*force_1[2] + dr_b[0]*force_3[2]);
  v_new[5] = -(dr_a[1]*force_1[2] + dr_b[1]*force_3[2]);

  // Compare
  for (int i = 0; i < 6; i++) {
    double scale = std::max(1.0, std::max(std::fabs(v_orig[i]), std::fabs(v_new[i])));
    double diff = std::fabs(v_orig[i] - v_new[i]);
    if (diff > tol * scale) {
      if (!mismatch) {
        error->warning(FLERR, "ev_tally verification FAILED at timestep {}",
                       update->ntimestep);
        mismatch = true;
      }
      fprintf(screen, "  VIRIAL[%d] mismatch: orig=%.15e  new=%.15e  diff=%.3e\n",
              i, v_orig[i], v_new[i], diff);
    }
  }

  // ===================================================================
  // 2. CENTROID VIRIAL: compare the two formulations
  // ===================================================================

  // --- Original LAMMPS centroid virial ---
  // Uses chain topology: r0 = (r1+r2+r3+r4)/4
  // a1 = r1 - r0 = ( 3*vb1 - 2*vb2 -   vb3)/4
  // a2 = r2 - r0 = (  -vb1 - 2*vb2 -   vb3)/4,  f2 = -f1-f3-f4
  // a3 = r3 - r0 = (  -vb1 + 2*vb2 -   vb3)/4
  // a4 = r4 - r0 = (  -vb1 + 2*vb2 + 3*vb3)/4

  double a_orig[4][3];
  for (int i = 0; i < 3; i++) {
    a_orig[0][i] = 0.25 * ( 3.0*vb1[i] - 2.0*vb2[i] -       vb3[i]);
    a_orig[1][i] = 0.25 * (    -vb1[i] - 2.0*vb2[i] -       vb3[i]);
    a_orig[2][i] = 0.25 * (    -vb1[i] + 2.0*vb2[i] -       vb3[i]);
    a_orig[3][i] = 0.25 * (    -vb1[i] + 2.0*vb2[i] + 3.0 * vb3[i]);
  }

  // f2 from original: reconstructed as -f1-f3-f4
  double f2_orig[3];
  for (int i = 0; i < 3; i++)
    f2_orig[i] = -force_1[i] - force_3[i] - force_4[i];

  // Build cv[atom][component] for original: (r_a - r0) ⊗ f_a
  // Components: 0=xx, 1=yy, 2=zz, 3=xy, 4=xz, 5=yz, 6=yx, 7=zx, 8=zy
  const double *f_orig[4] = {force_1, f2_orig, force_3, force_4};
  double cv_orig[4][9];
  for (int a = 0; a < 4; a++) {
    cv_orig[a][0] = a_orig[a][0] * f_orig[a][0];
    cv_orig[a][1] = a_orig[a][1] * f_orig[a][1];
    cv_orig[a][2] = a_orig[a][2] * f_orig[a][2];
    cv_orig[a][3] = a_orig[a][0] * f_orig[a][1];
    cv_orig[a][4] = a_orig[a][0] * f_orig[a][2];
    cv_orig[a][5] = a_orig[a][1] * f_orig[a][2];
    cv_orig[a][6] = a_orig[a][1] * f_orig[a][0];
    cv_orig[a][7] = a_orig[a][2] * f_orig[a][0];
    cv_orig[a][8] = a_orig[a][2] * f_orig[a][1];
  }

  // --- New centroid virial ---
  // All positions relative to r2 via minimum-image displacement vectors
  double dr12[3], dr42[3];
  for (int i = 0; i < 3; i++) {
    dr12[i] = -dr_a[i];                  // r1 - r2
    dr42[i] = dr32[i] + dr_b[i];         // r4 - r2 = (r3-r2) + (r4-r3)
  }

  // Centroid relative to r2
  double c[3];
  for (int i = 0; i < 3; i++)
    c[i] = 0.25 * (dr12[i] + dr32[i] + dr42[i]);

  // Offsets from centroid
  double a_new[4][3];
  for (int i = 0; i < 3; i++) {
    a_new[0][i] = dr12[i] - c[i];       // r1 - centroid
    a_new[1][i] =         - c[i];        // r2 - centroid
    a_new[2][i] = dr32[i] - c[i];       // r3 - centroid
    a_new[3][i] = dr42[i] - c[i];       // r4 - centroid
  }

  // Forces in the new convention: f2 = -force_1 = force_2, f4 = -force_3 = force_4
  // (note: force_2 = -force_1 and force_4 = -force_3 by construction)
  const double *f_new[4] = {force_1, force_2, force_3, force_4};
  double cv_new[4][9];
  for (int a = 0; a < 4; a++) {
    cv_new[a][0] = a_new[a][0] * f_new[a][0];
    cv_new[a][1] = a_new[a][1] * f_new[a][1];
    cv_new[a][2] = a_new[a][2] * f_new[a][2];
    cv_new[a][3] = a_new[a][0] * f_new[a][1];
    cv_new[a][4] = a_new[a][0] * f_new[a][2];
    cv_new[a][5] = a_new[a][1] * f_new[a][2];
    cv_new[a][6] = a_new[a][1] * f_new[a][0];
    cv_new[a][7] = a_new[a][2] * f_new[a][0];
    cv_new[a][8] = a_new[a][2] * f_new[a][1];
  }

  // Compare centroid virials per atom
  const char *atom_labels[4] = {"atom1", "atom2", "atom3", "atom4"};
  for (int a = 0; a < 4; a++) {
    for (int c = 0; c < 9; c++) {
      double scale = std::max(1.0, std::max(std::fabs(cv_orig[a][c]), std::fabs(cv_new[a][c])));
      double diff = std::fabs(cv_orig[a][c] - cv_new[a][c]);
      if (diff > tol * scale) {
        if (!mismatch) {
          error->warning(FLERR, "ev_tally verification FAILED at timestep {}",
                         update->ntimestep);
          mismatch = true;
        }
        fprintf(screen, "  CENTROID_VIRIAL[%s][%d] mismatch: orig=%.15e  new=%.15e  diff=%.3e\n",
                atom_labels[a], c, cv_orig[a][c], cv_new[a][c], diff);
      }
    }
  }

  // ===================================================================
  // 3. Also verify that the total centroid virial sums match the virial
  // ===================================================================
  // sum_a (r_a - r0) ⊗ f_a  should equal the total virial tensor
  // (only the symmetric part, components 0-5)

  double cv_sum_orig[6] = {0, 0, 0, 0, 0, 0};
  double cv_sum_new[6]  = {0, 0, 0, 0, 0, 0};
  for (int a = 0; a < 4; a++) {
    cv_sum_orig[0] += cv_orig[a][0];
    cv_sum_orig[1] += cv_orig[a][1];
    cv_sum_orig[2] += cv_orig[a][2];
    cv_sum_orig[3] += cv_orig[a][3];
    cv_sum_orig[4] += cv_orig[a][4];
    cv_sum_orig[5] += cv_orig[a][5];

    cv_sum_new[0] += cv_new[a][0];
    cv_sum_new[1] += cv_new[a][1];
    cv_sum_new[2] += cv_new[a][2];
    cv_sum_new[3] += cv_new[a][3];
    cv_sum_new[4] += cv_new[a][4];
    cv_sum_new[5] += cv_new[a][5];
  }

  for (int i = 0; i < 6; i++) {
    double scale = std::max(1.0, std::max(std::fabs(cv_sum_orig[i]), std::fabs(v_orig[i])));
    double diff_orig = std::fabs(cv_sum_orig[i] - v_orig[i]);
    double diff_new  = std::fabs(cv_sum_new[i]  - v_new[i]);
    if (diff_orig > tol * scale) {
      if (!mismatch) {
        error->warning(FLERR, "ev_tally verification FAILED at timestep {}",
                       update->ntimestep);
        mismatch = true;
      }
      fprintf(screen, "  ORIG sum(centroid_virial)[%d] != virial[%d]: sum=%.15e  vir=%.15e  diff=%.3e\n",
              i, i, cv_sum_orig[i], v_orig[i], diff_orig);
    }
    if (diff_new > tol * scale) {
      if (!mismatch) {
        error->warning(FLERR, "ev_tally verification FAILED at timestep {}",
                       update->ntimestep);
        mismatch = true;
      }
      fprintf(screen, "  NEW  sum(centroid_virial)[%d] != virial[%d]: sum=%.15e  vir=%.15e  diff=%.3e\n",
              i, i, cv_sum_new[i], v_new[i], diff_new);
    }
  }

  // ===================================================================
  // 4. Verify force balance: f1 + f2 + f3 + f4 = 0
  // ===================================================================
  for (int i = 0; i < 3; i++) {
    double fsum = force_1[i] + force_2[i] + force_3[i] + force_4[i];
    if (std::fabs(fsum) > tol) {
      if (!mismatch) {
        error->warning(FLERR, "ev_tally verification FAILED at timestep {}",
                       update->ntimestep);
        mismatch = true;
      }
      fprintf(screen, "  FORCE BALANCE[%d] violated: sum=%.15e\n", i, fsum);
    }
  }

  // ===================================================================
  // Summary
  // ===================================================================
  if (!mismatch) {
    // Only print on first call to avoid flooding the log
    static bool first_pass = true;
    if (first_pass) {
      fprintf(screen, "  ev_tally verification PASSED (timestep %ld)\n",
              (long)update->ntimestep);
      first_pass = false;
    }
  }
}