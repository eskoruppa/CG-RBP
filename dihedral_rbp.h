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


/* dihedral_rbp.h ----------------------------------------------- */
#ifdef DIHEDRAL_CLASS
// AngleStyle(rbp/gs,DihedralRBP);
DihedralStyle(rbp,DihedralRBP);
#else

#ifndef LMP_DIHEDRAL_RBP_H
#define LMP_DIHEDRAL_RBP_H

// clang-format off
#define RBP_DIHEDRAL_DEFAULT_SUBTRACT_GROUNDSTATE false
#define RBP_DIHEDRAL_USE_CUSTOM_EV_TALLY
#define DIHEDRAL_RBP_PRECOMPUTE_ACTIVE
// clang-format on

#include "dihedral.h"
#include "so3.h"
#include "parse_rbp.h"

namespace LAMMPS_NS {

static constexpr const char* RBPDIHEDRAL_STYLE = "rbp";

class DihedralRBP : public Dihedral {
 public:
  // DihedralRBP(class LAMMPS *, int, char **);
  DihedralRBP(LAMMPS *lmp);
  ~DihedralRBP() override;

  void compute(int, int) override;
  void coeff(int, char **) override;
  void init_style() override;

  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_data(FILE *) override;

 protected:
   struct RBPParams {
      double Ystatic1[6];     // Static wrench
      double Ystatic2[6];     // Static wrench
      double Smat1[3][3];     // Static rotation
      double Smat2[3][3];     // Static rotation
      double srot1[3];        // Static rotation
      double srot2[3];        // Static rotation
      double svec1[3];        // Static translation
      double svec2[3];        // Static translation
      double Mmat[6][6];      // Stiffness matrix
      double Mrr[3][3];       // rotational part of M
      double Mtt[3][3];       // translational part of M
      double Mrt[3][3];       // top-right cross terms
      double Mtr[3][3];       // bottom-left cross terms
      double Mrr_tp[3][3];    // transpose of Mrr
      double Mtt_tp[3][3];    // transpose of Mtt
      double Mrt_tp[3][3];    // transpose of Mrt
      double Mtr_tp[3][3];    // transpose of Mtr
      bool   subtract_groundstate;
      //  double equidist;
   };
   RBPParams *params;  // indexed by bond type
   void allocate();

   void assign_coeffs(int angle_type, const std::vector<double> &args, bool subtract_groundstate = false);
   void compute_derived_(int dihedral_type);
   void ev_tally_rbp(int i1, int i2, int i3, int i4, int nlocal, int newton_bond,
                double edihedral, const double *force_1, const double *force_3,
                const double *dr_a, const double *dr_b, const double *dr32);

   void verify_ev_tally(int i1, int i2, int i3, int i4, int nlocal, int newton_bond,
                        double edihedral, const double *force_1, const double *force_2,
                        const double *force_3, const double *force_4, const double *dr_a, 
                        const double *dr_b, const double *dr32);

#ifdef DIHEDRAL_RBP_PRECOMPUTE_ACTIVE
   class FixRBPLRF *fix_lrf;
#endif
};

}

#endif
#endif
