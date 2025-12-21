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


/* angle_rbp.h ----------------------------------------------- */
#ifdef ANGLE_CLASS
// AngleStyle(rbp/gs,AngleRBP);
AngleStyle(rbp,AngleRBP);
#else

#ifndef LMP_ANGLE_RBP_H
#define LMP_ANGLE_RBP_H

#define CHECK_TORQUE_BALANCE

#include "angle.h"
#include "so3.h"
#include "parse_rbp.h"

namespace LAMMPS_NS {

class AngleRBP : public Angle {
 public:
  // BondRBP(class LAMMPS *, int, char **);
  AngleRBP(LAMMPS *lmp);
  ~AngleRBP() override;
  void compute(int, int) override;
  void coeff(int, char **) override;
  void init_style() override;
  double equilibrium_angle(int) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_data(FILE *) override;
  double single(int, int, int, int) override;

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
    double Mmat[6][6];     // Stiffness matrix
    double Mrr[3][3];       // rotational part of M
    double Mtt[3][3];       // translational part of M
    double Mrt[3][3];   // bottom-left cross terms
    double Mtr[3][3];   // top-right cross terms
    double Mrr_tp[3][3];       // rotational part of M
    double Mtt_tp[3][3];       // translational part of M
    double Mrt_tp[3][3];   // bottom-left cross terms
    double Mtr_tp[3][3];   // top-right cross terms
   //  double equidist;
  };
  RBPParams *params;  // indexed by bond type
  void allocate();

  void assign_coeffs(int angle_type, const std::vector<double> &args);
};

}

#endif
#endif
