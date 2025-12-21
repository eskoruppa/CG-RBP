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


/* bond_rbp.h ----------------------------------------------- */
#ifdef BOND_CLASS
BondStyle(rbpfene,BondRBPFene);
#else

#ifndef LMP_BOND_RBPFENE_H
#define LMP_BOND_RBPFENE_H

#include "bond.h"
#include "so3.h"

namespace LAMMPS_NS {

class BondRBPFene : public Bond {
 public:
  // BondRBP(class LAMMPS *, int, char **);
  BondRBPFene(LAMMPS *lmp);
  ~BondRBPFene() override;
  void compute(int, int) override;
  void coeff(int, char **) override;
  void init_style() override;
  double equilibrium_distance(int) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_data(FILE *) override;
  double single(int, double, int, int, double &) override;

 protected:
  struct RBPParams {
    double Ystatic[6];     // Static wrench
    double Smat[3][3];     // Static rotation
    double srot[3];        // Static rotation
    double svec[3];        // Static translation
    double Mmat[6][6];     // Stiffness matrix
    double Mr[3][3];       // rotational part of M
    double Mt[3][3];       // translational part of M
    double Mtr_bl[3][3];   // bottom-left cross terms
    double Mtr_tr[3][3];   // top-right cross terms
    double equidist;
    double K;
    double r0;
    double epsilon;
    double sigma;
  };
  RBPParams *params;  // indexed by bond type
  void allocate();
};

}

#endif
#endif
