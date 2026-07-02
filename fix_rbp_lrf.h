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

#ifdef FIX_CLASS
// clang-format off
FixStyle(rbp/lrf,FixRBPLRF);
// clang-format on
#else

#ifndef LMP_FIX_RBP_LRF_H
#define LMP_FIX_RBP_LRF_H

#include "fix.h"
#include <string>

namespace LAMMPS_NS {

class FixRBPLRF : public Fix {
 public:
  FixRBPLRF(class LAMMPS *, int, char **);
  ~FixRBPLRF() override;
  int setmask() override;
  void min_pre_force(int) override;
  void setup_pre_force(int) override;
  void pre_force(int) override;

  double memory_usage() override;
  void grow_arrays(int) override;
  void copy_arrays(int, int, int) override;
  void set_arrays(int) override;
  int pack_exchange(int, double *) override;
  int unpack_exchange(int, double *) override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;

  // --- Junction registration (called from interaction style init_style) ---
  // bond_type : 1-based index
  // is_x      : true = X convention (subtract_groundstate), false = Y
  // srot/svec : ground-state rotation / translation vectors
  // caller    : name shown in error messages
  void register_bond_junction(int bond_type, bool is_x,
                              const double srot[3], const double svec[3],
                              const std::string &caller);
  void register_angle_junction(int angle_type, int sub_pair, bool is_x,
                               const double srot[3], const double svec[3],
                               const std::string &caller);
  void register_dihedral_junction(int dih_type, int sub_pair, bool is_x,
                                  const double srot[3], const double svec[3],
                                  const std::string &caller);

  // --- Per-atom arrays (read by interaction styles via pointer cast) ---
  // triads[i][0..8]      = row-major 3x3 triad from quaternion
  // fwd_euler[i][0..2]   = Om (X conv) or Phi_delta (Y conv)
  // fwd_Jinvtp[i][0..8]  = leftJacobianInverseTransposed of the euler vector
  double **triads;
  double **fwd_euler;
  double **fwd_Jinvtp;

 private:
  void compute_lrf();
  void validate_junctions();

  class AtomVecEllipsoid *avec;

  // forward-comm mode: 0 = lrf data (fwd_euler + fwd_Jinvtp), 1 = bond types
  // for one-time junction validation (see validate_junctions)
  int comm_mode;
  double *val_btype;    // transient per-atom bond type used during validation

  // Per-bond-type junction parameters (1-indexed, size nbond_specs+1)
  struct JunctionSpec {
    bool registered;
    bool is_x;
    double srot[3];
    double svec[3];
    double Smat[3][3];        // euler2rotmat(srot), for Y convention
    std::string first_caller;
  };
  JunctionSpec *bond_specs;
  int nbond_specs;

  // Pending angle / dihedral sub-junction specs for validation
  // Flattened: index = type * 2 + sub_pair  (sub_pair = 0 or 1)
  struct PendingSpec {
    bool registered;
    bool is_x;
    double srot[3];
    double svec[3];
    std::string caller;
  };
  PendingSpec *angle_specs;
  int nangle_specs;
  PendingSpec *dihedral_specs;
  int ndihedral_specs;

  bool validated;
};

}    // namespace LAMMPS_NS
#endif
#endif
