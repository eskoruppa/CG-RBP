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

#include "fix_rbp_lrf.h"

#include "atom.h"
#include "atom_vec_ellipsoid.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "math_extra.h"
#include "memory.h"
#include "neighbor.h"
#include "update.h"
#include "so3.h"
#include "lamath.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixRBPLRF::FixRBPLRF(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg),
    triads(nullptr), fwd_euler(nullptr), fwd_Jinvtp(nullptr),
    avec(nullptr),
    comm_mode(0), val_btype(nullptr),
    bond_specs(nullptr), nbond_specs(0),
    angle_specs(nullptr), nangle_specs(0),
    dihedral_specs(nullptr), ndihedral_specs(0),
    validated(false)
{
  // forward comm: fwd_euler(3) + fwd_Jinvtp(9) = 12
  comm_forward = 12;
  peratom_flag = 0;

  FixRBPLRF::grow_arrays(atom->nmax);
  atom->add_callback(Atom::GROW);

  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    memset(triads[i], 0, 9 * sizeof(double));
    memset(fwd_euler[i], 0, 3 * sizeof(double));
    memset(fwd_Jinvtp[i], 0, 9 * sizeof(double));
  }
}

/* ---------------------------------------------------------------------- */

FixRBPLRF::~FixRBPLRF()
{
  atom->delete_callback(id, Atom::GROW);
  memory->destroy(triads);
  memory->destroy(fwd_euler);
  memory->destroy(fwd_Jinvtp);

  delete[] bond_specs;
  delete[] angle_specs;
  delete[] dihedral_specs;
}

/* ---------------------------------------------------------------------- */

int FixRBPLRF::setmask()
{
  int mask = 0;
  mask |= MIN_PRE_FORCE;
  mask |= PRE_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixRBPLRF::min_pre_force(int /*vflag*/)
{
  compute_lrf();
}

/* ---------------------------------------------------------------------- */

void FixRBPLRF::setup_pre_force(int vflag)
{
  pre_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixRBPLRF::pre_force(int /*vflag*/)
{
  compute_lrf();
}

/* ---------------------------------------------------------------------- */

double FixRBPLRF::memory_usage()
{
  int nmax = atom->nmax;
  double bytes = 0.0;
  bytes += (double)nmax * 9 * sizeof(double);   // triads
  bytes += (double)nmax * 3 * sizeof(double);   // fwd_euler
  bytes += (double)nmax * 9 * sizeof(double);   // fwd_Jinvtp
  return bytes;
}

/* ---------------------------------------------------------------------- */

void FixRBPLRF::grow_arrays(int nmax)
{
  memory->grow(triads, nmax, 9, "fix_rbp_lrf:triads");
  memory->grow(fwd_euler, nmax, 3, "fix_rbp_lrf:fwd_euler");
  memory->grow(fwd_Jinvtp, nmax, 9, "fix_rbp_lrf:fwd_Jinvtp");
}

/* ---------------------------------------------------------------------- */

void FixRBPLRF::copy_arrays(int i, int j, int /*delflag*/)
{
  memcpy(triads[j], triads[i], 9 * sizeof(double));
  memcpy(fwd_euler[j], fwd_euler[i], 3 * sizeof(double));
  memcpy(fwd_Jinvtp[j], fwd_Jinvtp[i], 9 * sizeof(double));
}

/* ---------------------------------------------------------------------- */

void FixRBPLRF::set_arrays(int i)
{
  memset(triads[i], 0, 9 * sizeof(double));
  memset(fwd_euler[i], 0, 3 * sizeof(double));
  memset(fwd_Jinvtp[i], 0, 9 * sizeof(double));
}

/* ---------------------------------------------------------------------- */

int FixRBPLRF::pack_exchange(int i, double *buf)
{
  int m = 0;
  for (int k = 0; k < 9; k++) buf[m++] = triads[i][k];
  for (int k = 0; k < 3; k++) buf[m++] = fwd_euler[i][k];
  for (int k = 0; k < 9; k++) buf[m++] = fwd_Jinvtp[i][k];
  return m;  // 21
}

/* ---------------------------------------------------------------------- */

int FixRBPLRF::unpack_exchange(int nlocal, double *buf)
{
  int m = 0;
  for (int k = 0; k < 9; k++) triads[nlocal][k] = buf[m++];
  for (int k = 0; k < 3; k++) fwd_euler[nlocal][k] = buf[m++];
  for (int k = 0; k < 9; k++) fwd_Jinvtp[nlocal][k] = buf[m++];
  return m;  // 21
}

/* ---------------------------------------------------------------------- */

int FixRBPLRF::pack_forward_comm(int n, int *list, double *buf,
                                 int /*pbc_flag*/, int * /*pbc*/)
{
  int m = 0;
  if (comm_mode == 1) {
    // one-time validation comm: push owned bond types to ghosts
    for (int i = 0; i < n; i++) buf[m++] = val_btype[list[i]];
    return m;
  }
  for (int i = 0; i < n; i++) {
    int j = list[i];
    for (int k = 0; k < 3; k++) buf[m++] = fwd_euler[j][k];
    for (int k = 0; k < 9; k++) buf[m++] = fwd_Jinvtp[j][k];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void FixRBPLRF::unpack_forward_comm(int n, int first, double *buf)
{
  int m = 0;
  int last = first + n;
  if (comm_mode == 1) {
    for (int i = first; i < last; i++) val_btype[i] = buf[m++];
    return;
  }
  for (int i = first; i < last; i++) {
    for (int k = 0; k < 3; k++) fwd_euler[i][k] = buf[m++];
    for (int k = 0; k < 9; k++) fwd_Jinvtp[i][k] = buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

void FixRBPLRF::register_bond_junction(int bond_type, bool is_x,
                                       const double srot[3],
                                       const double svec[3],
                                       const std::string &caller)
{
  // Lazy allocation
  if (!bond_specs) {
    nbond_specs = atom->nbondtypes;
    bond_specs = new JunctionSpec[nbond_specs + 1];
    for (int i = 0; i <= nbond_specs; i++) bond_specs[i].registered = false;
  }

  if (bond_type < 1 || bond_type > nbond_specs)
    error->all(FLERR, "fix rbp/lrf: bond_type {} out of range in register_bond_junction "
               "from {}", bond_type, caller);

  JunctionSpec &js = bond_specs[bond_type];
  if (js.registered) {
    // Validate consistency
    if (js.is_x != is_x)
      error->all(FLERR, "fix rbp/lrf: convention mismatch for bond type {} — "
                 "{} says {} but {} says {}",
                 bond_type, js.first_caller, (js.is_x ? "X" : "Y"),
                 caller, (is_x ? "X" : "Y"));
    constexpr double tol = 1e-10;
    for (int k = 0; k < 3; k++) {
      if (std::fabs(js.srot[k] - srot[k]) > tol)
        error->all(FLERR, "fix rbp/lrf: srot mismatch for bond type {} between "
                   "{} and {}", bond_type, js.first_caller, caller);
      if (std::fabs(js.svec[k] - svec[k]) > tol)
        error->all(FLERR, "fix rbp/lrf: svec mismatch for bond type {} between "
                   "{} and {}", bond_type, js.first_caller, caller);
    }
    return;  // already registered with matching params
  }

  js.registered = true;
  js.is_x = is_x;
  for (int k = 0; k < 3; k++) { js.srot[k] = srot[k]; js.svec[k] = svec[k]; }
  so3::euler2rotmat(js.srot, js.Smat);
  js.first_caller = caller;
}

/* ---------------------------------------------------------------------- */

void FixRBPLRF::register_angle_junction(int angle_type, int sub_pair,
                                        bool is_x,
                                        const double srot[3],
                                        const double svec[3],
                                        const std::string &caller)
{
  if (!angle_specs) {
    nangle_specs = atom->nangletypes;
    angle_specs = new PendingSpec[(nangle_specs + 1) * 2];
    for (int i = 0; i < (nangle_specs + 1) * 2; i++)
      angle_specs[i].registered = false;
  }
  if (angle_type < 1 || angle_type > nangle_specs || sub_pair < 0 || sub_pair > 1)
    error->all(FLERR, "fix rbp/lrf: invalid angle_type/sub_pair in "
               "register_angle_junction from {}", caller);

  PendingSpec &ps = angle_specs[angle_type * 2 + sub_pair];
  ps.registered = true;
  ps.is_x = is_x;
  for (int k = 0; k < 3; k++) { ps.srot[k] = srot[k]; ps.svec[k] = svec[k]; }
  ps.caller = caller;
}

/* ---------------------------------------------------------------------- */

void FixRBPLRF::register_dihedral_junction(int dih_type, int sub_pair,
                                           bool is_x,
                                           const double srot[3],
                                           const double svec[3],
                                           const std::string &caller)
{
  if (!dihedral_specs) {
    ndihedral_specs = atom->ndihedraltypes;
    dihedral_specs = new PendingSpec[(ndihedral_specs + 1) * 2];
    for (int i = 0; i < (ndihedral_specs + 1) * 2; i++)
      dihedral_specs[i].registered = false;
  }
  if (dih_type < 1 || dih_type > ndihedral_specs || sub_pair < 0 || sub_pair > 1)
    error->all(FLERR, "fix rbp/lrf: invalid dih_type/sub_pair in "
               "register_dihedral_junction from {}", caller);

  PendingSpec &ps = dihedral_specs[dih_type * 2 + sub_pair];
  ps.registered = true;
  ps.is_x = is_x;
  for (int k = 0; k < 3; k++) { ps.srot[k] = srot[k]; ps.svec[k] = svec[k]; }
  ps.caller = caller;
}

/* ---------------------------------------------------------------------- */

void FixRBPLRF::validate_junctions()
{
  if (validated) return;
  validated = true;

  if (!bond_specs) {
    // No RBP bond junctions registered. Phase 2 of compute_lrf only fills
    // fwd_euler / fwd_Jinvtp from the bond list, so without any registered RBP
    // bond those arrays stay zero. If angle / dihedral junctions were
    // registered, they would silently read those zeros and produce wrong
    // forces and energies. Fail loudly instead.
    const char *offender = nullptr;
    if (angle_specs) {
      for (int i = 0; i < (nangle_specs + 1) * 2 && !offender; i++)
        if (angle_specs[i].registered) offender = angle_specs[i].caller.c_str();
    }
    if (dihedral_specs && !offender) {
      for (int i = 0; i < (ndihedral_specs + 1) * 2 && !offender; i++)
        if (dihedral_specs[i].registered) offender = dihedral_specs[i].caller.c_str();
    }
    if (offender)
      error->all(FLERR, "fix rbp/lrf: {} registered an RBP angle/dihedral "
                 "junction but no RBP bond style is present to precompute its "
                 "triad/euler data. Use an RBP bond style for these steps, or "
                 "disable the *_PRECOMPUTE_ACTIVE define for the angle/dihedral "
                 "style to fall back to the (bond-independent) local path.",
                 offender);
    // Otherwise nothing references the fix; harmless.
    return;
  }

  constexpr double tol = 1e-10;

  // -------------------------------------------------------------------
  // Build a ghost-aware atom -> bond_type map.
  //
  // Each bond is listed (bondlist[*][0] = local left atom) only on the proc
  // that owns its left atom, so a purely local map misses steps whose left
  // atom is a ghost here. An angle/dihedral, however, is owned by its central
  // atom's proc, where the step's left atom (id1 for sub 0) is typically a
  // ghost. We therefore fill the map for local left atoms and forward-comm it
  // to ghosts, so every junction can be validated on the proc that owns it.
  // -------------------------------------------------------------------

  int nall = atom->nlocal + atom->nghost;
  val_btype = new double[nall];
  for (int i = 0; i < nall; i++) val_btype[i] = 0.0;

  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;
  for (int bid = 0; bid < nbondlist; bid++)
    val_btype[bondlist[bid][0]] = (double) bondlist[bid][2];

  comm_mode = 1;
  comm->forward_comm(this);
  comm_mode = 0;

  // Helper lambda: compare pending spec against the bond spec at its step.
  // btype is the RBP bond type sitting at the junction's left atom (0 = none).
  auto check_match = [&](const PendingSpec &ps, int btype, int interaction_type,
                         int sub_pair, const char *style_name) {
    if (btype == 0)
      error->one(FLERR, "fix rbp/lrf: {} (type {} sub {}) sits on a step with no "
                 "RBP bond — the precompute fix has no triad/euler source for it. "
                 "Use an RBP bond style for these steps, or disable the "
                 "*_PRECOMPUTE_ACTIVE define to use the local path.",
                 style_name, interaction_type, sub_pair);
    if (btype < 1 || btype > nbond_specs)
      error->one(FLERR, "fix rbp/lrf: {} (type {} sub {}) maps to bond type {} "
                 "which is out of range",
                 style_name, interaction_type, sub_pair, btype);
    const JunctionSpec &bs = bond_specs[btype];
    if (!bs.registered)
      error->one(FLERR, "fix rbp/lrf: {} (type {} sub {}) references bond type {} "
                 "which has no registered RBP bond junction",
                 style_name, interaction_type, sub_pair, btype);

    if (ps.is_x != bs.is_x)
      error->one(FLERR, "fix rbp/lrf: convention mismatch for bond type {} — "
                 "{} says {} but {} (type {} sub {}) says {}",
                 btype, bs.first_caller, (bs.is_x ? "X" : "Y"),
                 style_name, interaction_type, sub_pair,
                 (ps.is_x ? "X" : "Y"));
    for (int k = 0; k < 3; k++) {
      if (std::fabs(ps.srot[k] - bs.srot[k]) > tol)
        error->one(FLERR, "fix rbp/lrf: srot[{}] mismatch for bond type {} — "
                   "{} has {:.6e} but {} (type {} sub {}) has {:.6e}",
                   k, btype, bs.first_caller, bs.srot[k],
                   style_name, interaction_type, sub_pair, ps.srot[k]);
      if (std::fabs(ps.svec[k] - bs.svec[k]) > tol)
        error->one(FLERR, "fix rbp/lrf: svec[{}] mismatch for bond type {} — "
                   "{} has {:.6e} but {} (type {} sub {}) has {:.6e}",
                   k, btype, bs.first_caller, bs.svec[k],
                   style_name, interaction_type, sub_pair, ps.svec[k]);
    }
  };

  // -------------------------------------------------------------------
  // Validate every LOCAL angle/dihedral instance against the bond at its step.
  //
  // Each instance is owned by exactly one proc (its central atom), so checking
  // local instances validates the whole system across procs with no double
  // counting. We do NOT require every registered *type* to be present locally:
  // under spatial decomposition each proc registers all types but owns only a
  // contiguous subset of instances, so a per-type "must be present" test would
  // wrongly fail (e.g. sequence-dependent models with one type per position).
  // -------------------------------------------------------------------

  // Validate angle sub-junctions
  if (angle_specs && neighbor->nanglelist > 0) {
    int **anglelist = neighbor->anglelist;
    int nanglelist = neighbor->nanglelist;

    for (int aid = 0; aid < nanglelist; aid++) {
      int id1 = anglelist[aid][0];   // left atom of sub-pair 0 step (id1,id2)
      int id2 = anglelist[aid][1];   // left atom of sub-pair 1 step (id2,id3)
      int atype = anglelist[aid][3];

      if (angle_specs[atype * 2].registered)
        check_match(angle_specs[atype * 2], (int) val_btype[id1],
                    atype, 0, "angle_rbp");
      if (angle_specs[atype * 2 + 1].registered)
        check_match(angle_specs[atype * 2 + 1], (int) val_btype[id2],
                    atype, 1, "angle_rbp");
    }
  }

  // Validate dihedral sub-junctions
  if (dihedral_specs && neighbor->ndihedrallist > 0) {
    int **dihedrallist = neighbor->dihedrallist;
    int ndihedrallist = neighbor->ndihedrallist;

    for (int did = 0; did < ndihedrallist; did++) {
      int id1 = dihedrallist[did][0];   // left atom of sub-pair 0 step (id1,id2)
      int id3 = dihedrallist[did][2];   // left atom of sub-pair 1 step (id3,id4)
      int dtype = dihedrallist[did][4];

      if (dihedral_specs[dtype * 2].registered)
        check_match(dihedral_specs[dtype * 2], (int) val_btype[id1],
                    dtype, 0, "dihedral_rbp");
      if (dihedral_specs[dtype * 2 + 1].registered)
        check_match(dihedral_specs[dtype * 2 + 1], (int) val_btype[id3],
                    dtype, 1, "dihedral_rbp");
    }
  }

  delete[] val_btype;
  val_btype = nullptr;
}

/* ---------------------------------------------------------------------- */

void FixRBPLRF::compute_lrf()
{
  avec = dynamic_cast<AtomVecEllipsoid *>(atom->style_match("ellipsoid"));
  AtomVecEllipsoid::Bonus *bonus = avec->bonus;
  int *ellipsoid = atom->ellipsoid;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;

  // -------------------------------------------------------------------
  // Phase 1: Compute triads for ALL atoms (local + ghost).
  //          Ghost quaternions are valid after the Verlet forward_comm.
  // -------------------------------------------------------------------

  for (int i = 0; i < nall; i++) {
    int n = ellipsoid[i];
    if (n < 0) continue;

    double *q = bonus[n].quat;
    double (*T)[3] = (double(*)[3]) triads[i];
    MathExtra::quat_to_mat(q, T);
  }

  // -------------------------------------------------------------------
  // Phase 2: Compute per-bond euler vector and Jinvtp.
  //          For each bond (id1, id2, bond_type):
  //            R = T1^T * T2
  //            X convention: Om = rotmat2euler(R)
  //                          Jinvtp = leftJacInvTp(Om)
  //            Y convention: D = Smat^T * R
  //                          Phi = rotmat2euler(D)
  //                          Jinvtp = leftJacInvTp(Phi)
  //
  // ORIENTATION REQUIREMENT (shared by the local-loop path as well):
  //   Every RBP step quantity is defined relative to the *left* triad of the
  //   step (w = T1^T (r2-r1), R = T1^T T2, Om = log(R)). Swapping the two base
  //   pairs of a step does NOT merely negate these, so step orientation is
  //   physically meaningful. Consequently the precomputed values are keyed to
  //   the bond's left atom id1 = bondlist[bid][0], and angle_rbp/dihedral_rbp
  //   read fwd_euler[id1]/[id2]/[id3] assuming the atom at the *start* of each
  //   of their sub-steps is the left atom of the matching bond. This holds only
  //   if bonds, angles and dihedrals are all declared with a consistent
  //   head->tail (5'->3') orientation. A reversed bond is caught as a hard
  //   error in validate_junctions ("no corresponding RBP bond found").
  //   It is also assumed that each atom is the left endpoint of at most one
  //   RBP bond (otherwise fwd_euler[id1] would be overwritten); this holds for
  //   a linear backbone.
  // -------------------------------------------------------------------

  if (bond_specs) {
    int **bondlist = neighbor->bondlist;
    int nbondlist = neighbor->nbondlist;

    for (int bid = 0; bid < nbondlist; bid++) {
      int id1 = bondlist[bid][0];
      int id2 = bondlist[bid][1];
      int btype = bondlist[bid][2];

      if (btype < 1 || btype > nbond_specs || !bond_specs[btype].registered)
        continue;

      // Triads as 3x3 (pointer cast, zero copy)
      double (*T1)[3] = (double(*)[3]) triads[id1];
      double (*T2)[3] = (double(*)[3]) triads[id2];

      // R = T1^T * T2  (local temporary)
      double R[3][3];
      lamath::mul_AtB(T1, T2, R);

      double euler[3];
      if (bond_specs[btype].is_x) {
        // X convention: Om = rotmat2euler(R)
        so3::rotmat2euler(R, euler);
      } else {
        // Y convention: D = Smat^T * R, Phi = rotmat2euler(D)
        double D[3][3];
        lamath::mul_AtB(bond_specs[btype].Smat, R, D);
        so3::rotmat2euler(D, euler);
      }

      fwd_euler[id1][0] = euler[0];
      fwd_euler[id1][1] = euler[1];
      fwd_euler[id1][2] = euler[2];

      // Jinvtp = leftJacobianInverseTransposed(euler)
      double (*Jinvtp)[3] = (double(*)[3]) fwd_Jinvtp[id1];
      so3::leftJacobianInverseTransposed(euler, Jinvtp);
    }
  }

  // -------------------------------------------------------------------
  // Phase 3: Validate junction parameter consistency (first call only).
  // -------------------------------------------------------------------

  validate_junctions();

  // -------------------------------------------------------------------
  // Phase 4: Communicate fwd_euler and fwd_Jinvtp to ghost atoms.
  // -------------------------------------------------------------------

  comm->forward_comm(this);
}
