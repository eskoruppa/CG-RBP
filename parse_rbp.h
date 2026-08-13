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

#ifndef LAMMPS_RBP_PARSE_DATABASE_H
#define LAMMPS_RBP_PARSE_DATABASE_H

#include <vector>
#include <string>

#include "lammps.h"

// clang-format off
// ============================================================
// Mandatory metadata fields.
// Comment out any define to make that field optional —
// no error will be raised if it is absent from the database file.
// ============================================================
#define RBP_META_REQUIRE_METADATA_SECTION   
#define RBP_META_REQUIRE_NUM_RIGID_BODIES
#define RBP_META_REQUIRE_COUPLING_RANGE
#define RBP_META_REQUIRE_NUM_BONDS
#define RBP_META_REQUIRE_NUM_ANGLES
#define RBP_META_REQUIRE_NUM_DIHEDRALS
#define RBP_META_REQUIRE_NUM_BOND_TYPES
#define RBP_META_REQUIRE_NUM_ANGLE_TYPES
#define RBP_META_REQUIRE_NUM_DIHEDRAL_TYPES
#define RBP_META_REQUIRE_CHARS_PER_ATOM
#define RBP_META_REQUIRE_BOND_STYLE
#define RBP_META_REQUIRE_ANGLE_STYLE
#define RBP_META_REQUIRE_DIHEDRAL_STYLE
#define RBP_META_REQUIRE_SUBTRACT_GROUNDSTATE
#define RBP_META_REQUIRE_SEQS_SET
#define RBP_META_REQUIRE_SEQS_CENTERED
#define RBP_META_REQUIRE_CLOSED
#define RBP_META_REQUIRE_UNIT_LENGTH
// #define RBP_META_REQUIRE_UNIT_ENERGY
// #define RBP_META_REQUIRE_ODD
// ============================================================
// Default values for optional metadata fields.
// Applied when a field is absent and not mandatory.
// ============================================================
#define RBP_META_DEFAULT_SUBTRACT_GROUNDSTATE false
#define RBP_META_DEFAULT_ODD                  false
#define RBP_META_DEFAULT_SEQS_SET             false
#define RBP_META_DEFAULT_SEQS_CENTERED        false
#define RBP_META_DEFAULT_CLOSED               false
#define RBP_META_DEFAULT_UNIT_LENGTH          1.0
#define RBP_META_DEFAULT_UNIT_ENERGY          1.0
// ============================================================
// Identifier strings for database file fields and sections. 
// These must match the strings used in the database files.
// ============================================================
// metadata field identifiers
#define RBP_META_IDENTIFIER_NUM_RIGID_BODIES "number of rigid bodies"
#define RBP_META_IDENTIFIER_COUPLING_RANGE "coupling range"
#define RBP_META_IDENTIFIER_NUM_BONDS "number of bonds"
#define RBP_META_IDENTIFIER_NUM_ANGLES "number of angles"
#define RBP_META_IDENTIFIER_NUM_DIHEDRALS "number of dihedrals"
#define RBP_META_IDENTIFIER_NUM_BOND_TYPES "number of bond types"
#define RBP_META_IDENTIFIER_NUM_ANGLE_TYPES "number of angle types"
#define RBP_META_IDENTIFIER_NUM_DIHEDRAL_TYPES "number of dihedral types"
#define RBP_META_IDENTIFIER_CHARS_PER_ATOM "chars per atom"
#define RBP_META_IDENTIFIER_BOND_STYLE "bond style"
#define RBP_META_IDENTIFIER_ANGLE_STYLE "angle style"
#define RBP_META_IDENTIFIER_DIHEDRAL_STYLE "dihedral style"
#define RBP_META_IDENTIFIER_SUBTRACT_GROUNDSTATE "subtract groundstate"
#define RBP_META_IDENTIFIER_ODD "odd elasticity"
#define RBP_META_IDENTIFIER_SEQS_SET "seqs set"
#define RBP_META_IDENTIFIER_SEQS_CENTERED "seqs centered"
#define RBP_META_IDENTIFIER_CLOSED "closed"
#define RBP_META_IDENTIFIER_UNIT_LENGTH "unit length"
#define RBP_META_IDENTIFIER_UNIT_ENERGY "unit energy"
// section header identifiers
#define RBP_SECTION_HEADER_BOND_COEFFS "Bond Coeffs"
#define RBP_SECTION_HEADER_ANGLE_COEFFS "Angle Coeffs"
#define RBP_SECTION_HEADER_DIHEDRAL_COEFFS "Dihedral Coeffs"
#define RBP_SECTION_HEADER_UNKNOWN "Unknown section"
// clang-format on

namespace LAMMPS_NS {
class Error;

struct RBPCoeffs {
  std::vector<double> coeffs;
};

struct RBPMetadata {
  // counts
  int num_rigid_bodies;
  int num_bonds;
  int num_angles;
  int num_dihedrals;
  int num_bond_types;
  int num_angle_types;
  int num_dihedral_types;
  int coupling_range;
  int chars_per_atom;

  // interaction styles
  std::string bond_style;
  std::string angle_style;
  std::string dihedral_style;

  // flags
  bool subtract_groundstate;
  bool odd;              // file declares non-symmetric (odd elastic) blocks
  bool seqs_set;
  bool seqs_centered;
  bool closed;

  // unit conversion factors
  double unit_length;
  double unit_energy;

  RBPMetadata()
    : num_rigid_bodies(-1), num_bonds(-1), num_angles(-1), num_dihedrals(-1),
      num_bond_types(-1), num_angle_types(-1), num_dihedral_types(-1),
      coupling_range(-1), chars_per_atom(-1),
      bond_style(""), angle_style(""), dihedral_style(""),
      subtract_groundstate(RBP_META_DEFAULT_SUBTRACT_GROUNDSTATE),
      odd(RBP_META_DEFAULT_ODD),
      seqs_set(RBP_META_DEFAULT_SEQS_SET),
      seqs_centered(RBP_META_DEFAULT_SEQS_CENTERED),
      closed(RBP_META_DEFAULT_CLOSED),
      unit_length(RBP_META_DEFAULT_UNIT_LENGTH),
      unit_energy(RBP_META_DEFAULT_UNIT_ENERGY) {}
};

class RBPDatabase {
 public:
  RBPDatabase(LAMMPS *lmp, Error *error);

  void read(const std::string &filename);

  const RBPCoeffs &bond(int id) const;
  const RBPCoeffs &angle(int id) const;
  const RBPCoeffs &dihedral(int id) const;
  
  const RBPMetadata &metadata() const { return db_metadata; }

  int max_bond_id() const { return bond_db.size() - 1; }
  int max_angle_id() const { return angle_db.size() - 1; }
  int max_dihedral_id() const { return dihedral_db.size() - 1; }

 private:
  LAMMPS *lmp;
  Error *error;
  FILE *screen;

  std::vector<RBPCoeffs> bond_db;
  std::vector<RBPCoeffs> angle_db;
  std::vector<RBPCoeffs> dihedral_db;
  RBPMetadata db_metadata;

  void parse_file(const std::string &filename);
};

} 

#endif