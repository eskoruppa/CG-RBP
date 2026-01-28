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


#ifndef LAMMPS_RBP_DATABASE_H
#define LAMMPS_RBP_DATABASE_H

#include <vector>
#include <string>

#include "lammps.h"

namespace LAMMPS_NS {
class Error;

struct RBPCoeffs {
  std::vector<double> coeffs;
};

struct RBPMetadata {
  int num_bond_types;
  int num_angle_types;
  int num_dihedral_types;
  int coupling_range;
  std::string bond_style;      // e.g., "rbp" or "rbpfene"
  std::string angle_style;     // e.g., "rbp"
  std::string dihedral_style;  // e.g., "rbp"
  
  RBPMetadata() 
    : num_bond_types(-1), num_angle_types(-1), num_dihedral_types(-1),
      coupling_range(-1), bond_style(""), angle_style(""), dihedral_style("") {}
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
