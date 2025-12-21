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


#include "parse_rbp.h"
#include "error.h"

#include "text_file_reader.h"
#include "tokenizer.h"  

using namespace LAMMPS_NS;

namespace {

enum Section { RBPDB_NONE, RBPDB_BOND, RBPDB_ANGLE, RBPDB_DIHEDRAL };

auto section_name = [](Section s) {
  switch (s) {
    case RBPDB_BOND:     return "Bond Coeffs";
    case RBPDB_ANGLE:    return "Angle Coeffs";
    case RBPDB_DIHEDRAL: return "Dihedral Coeffs";
    default:             return "Unknown section";
  }
};

inline std::string strip(const std::string &s)
{
  size_t start = 0;
  while (start < s.size() && std::isspace(static_cast<unsigned char>(s[start])))
    start++;

  size_t end = s.size();
  while (end > start && std::isspace(static_cast<unsigned char>(s[end-1])))
    end--;

  return s.substr(start, end - start);
}

} 

RBPDatabase::RBPDatabase(LAMMPS *lmp_in, Error *err)
  : lmp(lmp_in), error(err), screen(lmp_in->screen) {}

void RBPDatabase::read(const std::string &filename)
{
  parse_file(filename);
}

const RBPCoeffs &RBPDatabase::bond(int id) const
{
  if (id <= 0 || id >= (int)bond_db.size())
    error->all(FLERR, "Invalid bond database id");
  return bond_db[id];
}

const RBPCoeffs &RBPDatabase::angle(int id) const
{
  if (id <= 0 || id >= (int)angle_db.size())
    error->all(FLERR, "Invalid angle database id");
  return angle_db[id];
}

const RBPCoeffs &RBPDatabase::dihedral(int id) const
{
  if (id <= 0 || id >= (int)dihedral_db.size())
    error->all(FLERR, "Invalid dihedral database id");
  return dihedral_db[id];
}

void RBPDatabase::parse_file(const std::string &filename)
{

  TextFileReader reader(filename, "rbp database-file");
  reader.set_bufsize(2048);

  Section section = RBPDB_NONE;
  int max_bond = 0, max_angle = 0, max_dihedral = 0;

  // --------------------------------
  // First pass to determine max ids
  {
    char *cline = nullptr;
    while ((cline = reader.next_line(0))) {
      std::string line = strip(std::string(cline));
      if (line.empty()) continue;

      if (line == "Bond Coeffs")      { section = RBPDB_BOND; continue; }
      if (line == "Angle Coeffs")     { section = RBPDB_ANGLE; continue; }
      if (line == "Dihedral Coeffs")  { section = RBPDB_DIHEDRAL; continue; }

      int id = 0;
      try {
        ValueTokenizer values(line);
        id = values.next_int();
      } catch (TokenizerException &e) {
        std::string msg =
          "Malformed database entry in rbp database-file:\n " + filename +
          "\nSection: " + std::string(section_name(section)) +
          "\nLine expected to start with integer. Encountered: \n" + line + "\n";
        error->all(FLERR, msg.c_str());
      }

      if (id < 1) {
        std::string msg =
          "Database IDs must start at 1 in rbp database-file:\n " + filename +
          "\nSection: " + std::string(section_name(section)) +
          "\nEncountered id " + std::to_string(id) +
          " in line:\n" + line + "\n";
        error->all(FLERR, msg.c_str());
      }

      if (     section == RBPDB_BOND)     max_bond     = std::max(max_bond, id);
      else if (section == RBPDB_ANGLE)    max_angle    = std::max(max_angle, id);
      else if (section == RBPDB_DIHEDRAL) max_dihedral = std::max(max_dihedral, id);
    }
  }

  bond_db.assign(max_bond + 1, RBPCoeffs());
  angle_db.assign(max_angle + 1, RBPCoeffs());
  dihedral_db.assign(max_dihedral + 1, RBPCoeffs());

  // marker arrays to verify that all indices between 1 and max_x has been read
  std::vector<char> seen_bond(max_bond + 1, 0);
  std::vector<char> seen_angle(max_angle + 1, 0);
  std::vector<char> seen_dihedral(max_dihedral + 1, 0);

  // --------------------------------
  // Second pass to read coefficients
  reader.rewind();
  section = RBPDB_NONE;

  char *cline = nullptr;
  while ((cline = reader.next_line(0))) {
    std::string line = strip(std::string(cline));
    if (line.empty()) continue;

    if (line == "Bond Coeffs")      { section = RBPDB_BOND; continue; }
    if (line == "Angle Coeffs")     { section = RBPDB_ANGLE; continue; }
    if (line == "Dihedral Coeffs")  { section = RBPDB_DIHEDRAL; continue; }

    if (section == RBPDB_NONE) {
      std::string msg =
        "Coefficient line encountered in rbp database-file \"" + filename +
        "\" before any section header "
        "(Bond Coeffs / Angle Coeffs / Dihedral Coeffs)";
      error->all(FLERR, msg.c_str());
    }

    int expected = (section == RBPDB_BOND) ? 27 :
                   (section == RBPDB_ANGLE || section == RBPDB_DIHEDRAL) ? 48 : 0;

    int id = 0;
    RBPCoeffs *target = nullptr;

    try {
      ValueTokenizer values(line);
      id = values.next_int();

      if (section == RBPDB_BOND)      { target = &bond_db[id];     seen_bond[id] = 1; }
      if (section == RBPDB_ANGLE)     { target = &angle_db[id];    seen_angle[id] = 1; }
      if (section == RBPDB_DIHEDRAL)  { target = &dihedral_db[id]; seen_dihedral[id] = 1; }

      target->coeffs.resize(expected);
      for (int i = 0; i < expected; i++) {
        target->coeffs[i] = values.next_double();
      }

    } catch (TokenizerException &e) {
      error->all(FLERR, "Incorrect number of coefficients in database entry");
    }
  }

  // --------------------------------------------
  // Check contiguity of provided indices
  for (int i = 1; i <= max_bond; i++)
    if (!seen_bond[i])
      error->all(FLERR, "Missing bond database ID");

  for (int i = 1; i <= max_angle; i++)
    if (!seen_angle[i])
      error->all(FLERR, "Missing angle database ID");

  for (int i = 1; i <= max_dihedral; i++)
    if (!seen_dihedral[i])
      error->all(FLERR, "Missing dihedral database ID");
}