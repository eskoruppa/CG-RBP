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


#include "parse_rbp.h"
#include "error.h"

#include "text_file_reader.h"
#include "tokenizer.h"  

using namespace LAMMPS_NS;

namespace {

enum Section { RBPDB_NONE, RBPDB_BOND, RBPDB_ANGLE, RBPDB_DIHEDRAL };

auto section_name = [](Section s) {
  switch (s) {
    case RBPDB_BOND:     return RBP_SECTION_HEADER_BOND_COEFFS;
    case RBPDB_ANGLE:    return RBP_SECTION_HEADER_ANGLE_COEFFS;
    case RBPDB_DIHEDRAL: return RBP_SECTION_HEADER_DIHEDRAL_COEFFS;
    default:             return RBP_SECTION_HEADER_UNKNOWN;
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

// Helper function to check if a line is a metadata line (contains ':')
inline bool is_metadata_line(const std::string &line)
{
  return line.find(':') != std::string::npos;
}

// Helper function to parse metadata line and extract key and value
inline bool parse_metadata(const std::string &line, std::string &key, std::string &value)
{
  size_t colon_pos = line.find(':');
  if (colon_pos == std::string::npos) return false;
  
  key = strip(line.substr(0, colon_pos));
  value = strip(line.substr(colon_pos + 1));
  return true;
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
  bool found_metadata = false;
  db_metadata = RBPMetadata();

  // Presence flags for fields that have no integer sentinel
  bool seen_subtract_groundstate = false;
  bool seen_seqs_set             = false;
  bool seen_seqs_centered        = false;
  bool seen_closed               = false;
  bool seen_unit_length          = false;
  bool seen_unit_energy          = false;

  // --------------------------------
  // First pass: determine max ids and parse all metadata
  {
    char *cline = nullptr;
    while ((cline = reader.next_line(0))) {
      std::string line = strip(std::string(cline));
      if (line.empty()) continue;

      // Metadata lines contain ':' and appear before any section header
      if (is_metadata_line(line)) {
        found_metadata = true;
        std::string key, value;
        if (!parse_metadata(line, key, value)) continue;

        // Helpers scoped to (key, value) for clean error reporting
        auto parse_int = [&](int &dest) {
          try { dest = std::stoi(value); }
          catch (...) { error->warning(FLERR, ("Failed to parse metadata field '" + key + "'").c_str()); }
        };
        auto parse_bool = [&](bool &dest) {
          try { dest = std::stoi(value) != 0; }
          catch (...) { error->warning(FLERR, ("Failed to parse metadata field '" + key + "'").c_str()); }
        };
        auto parse_double = [&](double &dest) {
          try { dest = std::stod(value); }
          catch (...) { error->warning(FLERR, ("Failed to parse metadata field '" + key + "'").c_str()); }
        };

        if      (key == RBP_META_IDENTIFIER_NUM_RIGID_BODIES)       parse_int(db_metadata.num_rigid_bodies);
        else if (key == RBP_META_IDENTIFIER_COUPLING_RANGE)         parse_int(db_metadata.coupling_range);
        else if (key == RBP_META_IDENTIFIER_NUM_BONDS)              parse_int(db_metadata.num_bonds);
        else if (key == RBP_META_IDENTIFIER_NUM_ANGLES)             parse_int(db_metadata.num_angles);
        else if (key == RBP_META_IDENTIFIER_NUM_DIHEDRALS)          parse_int(db_metadata.num_dihedrals);
        else if (key == RBP_META_IDENTIFIER_NUM_BOND_TYPES)         parse_int(db_metadata.num_bond_types);
        else if (key == RBP_META_IDENTIFIER_NUM_ANGLE_TYPES)        parse_int(db_metadata.num_angle_types);
        else if (key == RBP_META_IDENTIFIER_NUM_DIHEDRAL_TYPES)     parse_int(db_metadata.num_dihedral_types);
        else if (key == RBP_META_IDENTIFIER_CHARS_PER_ATOM)         parse_int(db_metadata.chars_per_atom);
        else if (key == RBP_META_IDENTIFIER_BOND_STYLE)             db_metadata.bond_style = value;
        else if (key == RBP_META_IDENTIFIER_ANGLE_STYLE)            db_metadata.angle_style = value;
        else if (key == RBP_META_IDENTIFIER_DIHEDRAL_STYLE)         db_metadata.dihedral_style = value;
        else if (key == RBP_META_IDENTIFIER_SUBTRACT_GROUNDSTATE)   { parse_bool(db_metadata.subtract_groundstate); seen_subtract_groundstate = true; }
        else if (key == RBP_META_IDENTIFIER_SEQS_SET)               { parse_bool(db_metadata.seqs_set);             seen_seqs_set = true; }
        else if (key == RBP_META_IDENTIFIER_SEQS_CENTERED)          { parse_bool(db_metadata.seqs_centered);        seen_seqs_centered = true; }
        else if (key == RBP_META_IDENTIFIER_CLOSED)                 { parse_bool(db_metadata.closed);               seen_closed = true; }
        else if (key == RBP_META_IDENTIFIER_UNIT_LENGTH)            { parse_double(db_metadata.unit_length);         seen_unit_length = true; }
        else if (key == RBP_META_IDENTIFIER_UNIT_ENERGY)            { parse_double(db_metadata.unit_energy);         seen_unit_energy = true; }
        // unknown keys are silently ignored to allow forward-compatible database files
        continue;
      }

      if (line == RBP_SECTION_HEADER_BOND_COEFFS)      { section = RBPDB_BOND; continue; }
      if (line == RBP_SECTION_HEADER_ANGLE_COEFFS)     { section = RBPDB_ANGLE; continue; }
      if (line == RBP_SECTION_HEADER_DIHEDRAL_COEFFS)  { section = RBPDB_DIHEDRAL; continue; }

      int id = 0;
      try {
        ValueTokenizer values(line);
        id = values.next_int();
      } catch (TokenizerException &e) {
        // Skip lines that don't start with an integer (e.g., sequence or connectivity data)
        continue;
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

    // Skip metadata lines
    if (is_metadata_line(line)) continue;

    if (line == RBP_SECTION_HEADER_BOND_COEFFS)      { section = RBPDB_BOND; continue; }
    if (line == RBP_SECTION_HEADER_ANGLE_COEFFS)     { section = RBPDB_ANGLE; continue; }
    if (line == RBP_SECTION_HEADER_DIHEDRAL_COEFFS)  { section = RBPDB_DIHEDRAL; continue; }

    if (section == RBPDB_NONE) {
      continue;
    }

    int id = 0;
    RBPCoeffs *target = nullptr;

    // first token: the id
    ValueTokenizer values(line);
    try {
      id = values.next_int();
    } catch (TokenizerException &e) {
      // Skip lines that don't start with an integer (e.g., sequence or connectivity data)
      continue;
    }

    if (section == RBPDB_BOND)      { target = &bond_db[id];     seen_bond[id] = 1; }
    if (section == RBPDB_ANGLE)     { target = &angle_db[id];    seen_angle[id] = 1; }
    if (section == RBPDB_DIHEDRAL)  { target = &dihedral_db[id]; seen_dihedral[id] = 1; }

    // remaining tokens: arbitrary number of doubles
    while (values.has_next()) {
      std::string tok = values.next_string();
      try {
        double val = std::stod(tok);
        target->coeffs.push_back(val);
      } catch (...) {
        std::string msg =
          "Non-numeric coefficient '" + tok + "' in rbp database-file:\n " + filename +
          "\nSection: " + std::string(section_name(section)) +
          "\nLine: " + line + "\n";
        error->all(FLERR, msg.c_str());
      }
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
  
  // ============================================
  // VALIDATE METADATA
  // Each check is guarded by a preprocessor define in parse_rbp.h.
  // Comment out the respective define to make a field optional.
  // ============================================

  #ifdef RBP_META_REQUIRE_METADATA_SECTION
    if (!found_metadata)
      error->all(FLERR, "Metadata section is mandatory in rbp database-file:\n " + filename +
                "\nExpected metadata lines with format: key: value");
  #endif

    // Integer count fields
  #ifdef RBP_META_REQUIRE_NUM_RIGID_BODIES
    if (db_metadata.num_rigid_bodies < 0)
      error->all(FLERR, "Missing required metadata 'number of rigid bodies' in rbp database-file:\n " + filename);
  #endif
  #ifdef RBP_META_REQUIRE_COUPLING_RANGE
    if (db_metadata.coupling_range < 0)
      error->all(FLERR, "Missing required metadata 'coupling range' in rbp database-file:\n " + filename);
  #endif
  #ifdef RBP_META_REQUIRE_NUM_BONDS
    if (db_metadata.num_bonds < 0)
      error->all(FLERR, "Missing required metadata 'number of bonds' in rbp database-file:\n " + filename);
  #endif
  #ifdef RBP_META_REQUIRE_NUM_ANGLES
    if (db_metadata.num_angles < 0)
      error->all(FLERR, "Missing required metadata 'number of angles' in rbp database-file:\n " + filename);
  #endif
  #ifdef RBP_META_REQUIRE_NUM_DIHEDRALS
    if (db_metadata.num_dihedrals < 0)
      error->all(FLERR, "Missing required metadata 'number of dihedrals' in rbp database-file:\n " + filename);
  #endif
  #ifdef RBP_META_REQUIRE_NUM_BOND_TYPES
    if (db_metadata.num_bond_types < 0)
      error->all(FLERR, "Missing required metadata 'number of bond types' in rbp database-file:\n " + filename);
  #endif
  #ifdef RBP_META_REQUIRE_NUM_ANGLE_TYPES
    if (db_metadata.num_angle_types < 0)
      error->all(FLERR, "Missing required metadata 'number of angle types' in rbp database-file:\n " + filename);
  #endif
  #ifdef RBP_META_REQUIRE_NUM_DIHEDRAL_TYPES
    if (db_metadata.num_dihedral_types < 0)
      error->all(FLERR, "Missing required metadata 'number of dihedral types' in rbp database-file:\n " + filename);
  #endif
  #ifdef RBP_META_REQUIRE_CHARS_PER_ATOM
    if (db_metadata.chars_per_atom < 0)
      error->all(FLERR, "Missing required metadata 'chars per atom' in rbp database-file:\n " + filename);
  #endif

    // Style string fields
  #ifdef RBP_META_REQUIRE_BOND_STYLE
    if (db_metadata.bond_style.empty())
      error->all(FLERR, "Missing required metadata 'bond style' in rbp database-file:\n " + filename);
  #endif
  #ifdef RBP_META_REQUIRE_ANGLE_STYLE
    if (db_metadata.angle_style.empty())
      error->all(FLERR, "Missing required metadata 'angle style' in rbp database-file:\n " + filename);
  #endif
  #ifdef RBP_META_REQUIRE_DIHEDRAL_STYLE
    if (db_metadata.dihedral_style.empty())
      error->all(FLERR, "Missing required metadata 'dihedral style' in rbp database-file:\n " + filename);
  #endif

    // Boolean flag fields
  #ifdef RBP_META_REQUIRE_SUBTRACT_GROUNDSTATE
    if (!seen_subtract_groundstate)
      error->all(FLERR, "Missing required metadata 'subtract groundstate' in rbp database-file:\n " + filename);
  #else
    if (!seen_subtract_groundstate) db_metadata.subtract_groundstate = RBP_META_DEFAULT_SUBTRACT_GROUNDSTATE;
  #endif
  #ifdef RBP_META_REQUIRE_SEQS_SET
    if (!seen_seqs_set)
      error->all(FLERR, "Missing required metadata 'seqs set' in rbp database-file:\n " + filename);
  #else
    if (!seen_seqs_set) db_metadata.seqs_set = RBP_META_DEFAULT_SEQS_SET;
  #endif
  #ifdef RBP_META_REQUIRE_SEQS_CENTERED
    if (!seen_seqs_centered)
      error->all(FLERR, "Missing required metadata 'seqs centered' in rbp database-file:\n " + filename);
  #else
    if (!seen_seqs_centered) db_metadata.seqs_centered = RBP_META_DEFAULT_SEQS_CENTERED;
  #endif
  #ifdef RBP_META_REQUIRE_CLOSED
    if (!seen_closed)
      error->all(FLERR, "Missing required metadata 'closed' in rbp database-file:\n " + filename);
  #else
    if (!seen_closed) db_metadata.closed = RBP_META_DEFAULT_CLOSED;
  #endif

    // Double fields
  #ifdef RBP_META_REQUIRE_UNIT_LENGTH
    if (!seen_unit_length)
      error->all(FLERR, "Missing required metadata 'unit length' in rbp database-file:\n " + filename);
  #else
    if (!seen_unit_length) db_metadata.unit_length = RBP_META_DEFAULT_UNIT_LENGTH;
  #endif
  #ifdef RBP_META_REQUIRE_UNIT_ENERGY
    if (!seen_unit_energy)
      error->all(FLERR, "Missing required metadata 'unit energy' in rbp database-file:\n " + filename);
  #else
    if (!seen_unit_energy) db_metadata.unit_energy = RBP_META_DEFAULT_UNIT_ENERGY;
  #endif

  // Validate type counts against actual data (always run when the field was present)
  if (db_metadata.num_bond_types >= 0 && db_metadata.num_bond_types != max_bond) {
    std::string msg =
      "Inconsistent bond type count in rbp database-file:\n " + filename +
      "\nMetadata specifies " + std::to_string(db_metadata.num_bond_types) +
      " bond types, but found " + std::to_string(max_bond) + " in Bond Coeffs section.";
    error->all(FLERR, msg.c_str());
  }
  if (db_metadata.num_angle_types >= 0 && db_metadata.num_angle_types != max_angle) {
    std::string msg =
      "Inconsistent angle type count in rbp database-file:\n " + filename +
      "\nMetadata specifies " + std::to_string(db_metadata.num_angle_types) +
      " angle types, but found " + std::to_string(max_angle) + " in Angle Coeffs section.";
    error->all(FLERR, msg.c_str());
  }
  if (db_metadata.num_dihedral_types >= 0 && db_metadata.num_dihedral_types != max_dihedral) {
    std::string msg =
      "Inconsistent dihedral type count in rbp database-file:\n " + filename +
      "\nMetadata specifies " + std::to_string(db_metadata.num_dihedral_types) +
      " dihedral types, but found " + std::to_string(max_dihedral) + " in Dihedral Coeffs section.";
    error->all(FLERR, msg.c_str());
  }
}