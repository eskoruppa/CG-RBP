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
  
  // Variables to store metadata values (required)
  int metadata_num_bond_types = -1;
  int metadata_num_angle_types = -1;
  int metadata_num_dihedral_types = -1;
  int metadata_coupling_range = -1;
  std::string metadata_bond_style = "";
  std::string metadata_angle_style = "";
  std::string metadata_dihedral_style = "";
  
  bool found_metadata = false;

  // --------------------------------
  // First pass to determine max ids and extract metadata
  {
    char *cline = nullptr;
    while ((cline = reader.next_line(0))) {
      std::string line = strip(std::string(cline));
      if (line.empty()) continue;

      // Check for metadata lines
      if (is_metadata_line(line)) {
        found_metadata = true;
        std::string key, value;
        if (parse_metadata(line, key, value)) {
          // Extract required metadata fields
          if (key == "number of bond types") {
            try {
              metadata_num_bond_types = std::stoi(value);
            } catch (...) {
              error->warning(FLERR, "Failed to parse 'number of bond types' metadata");
            }
          }
          else if (key == "number of angle types") {
            try {
              metadata_num_angle_types = std::stoi(value);
            } catch (...) {
              error->warning(FLERR, "Failed to parse 'number of angle types' metadata");
            }
          }
          else if (key == "number of dihedral types") {
            try {
              metadata_num_dihedral_types = std::stoi(value);
            } catch (...) {
              error->warning(FLERR, "Failed to parse 'number of dihedral types' metadata");
            }
          }
          else if (key == "coupling range") {
            try {
              metadata_coupling_range = std::stoi(value);
            } catch (...) {
              error->warning(FLERR, "Failed to parse 'coupling range' metadata");
            }
          }
          else if (key == "bond style") {
            metadata_bond_style = value;
          }
          else if (key == "angle style") {
            metadata_angle_style = value;
          }
          else if (key == "dihedral style") {
            metadata_dihedral_style = value;
          }
          // Skip other metadata fields (seqs set, seqs centered, chars per atom)
        }
        // Skip all metadata lines
        continue;
      }

      if (line == "Bond Coeffs")      { section = RBPDB_BOND; continue; }
      if (line == "Angle Coeffs")     { section = RBPDB_ANGLE; continue; }
      if (line == "Dihedral Coeffs")  { section = RBPDB_DIHEDRAL; continue; }

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

    if (line == "Bond Coeffs")      { section = RBPDB_BOND; continue; }
    if (line == "Angle Coeffs")     { section = RBPDB_ANGLE; continue; }
    if (line == "Dihedral Coeffs")  { section = RBPDB_DIHEDRAL; continue; }

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
  // VALIDATE AND STORE METADATA (mandatory)
  // ============================================
  
  if (!found_metadata) {
    error->all(FLERR, "Metadata section is mandatory in rbp database-file:\n " + filename +
               "\nExpected metadata lines with format: key: value");
  }
  
  // Validate required metadata fields
  if (metadata_num_bond_types < 0) {
    error->all(FLERR, "Missing required metadata 'number of bond types' in rbp database-file:\n " + filename);
  }
  if (metadata_num_angle_types < 0) {
    error->all(FLERR, "Missing required metadata 'number of angle types' in rbp database-file:\n " + filename);
  }
  if (metadata_num_dihedral_types < 0) {
    error->all(FLERR, "Missing required metadata 'number of dihedral types' in rbp database-file:\n " + filename);
  }
  if (metadata_bond_style.empty()) {
    error->all(FLERR, "Missing required metadata 'bond style' in rbp database-file:\n " + filename);
  }
  if (metadata_angle_style.empty()) {
    error->all(FLERR, "Missing required metadata 'angle style' in rbp database-file:\n " + filename);
  }
  if (metadata_dihedral_style.empty()) {
    error->all(FLERR, "Missing required metadata 'dihedral style' in rbp database-file:\n " + filename);
  }
  
  // Validate type counts match actual data
  if (metadata_num_bond_types != max_bond) {
    std::string msg =
      "Inconsistent bond type count in rbp database-file:\n " + filename +
      "\nMetadata specifies " + std::to_string(metadata_num_bond_types) +
      " bond types, but found " + std::to_string(max_bond) + " in Bond Coeffs section.";
    error->all(FLERR, msg.c_str());
  }
  
  if (metadata_num_angle_types != max_angle) {
    std::string msg =
      "Inconsistent angle type count in rbp database-file:\n " + filename +
      "\nMetadata specifies " + std::to_string(metadata_num_angle_types) +
      " angle types, but found " + std::to_string(max_angle) + " in Angle Coeffs section.";
    error->all(FLERR, msg.c_str());
  }
  
  if (metadata_num_dihedral_types != max_dihedral) {
    std::string msg =
      "Inconsistent dihedral type count in rbp database-file:\n " + filename +
      "\nMetadata specifies " + std::to_string(metadata_num_dihedral_types) +
      " dihedral types, but found " + std::to_string(max_dihedral) + " in Dihedral Coeffs section.";
    error->all(FLERR, msg.c_str());
  }
  
  // Store metadata in the database
  db_metadata.num_bond_types = metadata_num_bond_types;
  db_metadata.num_angle_types = metadata_num_angle_types;
  db_metadata.num_dihedral_types = metadata_num_dihedral_types;
  db_metadata.coupling_range = metadata_coupling_range;
  db_metadata.bond_style = metadata_bond_style;
  db_metadata.angle_style = metadata_angle_style;
  db_metadata.dihedral_style = metadata_dihedral_style;
}