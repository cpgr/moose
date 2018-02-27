//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "GeochemicalDatabaseReader.h"

#include "MooseUtils.h"
#include "Conversion.h"

#include <fstream>

GeochemicalDatabaseReader::GeochemicalDatabaseReader(const std::string filename)
  : _filename(filename), _read(false)
{
}

void
GeochemicalDatabaseReader::read()
{
  MooseUtils::checkFileReadable(_filename);

  std::ifstream file_data(_filename.c_str());
  std::string line;
  std::vector<std::string> tokens;

  // Read the database
  while (std::getline(file_data, line))
  {
    // Skip lines starting with comment character #
    if (line.find('#') == 0)
      continue;

    // Skip empty lines
    if (line.empty())
      continue;

    // Skip lines containing whitespace only
    if (MooseUtils::trim(line).length() == 0)
      continue;

    // remove single quotes from strings in database
    line.erase(std::remove(line.begin(), line.end(), '\''), line.end());

    // Tokenize each line
    MooseUtils::tokenize(line, tokens, 1, " ");

    if (tokens[0] == "temperatures")
    {
      // Read the temperature data
      _temperature_points.resize(tokens.size() - 1);
      for (auto i = beginIndex(_temperature_points); i < _temperature_points.size(); ++i)
        _temperature_points[i] = std::stod(tokens[i + 1]);
    }

    else if (tokens[0] == "primary")
    {
      _section_enum = sectionEnum::PRIMARY;
      continue;
    }

    else if (tokens[0] == "secondary")
    {
      _section_enum = sectionEnum::SECONDARY;
      continue;
    }

    else if (tokens[0] == "minerals")
    {
      _section_enum = sectionEnum::MINERAL;
      continue;
    }
    else
    {
    }

    // Read the actual data
    switch (_section_enum)
    {
      case (sectionEnum::PRIMARY):
        // Read primary species if it is specified in the given list of species
        if (std::find(_ps_names.begin(), _ps_names.end(), tokens[0]) != _ps_names.end())
          readPrimarySpecies(tokens);

        break;

      case (sectionEnum::SECONDARY):
        // Read secondary equilibrium species if it is specified in the given list
        if (std::find(_es_names.begin(), _es_names.end(), tokens[0]) != _es_names.end())
          readEquilibriumSpecies(tokens);

        break;

      case (sectionEnum::MINERAL):
        // Read secondary mineral species if it is specified in the given list
        if (std::find(_ms_names.begin(), _ms_names.end(), tokens[0]) != _ms_names.end())
          readMineralSpecies(tokens);

        break;
    }
  }

  // Check to make sure that all primary species in the input file have
  // been read from the database
  std::vector<std::string> ps_db(_primary_species.size());
  for (auto i = beginIndex(ps_db); i < ps_db.size(); ++i)
    ps_db[i] = _primary_species[i].name;

  for (auto i = beginIndex(_ps_names); i < _ps_names.size(); ++i)
    if (std::find(ps_db.begin(), ps_db.end(), _ps_names[i]) == ps_db.end())
      mooseError("GeochemicalDatabaseReader: primary species ",
                 _ps_names[i],
                 " not found in database ",
                 _filename);

  // Check to make sure that all secondary equilibrium species in the input file have
  // been read from the database
  std::vector<std::string> es_db(_equilibrium_species.size());
  for (auto i = beginIndex(es_db); i < es_db.size(); ++i)
    es_db[i] = _equilibrium_species[i].name;

  for (auto i = beginIndex(_es_names); i < _es_names.size(); ++i)
    if (std::find(es_db.begin(), es_db.end(), _es_names[i]) == es_db.end())
      mooseError("GeochemicalDatabaseReader: secondary equilibrium species ",
                 _es_names[i],
                 " not found in database ",
                 _filename);

  // Check to make sure that all secondary mineral species in the input file have
  // been read from the database
  std::vector<std::string> ms_db(_mineral_species.size());
  for (auto i = beginIndex(ms_db); i < ms_db.size(); ++i)
    ms_db[i] = _mineral_species[i].name;

  for (auto i = beginIndex(_ms_names); i < _ms_names.size(); ++i)
    if (std::find(ms_db.begin(), ms_db.end(), _ms_names[i]) == ms_db.end())
      mooseError("GeochemicalDatabaseReader: secondary mineral species ",
                 _ms_names[i],
                 " not found in database ",
                 _filename);

  // Successfully read the database
  _read = true;
}

void
GeochemicalDatabaseReader::readPrimarySpecies(std::vector<std::string> & tokens)
{
  // Make sure that there are exactly 4 tokens
  if (tokens.size() != 4)
    mooseError(
        "GeochemicalDatabaseReader: incorrect number of primary species values encountered in ",
        _filename,
        " for species ",
        tokens[0]);

  // Save all values into a primary species struct
  GeochemicalDatabasePrimarySpecies pspecies;

  pspecies.name = tokens[0];
  pspecies.radius = std::stod(tokens[1]);
  pspecies.charge = std::stod(tokens[2]);
  pspecies.molar_mass = std::stod(tokens[3]);

  _primary_species.push_back(pspecies);
}

void
GeochemicalDatabaseReader::readEquilibriumSpecies(std::vector<std::string> & tokens)
{
  GeochemicalDatabaseEquilibriumSpecies sspecies;

  // The first four tokens are name, radius a0, charge and molar mass
  sspecies.name = tokens[0];
  sspecies.debye_a = std::stod(tokens[1]);
  sspecies.charge = std::stod(tokens[2]);
  sspecies.molar_mass = std::stod(tokens[3]);

  // The next element of tokens is the number of primary species involved
  // this secondary species equilibrium reaction
  unsigned int nprimary = std::stoi(tokens[4]);
  sspecies.nspecies = nprimary;

  // At this point, we can check the number of tokens is correct. There should be
  // values for name, the Debye a parameter, charge and molar mass, nprimary, 2 * nprimary
  // (corresponding to primary species), _temperature_points.size() equilibrium constant values
  const unsigned int num_tokens = 5 + 2 * nprimary + _temperature_points.size();
  if (tokens.size() != num_tokens)
    mooseError(
        "GeochemicalDatabaseReader::incorrect number of secondary species values encountered in ",
        _filename,
        " for species ",
        tokens[0]);

  // Read in the pairs of stoichiometric coefficients and primary variables for
  // the nprimary species involved
  std::vector<Real> stoichiometric_coeff;
  std::vector<std::string> primary_species;

  // The primary species begin at position 5
  const unsigned int element_pos = 5;

  for (unsigned int i = 0; i < nprimary; ++i)
  {
    stoichiometric_coeff.push_back(std::stod(tokens[element_pos + 2 * i]));
    primary_species.push_back(tokens[element_pos + 2 * i + 1]);
  }

  sspecies.stoichiometric_coeff = stoichiometric_coeff;
  sspecies.primary_species = primary_species;

  // Now read the array of equilibrium values (starting from position 5 + 2 * nprimary)
  const unsigned int eq_pos = 5 + 2 * nprimary;

  std::vector<Real> eq_const;
  for (unsigned int i = eq_pos; i < eq_pos + _temperature_points.size(); ++i)
    eq_const.push_back(std::stod(tokens[i]));

  sspecies.equilibrium_const = eq_const;

  _equilibrium_species.push_back(sspecies);
}

void
GeochemicalDatabaseReader::readMineralSpecies(std::vector<std::string> & tokens)
{
  GeochemicalDatabaseMineralSpecies mspecies;

  mspecies.name = tokens[0];
  mspecies.molar_mass = std::stod(tokens[1]);
  mspecies.molar_volume = std::stod(tokens[2]);

  // The next element of tokens is the number of primary species involved
  // this secondary species equilibrium reaction
  unsigned int nprimary = std::stoi(tokens[3]);
  mspecies.nspecies = nprimary;

  // At this point, we can check the number of tokens is correct. There should be
  // values for name, molar mass, molar volume, nprimary, 2 * nprimary
  // (corresponding to primary species), _temperature_points.size() equilibrium constant values
  const unsigned int num_tokens = 4 + 2 * nprimary + _temperature_points.size();
  if (tokens.size() != num_tokens)
    mooseError(
        "GeochemicalDatabaseReader::incorrect number of mineral species values encountered in ",
        _filename,
        " for species ",
        tokens[0]);

  // Read in the pairs of stoichiometric coefficients and primary variables for
  // the nprimary species involved
  std::vector<Real> stoichiometric_coeff;
  std::vector<std::string> primary_species;

  // The primary species begin at position 5
  const unsigned int element_pos = 4;

  for (unsigned int i = 0; i < nprimary; ++i)
  {
    stoichiometric_coeff.push_back(std::stod(tokens[element_pos + 2 * i]));
    primary_species.push_back(tokens[element_pos + 2 * i + 1]);
  }

  mspecies.stoichiometric_coeff = stoichiometric_coeff;
  mspecies.primary_species = primary_species;

  // Now read the array of equilibrium values (starting from position 4 + 2 * nprimary)
  const unsigned int eq_pos = 4 + 2 * nprimary;

  std::vector<Real> eq_const;
  for (unsigned int i = eq_pos; i < eq_pos + _temperature_points.size(); ++i)
    eq_const.push_back(std::stod(tokens[i]));

  mspecies.equilibrium_const = eq_const;

  _mineral_species.push_back(mspecies);
}

void
GeochemicalDatabaseReader::getPrimarySpecies(
    std::vector<GeochemicalDatabasePrimarySpecies> & primary_species) const
{
  if (!_read)
    mooseError(
        "GeochemicalDatabaseReader::read() must be called before calling getPrimarySpecies()");

  primary_species = _primary_species;
}

void
GeochemicalDatabaseReader::getEquilibriumSpecies(
    std::vector<GeochemicalDatabaseEquilibriumSpecies> & equilibrium_species) const
{
  if (!_read)
    mooseError(
        "GeochemicalDatabaseReader::read() must be called before calling getEquilibriumSpecies()");

  equilibrium_species = _equilibrium_species;
}

void
GeochemicalDatabaseReader::getMineralSpecies(
    std::vector<GeochemicalDatabaseMineralSpecies> & mineral_species) const
{
  if (!_read)
    mooseError(
        "GeochemicalDatabaseReader::read() must be called before calling getMineralSpecies()");

  mineral_species = _mineral_species;
}

void
GeochemicalDatabaseReader::equilibriumReactions(std::vector<std::string> & reactions)
{
  if (!_read)
    mooseError(
        "GeochemicalDatabaseReader::read() must be called before calling equilibriumReactions()");

  reactions.resize(_equilibrium_species.size());

  // Format the reactions for pretty printing
  for (auto i = beginIndex(_equilibrium_species); i < _equilibrium_species.size(); ++i)
  {
    std::string eq_reaction;

    // Add the primary species and stoichiometric coefficients for each reaction
    for (auto j = beginIndex(_equilibrium_species[i].primary_species);
         j < _equilibrium_species[i].primary_species.size();
         ++j)
    {
      Real sto = _equilibrium_species[i].stoichiometric_coeff[j];

      std::string sign = "+";
      if (sto < 0.0)
        sign = "-";

      if (j == 0 && sign == "-")
        eq_reaction += sign;

      if (j > 0)
        eq_reaction += sign;

      if (std::fabs(sto) != 1.0)
        eq_reaction += " " + Moose::stringify(std::fabs(sto));

      eq_reaction += " " + _equilibrium_species[i].primary_species[j] + " ";
    }

    // Add the equilibrium species name
    eq_reaction += "= " + _equilibrium_species[i].name;

    reactions[i] = eq_reaction;
  }
}

void
GeochemicalDatabaseReader::mineralReactions(std::vector<std::string> & reactions)
{
  if (!_read)
    mooseError(
        "GeochemicalDatabaseReader::read() must be called before calling mineralReactions()");

  reactions.resize(_mineral_species.size());

  // Format the reactions for pretty printing
  for (auto i = beginIndex(_mineral_species); i < _mineral_species.size(); ++i)
  {
    std::string mineral_reaction;

    // Add the primary species and stoichiometric coefficients for each reaction
    for (auto j = beginIndex(_mineral_species[i].primary_species);
         j < _mineral_species[i].primary_species.size();
         ++j)
    {
      Real sto = _mineral_species[i].stoichiometric_coeff[j];

      std::string sign = "+";
      if (sto < 0.0)
        sign = "-";

      if (j == 0 && sign == "-")
        mineral_reaction += sign;

      if (j > 0)
        mineral_reaction += sign;

      if (std::fabs(sto) != 1.0)
        mineral_reaction += " " + Moose::stringify(std::fabs(sto));

      mineral_reaction += " " + _mineral_species[i].primary_species[j] + " ";
    }

    // Add the mineral species name
    mineral_reaction += "= " + _mineral_species[i].name;

    reactions[i] = mineral_reaction;
  }
}
