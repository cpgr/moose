
//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef GEOCHEMICALDATABASEREADER_H
#define GEOCHEMICALDATABASEREADER_H

#include "MooseTypes.h"

/**
 * Data structure for primary species. Members are:
 * Species name
 * ionic radius (Angstrom)
 * charge
 * molar mass (g/mol)
 */
struct GeochemicalDatabasePrimarySpecies
{
  std::string name;
  Real radius;
  Real charge;
  Real molar_mass;
};

/**
 * Data structure for secondary equilibrium species. Members are:
 * Species name
 * Number of constituent species in the reaction
 * Names of primary species in the reaction
 * Stoichiometric coefficients of primarys species in the reaction
 * Equilibrium constant at given temperatures
 * Debye-Huckel parameter a
 * charge
 * molar mass (g/mol)
 */
struct GeochemicalDatabaseEquilibriumSpecies
{
  std::string name;
  unsigned int nspecies;
  std::vector<std::string> primary_species;
  std::vector<Real> stoichiometric_coeff;
  std::vector<Real> equilibrium_const;
  Real debye_a;
  Real charge;
  Real molar_mass;
};

/**
 * Data structure for mineral species. Members are:
 * Species name
 * Molar volume (mol/l)
 * Number of constituent species in the reaction
 * Names of primary species in the reaction
 * Stoichiometric coefficients of primarys species in the reaction
 * Equilibrium constant at given temperatures
 * molar mass (g/mol)
 */
struct GeochemicalDatabaseMineralSpecies
{
  std::string name;
  Real molar_volume;
  unsigned int nspecies;
  std::vector<std::string> primary_species;
  std::vector<Real> stoichiometric_coeff;
  std::vector<Real> equilibrium_const;
  Real molar_mass;
};

/**
 * Class for reading geochemical reactions from a geochemical database
 */
class GeochemicalDatabaseReader
{
public:
  GeochemicalDatabaseReader(const std::string filename);

  /**
   * Read the thermodynamic database
   */
  void read();

  /**
   * Set the list of primary species to read from database
   * @param names list of primary species names
   */
  void setPrimarySpeciesNames(std::vector<std::string> names) { _ps_names = names; };

  /**
   * Set the list of secondary equilibrium species to read from database
   * @param names list of equilibrium species names
   */
  void setEquilibriumSpeciesNames(std::vector<std::string> names) { _es_names = names; };

  /**
   * Set the list of secondary mineral species to read from database
   * @param names list of mineral species names
   */
  void setMineralSpeciesNames(std::vector<std::string> names) { _ms_names = names; };

  /**
   * Get the temperature points that the equilibrium constant is defined at.
   * Note: must be called after read()
   * @param[out] vector of temperature points (C)
   */
  void getTemperatures(std::vector<Real> & temperature) { temperature = _temperature_points; };

  /**
   * Get the primary species structure
   * Note: must be called after read()
   * @param[out] primary species structure
   */
  void getPrimarySpecies(std::vector<GeochemicalDatabasePrimarySpecies> & primary_species) const;

  /**
   * Get the secondary equilibrium species structure
   * Note: must be called after read()
   * @param[out] secondary species structure
   */
  void getEquilibriumSpecies(
      std::vector<GeochemicalDatabaseEquilibriumSpecies> & secondary_species) const;

  /**
   * Get the mineral species structure
   * Note: must be called after read()
   * @param[out] secondary species structure
   */
  void getMineralSpecies(std::vector<GeochemicalDatabaseMineralSpecies> & mineral_species) const;

  /**
   * Generates a formatted vector of strings representing all aqueous equilibrium
   * reactions
   * Note: must be called after read()
   * @param[out] reactions formatted equilibrium reactions
   */
  void equilibriumReactions(std::vector<std::string> & reactions);

  /**
   * Generates a formatted vector of strings representing all mineral reactions
   * Note: must be called after read()
   * @param[out] reactions formatted mineral reactions
   */
  void mineralReactions(std::vector<std::string> & reactions);

protected:
  /**
   * Reads primary species data from database
   * @param tokens tokenized line read from database
   */
  void readPrimarySpecies(std::vector<std::string> & tokens);

  /**
   * Reads secondary equilibrium species data from database
   * @param tokens tokenized line read from database
   */
  void readEquilibriumSpecies(std::vector<std::string> & tokens);

  /**
   * Reads mineral specise data from database
   * @param tokens tokenized line read from database
   */
  void readMineralSpecies(std::vector<std::string> & tokens);

  /// Database filename
  const FileName _filename;
  /// List of primary species names to read from database
  std::vector<std::string> _ps_names;
  /// List of secondary equilibrium species to read from database
  std::vector<std::string> _es_names;
  /// List of secondary mineral species to read from database
  std::vector<std::string> _ms_names;
  /// Flag to indicate that read() has been called
  bool _read;
  /// Temperature points in database
  std::vector<Real> _temperature_points;
  /// Primary species data read from the database
  std::vector<GeochemicalDatabasePrimarySpecies> _primary_species;
  /// Secondary equilibrium species data read from the database
  std::vector<GeochemicalDatabaseEquilibriumSpecies> _equilibrium_species;
  /// Mineral species data read from the database
  std::vector<GeochemicalDatabaseMineralSpecies> _mineral_species;
  /// enum of database sections
  enum class sectionEnum
  {
    PRIMARY,
    SECONDARY,
    MINERAL
  } _section_enum;
};

#endif // GEOCHEMICALDATABASEREADER_H
