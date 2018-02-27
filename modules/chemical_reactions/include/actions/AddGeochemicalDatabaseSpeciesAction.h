//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef ADDGEOCHEMICALDATABASESPECIESACTION_H
#define ADDGEOCHEMICALDATABASESPECIESACTION_H

#include "Action.h"
#include "GeochemicalDatabaseReader.h"

#include "libmesh/fe_type.h"

class AddGeochemicalDatabaseSpeciesAction;

template <>
InputParameters validParams<AddGeochemicalDatabaseSpeciesAction>();

class AddGeochemicalDatabaseSpeciesAction : public Action
{
public:
  AddGeochemicalDatabaseSpeciesAction(const InputParameters & params);

  virtual void act() override;

private:
  template <typename T>
  void removeH2O(std::vector<T> & species);

  template <typename T>
  void checkPrimarySpecies(std::vector<T> & species);

  /// Database filename
  const FileName _filename;
  /// List of primary species to read from database
  const std::vector<std::string> _primary_species_names;
  /// List of secondary equilibrium species to read from database
  const std::vector<std::string> _secondary_species_names;
  /// List of mineral species to read from database
  const std::vector<std::string> _mineral_species_names;
  /// Temperature (in K)
  const std::vector<VariableName> _temperature;
  /// Scaling of variables
  const Real _scaling;
  /// Finite element type
  const FEType _fe_type;
  /// String representing water species in database
  const std::string _h2o;
  /// Temperature points where the equilibrium constants are defined
  std::vector<Real> _temperature_points;
  /// Value of missing log(Keq) points
  const Real _invalid_logk;
  /// Primary species data read from the database
  std::vector<GeochemicalDatabasePrimarySpecies> _primary_species;
  /// Secondary equilibrium species data read from the database
  std::vector<GeochemicalDatabaseEquilibriumSpecies> _equilibrium_species;
  /// Mineral species data read from the database
  std::vector<GeochemicalDatabaseMineralSpecies> _mineral_species;
  /// Conversion from C to K
  const Real _T_c2k;
};

#endif // ADDGEOCHEMICALDATABASESPECIESACTION_H
