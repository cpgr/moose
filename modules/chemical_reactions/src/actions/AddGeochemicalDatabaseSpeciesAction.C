//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "AddGeochemicalDatabaseSpeciesAction.h"
#include "AddVariableAction.h"
#include "FEProblem.h"
#include "Factory.h"

#include "libmesh/string_to_enum.h"

template <>
InputParameters
validParams<AddGeochemicalDatabaseSpeciesAction>()
{
  InputParameters params = validParams<Action>();
  params.addRequiredParam<FileName>("filename", "The geochemical database");
  params.addRequiredParam<std::vector<std::string>>("primary_species",
                                                    "The list of primary variables to add");
  params.addParam<std::vector<std::string>>("secondary_species",
                                            "The list of secondary equilibrium species to add");
  params.addParam<std::vector<std::string>>("mineral_species",
                                            "The list of mineral species to add");
  params.addCoupledVar(
      "temperature", 298.15, "The temperature of the aqueous phase (K). Default is 298.15K");
  // Get MooseEnums for the possible order/family options for this variable
  MooseEnum families(AddVariableAction::getNonlinearVariableFamilies());
  MooseEnum orders(AddVariableAction::getNonlinearVariableOrders());
  params.addParam<MooseEnum>(
      "family", families, "Specifies the family of FE shape function to use");
  params.addParam<MooseEnum>(
      "order", orders, "Specifies the order of the FE shape function to use");
  params.addParam<Real>(
      "scaling", 1.0, "Specifies a scaling factor to apply to nonlinear variables");
  params.addParam<std::string>("water", "H2O", "Water species in database");
  params.addClassDescription("Adds all variables, auxilliary variables and kernels for given "
                             "species using the given database");
  return params;
}

AddGeochemicalDatabaseSpeciesAction::AddGeochemicalDatabaseSpeciesAction(
    const InputParameters & params)
  : Action(params),
    _filename(getParam<FileName>("filename")),
    _primary_species_names(getParam<std::vector<std::string>>("primary_species")),
    _secondary_species_names(getParam<std::vector<std::string>>("secondary_species")),
    _mineral_species_names(getParam<std::vector<std::string>>("mineral_species")),
    _temperature(getParam<std::vector<VariableName>>("temperature")),
    _scaling(getParam<Real>("scaling")),
    _fe_type(Utility::string_to_enum<Order>(getParam<MooseEnum>("order")),
             Utility::string_to_enum<FEFamily>(getParam<MooseEnum>("family"))),
    _h2o(getParam<std::string>("water")),
    _invalid_logk(500.0),
    _T_c2k(273.15)
{
  // Read database
  GeochemicalDatabaseReader database(_filename);
  database.setPrimarySpeciesNames(_primary_species_names);
  database.setEquilibriumSpeciesNames(_secondary_species_names);
  database.setMineralSpeciesNames(_mineral_species_names);

  database.read();

  database.getTemperatures(_temperature_points);
  database.getPrimarySpecies(_primary_species);
  database.getEquilibriumSpecies(_equilibrium_species);
  database.getMineralSpecies(_mineral_species);

  // Convert temperature from C to K
  for (auto & tp : _temperature_points)
    tp += _T_c2k;

  // Remove any water species as they are not used for any kernels etc
  removeH2O(_equilibrium_species);
  removeH2O(_mineral_species);

  // Check that all primary species in each secondary reaction (equilibrium and mineral)
  // read from the database have been included in the list of primary species
  checkPrimarySpecies(_equilibrium_species);
  checkPrimarySpecies(_mineral_species);

  // Print out details of the equilibrium reactions to the console
  if (!_equilibrium_species.empty())
  {
    std::vector<std::string> eq_reactions;
    database.equilibriumReactions(eq_reactions);

    _console << "Aqueous equilibrium reactions:\n";
    for (auto i = beginIndex(eq_reactions); i < eq_reactions.size(); ++i)
      _console << "  Reaction " << i + 1 << ": " << eq_reactions[i] << "\n";
    _console << "\n";
  }

  // Print out details of the mineral reactions to the console
  if (!_mineral_species.empty())
  {
    std::vector<std::string> mineral_reactions;
    database.mineralReactions(mineral_reactions);

    _console << "Mineral reactions:\n";
    for (auto i = beginIndex(mineral_reactions); i < mineral_reactions.size(); ++i)
      _console << "  Reaction " << i + 1 << ": " << mineral_reactions[i] << "\n";
    _console << "\n";
  }
}

void
AddGeochemicalDatabaseSpeciesAction::act()
{
  // Add nonlinear variables for each primary species
  if (_current_task == "add_variable")
    for (auto i = beginIndex(_primary_species); i < _primary_species.size(); ++i)
      _problem->addVariable(_primary_species[i].name, _fe_type, _scaling);

  // Add auxvariables for each secondary equilibrium and mineral species, as well
  // as equilibrium constants
  if (_current_task == "add_aux_variable")
  {
    if (!_equilibrium_species.empty())
      for (auto i = beginIndex(_equilibrium_species); i < _equilibrium_species.size(); ++i)
      {
        _problem->addAuxVariable(_equilibrium_species[i].name, _fe_type);
        _problem->addAuxVariable(_equilibrium_species[i].name + "_logk", _fe_type);
      }

    if (!_mineral_species.empty())
      for (auto i = beginIndex(_mineral_species); i < _mineral_species.size(); ++i)
      {
        _problem->addAuxVariable(_mineral_species[i].name, _fe_type);
        _problem->addAuxVariable(_mineral_species[i].name + "_logk", _fe_type);
      }
  }

  if (_current_task == "add_kernel")
  {
    // Add kernels for primary species
    for (auto i = beginIndex(_primary_species); i < _primary_species.size(); ++i)
    {
      InputParameters params_td = _factory.getValidParams("PrimaryTimeDerivative");
      params_td.set<NonlinearVariableName>("variable") = _primary_species[i].name;
      _problem->addKernel("PrimaryTimeDerivative", _primary_species[i].name + "_td", params_td);
    }

    // Add kernels for secondary equilibrium species
    if (!_equilibrium_species.empty())
    {
      for (auto i = beginIndex(_equilibrium_species); i < _equilibrium_species.size(); ++i)
        for (auto j = beginIndex(_equilibrium_species[i].primary_species);
             j < _equilibrium_species[i].primary_species.size();
             ++j)
        {
          // Temporary vector of stoichiometric coefficients and primary species with
          // the current primary species and water species removed
          std::vector<Real> sto_v = _equilibrium_species[i].stoichiometric_coeff;
          std::vector<VariableName> v_names(_equilibrium_species[i].primary_species.size());
          for (auto k = beginIndex(v_names); k < v_names.size(); ++k)
            v_names[k] = _equilibrium_species[i].primary_species[k];

          sto_v.erase(sto_v.begin() + j);
          v_names.erase(v_names.begin() + j);

          InputParameters params_sub = _factory.getValidParams("CoupledBEEquilibriumSub");
          params_sub.set<NonlinearVariableName>("variable") =
              _equilibrium_species[i].primary_species[j];
          params_sub.set<Real>("weight") = 1.0;
          params_sub.set<std::vector<VariableName>>("log_k") = {_equilibrium_species[i].name +
                                                                "_logk"};
          params_sub.set<Real>("sto_u") = _equilibrium_species[i].stoichiometric_coeff[j];
          params_sub.set<std::vector<Real>>("sto_v") = sto_v;
          params_sub.set<std::vector<VariableName>>("v") = v_names;
          _problem->addKernel("CoupledBEEquilibriumSub",
                              _equilibrium_species[i].primary_species[j] + "_" +
                                  _equilibrium_species[i].name + "_sub",
                              params_sub);
        }
    }
  }

  if (_current_task == "add_aux_kernel")
  {
    // Add EquilibriumConstantAux AuxKernels for each equilibrium species
    if (!_equilibrium_species.empty())
    {
      for (auto i = beginIndex(_equilibrium_species); i < _equilibrium_species.size(); ++i)
      {
        // Temporary vectors of temperatures and logk with invalid data removed
        std::vector<Real> tp;
        std::vector<Real> logk;

        for (auto j = beginIndex(_temperature_points); j < _temperature_points.size(); ++j)
          if (_equilibrium_species[i].equilibrium_const[j] != _invalid_logk)
          {
            tp.push_back(_temperature_points[j]);
            logk.push_back(_equilibrium_species[i].equilibrium_const[j]);
          }

        InputParameters params_eqaux = _factory.getValidParams("EquilibriumConstantAux");
        params_eqaux.set<AuxVariableName>("variable") = _equilibrium_species[i].name + "_logk";
        params_eqaux.applySpecificParameters(parameters(), {"temperature"});
        params_eqaux.set<std::vector<Real>>("temperature_points") = tp;
        params_eqaux.set<std::vector<Real>>("logk_points") = logk;
        _problem->addAuxKernel("EquilibriumConstantAux",
                               "aux_" + _equilibrium_species[i].name + "_logk",
                               params_eqaux);
      }

      // Add AqueousEquilibriumRxnAux AuxKernels for equilibrium species
      for (auto i = beginIndex(_equilibrium_species); i < _equilibrium_species.size(); ++i)
      {
        // Temporary vector of variable names
        std::vector<VariableName> pspecies(_equilibrium_species[i].primary_species.size());
        for (auto j = beginIndex(pspecies); j < pspecies.size(); ++j)
          pspecies[j] = _equilibrium_species[i].primary_species[j];

        InputParameters params_eq = _factory.getValidParams("AqueousEquilibriumRxnAux");
        params_eq.set<AuxVariableName>("variable") = _equilibrium_species[i].name;
        params_eq.set<std::vector<VariableName>>("log_k") = {_equilibrium_species[i].name +
                                                             "_logk"};
        params_eq.set<std::vector<Real>>("sto_v") = _equilibrium_species[i].stoichiometric_coeff;
        params_eq.set<std::vector<VariableName>>("v") = pspecies;
        _problem->addAuxKernel(
            "AqueousEquilibriumRxnAux", "aux_" + _equilibrium_species[i].name, params_eq);
      }
    }
  }
}

template <typename T>
void
AddGeochemicalDatabaseSpeciesAction::removeH2O(std::vector<T> & species)
{
  // Remove any water species as they are not used for any kernels etc
  for (auto i = beginIndex(species); i < species.size(); ++i)
  {
    auto & ps = species[i].primary_species;
    auto & sc = species[i].stoichiometric_coeff;

    for (auto j = beginIndex(ps); j < ps.size(); ++j)
    {
      auto it = std::find(ps.begin(), ps.end(), _h2o);

      if (it != ps.end())
      {
        std::size_t idx = std::distance(ps.begin(), it);
        sc.erase(sc.begin() + idx);
        ps.erase(ps.begin() + idx);
      }
    }
  }
}

template <typename T>
void
AddGeochemicalDatabaseSpeciesAction::checkPrimarySpecies(std::vector<T> & species)
{
  for (auto i = beginIndex(species); i < species.size(); ++i)
    for (auto j = beginIndex(species[i].primary_species); j < species[i].primary_species.size();
         ++j)
    {
      if (std::find(_primary_species_names.begin(),
                    _primary_species_names.end(),
                    species[i].primary_species[j]) == _primary_species_names.end())
        mooseError("Primary species ",
                   species[i].primary_species[j],
                   " is required for ",
                   species[i].name,
                   " but has not been included in the list of primary species");
    }
}
