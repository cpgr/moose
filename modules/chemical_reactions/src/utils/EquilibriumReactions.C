//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "EquilibriumReactions.h"
#include "MooseError.h"
#include "MooseArray.h"

namespace EquilibriumReactions
{
Real
equilibriumSpecies(std::vector<Real> & c_j,
                   std::vector<Real> & gamma_j,
                   std::vector<Real> & sto_j,
                   Real gamma_i,
                   Real log_Ki)
{
  return equilibriumSpeciesTempl<Real>(c_j, gamma_j, sto_j, gamma_i, log_Ki);
}

DualReal
equilibriumSpecies(std::vector<DualReal> & c_j,
                   std::vector<DualReal> & gamma_j,
                   std::vector<DualReal> & sto_j,
                   DualReal gamma_i,
                   DualReal log_Ki)
{
  return equilibriumSpeciesTempl<DualReal>(c_j, gamma_j, sto_j, gamma_i, log_Ki);
}

Real
equilibriumSpecies(std::vector<const VariableValue *> & c_j,
                   std::vector<const VariableValue *> & gamma_j,
                   const std::vector<Real> & sto_j,
                   const VariableValue & gamma_i,
                   const VariableValue & log_Ki,
                   unsigned int qp)
{
  Real c_i = 1.0;

  for (auto i = beginIndex(c_j); i < c_j.size(); ++i)
    c_i *= std::pow((*gamma_j[i])[qp] * (*c_j[i])[qp], sto_j[i]);

  return std::pow(10.0, log_Ki[qp]) * c_i / gamma_i[qp];
}

template <typename T>
T
equilibriumSpeciesTempl(
    std::vector<T> & c_j, std::vector<T> & gamma_j, std::vector<T> & sto_j, T gamma_i, T log_Ki)
{
  mooseAssert(c_j.size() == gamma_j.size(),
              "Number of primary species must be equal to number of activity coefficients");
  mooseAssert(c_j.size() == sto_j.size(),
              "Number of primary species must be equal to number of stoichiometric coefficients");
  mooseAssert(gamma_i > 0.0, "Activity coefficient must be greater than zero");

  T c_i = 1.0;

  for (std::size_t i = 0; i < c_j.size(); ++i)
    c_i *= std::pow(gamma_j[i] * c_j[i], sto_j[i]);

  return std::pow(10.0, log_Ki) * c_i / gamma_i;
}
} // namespace EquilibriumReactions
