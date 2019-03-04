//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef ADACTIVITYCOEFFICIENT_H
#define ADACTIVITYCOEFFICIENT_H

#include "ADMaterial.h"

template <ComputeStage>
class ADActivityCoefficient;

declareADValidParams(ADActivityCoefficient);

template <ComputeStage compute_stage>
class ADActivityCoefficient : public ADMaterial<compute_stage>
{
public:
  ADActivityCoefficient(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  /// Primary species concentrations (Nonlinear variables)
  std::vector<const ADVariableValue *> _primary_species;
  /// Number of primary species
  const unsigned int _n_primary;
  /// Secondary (equilibrium) species present (AuxVariables)
  std::vector<const VariableValue *> _secondary_species;
  /// Number of secondary (equilibrium) species
  const unsigned int _n_secondary;
  /// Activity coefficients for each primary species
  std::vector<ADMaterialProperty(Real) *> _gamma_primary;
  /// Activity coefficients for each secondary species
  std::vector<ADMaterialProperty(Real) *> _gamma_secondary;

  usingMaterialMembers;
};

#endif // ADACTIVITYCOEFFICIENT_H
