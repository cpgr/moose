//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef ADSECONDARYSPECIESDIFFUSION_H
#define ADSECONDARYSPECIESDIFFUSION_H

#include "ADKernelGrad.h"

template <ComputeStage compute_stage>
class ADSecondarySpeciesDiffusion;

declareADValidParams(ADPrimarySpeciesDiffusion);

template <ComputeStage compute_stage>
class ADSecondarySpeciesDiffusion : public ADKernelGrad<compute_stage>
{
public:
  ADSecondarySpeciesDiffusion(const InputParameters & parameters);

protected:
  virtual ADVectorResidual precomputeQpResidual() override;

  /// Coupled primary species concentrations
  std::vector<const ADVariableValue *> _primary_species;
  /// Gradient of coupled primary species concentrations
  std::vector<const ADVariableGradient *> _grad_primary_species;
  /// Number of coupled primary species
  const unsigned int _n_primary;
  /// Name of the secondary (equilibrium) species in this kernel
  const AuxVariableName _secondary_name;
  /// Diffusion coefficient
  const MaterialProperty<Real> & _diffusivity;
  /// Weight of the equilibrium species
  const Real _weight;
  /// Equilibrium constant for the equilibrium species
  const ADVariableValue & _log_k;
  /// Stoichiometric coefficient of the primary species that this kernel acts on in the equilibrium species
  const Real _sto;
  /// Stoichiometric coefficients of the additional coupled primary species in the equilibrium species
  const std::vector<Real> _sto_coupled;
  /// Activity coefficients for this primary species
  const ADMaterialProperty(Real) & _gamma_primary;
  /// Activity coefficients for each coupled primary species
  std::vector<const ADMaterialProperty(Real) *> _gamma_coupled;
  /// Activity coefficients for this secondary species
  const ADMaterialProperty(Real) & _gamma_secondary;

  usingKernelMembers;
};

#endif // ADSECONDARYSPECIESDIFFUSION_H
