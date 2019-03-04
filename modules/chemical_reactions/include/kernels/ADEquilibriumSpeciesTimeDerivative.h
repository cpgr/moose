//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef ADEQUILIBRIUMSPECIESTIMEDERIVATIVE_H
#define ADEQUILIBRIUMSPECIESTIMEDERIVATIVE_H

#include "ADTimeDerivative.h"

template <ComputeStage compute_stage>
class ADEquilibriumSpeciesTimeDerivative;

declareADValidParams(ADEquilibriumSpeciesTimeDerivative);

template <ComputeStage compute_stage>
class ADEquilibriumSpeciesTimeDerivative : public ADTimeDerivative<compute_stage>
{
public:
  ADEquilibriumSpeciesTimeDerivative(const InputParameters & parameters);

protected:
  virtual ADResidual precomputeQpResidual() override;

  /// Weight of the equilibrium species in the total primary species
  const Real _weight;
  /// Equilibrium constant for the equilibrium species
  const ADVariableValue & _log_k;
  /// Stoichiometric coefficient of the primary species in the equilibrium species
  const Real _sto_u;
  /// Stoichiometric coefficients of the coupled primary species in the equilibrium species
  const std::vector<Real> _sto_v;
  // /// Activity coefficient of primary species in the equilibrium species
  // const ADMaterialProperty(Real) & _gamma_u;
  // /// Old activity coefficient of primary species in the equilibrium species
  // const ADMaterialProperty(Real) & _gamma_u_old;
  // /// Activity coefficients of coupled primary species in the equilibrium species
  // std::vector<const ADMaterialProperty(Real) *> _gamma_v;
  // /// Old activity coefficients of coupled primary species in the equilibrium species
  // std::vector<const ADMaterialProperty(Real) *> _gamma_v_old;
  // /// Activity coefficient of equilibrium species
  // const ADMaterialProperty(Real) & _gamma_eq;
  // /// Old activity coefficient of equilibrium species
  // const ADMaterialProperty(Real) & _gamma_eq_old;
  /// Porosity
  const MaterialProperty<Real> & _porosity;
  /// Old porosity
  const MaterialProperty<Real> & _porosity_old;
  // /// Coupled primary species variable numbers
  // std::vector<unsigned int> _vars;
  // /// Coupled primary species concentrations
  // std::vector<const ADVariableValue *> _v_vals;
  // /// Old values of coupled primary species concentrations
  // std::vector<const ADVariableValue *> _v_vals_old;
  /// Old value of the primary species concentration.
  // const ADVariableValue & _u_old;

  usingTimeKernelMembers;
};

#endif // ADEQUILIBRIUMSPECIESTIMEDERIVATIVE_H
