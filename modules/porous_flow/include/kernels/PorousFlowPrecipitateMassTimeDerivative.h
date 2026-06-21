//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "TimeKernel.h"
#include "PorousFlowDictator.h"

/**
 * Time derivative of the precipitated (solid) halite mass, applied to the salt-component equation.
 *
 * Computes d/dt(rho_halite * c_halite), where c_halite is the halite volume fraction
 * (PorousFlow_halite_volume_fraction, produced by PorousFlowHaliteVolumeFraction). Together with
 * PorousFlowMassTimeDerivative on the dissolved salt, this conserves total salt as the brine dries
 * out. The mass is lumped to the nodes, consistent with PorousFlowMassTimeDerivative.
 */
class PorousFlowPrecipitateMassTimeDerivative : public TimeKernel
{
public:
  static InputParameters validParams();

  PorousFlowPrecipitateMassTimeDerivative(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  /// Derivative of the residual wrt the PorousFlow variable pvar
  Real computeQpJac(unsigned int pvar) const;

  /// PorousFlow UserObject holding the list of PorousFlow variables
  const PorousFlowDictator & _dictator;

  /// Whether the kernel variable is a PorousFlow variable
  const bool _var_is_porflow_var;

  /// Density of solid halite (kg/m^3)
  const Real _halite_density;

  /// Halite volume fraction (m^3 halite / m^3 porous medium), lumped to nodes
  const MaterialProperty<Real> & _halite_volume_fraction;

  /// Old value of the halite volume fraction
  const MaterialProperty<Real> & _halite_volume_fraction_old;

  /// Derivative of the halite volume fraction wrt the PorousFlow variables
  const MaterialProperty<std::vector<Real>> & _dhalite_volume_fraction_dvar;
};
