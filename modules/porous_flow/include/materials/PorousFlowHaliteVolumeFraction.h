//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "PorousFlowMaterialVectorBase.h"

/**
 * Material that computes the solid halite volume fraction (m^3 halite / m^3 porous medium) from
 * the local-equilibrium precipitated-salt mass fraction reported by a salt-precipitating fluid
 * state (e.g. PorousFlowBrineCO2 with precipitate_salt = true).
 *
 * The fluid state reports precipitated_salt = m_h / (m_g + m_l), the solid halite mass per unit
 * fluid mass. Multiplying by the fluid mass per unit medium volume, phi * Sum_ph(S_ph rho_ph),
 * and dividing by the halite density gives the volume fraction
 *
 *   c_halite = precipitated_salt * phi * Sum_ph(S_ph rho_ph) / rho_halite.
 *
 * The OLD porosity is used to break the cyclic dependence between porosity and halite volume
 * fraction, exactly as PorousFlowAqueousPreDisMineral does for kinetic minerals.
 */
class PorousFlowHaliteVolumeFraction : public PorousFlowMaterialVectorBase
{
public:
  static InputParameters validParams();

  PorousFlowHaliteVolumeFraction(const InputParameters & parameters);

protected:
  void initQpStatefulProperties() override;
  void computeQpProperties() override;

  /// Density of solid halite (kg/m^3)
  const Real _halite_density;

  /// Precipitated salt mass fraction (kg solid halite / kg fluid) from the fluid state
  const MaterialProperty<Real> & _precipitated_salt;

  /// Derivative of the precipitated salt mass fraction wrt the PorousFlow variables
  const MaterialProperty<std::vector<Real>> & _dprecipitated_salt_dvar;

  /// Phase saturations
  const MaterialProperty<std::vector<Real>> & _saturation;

  /// Derivative of the phase saturations wrt the PorousFlow variables
  const MaterialProperty<std::vector<std::vector<Real>>> & _dsaturation_dvar;

  /// Phase densities
  const MaterialProperty<std::vector<Real>> & _fluid_density;

  /// Derivative of the phase densities wrt the PorousFlow variables
  const MaterialProperty<std::vector<std::vector<Real>>> & _dfluid_density_dvar;

  /// Old porosity (used to break the porosity <-> halite-concentration cyclic dependency)
  const MaterialProperty<Real> & _porosity_old;

  /// Computed halite volume fraction (m^3 halite / m^3 porous medium)
  MaterialProperty<Real> & _halite_volume_fraction;

  /// Derivative of the halite volume fraction wrt the PorousFlow variables
  MaterialProperty<std::vector<Real>> & _dhalite_volume_fraction_dvar;
};
