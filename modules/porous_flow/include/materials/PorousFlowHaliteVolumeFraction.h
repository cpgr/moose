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
 *
 * Templated on is_ad so an AD version (ADPorousFlowHaliteVolumeFraction) is available for the
 * (AD-only) finite-volume models. The AD path propagates derivatives automatically through the
 * generic precipitated-salt, saturation and density properties; the old porosity is deliberately
 * read as a plain (non-AD) coefficient so it contributes no current-step Jacobian term.
 */
template <bool is_ad>
class PorousFlowHaliteVolumeFractionTempl : public PorousFlowMaterialVectorBase
{
public:
  static InputParameters validParams();

  PorousFlowHaliteVolumeFractionTempl(const InputParameters & parameters);

protected:
  void initQpStatefulProperties() override;
  void computeQpProperties() override;

  /// Density of solid halite (kg/m^3)
  const Real _halite_density;

  /// Precipitated salt mass fraction (kg solid halite / kg fluid) from the fluid state
  const GenericMaterialProperty<Real, is_ad> & _precipitated_salt;

  /// Derivative of the precipitated salt mass fraction wrt the PorousFlow variables (non-AD only)
  const MaterialProperty<std::vector<Real>> * const _dprecipitated_salt_dvar;

  /// Phase saturations
  const GenericMaterialProperty<std::vector<Real>, is_ad> & _saturation;

  /// Derivative of the phase saturations wrt the PorousFlow variables (non-AD only)
  const MaterialProperty<std::vector<std::vector<Real>>> * const _dsaturation_dvar;

  /// Phase densities
  const GenericMaterialProperty<std::vector<Real>, is_ad> & _fluid_density;

  /// Derivative of the phase densities wrt the PorousFlow variables (non-AD only)
  const MaterialProperty<std::vector<std::vector<Real>>> * const _dfluid_density_dvar;

  /// Old porosity (used to break the porosity <-> halite-concentration cyclic dependency). Read as a
  /// plain Real even in the AD path: the old value is fixed data carrying no current-step derivative.
  const MaterialProperty<Real> & _porosity_old;

  /// Current porosity (used only at t = 0, where the old porosity is unavailable, to seed an
  /// oversaturated initial condition consistently)
  const GenericMaterialProperty<Real, is_ad> & _porosity;

  /// Computed halite volume fraction (m^3 halite / m^3 porous medium)
  GenericMaterialProperty<Real, is_ad> & _halite_volume_fraction;

  /// Derivative of the halite volume fraction wrt the PorousFlow variables (non-AD only)
  MaterialProperty<std::vector<Real>> * const _dhalite_volume_fraction_dvar;
};

typedef PorousFlowHaliteVolumeFractionTempl<false> PorousFlowHaliteVolumeFraction;
typedef PorousFlowHaliteVolumeFractionTempl<true> ADPorousFlowHaliteVolumeFraction;
