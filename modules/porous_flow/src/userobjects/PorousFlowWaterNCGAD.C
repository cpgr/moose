//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PorousFlowWaterNCGAD.h"
#include "SinglePhaseFluidProperties.h"
#include "Conversion.h"

registerMooseObject("PorousFlowApp", PorousFlowWaterNCGAD);

template <>
InputParameters
validParams<PorousFlowWaterNCGAD>()
{
  InputParameters params = validParams<PorousFlowFluidStateBase>();
  params.addRequiredParam<UserObjectName>("water_fp", "The name of the user object for water");
  params.addRequiredParam<UserObjectName>(
      "gas_fp", "The name of the user object for the non-condensable gas");
  params.addClassDescription("Fluid state class for water and non-condensable gas");
  return params;
}

PorousFlowWaterNCGAD::PorousFlowWaterNCGAD(const InputParameters & parameters)
  : PorousFlowFluidStateBase(parameters),
    _water_fp(getUserObject<SinglePhaseFluidProperties>("water_fp")),
    _ncg_fp(getUserObject<SinglePhaseFluidProperties>("gas_fp")),
    _Mh2o(_water_fp.molarMass()),
    _Mncg(_ncg_fp.molarMass()),
    _water_triple_temperature(_water_fp.triplePointTemperature()),
    _water_critical_temperature(_water_fp.criticalTemperature())
{
  // Check that the correct FluidProperties UserObjects have been provided
  if (_water_fp.fluidName() != "water")
    paramError("water_fp", "A valid water FluidProperties UserObject must be provided in water_fp");

  // Set the number of phases and components, and their indexes
  _num_phases = 2;
  _num_components = 2;
  _gas_phase_number = 1 - _aqueous_phase_number;
  _gas_fluid_component = 1 - _aqueous_fluid_component;

  // Check that _aqueous_phase_number is <= total number of phases
  if (_aqueous_phase_number >= _num_phases)
    paramError("liquid_phase_number",
               "This value is larger than the possible number of phases ",
               _num_phases);

  // Check that _aqueous_fluid_component is <= total number of fluid components
  if (_aqueous_fluid_component >= _num_components)
    paramError("liquid_fluid_component",
               "This value is larger than the possible number of fluid components",
               _num_components);

  // Set the size of the FluidStateProperties vector
  _ad_fsp.resize(_num_phases, ADFluidStateProperties(_num_components));
}

std::string
PorousFlowWaterNCGAD::fluidStateName() const
{
  return "water-ncg";
}

void
PorousFlowWaterNCGAD::thermophysicalProperties(Real pressure,
                                               Real temperature,
                                               Real /* Xnacl */,
                                               Real Z,
                                               unsigned int qp,
                                               std::vector<FluidStateProperties> & fsp) const
{
  DualReal p = pressure;
  p.derivatives()[0] = 1.0;

  DualReal T = temperature;
  T.derivatives()[1] = 1.0;

  DualReal z = Z;
  z.derivatives()[2] = 1.0;

  // Check whether the input temperature is within the region of validity
  checkVariables(temperature);

  // Clear all of the FluidStateProperties data
  clearFluidStateProperties(fsp);

  std::vector<ADFluidStateProperties> adfsp(_num_phases, ADFluidStateProperties(_num_components));
  auto & liquid = adfsp[_aqueous_phase_number];
  auto & gas = adfsp[_gas_phase_number];

  FluidStatePhaseEnum phase_state;
  massFractions(p, T, z, phase_state, adfsp);

  switch (phase_state)
  {
    case FluidStatePhaseEnum::GAS:
    {
      // Set the gas saturations
      gas.saturation = 1.0;

      // Calculate gas properties
      gasProperties(p, T, adfsp);

      break;
    }

    case FluidStatePhaseEnum::LIQUID:
    {
      // Calculate the liquid properties
      liquidProperties(p, T, qp, adfsp);

      break;
    }

    case FluidStatePhaseEnum::TWOPHASE:
    {
      // Calculate the gas and liquid properties in the two phase region
      twoPhaseProperties(p, T, z, qp, adfsp);

      break;
    }
  }

  // Liquid saturations can now be set
  liquid.saturation = 1.0 - gas.saturation;

  // Save pressures to FluidStateProperties object
  gas.pressure = pressure;
  liquid.pressure = pressure - _pc.capillaryPressure(liquid.saturation, qp);

  // The last thing is to fill in the fsp
  fsp[_aqueous_phase_number].pressure = liquid.pressure.value();

  fsp[_aqueous_phase_number].density = liquid.density.value();
  fsp[_aqueous_phase_number].ddensity_dp = liquid.density.derivatives()[0];
  fsp[_aqueous_phase_number].ddensity_dT = liquid.density.derivatives()[1];
  fsp[_aqueous_phase_number].ddensity_dZ = liquid.density.derivatives()[2];

  fsp[_aqueous_phase_number].viscosity = liquid.viscosity.value();
  fsp[_aqueous_phase_number].dviscosity_dp = liquid.viscosity.derivatives()[0];
  fsp[_aqueous_phase_number].dviscosity_dT = liquid.viscosity.derivatives()[1];
  fsp[_aqueous_phase_number].dviscosity_dZ = liquid.viscosity.derivatives()[2];

  fsp[_aqueous_phase_number].enthalpy = liquid.enthalpy.value();
  fsp[_aqueous_phase_number].denthalpy_dp = liquid.enthalpy.derivatives()[0];
  fsp[_aqueous_phase_number].denthalpy_dT = liquid.enthalpy.derivatives()[1];
  fsp[_aqueous_phase_number].denthalpy_dZ = liquid.enthalpy.derivatives()[2];

  fsp[_aqueous_phase_number].mass_fraction[_aqueous_fluid_component] =
      liquid.mass_fraction[_aqueous_fluid_component].value();
  fsp[_aqueous_phase_number].dmass_fraction_dp[_aqueous_fluid_component] =
      liquid.mass_fraction[_aqueous_fluid_component].derivatives()[0];
  fsp[_aqueous_phase_number].dmass_fraction_dT[_aqueous_fluid_component] =
      liquid.mass_fraction[_aqueous_fluid_component].derivatives()[1];
  fsp[_aqueous_phase_number].dmass_fraction_dZ[_aqueous_fluid_component] =
      liquid.mass_fraction[_aqueous_fluid_component].derivatives()[2];

  fsp[_aqueous_phase_number].mass_fraction[_gas_fluid_component] =
      liquid.mass_fraction[_gas_fluid_component].value();
  fsp[_aqueous_phase_number].dmass_fraction_dp[_gas_fluid_component] =
      liquid.mass_fraction[_gas_fluid_component].derivatives()[0];
  fsp[_aqueous_phase_number].dmass_fraction_dT[_gas_fluid_component] =
      liquid.mass_fraction[_gas_fluid_component].derivatives()[1];
  fsp[_aqueous_phase_number].dmass_fraction_dZ[_gas_fluid_component] =
      liquid.mass_fraction[_gas_fluid_component].derivatives()[2];

  fsp[_aqueous_phase_number].saturation = liquid.saturation.value();
  fsp[_aqueous_phase_number].dsaturation_dp = liquid.saturation.derivatives()[0];
  fsp[_aqueous_phase_number].dsaturation_dT = liquid.saturation.derivatives()[1];
  fsp[_aqueous_phase_number].dsaturation_dZ = liquid.saturation.derivatives()[2];

  // Gas
  fsp[_gas_phase_number].pressure = gas.pressure.value();

  fsp[_gas_phase_number].density = gas.density.value();
  fsp[_gas_phase_number].ddensity_dp = gas.density.derivatives()[0];
  fsp[_gas_phase_number].ddensity_dT = gas.density.derivatives()[1];
  fsp[_gas_phase_number].ddensity_dZ = gas.density.derivatives()[2];

  fsp[_gas_phase_number].viscosity = gas.viscosity.value();
  fsp[_gas_phase_number].dviscosity_dp = gas.viscosity.derivatives()[0];
  fsp[_gas_phase_number].dviscosity_dT = gas.viscosity.derivatives()[1];
  fsp[_gas_phase_number].dviscosity_dZ = gas.viscosity.derivatives()[2];

  fsp[_gas_phase_number].enthalpy = gas.enthalpy.value();
  fsp[_gas_phase_number].denthalpy_dp = gas.enthalpy.derivatives()[0];
  fsp[_gas_phase_number].denthalpy_dT = gas.enthalpy.derivatives()[1];
  fsp[_gas_phase_number].denthalpy_dZ = gas.enthalpy.derivatives()[2];

  fsp[_gas_phase_number].mass_fraction[_aqueous_fluid_component] =
      gas.mass_fraction[_aqueous_fluid_component].value();
  fsp[_gas_phase_number].dmass_fraction_dp[_aqueous_fluid_component] =
      gas.mass_fraction[_aqueous_fluid_component].derivatives()[0];
  fsp[_gas_phase_number].dmass_fraction_dT[_aqueous_fluid_component] =
      gas.mass_fraction[_aqueous_fluid_component].derivatives()[1];
  fsp[_gas_phase_number].dmass_fraction_dZ[_aqueous_fluid_component] =
      gas.mass_fraction[_aqueous_fluid_component].derivatives()[2];

  fsp[_gas_phase_number].mass_fraction[_gas_fluid_component] =
      gas.mass_fraction[_gas_fluid_component].value();
  fsp[_gas_phase_number].dmass_fraction_dp[_gas_fluid_component] =
      gas.mass_fraction[_gas_fluid_component].derivatives()[0];
  fsp[_gas_phase_number].dmass_fraction_dT[_gas_fluid_component] =
      gas.mass_fraction[_gas_fluid_component].derivatives()[1];
  fsp[_gas_phase_number].dmass_fraction_dZ[_gas_fluid_component] =
      gas.mass_fraction[_gas_fluid_component].derivatives()[2];

  fsp[_gas_phase_number].saturation = gas.saturation.value();
  fsp[_gas_phase_number].dsaturation_dp = gas.saturation.derivatives()[0];
  fsp[_gas_phase_number].dsaturation_dT = gas.saturation.derivatives()[1];
  fsp[_gas_phase_number].dsaturation_dZ = gas.saturation.derivatives()[2];

  fsp[_gas_phase_number].d2saturation_dp2 = gas.d2saturation_dp2;
  fsp[_gas_phase_number].d2saturation_dpZ = gas.d2saturation_dpZ;
  fsp[_gas_phase_number].d2saturation_dZ2 = gas.d2saturation_dZ2;
  fsp[_aqueous_phase_number].d2saturation_dp2 = -gas.d2saturation_dp2;
  fsp[_aqueous_phase_number].d2saturation_dpZ = -gas.d2saturation_dpZ;
  fsp[_aqueous_phase_number].d2saturation_dZ2 = -gas.d2saturation_dZ2;
}

void
PorousFlowWaterNCGAD::massFractions(const DualReal & pressure,
                                    const DualReal & temperature,
                                    const DualReal & Z,
                                    FluidStatePhaseEnum & phase_state,
                                    std::vector<ADFluidStateProperties> & adfsp) const
{
  auto & liquid = adfsp[_aqueous_phase_number];
  auto & gas = adfsp[_gas_phase_number];

  // Equilibrium mass fraction of NCG in liquid and H2O in gas phases
  DualReal Xncg, Yh2o;
  equilibriumMassFractions(pressure, temperature, Xncg, Yh2o);

  DualReal Yncg = 1.0 - Yh2o;

  // Determine which phases are present based on the value of Z
  phaseState(Z.value(), Xncg.value(), Yncg.value(), phase_state);

  // The equilibrium mass fractions calculated above are only correct in the two phase
  // state. If only liquid or gas phases are present, the mass fractions are given by
  // the total mass fraction Z.
  DualReal Xh2o = 0.0;

  switch (phase_state)
  {
    case FluidStatePhaseEnum::LIQUID:
    {
      Xncg = Z;
      Yncg = 0.0;
      Xh2o = 1.0 - Z;
      Yh2o = 0.0;
      Xncg.derivatives()[0] = 0.0;
      Xncg.derivatives()[1] = 0.0;
      Xncg.derivatives()[2] = 1.0;
      break;
    }

    case FluidStatePhaseEnum::GAS:
    {
      Xncg = 0.0;
      Yncg = Z;
      Yh2o = 1.0 - Z;
      Yncg.derivatives()[0] = 0.0;
      Yncg.derivatives()[1] = 0.0;
      Yncg.derivatives()[2] = 1.0;
      break;
    }

    case FluidStatePhaseEnum::TWOPHASE:
    {
      // Keep equilibrium mass fractions
      Xh2o = 1.0 - Xncg;
      break;
    }
  }

  // Save the mass fractions in the FluidStateMassFractions object
  liquid.mass_fraction[_aqueous_fluid_component] = Xh2o;
  liquid.mass_fraction[_gas_fluid_component] = Xncg;
  gas.mass_fraction[_aqueous_fluid_component] = Yh2o;
  gas.mass_fraction[_gas_fluid_component] = Yncg;
}

void
PorousFlowWaterNCGAD::gasProperties(DualReal pressure,
                                    DualReal temperature,
                                    std::vector<ADFluidStateProperties> & adfsp) const
{
  auto & liquid = adfsp[_aqueous_phase_number];
  auto & gas = adfsp[_gas_phase_number];

  DualReal psat = _water_fp.vaporPressure(temperature);

  const DualReal Yncg = gas.mass_fraction[_gas_fluid_component];
  const DualReal Xncg = liquid.mass_fraction[_gas_fluid_component];

  // Note: take derivative wrt (Yncg * p) etc, hence the dp0
  DualReal ncg_density, ncg_viscosity, ncg_enthalpy;
  DualReal vapor_density, vapor_viscosity, vapor_enthalpy;

  // NCG density, viscosity and enthalpy calculated using partial pressure
  // Yncg * gas_poreressure (Dalton's law)
  _ncg_fp.rho_mu_from_p_T(Yncg * pressure, temperature, ncg_density, ncg_viscosity);
  ncg_enthalpy = _ncg_fp.h_from_p_T(Yncg * pressure, temperature);

  // Vapor density, viscosity and enthalpy calculated using partial pressure
  // X1 * psat (Raoult's law)
  _water_fp.rho_mu_from_p_T((1.0 - Xncg) * psat, temperature, vapor_density, vapor_viscosity);

  vapor_enthalpy = _water_fp.h_from_p_T((1.0 - Xncg) * psat, temperature);

  // Density is just the sum of individual component densities
  gas.density = ncg_density + vapor_density;

  // Viscosity of the gas phase is a weighted sum of the individual viscosities
  gas.viscosity = Yncg * ncg_viscosity + (1.0 - Yncg) * vapor_viscosity;

  // Enthalpy of the gas phase is a weighted sum of the individual enthalpies
  gas.enthalpy = Yncg * ncg_enthalpy + (1.0 - Yncg) * vapor_enthalpy;
}

void
PorousFlowWaterNCGAD::liquidProperties(DualReal pressure,
                                       DualReal temperature,
                                       unsigned int qp,
                                       std::vector<ADFluidStateProperties> & adfsp) const
{
  auto & liquid = adfsp[_aqueous_phase_number];
  auto & gas = adfsp[_gas_phase_number];

  const DualReal liquid_saturation = 1.0 - gas.saturation;
  DualReal liquid_pressure = pressure - _pc.capillaryPressure(liquid_saturation, qp);

  // Calculate liquid density and viscosity if in the two phase or single phase
  // liquid region, assuming they are not affected by the presence of dissolved
  // NCG. Note: the (small) contribution due to derivative of capillary pressure
  // wrt pressure (using the chain rule) is not implemented.
  DualReal liquid_density, liquid_viscosity;
  _water_fp.rho_mu_from_p_T(liquid_pressure, temperature, liquid_density, liquid_viscosity);

  liquid.density = liquid_density;
  liquid.viscosity = liquid_viscosity;

  // Enthalpy does include a contribution due to the enthalpy of dissolution
  DualReal hdis = enthalpyOfDissolution(temperature);

  DualReal water_enthalpy = _water_fp.h_from_p_T(liquid_pressure, temperature);
  DualReal ncg_enthalpy = _ncg_fp.h_from_p_T(pressure, temperature);

  const DualReal Xncg = liquid.mass_fraction[_gas_fluid_component];
  liquid.enthalpy = (1.0 - Xncg) * water_enthalpy + Xncg * (ncg_enthalpy + hdis);
}

DualReal
PorousFlowWaterNCGAD::liquidDensity(DualReal pressure, DualReal temperature) const
{
  return _water_fp.rho_from_p_T(pressure, temperature);
}

DualReal
PorousFlowWaterNCGAD::gasDensity(DualReal pressure,
                                 DualReal temperature,
                                 std::vector<ADFluidStateProperties> & adfsp) const
{
  auto & liquid = adfsp[_aqueous_phase_number];
  auto & gas = adfsp[_gas_phase_number];

  DualReal psat = _water_fp.vaporPressure(temperature);

  const DualReal Yncg = gas.mass_fraction[_gas_fluid_component];
  const DualReal Xncg = liquid.mass_fraction[_gas_fluid_component];

  // Note: take derivative wrt (Yncg * p) etc, hence the dp0
  DualReal ncg_density = _ncg_fp.rho_from_p_T(Yncg * pressure, temperature);
  DualReal vapor_density = _water_fp.rho_from_p_T((1.0 - Xncg) * psat, temperature);

  // Density is just the sum of individual component densities
  return ncg_density + vapor_density;
}

DualReal
PorousFlowWaterNCGAD::saturation(DualReal pressure,
                                 DualReal temperature,
                                 DualReal Z,
                                 std::vector<ADFluidStateProperties> & adfsp) const
{
  auto & gas = adfsp[_gas_phase_number];
  auto & liquid = adfsp[_aqueous_fluid_component];

  // Approximate liquid density as saturation isn't known yet, by using the gas
  // pressure rather than the liquid pressure. This does result in a small error
  // in the calculated saturation, but this is below the error associated with
  // the correlations. A more accurate saturation could be found iteraviely,
  // at the cost of increased computational expense

  // The gas and liquid densities
  const DualReal gas_density = gasDensity(pressure, temperature, adfsp);
  const DualReal liquid_density = liquidDensity(pressure, temperature);

  // Set mass equilibrium constants used in the calculation of vapor mass fraction
  const DualReal Xncg = liquid.mass_fraction[_gas_fluid_component];
  const DualReal Yncg = gas.mass_fraction[_gas_fluid_component];

  const DualReal K0 = Yncg / Xncg;
  const DualReal K1 = (1.0 - Yncg) / (1.0 - Xncg);
  const DualReal vapor_mass_fraction = vaporMassFraction(Z, K0, K1);

  // The gas saturation in the two phase case
  const DualReal saturation = vapor_mass_fraction * liquid_density /
                              (gas_density + vapor_mass_fraction * (liquid_density - gas_density));

  // std::cout << "gas density " << gas_density << std::endl;
  // std::cout << "liquid density " << liquid_density << std::endl;
  //
  // const DualReal saturation = Z * Z;
  return saturation;
}

void
PorousFlowWaterNCGAD::twoPhaseProperties(DualReal pressure,
                                         DualReal temperature,
                                         DualReal Z,
                                         unsigned int qp,
                                         std::vector<ADFluidStateProperties> & adfsp) const
{
  auto & gas = adfsp[_gas_phase_number];
  auto & liquid = adfsp[_aqueous_fluid_component];

  const DualReal Xncg = liquid.mass_fraction[_gas_fluid_component];
  const DualReal Yncg = gas.mass_fraction[_gas_fluid_component];

  // Calculate all of the gas phase properties, as these don't depend on saturation
  gasProperties(pressure, temperature, adfsp);

  // The gas saturation in the two phase case
  gas.saturation = saturation(pressure, temperature, Z, adfsp);

  // The liquid pressure and properties can now be calculated
  liquidProperties(pressure, temperature, qp, adfsp);

  // In the two-phase region, we need second derivatives of saturation as well. We can
  // calculate these using finite differences
  FluidStatePhaseEnum phase_state;
  std::vector<ADFluidStateProperties> fsptmp(2, ADFluidStateProperties(2));

  const Real dp = 1.0e1;
  massFractions(pressure + dp, temperature, Z, phase_state, fsptmp);
  DualReal s1 = saturation(pressure + dp, temperature, Z, fsptmp);

  clearFluidStateProperties(fsptmp);
  massFractions(pressure - dp, temperature, Z, phase_state, fsptmp);
  DualReal s2 = saturation(pressure - dp, temperature, Z, fsptmp);

  gas.d2saturation_dp2 = (s1.value() - 2.0 * gas.saturation.value() + s2.value()) / dp / dp;

  const Real dZ = 1.0e-6;
  clearFluidStateProperties(fsptmp);
  massFractions(pressure, temperature, Z + dZ, phase_state, fsptmp);
  s1 = saturation(pressure, temperature, Z + dZ, fsptmp);

  clearFluidStateProperties(fsptmp);
  massFractions(pressure, temperature, Z - dZ, phase_state, fsptmp);
  s2 = saturation(pressure, temperature, Z - dZ, fsptmp);

  gas.d2saturation_dZ2 = (s1.value() - 2.0 * gas.saturation.value() + s2.value()) / dZ / dZ;

  clearFluidStateProperties(fsptmp);
  massFractions(pressure + dp, temperature, Z + dZ, phase_state, fsptmp);
  s1 = saturation(pressure + dp, temperature, Z + dZ, fsptmp);

  clearFluidStateProperties(fsptmp);
  massFractions(pressure + dp, temperature, Z - dZ, phase_state, fsptmp);
  s2 = saturation(pressure + dp, temperature, Z - dZ, fsptmp);

  clearFluidStateProperties(fsptmp);
  massFractions(pressure - dp, temperature, Z + dZ, phase_state, fsptmp);
  DualReal s3 = saturation(pressure - dp, temperature, Z + dZ, fsptmp);

  clearFluidStateProperties(fsptmp);
  massFractions(pressure - dp, temperature, Z - dZ, phase_state, fsptmp);
  DualReal s4 = saturation(pressure - dp, temperature, Z - dZ, fsptmp);

  gas.d2saturation_dpZ = (s1.value() - s2.value() - s3.value() + s4.value()) / (4.0 * dp * dZ);
  // gas.d2saturation_dpZ = 0.0;
  // // Derivatives of saturation wrt primary variables are approximated using finite
  // // differences. Note: must calculate change in mass fraction due to increased
  // // primary variable as well
  // const Real eps = 1.0e-6;
  // const Real dp = 1.0e0;
  // const Real dXncg_dp = liquid.dmass_fraction_dp[_gas_fluid_component];
  // const Real dYncg_dp = gas.dmass_fraction_dp[_gas_fluid_component];
  // Real s1 = saturation(pressure - dp, temperature, Z, Xncg - dp * dXncg_dp, Yncg - dp *
  // dYncg_dp); Real s2 = saturation(pressure + dp, temperature, Z, Xncg + dp * dXncg_dp, Yncg + dp
  // * dYncg_dp); gas.dsaturation_dp = (s2 - s1) / (2.0 * dp);
  // // gas.d2saturation_dp2 = (s2 - 2.0 * gas.saturation + s1) / dp / dp;
  //
  // const Real dT = temperature * eps;
  // const Real dXncg_dT = liquid.dmass_fraction_dT[_gas_fluid_component];
  // const Real dYncg_dT = gas.dmass_fraction_dT[_gas_fluid_component];
  // s2 = saturation(pressure, temperature + dT, Z, Xncg + dT * dXncg_dT, Yncg + dT * dYncg_dT);
  // gas.dsaturation_dT = (s2 - gas.saturation) / dT;
  //
  // const Real dZ = 1.0e-6;
  // // Xncg and Yncg don't depend on Z in the two phase region
  // s1 = saturation(pressure, temperature, Z - dZ, Xncg, Yncg);
  // s2 = saturation(pressure, temperature, Z + dZ, Xncg, Yncg);
  // gas.dsaturation_dZ = (s2 - s1) / (2.0 * dZ);
  // gas.d2saturation_dZ2 = (s2 - 2.0 * gas.saturation + s1) / dZ / dZ;

  // Mixed derivative d2S/dpdZ
  // Real sm0 =
  //     saturation(pressure + dp, temperature, Z + dZ, Xncg + dp * dXncg_dp, Yncg + dp * dYncg_dp);
  // Real sm1 =
  //     saturation(pressure + dp, temperature, Z - dZ, Xncg + dp * dXncg_dp, Yncg + dp * dYncg_dp);
  // Real sm2 =
  //     saturation(pressure - dp, temperature, Z + dZ, Xncg - dp * dXncg_dp, Yncg - dp * dYncg_dp);
  // Real sm3 =
  //     saturation(pressure - dp, temperature, Z - dZ, Xncg - dp * dXncg_dp, Yncg - dp * dYncg_dp);
  //
  // gas.d2saturation_dpZ = (sm0 - sm1 - sm2 + sm3) / 4.0 / dp / dZ;
}

void
PorousFlowWaterNCGAD::equilibriumMassFractions(const DualReal & pressure,
                                               const DualReal & temperature,
                                               DualReal & Xncg,
                                               DualReal & Yh2o) const
{
  // Equilibrium constants for each component (Henry's law for the NCG
  // component, and Raoult's law for water).
  const DualReal Kh = _ncg_fp.henryConstant(temperature);
  const DualReal psat = _water_fp.vaporPressure(temperature);

  const DualReal Kncg = Kh / pressure;
  const DualReal Kh2o = psat / pressure;

  // The mole fractions for the NCG component in the two component
  // case can be expressed in terms of the equilibrium constants only
  const DualReal xncg = (1.0 - Kh2o) / (Kncg - Kh2o);
  const DualReal yncg = Kncg * xncg;

  // Convert mole fractions to mass fractions
  Xncg = moleFractionToMassFraction(xncg);
  Yh2o = 1.0 - moleFractionToMassFraction(yncg);
}

DualReal
PorousFlowWaterNCGAD::moleFractionToMassFraction(const DualReal & xmol) const
{
  return xmol * _Mncg / (xmol * _Mncg + (1.0 - xmol) * _Mh2o);
}

void
PorousFlowWaterNCGAD::checkVariables(Real temperature) const
{
  // Check whether the input temperature is within the region of validity of this equation
  // of state (T_triple <= T <= T_critical)
  if (temperature < _water_triple_temperature || temperature > _water_critical_temperature)
    mooseException(name() + ": temperature " + Moose::stringify(temperature) +
                   " is outside range 273.16 K <= T <= 647.096 K");
}

DualReal
PorousFlowWaterNCGAD::enthalpyOfDissolution(DualReal temperature) const
{
  // Henry's constant
  DualReal Kh = _ncg_fp.henryConstant(temperature);

  DualReal hdis = -_R * temperature * temperature * Kh.derivatives()[1] / Kh / _Mncg;

  // Derivative of enthalpy of dissolution wrt temperature requires the second derivative of
  // Henry's constant wrt temperature. For simplicity, approximate this numerically
  const Real dT = temperature.value() * 1.0e-8;
  const DualReal t2 = temperature + dT;
  DualReal Kh2 = _ncg_fp.henryConstant(t2);

  const Real dhdis_dT = (-_R * t2 * t2 * Kh2.derivatives()[1] / Kh2 / _Mncg - hdis).value() / dT;

  for (std::size_t i = 0; i < temperature.derivatives().size(); ++i)
    hdis.derivatives()[i] = temperature.derivatives()[i] * dhdis_dT;

  return hdis;
}

Real
PorousFlowWaterNCGAD::totalMassFraction(
    Real pressure, Real temperature, Real /* Xnacl */, Real saturation, unsigned int qp) const
{
  // Check whether the input temperature is within the region of validity
  checkVariables(temperature);

  // As we do not require derivatives, we can simply ignore their initialisation
  const DualReal p = pressure;
  const DualReal T = temperature;

  // FluidStateProperties data structure
  std::vector<ADFluidStateProperties> fsp(_num_phases, ADFluidStateProperties(_num_components));
  auto & liquid = fsp[_aqueous_phase_number];
  auto & gas = fsp[_gas_phase_number];

  // Calculate equilibrium mass fractions in the two-phase state
  DualReal Xncg, Yh2o;
  equilibriumMassFractions(p, T, Xncg, Yh2o);

  // Save the mass fractions in the FluidStateMassFractions object
  const DualReal Yncg = 1.0 - Yh2o;
  liquid.mass_fraction[_aqueous_fluid_component] = 1.0 - Xncg;
  liquid.mass_fraction[_gas_fluid_component] = Xncg;
  gas.mass_fraction[_aqueous_fluid_component] = Yh2o;
  gas.mass_fraction[_gas_fluid_component] = Yncg;

  // Gas properties
  gasProperties(p, T, fsp);

  // Saturation is given
  gas.saturation = saturation;

  // Liquid properties
  liquidProperties(pressure, temperature, qp, fsp);

  // The total mass fraction of ncg (Z) can now be calculated
  const DualReal Z =
      (saturation * gas.density * Yncg + (1.0 - saturation) * liquid.density * Xncg) /
      (saturation * gas.density + (1.0 - saturation) * liquid.density);
  return Z.value();
}
