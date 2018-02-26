//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PorousFlowBrineMethane.h"
#include "BrineFluidProperties.h"
#include "SinglePhaseFluidPropertiesPT.h"

template <>
InputParameters
validParams<PorousFlowBrineMethane>()
{
  InputParameters params = validParams<PorousFlowFluidStateBase>();
  params.addRequiredParam<UserObjectName>("brine_fp", "The name of the user object for brine");
  params.addRequiredParam<UserObjectName>("methane_fp", "The name of the user object for methane");
  params.addClassDescription("Fluid state class for brine and methane (CH4)");
  return params;
}

PorousFlowBrineMethane::PorousFlowBrineMethane(const InputParameters & parameters)
  : PorousFlowFluidStateBase(parameters),
    _brine_fp(getUserObject<BrineFluidProperties>("brine_fp")),
    _ch4_fp(getUserObject<SinglePhaseFluidPropertiesPT>("methane_fp")),
    _water_fp(_brine_fp.getComponent(BrineFluidProperties::WATER)),
    _Mh2o(_brine_fp.molarMassH2O()),
    _invMh2o(1.0 / _Mh2o),
    _Mch4(_ch4_fp.molarMass()),
    _Mnacl(_brine_fp.molarMassNaCl())
{
  // Check that the correct FluidProperties UserObjects have been provided
  if (_ch4_fp.fluidName() != "methane")
    mooseError("Only a valid methane FluidProperties UserObject can be provided in methane_fp");
  if (_brine_fp.fluidName() != "brine")
    mooseError("Only a valid brine FluidProperties UserObject can be provided in brine_fp");
}

std::string
PorousFlowBrineMethane::fluidStateName() const
{
  return "brine-methane";
}

void
PorousFlowBrineMethane::thermophysicalProperties(Real pressure,
                                                 Real temperature,
                                                 Real xnacl,
                                                 Real z,
                                                 std::vector<FluidStateProperties> & fsp) const
{
  FluidStateProperties & liquid = fsp[_aqueous_phase_number];
  FluidStateProperties & gas = fsp[_gas_phase_number];

  // Check whether the input temperature is within the region of validity
  checkVariables(pressure, temperature);

  // Clear all of the FluidStateProperties data
  clearFluidStateProperties(fsp);

  FluidStatePhaseEnum phase_state;
  massFractions(pressure, temperature, xnacl, z, phase_state, fsp);

  switch (phase_state)
  {
    case FluidStatePhaseEnum::GAS:
    {
      // Set the gas saturations
      gas.saturation = 1.0;

      // Calculate gas properties
      gasProperties(pressure, temperature, fsp);

      break;
    }

    case FluidStatePhaseEnum::LIQUID:
    {
      // Calculate the liquid properties
      Real liquid_pressure = pressure - _pc_uo.capillaryPressure(1.0);
      liquidProperties(liquid_pressure, temperature, xnacl, fsp);

      break;
    }

    case FluidStatePhaseEnum::TWOPHASE:
    {
      // Calculate the gas properties
      gasProperties(pressure, temperature, fsp);

      // Calculate the saturation
      saturationTwoPhase(pressure, temperature, xnacl, z, fsp);

      // Calculate the liquid properties
      Real liquid_pressure = pressure - _pc_uo.capillaryPressure(1.0 - gas.saturation);
      liquidProperties(liquid_pressure, temperature, xnacl, fsp);

      break;
    }
  }

  // Liquid saturations can now be set
  liquid.saturation = 1.0 - gas.saturation;
  liquid.dsaturation_dp = -gas.dsaturation_dp;
  liquid.dsaturation_dT = -gas.dsaturation_dT;
  liquid.dsaturation_dz = -gas.dsaturation_dz;

  // Save pressures to FluidStateProperties object
  gas.pressure = pressure;
  liquid.pressure = pressure - _pc_uo.capillaryPressure(liquid.saturation);
}

void
PorousFlowBrineMethane::massFractions(Real pressure,
                                      Real temperature,
                                      Real xnacl,
                                      Real z,
                                      FluidStatePhaseEnum & phase_state,
                                      std::vector<FluidStateProperties> & fsp) const
{
  FluidStateProperties & liquid = fsp[_aqueous_phase_number];
  FluidStateProperties & gas = fsp[_gas_phase_number];

  // Equilibrium mass fraction of Methane in liquid and H2O in gas phases
  Real Xch4, dXch4_dp, dXch4_dT, Yh2o, dYh2o_dp, dYh2o_dT;
  equilibriumMassFractions(
      pressure, temperature, xnacl, Xch4, dXch4_dp, dXch4_dT, Yh2o, dYh2o_dp, dYh2o_dT);

  Real Ych4 = 1.0 - Yh2o;
  Real dYch4_dp = -dYh2o_dp;
  Real dYch4_dT = -dYh2o_dT;

  // Determine which phases are present based on the value of z
  phaseState(z, Xch4, Ych4, phase_state);

  // The equilibrium mass fractions calculated above are only correct in the two phase
  // state. If only liquid or gas phases are present, the mass fractions are given by
  // the total mass fraction z
  Real Xh2o = 0.0;
  Real dXch4_dz = 0.0, dYch4_dz = 0.0;

  switch (phase_state)
  {
    case FluidStatePhaseEnum::LIQUID:
    {
      Xch4 = z;
      Ych4 = 0.0;
      Xh2o = 1.0 - z;
      Yh2o = 0.0;
      dXch4_dp = 0.0;
      dXch4_dT = 0.0;
      dXch4_dz = 1.0;
      dYch4_dp = 0.0;
      dYch4_dT = 0.0;
      break;
    }

    case FluidStatePhaseEnum::GAS:
    {
      Xch4 = 0.0;
      Ych4 = z;
      Yh2o = 1.0 - z;
      dXch4_dp = 0.0;
      dXch4_dT = 0.0;
      dYch4_dz = 1.0;
      dYch4_dp = 0.0;
      dYch4_dT = 0.0;
      break;
    }

    case FluidStatePhaseEnum::TWOPHASE:
    {
      // Keep equilibrium mass fractions
      Xh2o = 1.0 - Xch4;
      break;
    }
  }

  // Save the mass fractions in the FluidStateMassFractions object
  liquid.mass_fraction[_aqueous_fluid_component] = Xh2o;
  liquid.mass_fraction[_gas_fluid_component] = Xch4;
  gas.mass_fraction[_aqueous_fluid_component] = Yh2o;
  gas.mass_fraction[_gas_fluid_component] = Ych4;

  // Save the derivatives wrt PorousFlow variables
  liquid.dmass_fraction_dp[_aqueous_fluid_component] = -dXch4_dp;
  liquid.dmass_fraction_dp[_gas_fluid_component] = dXch4_dp;
  liquid.dmass_fraction_dT[_aqueous_fluid_component] = -dXch4_dT;
  liquid.dmass_fraction_dT[_gas_fluid_component] = dXch4_dT;
  liquid.dmass_fraction_dz[_aqueous_fluid_component] = -dXch4_dz;
  liquid.dmass_fraction_dz[_gas_fluid_component] = dXch4_dz;

  gas.dmass_fraction_dp[_aqueous_fluid_component] = -dYch4_dp;
  gas.dmass_fraction_dp[_gas_fluid_component] = dYch4_dp;
  gas.dmass_fraction_dT[_aqueous_fluid_component] = -dYch4_dT;
  gas.dmass_fraction_dT[_gas_fluid_component] = dYch4_dT;
  gas.dmass_fraction_dz[_aqueous_fluid_component] = -dYch4_dz;
  gas.dmass_fraction_dz[_gas_fluid_component] = dYch4_dz;
}

void
PorousFlowBrineMethane::gasProperties(Real pressure,
                                      Real temperature,
                                      std::vector<FluidStateProperties> & fsp) const
{
  FluidStateProperties & gas = fsp[_gas_phase_number];

  // Gas density and viscosity are approximated with pure Methane - no correction due
  // to the small amount of water vapor is made
  Real ch4_density, dch4_density_dp, dch4_density_dT;
  Real ch4_viscosity, dch4_viscosity_drho, dch4_viscosity_dT;
  _ch4_fp.rho_dpT(pressure, temperature, ch4_density, dch4_density_dp, dch4_density_dT);
  _ch4_fp.mu_drhoT_from_rho_T(ch4_density,
                              temperature,
                              dch4_density_dT,
                              ch4_viscosity,
                              dch4_viscosity_drho,
                              dch4_viscosity_dT);

  // Save the values to the FluidStateProperties object. Note that derivatives wrt z are 0
  gas.density = ch4_density;
  gas.ddensity_dp = dch4_density_dp;
  gas.ddensity_dT = dch4_density_dT;
  gas.ddensity_dz = 0.0;

  gas.viscosity = ch4_viscosity;
  gas.dviscosity_dp = dch4_viscosity_drho * dch4_density_dp;
  gas.dviscosity_dT = dch4_viscosity_dT;
  gas.dviscosity_dz = 0.0;
}

void
PorousFlowBrineMethane::liquidProperties(Real pressure,
                                         Real temperature,
                                         Real xnacl,
                                         std::vector<FluidStateProperties> & fsp) const
{
  FluidStateProperties & liquid = fsp[_aqueous_phase_number];

  // The liquid density and viscosity should be calcualted here
  Real liquid_density;
  Real dliquid_density_dp;
  Real dliquid_density_dT;
  Real dliquid_density_dz;

  Real liquid_viscosity;
  Real dliquid_viscosity_drho;
  Real dliquid_viscosity_dT;
  Real dliquid_viscosity_dz;

  // Note: brine viscosity (and derivatives) requires water density (and derivatives)
  Real water_density, dwater_density_dp, dwater_density_dT;
  _water_fp.rho_dpT(pressure, temperature, water_density, dwater_density_dp, dwater_density_dT);

  // Save the values to the FluidStateProperties object
  liquid.density = liquid_density;
  liquid.ddensity_dp = dliquid_density_dp;
  liquid.ddensity_dT = dliquid_density_dT;
  liquid.ddensity_dz = dliquid_density_dz;

  liquid.viscosity = liquid_viscosity;
  liquid.dviscosity_dp = dliquid_viscosity_drho * dwater_density_dp;
  liquid.dviscosity_dT = dliquid_viscosity_dT;
  liquid.dviscosity_dz = dliquid_viscosity_dz;
}

void
PorousFlowBrineMethane::saturationTwoPhase(Real pressure,
                                           Real temperature,
                                           Real xnacl,
                                           Real z,
                                           std::vector<FluidStateProperties> & fsp) const
{
  FluidStateProperties & liquid = fsp[_aqueous_phase_number];
  FluidStateProperties & gas = fsp[_gas_phase_number];

  // Approximate liquid density as saturation isn't known yet
  Real brine_density, dbrine_density_dp, dbrine_density_dT, dbrine_density_dx;
  _brine_fp.rho_dpTx(pressure,
                     temperature,
                     xnacl,
                     brine_density,
                     dbrine_density_dp,
                     dbrine_density_dT,
                     dbrine_density_dx);

  // Mass fraction of Methane in liquid phase
  Real Xch4 = liquid.mass_fraction[_gas_fluid_component];
  Real dXch4_dp = liquid.dmass_fraction_dp[_gas_fluid_component];
  Real dXch4_dT = liquid.dmass_fraction_dT[_gas_fluid_component];

  // The liquid density should be calcualted here
  Real liquid_density;

  Real dliquid_density_dp;

  Real dliquid_density_dT;

  Real Ych4 = gas.mass_fraction[_gas_fluid_component];
  Real dYch4_dp = gas.dmass_fraction_dp[_gas_fluid_component];
  Real dYch4_dT = gas.dmass_fraction_dT[_gas_fluid_component];

  // Set mass equilibrium constants used in the calculation of vapor mass fraction
  Real K0 = Ych4 / Xch4;
  Real K1 = (1.0 - Ych4) / (1.0 - Xch4);
  Real vapor_mass_fraction = vaporMassFraction(z, K0, K1);

  // The gas saturation in the two phase case
  gas.saturation = vapor_mass_fraction * liquid_density /
                   (gas.density + vapor_mass_fraction * (liquid_density - gas.density));

  Real dv_dz = (K1 - K0) / ((K0 - 1.0) * (K1 - 1.0));
  Real denominator = (gas.density + vapor_mass_fraction * (liquid_density - gas.density)) *
                     (gas.density + vapor_mass_fraction * (liquid_density - gas.density));

  Real ds_dz = gas.density * liquid_density * dv_dz / denominator;

  Real dK0_dp = (Xch4 * dYch4_dp - Ych4 * dXch4_dp) / Xch4 / Xch4;
  Real dK0_dT = (Xch4 * dYch4_dT - Ych4 * dXch4_dT) / Xch4 / Xch4;

  Real dK1_dp = ((1.0 - Ych4) * dXch4_dp - (1.0 - Xch4) * dYch4_dp) / (1.0 - Xch4) / (1.0 - Xch4);
  Real dK1_dT = ((1.0 - Ych4) * dXch4_dT - (1.0 - Xch4) * dYch4_dT) / (1.0 - Xch4) / (1.0 - Xch4);

  Real dv_dp = z * dK1_dp / (K1 - 1.0) / (K1 - 1.0) + (1.0 - z) * dK0_dp / (K0 - 1.0) / (K0 - 1.0);

  Real ds_dp = gas.density * liquid_density * dv_dp +
               vapor_mass_fraction * (1.0 - vapor_mass_fraction) *
                   (gas.density * dliquid_density_dp - gas.ddensity_dp * liquid_density);
  ds_dp /= denominator;

  Real dv_dT = z * dK1_dT / (K1 - 1.0) / (K1 - 1.0) + (1.0 - z) * dK0_dT / (K0 - 1.0) / (K0 - 1.0);

  Real ds_dT = gas.density * liquid_density * dv_dT +
               vapor_mass_fraction * (1.0 - vapor_mass_fraction) *
                   (gas.density * dliquid_density_dT - gas.ddensity_dT * liquid_density);
  ds_dT /= denominator;

  gas.dsaturation_dp = ds_dp;
  gas.dsaturation_dT = ds_dT;
  gas.dsaturation_dz = ds_dz;
}

void
PorousFlowBrineMethane::equilibriumMassFractions(Real pressure,
                                                 Real temperature,
                                                 Real xnacl,
                                                 Real & Xch4,
                                                 Real & dXch4_dp,
                                                 Real & dXch4_dT,
                                                 Real & Yh2o,
                                                 Real & dYh2o_dp,
                                                 Real & dYh2o_dT) const
{
  // Mass fraction of water in gas: Yh2o given by Eq 5 from Duan and Mao

  // Mass fraction of methane in liquid: Xch4 calculated from Eq 3 by solving
  // for mch4 (moles/kg of CH4 in liquid) and dividing by _Mch4 (molar mass of methane)
}

void
PorousFlowBrineMethane::fugacityCoefficientMethane(
    Real pressure, Real temperature, Real & fch4, Real & dfch4_dp, Real & dfch4_dT) const
{
  // Eq. A4 from Duan and Mao
}

void
PorousFlowBrineMethane::fugacityCoefficientH2O(
    Real pressure, Real temperature, Real & fh2o, Real & dfh2o_dp, Real & dfh2o_dT) const
{
  // Eq. 6 from Duan and Mao
}

void
PorousFlowBrineMethane::activityCoefficient(Real pressure,
                                            Real temperature,
                                            Real xnacl,
                                            Real & gamma,
                                            Real & dgamma_dp,
                                            Real & dgamma_dT) const
{
  // I think Eq 10 from Duan et al, 199b reference of Duan and Mao
}

void
PorousFlowBrineMethane::equilibriumConstantH2O(Real temperature, Real & kh2o, Real & dkh2o_dT) const
{
  // Might not be needed?
}

void
PorousFlowBrineMethane::equilibriumConstantMethane(Real temperature,
                                                   Real & kch4,
                                                   Real & dkch4_dT) const
{
  // Might not be needed?
}

Real
PorousFlowBrineMethane::totalMassFraction(Real pressure,
                                          Real temperature,
                                          Real xnacl,
                                          Real saturation) const
{
  // Check whether the input pressure and temperature are within the region of validity
  checkVariables(pressure, temperature);

  // FluidStateProperties data structure
  std::vector<FluidStateProperties> fsp(_num_phases, FluidStateProperties(_num_components));
  FluidStateProperties & liquid = fsp[_aqueous_phase_number];
  FluidStateProperties & gas = fsp[_gas_phase_number];

  // Calculate equilibrium mass fractions in the two-phase state
  Real Xch4, dXch4_dp, dXch4_dT, Yh2o, dYh2o_dp, dYh2o_dT;
  equilibriumMassFractions(
      pressure, temperature, xnacl, Xch4, dXch4_dp, dXch4_dT, Yh2o, dYh2o_dp, dYh2o_dT);

  // Save the mass fractions in the FluidStateMassFractions object
  Real Ych4 = 1.0 - Yh2o;
  liquid.mass_fraction[_aqueous_fluid_component] = 1.0 - Xch4;
  liquid.mass_fraction[_gas_fluid_component] = Xch4;
  gas.mass_fraction[_aqueous_fluid_component] = Yh2o;
  gas.mass_fraction[_gas_fluid_component] = Ych4;

  // Gas properties
  gasProperties(pressure, temperature, fsp);

  // Liquid properties
  Real liquid_saturation = 1.0 - saturation;
  Real liquid_pressure = pressure - _pc_uo.capillaryPressure(liquid_saturation);
  liquidProperties(liquid_pressure, temperature, xnacl, fsp);

  // The total mass fraction of ncg (z) can now be calculated
  Real z = (saturation * gas.density * Ych4 + liquid_saturation * liquid.density * Xch4) /
           (saturation * gas.density + liquid_saturation * liquid.density);

  return z;
}

void
PorousFlowBrineMethane::checkVariables(Real pressure, Real temperature) const
{
  // The calculation of mass fractions is valid from 273K <= T <= 523K, and
  // pressure less than 60 MPa
  if (temperature < 273.0 || temperature > 523.0)
    mooseError("PorousFlowBrineMethane: Temperature is outside range 285.15 K <= T "
               "<= 373.15 K");

  if (pressure > 2.0e8)
    mooseError("PorousFlowBrineMethane: Pressure must be less than 200 MPa");
}
