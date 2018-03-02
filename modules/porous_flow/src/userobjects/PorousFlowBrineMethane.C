//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
//* 
//* References:
//* [1] Duan, Moller and Weare, "An equation of state for the CH4-CO2-H2O system:
//*     I. Pure systems from 0 to 1000 degree C and 0 to 8000 bar", 
//*     Geochimica et Cosmochimica Acta Vol. 56, pp. 2605-2617, 1992
//* [2] Duan and Mao, "A thermodynamic model for calculating methane solubility, density
//*     and gas phase composition of methane-bearing aqueous fluids from 273 to 523 K
//*     and from 1 to 2000 bar", Geochimica et Cosmochimica Acta, 70, 3369-3386 (2006)

#include "PorousFlowBrineMethane.h"
#include "BrineFluidProperties.h"
#include "SinglePhaseFluidPropertiesPT.h"

// EOS parameters for methane. [Ref 1] Table 3
const Real CH4a1 = 8.72553928e-02; 
const Real CH4a2 = -7.52599476e-01; 
const Real CH4a3 = 3.75419887e-01;
const Real CH4a4 = 1.07291342e-02;
const Real CH4a5 = 5.49626360e-03;
const Real CH4a6 = -1.84772802e-02;
const Real CH4a7 = 3.18993183e-04;
const Real CH4a8 = 2.11079375e-04;
const Real CH4a9 = 2.01682801e-05;
const Real CH4a10 = -1.65606189e-05; 
const Real CH4a11 = 1.19614546e-04;
const Real CH4a12 = -1.08087289e-04; 
const Real CH4alpha = 4.48262295e-02;
const Real CH4beta = 7.5397e-01;
const Real CH4gamma = 7.7167e-02;


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

  // The liquid density should be calculated here
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
PorousFlowBrineMethane::EOSTempDependentCoefsMethane(
  Real Tr, Real & B, Real & C, Real & D, Real & E, Real & F) const
{
  // Compute temperature-dependent coefficients of EOS for methane.
  // Input is reduced temperature.

  // To avoid expensive recalculation
  Real inv_Tr_squared = 1.0 / (Tr*Tr);
  Real inv_Tr_cubed = 1.0 / std::pow(Tr,3.0);

  Real B = CH4a1 + CH4a2*inv_Tr_squared + CH4a3*inv_Tr_cubed;
  Real C = CH4a4 + CH4a5*inv_Tr_squared + CH4a6*inv_Tr_cubed;
  Real D = CH4a7 + CH4a8*inv_Tr_squared + CH4a9*inv_Tr_cubed;
  Real E = CH4a10 + CH4a11*inv_Tr_squared + CH4a12*inv_Tr_cubed;
  Real F = CH4alpha * inv_Tr_cubed;
}

Real 
PorousFlowBrineMethane::EOSBalanceMethane(
  Real Pr, Real Tr, Real Vr) const
{
  // EOS for methane, rearranged so that the output is zero when it is perfectly balanced.
  // All input are reduced properties.

  // To avoid expensive recalculation
  Real inv_Vr_squared = 1.0 / (Vr*Vr);

  EOSTempDependentCoefsMethane(Tr, B, C, D, E, F);

  return(1.0 + B/Vr + C*inv_Vr_squared + D/std::pow(Vr,4.0) + E/std::pow(Vr,5.0) +
      F * inv_Vr_squared * (CH4beta+CH4gamma*inv_Vr_squared) * std::exp(-CH4gamma*inv_Vr_squared) -
      Pr * Vr / Tr);
}

Real
PorousFlowBrineMethane::solveMolarVolumeMethane(Real Pr, Real Tr) const
{
  // Based on the Matlab code.
  // Given reduced pressure and temperature, return reduced molar volume.
  // Tolerance for Vr. 
  // TODO: Check for input valid range to avoid converging to wrong root.
  // QUESTION: Given the order of parameters, do we need to be this accurate?
  const Real tolerance = 1.0e-12;  
  Real Vr0 = 0.02;
  Real Vr1 = 110;
  Real balance0 = EOSBalanceMethane(Pr, Tr, Vr0);
  Real balance1 = EOSBalanceMethane(Pr, Tr, Vr1);
  if (balance0 * balance1 > 0.0) {
    // The method works only when balance0 * balance1 <= 0.0
    // TODO: capture this error
    return 0.0;
  }
  for (int i=0; i<100; ++i) {
    Vr = (Vr0+Vr1) * 0.5;
    balance = EOSBalanceMethane(Pr, Tr, Vr);
    if (abs(Vr0-Vr1) < tolerance) {
      return Vr;
    } else {
      if (balance * balance0 > 0.0) {
        Vr0 = Vr;
      } else {
        Vr1 = Vr;
      }
    }
  }
  // Finished without sufficiently converged.
  // TODO: capture this error
}

void
PorousFlowBrineMethane::fugacityCoefficientMethane(
    Real pressure, Real temperature, Real & fch4, Real & dfch4_dp, Real & dfch4_dT) const
{
  // Eq. A4 from Duan and Mao, which refers to [Ref. 1].
  //
  // Convert pressure from Pa to bar
  Real pbar = pressure * 1.0e-5;

  // Universal gas constant
  // TO DO: Maybe it already exists somewhere? Check
  const Real R = 0.08314472;  // bar.l/(mol.K)

  // Critical pressure and temperature of methane.
  // Note that Vc is the molar volume at critical pressure and temperature, 
  // and not critical molar volume.
  const Real Pc = 46.41;            // Pressure in bar
  const Real Tc = -82.55 + 273.15;  // Temperature in Kelvin
  const Real Vc = R * Tc / Pc;      // Molar volume in l/mol

  // Reduced properties
  Real Pr = pbar / Pc;
  Real Tr = temperature / Tc;
  Real Vr = solveMolarVolumeMethane(Pr, Tr);

  Real Z = Pr * Vr / Tr;
  EOSTempDependentCoefsMethane(Tr, B, C, D, E, F);
  Real gamma_by_square_Vr = CH4gamma / (Vr * Vr);
  Real G = 0.5 * F / CH4gamma * 
      (CH4beta + 1.0 - (beta + 1.0 + gamma_by_square_Vr) * std::exp(-gamma_by_square_Vr));

  // Fugacity coefficient
  fch4 = std::exp(Z - 1.0 - std::log(Z) + 
      B/Vr + C/(2.0*Vr*Vr) + D/(4.0*std::pow(Vr,4.0)) + E/(5.0*std::pow(Vr,5.0)) + G);

  // The derivative of fugacity coefficient wrt pressure
  dfch4_dp = ;

  // The derivative of fugacity coefficient wrt temperature
  dfch4_dT = ;
}

void
PorousFlowBrineMethane::fugacityCoefficientH2O(
    Real pressure, Real temperature, Real & fh2o, Real & dfh2o_dp, Real & dfh2o_dT) const
{
  // Eq. 6 from Duan and Mao
  // Convert pressure from Pa to bar
  Real pbar = pressure * 1.0e-5;

  // Parameters for eq. 6
  const Real a1 = -1.42006707e-2;
  const Real a2 = 1.08369910e-2;
  const Real a3 = -1.59213160e-6;
  const Real a4 = -1.10804676e-5;
  const Real a5 = -3.14287155;
  const Real a6 = 1.06338095e-3;

  // Fugacity coefficient
  fh2o = std::exp(a1 + a2*p + a3*pbar*pbar + a4*pbar*temperature + a5*pbar/temperature + a6*pbar*pbar/temperature);

  // The derivative of fugacity coefficient wrt pressure
  dfh2o_dp = fh2o * (a2 + 2*a3*pbar + a4*temperature + a5/temperature + 2*a6*pbar/temperature);

  // The derivative of fugacity coefficient wrt temperature
  inv_temp_squared = 1.0 / (temperature*temperature);
  dfh2o_dT = fh2o * (a4*pbar - a5*pbar*inv_temp_squared - a6*pbar*pbar*inv_temp_squared);
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
