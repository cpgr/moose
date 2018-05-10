//* This module contains a thermodynamic model for methane in water (H2O-CH4 system) 
//* or brine (H2O-CH4-NaCl system).
//* Note:
//* - The liquid is mainly water, with some methane and perhaps NaCl dissolved in it
//* - The gas is mainly methane, with some water vapour in it
//* 
//* This file can be used as a standalone program or part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
//* 
//* References:
//* [1] A thermodynamics model for calculating methane solubility, density and gas phase
//* composition of methane-bearing aqueous fluids from 273 to 523 K and from 1 to 2000 bar,
//* Z. Duan and S. Mao, 
//* Geochimica et Cosmochimica Acta 70, 2006, pp. 3369-3386.
//* 
//* [2] An equation of state for the CH4-CO2-H2O system: I. Pure systems from 
//* 0 to 1000 degree C and 0 to 8000 bar, 
//* Z. Duan, N, Moller, JH. Weare,
//* Geochimica et Cosmochimica Acta Vol. 56, 1992, pp. 2605-2617.
//* 
//* [3] International equations for the saturation properties of ordinary water substance.
//* Revised according to the International Temperature Scale of 1990.
//* Addendum to J. Phys. Chem. Ref. Data 16, 893 (1987),
//* W. Wagner and A. Pruss, 
//* J, Phys. Chem. Ref. Data, Vol. 22, No. 3, 1993, pp.783-787.
//* 
//* Authors: Chris Green and Paulus Lahur


// Comment this out when integrating with MOOSE.
// Uncomment to use this file by itself
// #define STANDALONE

//#define VERBOSE

#ifdef STANDALONE

#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <vector>
#include <iomanip>
using std::cout;
using std::fixed;
using std::vector;
typedef double Real;

#else // MOOSE

#include "PorousFlowBrineMethane.h"
#include "BrineFluidProperties.h"
#include "SinglePhaseFluidPropertiesPT.h"

#endif // MOOSE


// Universal gas constant
// The unit is (pressure).litre/mol.Kelvin, where the only difference is pressure
const Real R_bar = 0.083144598;  // bar.l/(mol.K)
const Real R_kPa = 8.3144598;    // kPa.l/(mol.K) = joule/(mol.K)
const Real R_atm = 0.082057338;  // atm.l/(mol.K)

struct molecule_s {
  const Real molarMass;  // kg/mol
  const Real TCCelsius;  // Critical temperature (Celsius)
  const Real PCBar;      // Critical pressure (bar)
  const Real a[15];      // EOS parameters ([Ref. 1], Table 3)

  Real criticalTemperatureCelsius() const {
    return TCCelsius;
  }

  Real criticalTemperatureKelvin() const {
    return TCCelsius + 273.15;
  }

  Real criticalPressureBar() const {
    return PCBar;
  }

  Real criticalPressurePascal() const {
    return PCBar * 1.0e5;
  }
};

struct bracket_s {
  Real min;
  Real max;
};

molecule_s CH4 = {
  16.04e-3,  // molar mass
  -82.55,    // Tc
  46.41,     // Pc
  {          // EOS parameters
    8.72553928e-02,
    -7.52599476e-01,
    3.75419887e-01,
    1.07291342e-02,
    5.49626360e-03,
    -1.84772802e-02,
    3.18993183e-04,
    2.11079375e-04,
    2.01682801e-05,
    -1.65606189e-05,
    1.19614546e-04,
    -1.08087289e-04,
    4.48262295e-02,
    7.5397e-01,
    7.7167e-02
  }
};
// Note: triple point:
// Tt = -182.48;   // Celsius
// Pt = 0.117; // bar

molecule_s CO2 = {
  44.01e-3,  // molar mass
  31.05,     // Tc
  73.825,    // Pc
  {          // EOS parameters
    8.99288497E-02,
    -4.94783127E-01,
    4.77922245E-02,
    1.03808883E-02,
    -2.82516861E-02,
    9.49887563E-02,
    5.20600880E-04,
    -2.93540971E-04,
    -1.77265112E-03,
    -2.51101973E-05,
    8.93353441E-05,
    7.88998563E-05,
    -1.66727022E-02,
    1.39800000E+00,
    2.96000000E-02
  }
};
// Note: triple point:
// Tt = -56.15;   // Celsius
// Pt = 5.1; // bar

molecule_s H2O = {
  18.01528e-3, // molar mass
  374.1, //373.946, //374.1,       // Tc
  221.19, //220.64, //221.19,      // Pc
  {            // EOS parameters
    8.64449220E-02,
    -3.96918955E-01,
    -5.73334886E-02,
    -2.93893000E-04,
    -4.15775512E-03,
    1.99496791E-02,
    1.18901426E-04,
    1.55212063E-04,
    -1.06855859E-04,
    -4.93197687E-06,
    -2.73739155E-06,
    2.65571238E-06,
    8.96079018E-03,
    4.02000000E+00,
    2.57000000E-02
  }
};
// Note: triple point:
// Tt = 0.01;       // Celsius
// Pt = 0.00611657; // bar


molecule_s NaCl = {
  58.4e-3,
  0.0,
  0.0,
  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}
};


Real 
pressureBarToPascal (const Real PBar) {
  return PBar * 1.0e5;
}

Real 
pressurePascalToBar (const Real PPascal) {
  return PPascal * 1.0e-5;
}

Real 
temperatureCelsiusToKelvin (const Real tCelsius) {
  return tCelsius + 273.15;
}

Real 
temperatureKelvinToCelsius (const Real TKelvin) {
  return TKelvin - 273.15;
}


#ifdef STANDALONE
#else // MOOSE
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
#endif


#ifdef STANDALONE
#else // MOOSE
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
#endif


#ifdef STANDALONE
#else // MOOSE
std::string
PorousFlowBrineMethane::fluidStateName() const
{
  return "brine-methane";
}
#endif


#ifdef STANDALONE
#else // MOOSE
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
#endif


#ifdef STANDALONE
#else // MOOSE
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
#endif


#ifdef STANDALONE
#else // MOOSE
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
#endif


#ifdef STANDALONE
#else // MOOSE
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
#endif


#ifdef STANDALONE
#else // MOOSE
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
#endif // MOOSE


Real 
waterVapourReducedPressureAtSaturation(const Real Tr) {
  // Compute water vapour pressure at saturation.
  // Ref 3, Eq. 1
  // Input : Reduced temperature
  // Output: Reduced pressure (of water vapour at saturation)

  if (Tr < 0.0) {
#ifdef VERBOSE
    cout << "ERROR in waterVapourReducedPressureAtSaturation: reduced temperature is too low: "
        << Tr << "\n";
#endif
    return 0.0;
  }
  if (Tr > 1.0) {
#ifdef VERBOSE
    cout << "ERROR in waterVapourReducedPressureAtSaturation: reduced temperature is too high: "
        << Tr << "\n";
#endif
    return 0.0;
  }
  
  const Real a1 = -7.85951783;
  const Real a2 = 1.84408259;
  const Real a3 = -11.7866497;
  const Real a4 = 22.6807411;
  const Real a5 = -15.9618719;
  const Real a6 = 1.80122502;

  Real tau = 1.0 - Tr;
  return (exp((a1*tau + a2*pow(tau,1.5) + a3*pow(tau,3.0) + a4*pow(tau,3.5) + 
      a5*pow(tau,4.0) + a6*pow(tau,7.5)) / Tr));
}


Real
waterVapourPressureAtSaturation(const Real TKelvin) {
  // Compute water vapour pressure at saturation.
  // Input: Temperature in Kelvin
  // Output: Pressure in bar
  Real Tr = TKelvin / H2O.criticalTemperatureKelvin();
  return (waterVapourReducedPressureAtSaturation(Tr) * H2O.criticalPressureBar());
}


Real 
waterliquidReducedDensityAtSaturation(const Real Tr) {
  // Compute water liquid density at saturation.
  // Ref 3, Eq. 2
  // Input : Reduced temperature
  // Output: Reduced density
  // Note that when output is negative, water is not liquid

  if (Tr < 0.0) {
#ifdef VERBOSE
    cout << "ERROR in waterliquidReducedDensityAtSaturation: reduced temperature is too low: "
        << Tr << "\n";
#endif
    return 0.0;
  }
  if (Tr > 1.0) {
#ifdef VERBOSE
    cout << "ERROR in waterliquidReducedDensityAtSaturation: reduced temperature is too high: "
        << Tr << "\n";
#endif
    return 0.0;
  }
  
  const Real a1 = 1.99274064;
  const Real a2 = 1.09965342;
  const Real a3 = -0.510839303;
  const Real a4 = -1.75493479;
  const Real a5 = -45.5170352;
  const Real a6 = -6.74694450e5;

  Real tau = 1.0 - Tr;
  return (1.0 + a1*pow(tau,1.0/3.0) + a2*pow(tau,2.0/3.0) + a3*pow(tau,5.0/3.0) + 
      a4*pow(tau,16.0/3.0) + a5*pow(tau,43.0/3.0) + a6*pow(tau,110.0/3.0));
}


void
EOSTempDependentCoefs(
  const molecule_s m, const Real Tr, 
  Real & B, Real & C, Real & D, Real & E, Real & F)
{
  // Compute temperature-dependent coefficients of EOS of the given molecule.
  // Input is reduced temperature.

  // To avoid expensive recalculation
  Real inv_Tr_squared = 1.0 / (Tr*Tr);
  Real inv_Tr_cubed = 1.0 / pow(Tr,3.0);

  B = m.a[0] + m.a[1]*inv_Tr_squared + m.a[2]*inv_Tr_cubed;
  C = m.a[3] + m.a[4]*inv_Tr_squared + m.a[5]*inv_Tr_cubed;
  D = m.a[6] + m.a[7]*inv_Tr_squared + m.a[8]*inv_Tr_cubed;
  E = m.a[9] + m.a[10]*inv_Tr_squared + m.a[11]*inv_Tr_cubed;
  F = m.a[12] * inv_Tr_cubed;
}


Real 
EOSBalance(
  const molecule_s m, const Real Pr, const Real Tr, const Real Vr)
{
  // EOS for methane, rearranged so that the output is zero when it is perfectly balanced.
  // All input are reduced properties.

  // To avoid expensive recalculation
  Real inv_Vr_squared = 1.0 / (Vr*Vr);

  Real B, C, D, E, F;
  EOSTempDependentCoefs(m, Tr, B, C, D, E, F);

  return(1.0 + B/Vr + C*inv_Vr_squared + D/pow(Vr,4.0) + E/pow(Vr,5.0) +
      F * inv_Vr_squared * (m.a[13]+m.a[14]*inv_Vr_squared) * exp(-m.a[14]*inv_Vr_squared) -
      Pr * Vr / Tr);
}


Real
solveMolarVolume(const molecule_s m, const Real Pr, const Real Tr)
{  
  const Real tolerance = 1.0e-6;
  const int maxIteration = 100;

  Real balance, balance0, balance1;
  std::vector<bracket_s> brackets;

  // Find the bracket that contains the smallest root
  // cout << "Searching for bracket ...\n";
  Real Vrguess = Tr / Pr;
  Real Vrmin = 0.0001 * Vrguess;
  Real Vrmax = 100.0 * Vrguess;
  Real increment = 0.0005 * (Vrmax-Vrmin);  // molar volume increment
  bracket_s bracket;
  bracket.min = Vrmin;
  balance0 = EOSBalance(m, Pr, Tr, bracket.min);
  bracket.max = bracket.min + increment;
  bool bracketFound = false;
  while (bracket.max < Vrmax) {
    balance1 = EOSBalance(m, Pr, Tr, bracket.max);
    if (balance0 * balance1 <= 0.0) {
      bracketFound = true;
      brackets.push_back(bracket);
    }
    bracket.min = bracket.max;
    balance0 = balance1;
    bracket.max = bracket.min + increment;
  }
  if (!bracketFound) {
    // cout << "ERROR in solveMolarVolume: bracket is not found\n";
    return -1.0;
  }
  // cout << "Bracket count: " << brackets.size() << "\n";

  // cout << "Searching for solution inside bracket ...\n";
  Real Vr, Vr0, Vr1;
  std::vector<Real> Vrlist;
  // cout << "Root count: " << brackets.size() << "\n";
  // cout << "Root(s): ";
  for (int i = 0; i < brackets.size(); ++i) {
    Vr0 = brackets[i].min;
    Vr1 = brackets[i].max;
    balance0 = EOSBalance(m, Pr, Tr, Vr0);
    balance1 = EOSBalance(m, Pr, Tr, Vr1);
    // cout << "Balances: " << balance0 << ", " << balance1 << "\n";
    if (balance0 * balance1 <= 0.0) {
      // There must be at least 1 root inside this bracket
      for (int j = 0; j < maxIteration; ++j) {
        Vr = (Vr0+Vr1) * 0.5;
        balance = EOSBalance(m, Pr, Tr, Vr);
        if (abs(Vr0-Vr1) < tolerance) {
#ifdef VERBOSE
          cout << "solveMolarVolume: balance achieved at " << i << " iterations" << "\n";
          cout << Vr ;
#endif
          // cout << Vr << ", ";
          Vrlist.push_back(Vr);
          break;
        } else {
          if (balance * balance0 > 0.0) {
            Vr0 = Vr;
          } else {
            Vr1 = Vr;
          }
        }
      } // next iteration
    } else {
#ifdef VERBOSE
      cout << "ERROR in solveMolarVolume: illegal bracket\n";
#endif
      return -1.0;
    }
  } // next bracket
  // cout << "\n";

  if (Vrlist.size() == 1) {
    return Vrlist[0];
  } else {
    // cout << "Multiple roots: " << Vrlist.size() << "\n";
    // Multiple solutions. Resolve this!
    // Pick the molar volume that is the furthest from the guessed value
    Real difference;
    Real maxDifference = -9.9e99;
    Vr = -1.0;
    // cout << "Vrlist: ";
    for (int i = 0; i < Vrlist.size(); ++i) {
      // cout << Vrlist[i] << ", ";
      difference = Vrlist[i] - Vrguess;
      if (difference > maxDifference) {
        Vr = Vrlist[i];
        maxDifference = difference;
      }
    }
    // cout << "Vrguess: " << Vrguess << "\n";
    return Vr;
  }
}


Real
pureFugacityCoefficient(
    const molecule_s m, const Real PPascal, const Real TKelvin)
{
#ifdef VERBOSE
  cout << "pureFugacityCoefficient\n";
#endif
  // Compute fugacity coefficient for CH4, CO2 and H20 in a generic way.
  // Ref. 1, Eq. A4, which is from Ref. 2 Eq. 1-6
  // Input: pressure (Pascal), temperature (Kelvin)
  // Output: fugacity
  // WARNING: 
  // For the present purpose, we care only for fugacity of CH4, 
  // where the results reproduce the values in Table 4 quite closely.
  // On the other hand, results for H2O and CO2 are a bit off in various places,
  // due to less-than-optimum root-finding routine.
  // TODO: Find a better way for root finding (ie. based in physics)

  Real PBar = pressurePascalToBar(PPascal);
  // Check valid range:
  // Note that the pressure range indicated by the title of Ref. 1 is 0~8000 bar,
  // but Table 4 shows 1~10000.
  // The range here takes into account small perturbation
  if ((PBar < 1.0) || (PBar > 10000.0001)) {
#ifdef VERBOSE
    cout << "ERROR in pureFugacityCoefficient: Illegal pressure range";
#endif
    return 0.0;
  }
  // Note that the temperature range indicated by the title of Ref. 1 is 0~1000 deg Celsius,
  // but Table 4 shows 0~1200.
  // The range here takes into account small perturbation
  if ((TKelvin < 273.15) || (TKelvin > 1474.15)) {
#ifdef VERBOSE
    cout << "ERROR in pureFugacityCoefficient: Illegal temperature range";
#endif
    return 0.0;
  }
#ifdef VERBOSE
  cout << "PBar: " << PBar << ", temp: " << temperature << "\n";
#endif

  // Critical pressure and temperature.
  // Note that Vc is the molar volume at critical pressure and temperature, 
  // and not critical molar volume.
  const Real Pc = m.criticalPressureBar();
  const Real Tc = m.criticalTemperatureKelvin();
  const Real Vc = R_bar * Tc / Pc;    // Molar volume in l/mol

  // Reduced properties
  Real Pr = PBar / Pc;
  Real Tr = TKelvin / Tc;
  Real Vr = solveMolarVolume(m, Pr, Tr);
  // cout << "<" << Vr/(Tr/Pr) << ">";
  if (Vr <= 0.0) { // error
#ifdef VERBOSE
    cout << "ERROR in pureFugacityCoefficient: illegal Vr: " << Vr << "\n";
#endif
    return 0.0;
  }
#ifdef VERBOSE
  cout << "Pr: " << Pr << ", Vr: " << Vr << ", Tr: " << Tr << "\n";
#endif
  Real Z = Pr * Vr / Tr;

  Real B, C, D, E, F;
  EOSTempDependentCoefs(m, Tr, B, C, D, E, F);
  Real a15_by_square_Vr = m.a[14] / (Vr * Vr);
  Real G = 0.5 * F / m.a[14] * 
      (m.a[13] + 1.0 - (m.a[13] + 1.0 + a15_by_square_Vr) * exp(-a15_by_square_Vr));

  // Fugacity coefficient
  return (exp(Z - 1.0 - log(Z) + 
      B/Vr + C/(2.0*Vr*Vr) + D/(4.0*pow(Vr,4.0)) + E/(5.0*pow(Vr,5.0)) + G));
}


Real
fugacityCoefficientH2O(
    const Real PPascal, const Real TKelvin)
{
  // Fugacity coefficient in the gas phase for CH4-H2O system.
  // Ref. 1, Eq. 6
  // Input: pressure in Pascal, temperature in Kelvin

  // Check valid input range
  if ((TKelvin < 273.0) || (TKelvin > 574.16)) {
#ifdef VERBOSE
    cout << "ERROR in fugacityCoefficientH2O: temperature is beyond valid range\n";
#endif
    return 0.0;
  }
  Real PBar = pressurePascalToBar(PPascal);
  if ((PBar < 0.999) || (PBar > 2000.01)) {
#ifdef VERBOSE
    cout << "ERROR in fugacityCoefficientH2O: pressure is beyond valid range\n";
#endif
    return 0.0;
  }

  // Parameters for eq. 6
  const Real a1 = -1.42006707e-2;
  const Real a2 = 1.08369910e-2;
  const Real a3 = -1.59213160e-6;
  const Real a4 = -1.10804676e-5;
  const Real a5 = -3.14287155;
  const Real a6 = 1.06338095e-3;

  // Fugacity coefficient
  return (exp(a1 + a2*PBar + a3*PBar*PBar + a4*PBar*TKelvin + a5*PBar/TKelvin + 
      a6*PBar*PBar/TKelvin));
}


Real
molFractionOfWaterInGas(const Real PPascal, const Real TKelvin, const Real Xh2o)
{
  // Compute mass fraction of water in gas, as in Ref. 1, Eq 5
  // Input: pressure in Pascal, temperature in Kelvin, mol fraction of H2O in liquid.
  // Output: mol fraction of water in gas phase.

  const Real criticalDensityH2O = 322.0;  // kg/m^3

  // Reduced temperature
  Real Tr = TKelvin / H2O.criticalTemperatureKelvin();
  // Reduced pressure of water vapour at saturation
  Real Psh2oBar = waterVapourReducedPressureAtSaturation(Tr) * H2O.criticalPressureBar();
#ifdef VERBOSE
  cout << "Water vapour reduced pressure at saturation: " << Psh2oBar << "\n";
#endif
  Real fh2o = fugacityCoefficientH2O(PPascal, TKelvin);
#ifdef VERBOSE
  cout << "Fugacity coefficient: " << fh2o << "\n";
#endif
  // Density of water
  Real densityh2o = waterliquidReducedDensityAtSaturation(Tr) * criticalDensityH2O;  // kg/m^3
#ifdef VERBOSE
  cout << "water density: " << densityh2o << "\n";
#endif
  // Molar volume of liquid water (dm^3/mole)
  Real Vlh2o = H2O.molarMass / densityh2o * 1.0e3;  // dm^3/mole
  Real PBar = pressurePascalToBar(PPascal);

  Real Yh2o = Xh2o * Psh2oBar / (fh2o * PBar) * exp(Vlh2o * (PBar - Psh2oBar) / (R_bar * TKelvin));
  return (std::min(Yh2o, 1.0));  // Filter out physically impossible fraction
}


Real 
genericEquation(const Real c[], const Real PBar, const Real TKelvin) 
{
  // Equation 9 of Reference 1
  // This is the generic form of the following equations:
  // - Chemical potential of standard liquid methane, divided by R and temperature.
  // - Lambda of CH4-Na
  // - Xi of CH4-Na-C
  // Input: a set of parameters specific to the equation above, pressure & temperature
  // Output: depends on specific equation

  if (PBar <= 0.0) {
#ifdef VERBOSE
    cout << "ERROR in genericEquation: illegal value of PBar: ", PBar, "\n";
#endif
    return 0.0;
  }
  if (TKelvin <= 0.0) {
#ifdef VERBOSE
    cout << "ERROR in genericEquation: illegal value of TKelvin: ", TKelvin, "\n";
#endif
    return 0.0;
  }
  // cout << "(" << PBar << ", " << TKelvin << ")";
  // cout << "(" << c[0] << ")";
  Real squareT = TKelvin * TKelvin;
  return (c[0] + c[1]*TKelvin + c[2]/TKelvin + c[3]*squareT + c[4]/squareT + 
      PBar * (c[5] + c[6]*TKelvin + c[7]/TKelvin + c[8]/squareT + c[9]*PBar*TKelvin));
}


Real
muCH4StandardLiquidByRT(const Real PBar, const Real TKelvin)
{
  // Chemical potential of standard liquid methane, divided by R and temperature.
  // Table 3 of Reference 1
  const Real c[10] = 
  {
    0.83143711e1,
    -0.72772168e-3,
    0.21489858e4,
    -0.14019672e-4,
    -0.66743449e6,
    0.76985890e-2,
    -0.50253331e-5,
    -0.30092013e1,
    0.48468502e3,
    0.0
  };
  // cout << "(" << PBar << ", " << TKelvin << ")";
  return (genericEquation(c, PBar, TKelvin));
}


Real 
lambdaCH4_Na(const Real PBar, const Real TKelvin)
{
  // Table 3 of Reference 1
  const Real c[10] =
  {
    -0.81222036,
    0.10635172e-2,
    0.18894036e3,
    0.0,
    0.0,
    0.44105635e-4,
    0.0,
    0.0,
    0.0,
    -0.46797718e-10
  };
  return (genericEquation(c, PBar, TKelvin));
}


Real
xiCH4_Na_Cl(const Real PBar, const Real TKelvin)
{
  // Table 3 of Reference 1
  // Actually, the result is just constant wrt pressure & temperature.
  const Real c[10] = 
  {
    -0.29903571e-2,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0
  };
  return (genericEquation(c, PBar, TKelvin));
}


#ifdef STANDALONE
Real
activityCoef(const Real pressure, const Real temperature, const Real Xnacl)
#else
Real
activityCoef(Real pressure, Real temperature, Real Xnacl)
#endif
{
  // Compute activity coefficient of methane in brine.
  // Ref. 1 Eq. 7, with the term for CH4-Cl set to zero.
  // Input: pressure in Pascal, temperature in Kelvin, mol fraction of NaCl
  // Output: activity coefficient

  Real PBar = pressurePascalToBar(pressure);
  // Molality of NaCl in liquid 
  // Note the constant 2 is due to NaCl dissociation
  Real mnacl = Xnacl / (2.0 * H2O.molarMass * (1.0 - Xnacl));

  return std::exp(2.0 * lambdaCH4_Na(PBar, temperature) * mnacl
      + xiCH4_Na_Cl(PBar, temperature) * mnacl * mnacl);
}


#ifdef STANDALONE
void
activityCoefficient(const Real pressure, const Real temperature, const Real xnacl,
    Real & gamma, Real & dgamma_dp, Real & dgamma_dT)
#else // MOOSE
void
PorousFlowBrineMethane::activityCoefficient(Real pressure,
                                            Real temperature,
                                            Real xnacl,
                                            Real & gamma,
                                            Real & dgamma_dp,
                                            Real & dgamma_dT) const
#endif
{
  gamma = activityCoef(pressure, temperature, xnacl);
  dgamma_dp = activityCoef(pressure + 1.0, temperature, xnacl) - gamma;
  dgamma_dT = activityCoef(pressure, temperature + 1.0, xnacl) - gamma;
}


#ifdef STANDALONE
Real
methaneSolubilityInLiquid(const Real pressure, const Real temperature, const Real Xnacl)
#else // MOOSE
Real
PorousFlowBrineMethane::methaneSolubilityInLiquid(
    Real pressure, Real temperature, Real Xnacl) const
#endif
{
  // Compute methane solubility in liquid, based on Ref. 1, Eq. 8.
  // Input: pressure in Pascal, temperature in Kelvin, mol fraction of NaCl
  // Output: methane solubility in the liquid (mol/kg).

  Real PBar = pressurePascalToBar(pressure);
  // Approximate mol fraction of water in liquid for CH4-H2O-NaCl system.
  // Note:
  // - For the case of CH4-H2O system (ie. no NaCl), it reduces to just 1.0
  // - The fraction of NaCl is multiplied by 2, because NaCl dissociates into Na+ and Cl-
  Real Xh2o = 1.0 - 2.0*Xnacl;  
  // Mol fraction of water in gas: Yh2o given by Ref. 1, Eq 5
  Real Yh2o = molFractionOfWaterInGas(pressure, temperature, Xh2o);  
  // Mol fraction of CH4 in gas
  Real Ych4 = 1.0 - Yh2o;
  // Solubility of CH4 in liquid (mol/kg)
  Real RHS = muCH4StandardLiquidByRT(PBar, temperature) 
      - log(pureFugacityCoefficient(CH4, pressure, temperature))
      + log(activityCoef(pressure, temperature, Xnacl));
  Real mch4 = std::max(0.0, Ych4 * PBar / exp(RHS));  // Filter out negative result

  return mch4;
}


#ifdef STANDALONE
void
equilibriumMassFractions(const Real pressure, const Real temperature, const Real Xnacl,
    Real & wch4, Real & wh2o)
#else
void
PorousFlowBrineMethane::equilibriumMassFractions(Real pressure, Real temperature, Real Xnacl,
    Real & wch4, Real & wh2o) const
#endif
{
  // Compute mol fraction of CH4 in liquid phase and H2O in gas phase.
  // Ref. 1, Eq. 8 (derived from Eq. 3 & 7).
  // Input: pressure in Pascal, temperature in Kelvin, mol fraction of NaCl
  // Output: mol fraction of CH4 in liquid phase and H2O in gas phase.

  Real PBar = pressurePascalToBar(pressure);
  // Approximate mol fraction of water in liquid for CH4-H2O-NaCl system.
  // Note that for the case of CH4-H2O system (ie. no NaCl), it reduces to just 1.0.
  Real Xh2o = 1.0 - 2.0*Xnacl;  
  // Mol fraction of water in gas: Yh2o given by Ref. 1, Eq 5
  Real Yh2o = molFractionOfWaterInGas(pressure, temperature, Xh2o);
  // Mol fraction of CH4 in gas
  Real Ych4 = 1.0 - Yh2o;
  // Total molar mass of gas: CH4, H2O (gas)
  Real totalGasMolarMass = Ych4 * CH4.molarMass + Yh2o * H2O.molarMass;
  // Mass fraction of water in gas
  wh2o = Yh2o * H2O.molarMass / totalGasMolarMass;

  // Solubility of CH4 in liquid (mol/kg)
  // Real RHS = muCH4StandardLiquidByRT(PBar, temperature) 
  //     - log(pureFugacityCoefficient(CH4, pressure, temperature))
  //     + log(activityCoef(pressure, temperature, Xnacl));
  // Real mch4 = std::max(0.0, Ych4 * PBar/ exp(RHS));  // Filter out negative result
  Real mch4 = methaneSolubilityInLiquid(pressure, temperature, Xnacl);

  Real massh2o = Xh2o * H2O.molarMass;
  Real massnacl = Xnacl * NaCl.molarMass;
  Real totalMass = massh2o + massnacl;  // Assume mass of CH4 in liquid is negligible

  // Mol fraction of CH4 in liquid
  Real Xch4 = mch4 * totalMass;  
  // Total molar mass of liquid: H2O, NaCl, CH4 (gas)
  Real totalLiquidMolarMass = Xh2o * H2O.molarMass + Xnacl * NaCl.molarMass + Xch4 * CH4.molarMass;
  // Mass fraction of CH4 in liquid
  wch4 = Xch4 * CH4.molarMass / totalLiquidMolarMass;
}


#ifdef STANDALONE
void
equilibriumMassFractions(
    const Real pressure, const Real temperature, const Real Xnacl,
    Real & Xch4, Real & dXch4_dp, Real & dXch4_dT,
    Real & Yh2o, Real & dYh2o_dp, Real & dYh2o_dT)

#else // MOOSE
void
PorousFlowBrineMethane::equilibriumMassFractions(Real pressure,
                                                 Real temperature,
                                                 Real Xnacl,
                                                 Real & Xch4,
                                                 Real & dXch4_dp,
                                                 Real & dXch4_dT,
                                                 Real & Yh2o,
                                                 Real & dYh2o_dp,
                                                 Real & dYh2o_dT) const
#endif // MOOSE
{
  // Mass fraction of water in gas: Yh2o given by Eq 5 from Duan and Mao
  // Mass fraction of methane in liquid: Xch4 calculated from Eq 3 by solving
  // for mch4 (moles/kg of CH4 in liquid) and dividing by _Mch4 (molar mass of methane)
  // Perturbed values:
  Real Xch4p, Yh2op;

  equilibriumMassFractions(pressure, temperature, Xnacl, Xch4, Yh2o);

  equilibriumMassFractions(pressure + 1.0, temperature, Xnacl, Xch4p, Yh2op);
  dXch4_dp = Xch4p - Xch4;
  dYh2o_dp = Yh2op - Yh2o;

  equilibriumMassFractions(pressure, temperature + 1.0, Xnacl, Xch4p, Yh2op);
  dXch4_dT = Xch4p - Xch4;
  dYh2o_dT = Yh2op - Yh2o;
}


#ifdef STANDALONE
void
equilibriumMolFractions(const Real pressure, const Real temperature, const Real Xnacl,
    Real & Xch4, Real & Yh2o)
{
  // Compute mol fraction of CH4 in liquid phase and H2O in gas phase.
  // Ref. 1, Eq. 8 (derived from Eq. 3 & 7).
  // Input: pressure in Pascal, temperature in Kelvin, mol fraction of NaCl
  // Output: mol fraction of CH4 in liquid phase and H2O in gas phase.

  Real PBar = pressurePascalToBar(pressure);
  // Approximate mol fraction of water in liquid for CH4-H2O-NaCl system.
  // Note that for the case of CH4-H2O system (ie. no NaCl), it reduces to just 1.0.
  Real Xh2o = 1.0 - 2.0*Xnacl;  
  // Mol fraction of water in gas: Yh2o given by Ref. 1, Eq 5
  Yh2o = molFractionOfWaterInGas(pressure, temperature, Xh2o);

  // const Real molarMassNa = 22.9898e-3;  // kg/mol
  // const Real molarMassCl = 35.453e-3;
  
  // Mol fraction of CH4 in gas
  Real Ych4 = 1.0 - Yh2o;
  // Molality of NaCl in liquid
  // Real mnacl = Xnacl / (1.0 - Xnacl) / H2O.molarMass;

  // Solubility of CH4 in liquid (mol/kg)
  // Real RHS = muCH4StandardLiquidByRT(PBar, temperature) 
  //     - log(pureFugacityCoefficient(CH4, pressure, temperature))
  //     + log(activityCoef(pressure, temperature, Xnacl));
  // Real mch4 = std::max(0.0, Ych4 * PBar/ exp(RHS));  // Filter out negative result
  Real mch4 = methaneSolubilityInLiquid(pressure, temperature, Xnacl);

  Real massh2o = Xh2o * H2O.molarMass;
  Real massnacl = Xnacl * NaCl.molarMass;
  Real totalMass = massh2o + massnacl;  // Assume mass of CH4 in liquid is negligible

  Xch4 = mch4 * totalMass;  // mol fraction of CH4 in liquid
}
#endif


#ifdef STANDALONE
void
equilibriumMolFractionsAndDerivatives(
    const Real PPascal, const Real TKelvin, const Real Xnacl,
    Real & Xch4, Real & dXch4_dp, Real & dXch4_dT,
    Real & Yh2o, Real & dYh2o_dp, Real & dYh2o_dT)
{
  // Perturbed values:
  Real Xch4p, Yh2op;

  equilibriumMolFractions(PPascal, TKelvin, Xnacl, Xch4, Yh2o);

  equilibriumMolFractions(PPascal + 1.0, TKelvin, Xnacl, Xch4p, Yh2op);
  dXch4_dp = Xch4p - Xch4;
  dYh2o_dp = Yh2op - Yh2o;

  equilibriumMolFractions(PPascal, TKelvin + 1.0, Xnacl, Xch4p, Yh2op);
  dXch4_dT = Xch4p - Xch4;
  dYh2o_dT = Yh2op - Yh2o;
}
#endif


#ifdef STANDALONE
void
fugacityCoefficientMethane(
    const Real pressure, const Real temperature, Real & fch4, Real & dfch4_dp, Real & dfch4_dT)
#else
void
PorousFlowBrineMethane::fugacityCoefficientMethane(
    Real pressure, Real temperature, Real & fch4, Real & dfch4_dp, Real & dfch4_dT) const
#endif
{
  // Ref. 1, Eq. A4, which is from Ref. 2 Eq. 1-6
  // Input: pressure (Pascal), temperature (Kelvin)
  // Output: fugacity, its derivative wrt pressure & temperature

  fch4 = pureFugacityCoefficient(CH4, pressure, temperature);
  dfch4_dp = pureFugacityCoefficient(CH4, pressure + 1.0, temperature) - fch4;
  dfch4_dT = pureFugacityCoefficient(CH4, pressure, temperature + 1.0) - fch4;
}


#ifdef STANDALONE
void
fugacityCoefficientH2O(
    const Real pressure, const Real temperature, Real & fh2o, Real & dfh2o_dp, Real & dfh2o_dT)
#else // MOOSE
void
PorousFlowBrineMethane::fugacityCoefficientH2O(
    Real pressure, Real temperature, Real & fh2o, Real & dfh2o_dp, Real & dfh2o_dT) const
#endif
{
  // Fugacity coefficient in the gas phase.
  // Ref. 1, Eq. 6
  // Input: pressure in Pascal, temperature in Kelvin
  fh2o = 0.0;
  dfh2o_dp = 0.0;
  dfh2o_dT = 0.0;

  // Check valid input range
  if ((temperature < 273.0) || (temperature > 523.16)) {
#ifdef VERBOSE
    cout << "ERROR in fugacityCoefficientH2O: temperature is beyond valid range\n";
#endif
  }
  Real PBar = pressurePascalToBar(pressure);
  if ((PBar < 0.999) || (PBar > 2000.01)) {
#ifdef VERBOSE
    cout << "ERROR in fugacityCoefficientH2O: pressure is beyond valid range\n";
#endif
  }

  // Parameters for eq. 6
  const Real a1 = -1.42006707e-2;
  const Real a2 = 1.08369910e-2;
  const Real a3 = -1.59213160e-6;
  const Real a4 = -1.10804676e-5;
  const Real a5 = -3.14287155;
  const Real a6 = 1.06338095e-3;

  // Fugacity coefficient
  fh2o = exp(a1 + a2*PBar + a3*PBar*PBar + a4*PBar*temperature + a5*PBar/temperature + 
      a6*PBar*PBar/temperature);

  // The derivative of fugacity coefficient wrt pressure
  dfh2o_dp = fh2o * (a2 + 2.0*a3*PBar + a4*temperature + a5/temperature + 2.0*a6*PBar/temperature);

  // The derivative of fugacity coefficient wrt temperature
  Real inv_temp_squared = 1.0 / (temperature*temperature);
  dfh2o_dT = fh2o * (a4*PBar - a5*PBar*inv_temp_squared - a6*PBar*PBar*inv_temp_squared);
}


#ifdef STANDALONE
#else // MOOSE
void
PorousFlowBrineMethane::equilibriumConstantH2O(Real temperature, Real & kh2o, Real & dkh2o_dT) const
{
  // Might not be needed?
}
#endif


#ifdef STANDALONE
#else // MOOSE
void
PorousFlowBrineMethane::equilibriumConstantMethane(Real temperature,
                                                   Real & kch4,
                                                   Real & dkch4_dT) const
{
  // Might not be needed?
}
#endif


#ifdef STANDALONE
#else // MOOSE
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
#endif


#ifdef STANDALONE
#else // MOOSE
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
#endif


#ifdef STANDALONE

//------------------------------------------------------------------------------------
// Test all functions in stand-alone mode

int
test_waterliquidReducedDensityAtSaturation() 
{
  cout << "test_waterliquidReducedDensityAtSaturation: \n";

  const Real criticalDensity = 322.0;
  Real Tr;
  Real reducedDensity, density;
  cout << "Tr range from 0.0 to 1.0\n";
  cout << "Reduced temperature (Temperature in Kelvin): Reduced density (Density in kg/m^3)\n";
  for (int i = 0; i < 11; ++i) {
    Real Tr = 0.1 * i;
    reducedDensity = waterliquidReducedDensityAtSaturation(Tr);
    density = reducedDensity * criticalDensity;
    cout << Tr << " (" << Tr * H2O.criticalTemperatureKelvin() << "): " << 
        reducedDensity << " (" << density << ")\n";
  }
  cout << "Temperature range from freezing point to boiling point\n";
  // cout << "At freezing point (0 degree Celsius):\n";
  cout << "Temperature in Celsius (Reduced temperature): Reduced density (Density in kg/m^3)\n";
  for (int i = 0; i < 11; ++i) {
    Real t = 10.0 * i;
    Tr = temperatureCelsiusToKelvin(t) / H2O.criticalTemperatureKelvin();
    reducedDensity = waterliquidReducedDensityAtSaturation(Tr);
    density = reducedDensity * criticalDensity;
    cout << t << " (" << Tr << "): " << 
      reducedDensity << " (" << density << ")\n";
  }
  cout << "\n";
  
  return 0;
}


int
test_waterVapourReducedPressureAtSaturation() 
{
  cout << "test_waterVapourReducedPressureAtSaturation: \n";

  const Real criticalDensity = 322.0;
  Real Tr;
  Real reducedPressure, PBar;
  cout << "Tr range from 0.0 to 1.0\n";
  cout << "Reduced temperature (Temperature in Kelvin): Reduced pressure (Pressure in bar)\n";
  for (int i = 0; i < 11; ++i) {
    Real Tr = 0.1 * i;
    reducedPressure = waterVapourReducedPressureAtSaturation(Tr);
    PBar = reducedPressure * H2O.criticalPressureBar();
    cout << Tr << " (" << Tr * H2O.criticalTemperatureKelvin() << "): " << 
        reducedPressure << " (" << PBar << ")\n";
  }
  cout << "Temperature range from freezing point to boiling point\n";
  // cout << "At freezing point (0 degree Celsius):\n";
  cout << "Temperature in Celsius (Reduced temperature): Reduced pressure (Pressure in bar)\n";
  for (int i = 0; i < 11; ++i) {
    Real t = 10.0 * i;
    Tr = temperatureCelsiusToKelvin(t) / H2O.criticalTemperatureKelvin();
    reducedPressure = waterliquidReducedDensityAtSaturation(Tr);
    PBar = reducedPressure * H2O.criticalPressureBar();
    cout << t << " (" << Tr << "): " << 
      reducedPressure << " (" << PBar << ")\n";
  }
  cout << "\n";

  return 0;
}


int
testSetUpPressure(const int max, vector<Real> & plist)
{
  // Set up the list of input pressure (bar)
  plist.push_back(1.0);
  if (max == 2000) {
    // Range: 1~2000
    plist.push_back(50.0);
    plist.push_back(100.0);
    plist.push_back(150.0);
    for (int i = 2; i < 21; ++i) {
      plist.push_back(i * 100.0);
    }
  } else if (max == 10000) {
    // Range: 1~10000
    plist.push_back(20.0);
    plist.push_back(50.0);
    for (int i = 1; i < 10; ++i) {
      plist.push_back(i*100.0);
    }
    for (int i = 0; i < 20; ++i) {
      plist.push_back(1000.0 + i*200.0);
    }
    for (int i = 0; i < 11; ++i) {
      plist.push_back(5000.0 + i*500.0);
    }
  } else {
    cout << "ERROR in testSetUpPressure: Illegal max value: " << max << "\n";
    return -1;
  }
  // cout << "Pressure entries: " << plist.size() << "\n";
  return 0;
}


int
testSetUpTemperature(const int max, vector<Real> & tlist)
{
  // Set up the list of input temperature (celsius)
  if (max == 300) {
    for (int i = 0; i < 11; ++i) {
      tlist.push_back(i * 30.0);
    }
  } else if (max == 1200) {
    for (int i = 0; i < 7; ++i) {
      tlist.push_back(i * 100.0);
    }
    for (int i = 0; i < 3; ++i) {
      tlist.push_back(800 + i*200.0);
    }
  } else {
    cout << "ERROR in testSetUpTemperature: Illegal max value: " << max << "\n";
    return -1;
  }
  // cout << "Temperature entries: " << tlist.size() << "\n";
  return 0;
}


int
test_pureFugacityCoefficient(const molecule_s & m, const std::string name)
{
  cout << "test_fugacityCoefficient: " << name << "\n";

  vector<Real> plist;
  vector<Real> tlist;
  Real PPascal;
  Real TKelvin;
  testSetUpPressure(10000, plist);
  testSetUpTemperature(1200, tlist);
  Real f[plist.size()][tlist.size()];
  for (int i = 0; i < plist.size(); ++i) {
    for (int j = 0; j < tlist.size(); ++j) {
      PPascal = pressureBarToPascal(plist[i]);
      TKelvin = temperatureCelsiusToKelvin(tlist[j]);
      f[i][j] = pureFugacityCoefficient(m, PPascal, TKelvin);
    }
  }
  cout << "\n";

  cout << "Fugacity coefficient\n";
  cout << std::setprecision(7);
  cout << std::setw(5) << "P";
  for (int i=0; i<tlist.size(); ++i) {
    cout << "," << std::setw(11) << tlist[i];
  }
  cout << "\n";
  for (int i = 0; i < plist.size(); ++i) {
    cout << std::setw(5) << plist[i];
    for (int j = 0; j < tlist.size(); ++j) {
      cout << "," << std::setw(11) << f[i][j]; 
    }
    cout << "\n";
  }
  cout << "\n";
  return 0;
}


int
test_molFractionOfWaterInGas()
{
  cout << "test_molFractionOfWaterInGas:\n";

  // Set up pressure list
  vector<Real> plist;
  plist.push_back(1.0);
  plist.push_back(20.0);
  plist.push_back(50.0);
  for (int i=1; i<10; ++i) {
    plist.push_back(i*100.0);
  }
  for (int i=0; i<6; ++i) {
    plist.push_back(1000.0 + i*200.0);
  }

  // Set up temperature list to match experimental data shown in Ref. 1, Fig. 1
  vector<Real> Tlist;
  Real molalityNaCl;

  for (int k = 0; k < 2; ++k) {
    Tlist.clear();
    if (k == 0) {
      // No NaCl
      molalityNaCl = 0.0;  // mol/kg
      Tlist.push_back(278.0);  // Kelvin
      Tlist.push_back(293.0);
      Tlist.push_back(298.0);
      Tlist.push_back(318.12);
      Tlist.push_back(377.59);
      Tlist.push_back(410.93);
      Tlist.push_back(510.93);
      Tlist.push_back(573.2);
    } else {
      // With NaCl
      molalityNaCl = 3.0;
      Tlist.push_back(298.0);
      Tlist.push_back(511.0);
    }

    Real PPascal;
    Real TKelvin;
    Real Xh2oGas[plist.size()][Tlist.size()];   // Fraction of water in gas phase
    // Molar fraction of NaCl in gas phase
    Real XNaCl = molalityNaCl / (1.0/H2O.molarMass + molalityNaCl);
    cout << "NaCl: molality: " << molalityNaCl << " mol/kg, mol fraction: " << XNaCl << "\n";
    Real Xh2o = 1.0 - 2*XNaCl;  // Fraction of water in liquid
    for (int i = 0; i < plist.size(); ++i) {
      for (int j = 0; j < Tlist.size(); ++j) {
        PPascal = pressureBarToPascal(plist[i]);
        // TKelvin = temperatureCelsiusToKelvin(tlist[j]);
        TKelvin = Tlist[j];
        Xh2oGas[i][j] = molFractionOfWaterInGas(PPascal, TKelvin, Xh2o);
      }
    }

    cout << "Fraction of water in gas for various pressure (bar) and temperature (K)\n";
    cout << std::setprecision(5);
    cout << std::setw(5) << "P";
    for (int i = 0; i < Tlist.size(); ++i) {
      cout << "," << std::setw(11) << Tlist[i];
    }
    cout << "\n";
    for (int i = 0; i < plist.size(); ++i) {
      cout << std::setw(5) << plist[i];
      for (int j = 0; j < Tlist.size(); ++j) {
        cout << "," << std::setw(11) << Xh2oGas[i][j]; 
      }
      cout << "\n";
    }
    cout << "\n";
  }
  return 0;
}


int
test_muCH4StandardLiquidByRT()
{
  cout << "test_muCH4StandardLiquidByRT\n";

  vector<Real> plist;
  testSetUpPressure(2000, plist);
  vector<Real> tlist;
  testSetUpTemperature(300, tlist);
  Real mu;
  cout << std::setprecision(5);
  cout << std::setw(5) << "P";
  for (int i = 0; i < tlist.size(); ++i) {
    cout << "," << std::setw(9) << tlist[i];
  }
  cout << "\n";
  for (int i = 0; i < plist.size(); ++i) {
    cout << std::setw(5) << plist[i];
    for (int j = 0; j < tlist.size(); ++j) {
      mu = muCH4StandardLiquidByRT(plist[i], temperatureCelsiusToKelvin(tlist[j]));
      cout << "," << std::setw(9) << mu;
    }
    cout << "\n";
  }
  cout << "\n";
}


int
test_lambdaCH4_Na()
{
  cout << "test_lambdaCH4_Na\n";

  vector<Real> plist;
  testSetUpPressure(2000, plist);
  vector<Real> tlist;
  testSetUpTemperature(300, tlist);
  Real mu;
  cout << std::setprecision(5);
  cout << std::setw(5) << "P";
  for (int i = 0; i < tlist.size(); ++i) {
    cout << "," << std::setw(9) << tlist[i];
  }
  cout << "\n";
  for (int i = 0; i < plist.size(); ++i) {
    cout << std::setw(5) << plist[i];
    for (int j = 0; j < tlist.size(); ++j) {
      mu = lambdaCH4_Na(plist[i], temperatureCelsiusToKelvin(tlist[j]));
      cout << "," << std::setw(9) << mu;
    }
    cout << "\n";
  }
  cout << "\n";
}


int
test_xiCH4_Na_Cl()
{
  cout << "test_xiCH4_Na_Cl\n";

  vector<Real> plist;
  testSetUpPressure(2000, plist);
  vector<Real> tlist;
  testSetUpTemperature(300, tlist);
  Real mu;
  cout << std::setprecision(5);
  cout << std::setw(5) << "P";
  for (int i = 0; i < tlist.size(); ++i) {
    cout << "," << std::setw(9) << tlist[i];
  }
  cout << "\n";
  for (int i = 0; i < plist.size(); ++i) {
    cout << std::setw(5) << plist[i];
    for (int j = 0; j < tlist.size(); ++j) {
      mu = xiCH4_Na_Cl(plist[i], temperatureCelsiusToKelvin(tlist[j]));
      cout << "," << std::setw(9) << mu;
    }
    cout << "\n";
  }
  cout << "\n";
}


// int test_molFractionOfMethaneInLiquid() 
// {
//   cout << "test_molFractionOfMethaneInLiquid\n";

//   vector<Real> plist;
//   testSetUpPressure(2000, plist);
//   vector<Real> tlist;
//   testSetUpTemperature(300, tlist);
//   Real XNaCl = 0.0;
//   Real YCH4 = 1.0;
//   Real PPascal;
//   Real XCH4;
//   cout << std::setprecision(5);
//   cout << std::setw(5) << "P";
//   for (int i = 0; i < tlist.size(); ++i) {
//     cout << "," << std::setw(9) << tlist[i];
//   }
//   cout << "\n";
//   for (int i = 0; i < plist.size(); ++i) {
//     cout << std::setw(5) << plist[i];
//     for (int j = 0; j < tlist.size(); ++j) {
//       PPascal = pressureBarToPascal(plist[i]);
//       XCH4 = molFractionOfMethaneInLiquid(PPascal, temperatureCelsiusToKelvin(tlist[j]),
//           XNaCl, YCH4);
//       cout << "," << std::setw(9) << XCH4;
//     }
//     cout << "\n";
//   }
//   cout << "\n";
// }


int
test_methaneSolubilityInLiquid()
{
  // The results must match Ref. 1 Table 4-8
  cout << "test_methaneSolubilityInLiquid\n";

  // Set the NaCl concentration to test (mol/kg of SOLVENT)
  vector<Real>mnaclList;
  mnaclList.push_back(0.0);  // pure water
  mnaclList.push_back(1.0);
  mnaclList.push_back(2.0);
  mnaclList.push_back(4.0);
  mnaclList.push_back(6.0);

  vector<Real> plist;
  testSetUpPressure(2000, plist);
  vector<Real> tlist;
  testSetUpTemperature(300, tlist);

  Real Xnacl;
  Real PPascal;
  Real TKelvin;
  Real mch4;

  for (int k = 0; k < mnaclList.size(); ++k) {
    // Compute mol fraction of NaCl
    // Note:
    // - NaCl is expressed as molality (mol/kg of SOLVENT)
    // - The concentration is multiplied by 2 because NaCl dissociate into Na+ and Cl- 
    Xnacl = 2.0 * mnaclList[k] * H2O.molarMass / (1.0 + 2 * mnaclList[k] * H2O.molarMass);
    cout << "NaCl concentration: " << mnaclList[k] << "mol/kg = " << Xnacl << "\n";
    cout << std::setprecision(6);
    cout << std::setw(5) << "P";
    for (int i = 0; i < tlist.size(); ++i) {
      cout << "," << std::setw(9) << tlist[i];
    }
    cout << "\n";
    for (int i = 0; i < plist.size(); ++i) {
      cout << std::setw(5) << plist[i];
      for (int j = 0; j < tlist.size(); ++j) {
        PPascal = pressureBarToPascal(plist[i]);
        TKelvin = temperatureCelsiusToKelvin(tlist[j]);
        mch4 = methaneSolubilityInLiquid(PPascal, TKelvin, Xnacl);
        cout << "," << std::setw(9) << mch4;
      }
      cout << "\n";
    }
    cout << "\n";
  }
}


int
test_methaneSolubilityInLiquid2()
{
  // The results must match Ref. 1 Table 4-8
  cout << "test_methaneSolubilityInLiquid2\n";

  const double maxRelDiff = 0.02;  // max relative difference
  const double maxDiff = 1.0e-7;    // max difference
  
  std::string fileNamePre = "CH4SolubilityInBrine";
  std::string fileNamePost = ".csv";
  std::string fileName;

  // Set the NaCl concentration to test (mol/kg of SOLVENT)
  vector<Real>mnaclList;
  mnaclList.push_back(0.0);  // pure water
  mnaclList.push_back(1.0);
  mnaclList.push_back(2.0);
  mnaclList.push_back(4.0);
  mnaclList.push_back(6.0);

  vector<Real> plist;
  testSetUpPressure(2000, plist);
  vector<Real> tlist;
  testSetUpTemperature(300, tlist);

  Real Xnacl;
  Real PPascal;
  Real TKelvin;
  Real mch4;

  for (int k = 0; k < mnaclList.size(); ++k) {
    // Compute mol fraction of NaCl
    // Note:
    // - NaCl is expressed as molality (mol/kg of SOLVENT)
    // - The concentration is multiplied by 2 because NaCl dissociate into Na+ and Cl- 
    Xnacl = 2.0 * mnaclList[k] * H2O.molarMass / (1.0 + 2 * mnaclList[k] * H2O.molarMass);
    cout << "NaCl concentration: " << mnaclList[k] << "mol/kg = " << Xnacl << "\n";
    cout << std::setprecision(6);
    cout << std::setw(5) << "P";
    for (int i = 0; i < tlist.size(); ++i) {
      cout << "," << std::setw(9) << tlist[i];
    }
    cout << "\n";

    std::ostringstream convert;
    convert << k;
    fileName = fileNamePre + convert.str() + fileNamePost;
    cout << "Reference data file: " << fileName << "\n";
    std::ifstream file(fileName);
    if (!file.good()) {
      cout << "ERROR in test_methaneSolubilityInLiquid2: cannot open file\n";
      return -1;
    }
    std::string line;
    std::getline(file, line);  // Read the header row
    if (!file.good()) {
      cout << "ERROR in test_methaneSolubilityInLiquid2: cannot read the row header\n";
      return -2;
    }
    for (int i = 0; i < plist.size(); ++i) { // For each row
      cout << std::setw(5) << plist[i];
      std::getline(file, line);  // Read a row
      if (!file.good()) {
        cout << "ERROR in test_methaneSolubilityInLiquid2: cannot read row: " << i << "\n";
        return -3;
      }
      std::stringstream iss(line);
      std::string val;
      std::getline(iss, val, ',');  // Read the header column
      if (!iss.good()) {
        cout << "ERROR in test_methaneSolubilityInLiquid2: cannot get the header column\n";
        return -4;
      }

      for (int j = 0; j < tlist.size(); ++j) { // For each column
        PPascal = pressureBarToPascal(plist[i]);
        TKelvin = temperatureCelsiusToKelvin(tlist[j]);
        mch4 = methaneSolubilityInLiquid(PPascal, TKelvin, Xnacl);
        cout << "," << std::setw(9) << mch4;
        std::getline(iss, val, ',');
        // cout << "val: " << val << "\n";
        if (val != "") {
          std::stringstream convertor(val);
          double refData;
          convertor >> refData;
          double diff = abs(mch4 - refData);
          if (diff > 0.0) {
            if (abs(refData) > 0.0) {
              // Compute relative difference
              double relDiff = diff / abs(refData);
              if (relDiff > maxRelDiff) {
                cout << "ERROR in test_methaneSolubilityInLiquid2: value deviates too much from reference\n";
                cout << mch4 << " != " << refData << " (ref)\n";
                cout << "Relative difference: " << relDiff << "\n";
                return -5;
              } // else small deviation
            } else { // reference data is zero. Use absolute difference instead
              if (diff > maxDiff) {
                cout << "ERROR in test_methaneSolubilityInLiquid2: value deviates too much from reference\n";
                cout << mch4 << " != " << refData << " (ref)\n";
                return -6;
              } // else small deviation
            }
          } // else exactly the same
        } // else reference data is empty
      } // next column
      cout << "\n";
    } // next row
    cout << "\n";
  } // next 
  return 0;
}


int
test_equilibriumMolFractions()
{
  cout << "test_equilibriumMolFractions\n";

  // Set the NaCl concentration to test
  vector<Real>mnaclList;
  mnaclList.push_back(0.0);
  // mnaclList.push_back(1.0);
  // mnaclList.push_back(2.0);
  // mnaclList.push_back(4.0);
  mnaclList.push_back(6.0);

  vector<Real> plist;
  testSetUpPressure(2000, plist);
  vector<Real> tlist;
  testSetUpTemperature(300, tlist);

  Real Xnacl;
  Real PPascal;
  Real TKelvin;

  Real Xch4[plist.size()][tlist.size()];
  Real Yh2o[plist.size()][tlist.size()];
  
  // For each value of NaCl concentration
  for (int k = 0; k < mnaclList.size(); ++k) {
    Xnacl = mnaclList[k] / (55.56 + mnaclList[k]);  // Na and Cl
    for (int i = 0; i < plist.size(); ++i) {
      for (int j = 0; j < tlist.size(); ++j) {
        PPascal = pressureBarToPascal(plist[i]);
        TKelvin = temperatureCelsiusToKelvin(tlist[j]);
        equilibriumMolFractions(PPascal, TKelvin, Xnacl, Xch4[i][j], Yh2o[i][j]);
      }
    }

    // Display the results
    cout << "Mol fraction of methane in liquid\n";
    cout << "NaCl concentration: " << mnaclList[k] << "mol/kg = " << Xnacl << "\n";
    cout << std::setprecision(4);
    cout << std::setw(5) << "P";
    for (int i = 0; i < tlist.size(); ++i) {
      cout << "," << std::setw(9) << tlist[i];
    }
    cout << "\n";
    for (int i = 0; i < plist.size(); ++i) {
      cout << std::setw(5) << plist[i];
      for (int j = 0; j < tlist.size(); ++j) {
        cout << "," << std::setw(9) << Xch4[i][j];
      }
      cout << "\n";
    }
    cout << "\n";
    cout << "Mol fraction of water in gas\n";
    cout << "NaCl concentration: " << mnaclList[k] << "mol/kg = " << Xnacl << "\n";
    cout << std::setprecision(4);
    cout << std::setw(5) << "P";
    for (int i = 0; i < tlist.size(); ++i) {
      cout << "," << std::setw(9) << tlist[i];
    }
    cout << "\n";
    for (int i = 0; i < plist.size(); ++i) {
      cout << std::setw(5) << plist[i];
      for (int j = 0; j < tlist.size(); ++j) {
        cout << "," << std::setw(9) << Yh2o[i][j];
      }
      cout << "\n";
    }
    cout << "\n";
  }
}


int
test_equilibriumMassFrac()
{
  cout << "test_equilibriumMassFrac\n";

  // Set the NaCl concentration to test
  vector<Real>mnaclList;
  mnaclList.push_back(0.0);
  // mnaclList.push_back(1.0);
  // mnaclList.push_back(2.0);
  // mnaclList.push_back(4.0);
  mnaclList.push_back(6.0);

  vector<Real> plist;
  testSetUpPressure(2000, plist);
  vector<Real> tlist;
  testSetUpTemperature(300, tlist);

  Real Xnacl;
  Real PPascal;
  Real TKelvin;

  Real Xch4[plist.size()][tlist.size()];
  Real Yh2o[plist.size()][tlist.size()];
  
  // For each value of NaCl concentration
  for (int k = 0; k < mnaclList.size(); ++k) {
    Xnacl = mnaclList[k] / (55.56 + mnaclList[k]);  // Na and Cl
    for (int i = 0; i < plist.size(); ++i) {
      for (int j = 0; j < tlist.size(); ++j) {
        PPascal = pressureBarToPascal(plist[i]);
        TKelvin = temperatureCelsiusToKelvin(tlist[j]);
        equilibriumMassFrac(PPascal, TKelvin, Xnacl, Xch4[i][j], Yh2o[i][j]);
      }
    }

    // Display the results
    cout << "Mass fraction of methane in liquid\n";
    cout << "NaCl concentration: " << mnaclList[k] << "mol/kg = " << Xnacl << "\n";
    cout << std::setprecision(4);
    cout << std::setw(5) << "P";
    for (int i = 0; i < tlist.size(); ++i) {
      cout << "," << std::setw(9) << tlist[i];
    }
    cout << "\n";
    for (int i = 0; i < plist.size(); ++i) {
      cout << std::setw(5) << plist[i];
      for (int j = 0; j < tlist.size(); ++j) {
        cout << "," << std::setw(9) << Xch4[i][j];
      }
      cout << "\n";
    }
    cout << "\n";
    cout << "Mass fraction of water in gas\n";
    cout << "NaCl concentration: " << mnaclList[k] << "mol/kg = " << Xnacl << "\n";
    cout << std::setprecision(4);
    cout << std::setw(5) << "P";
    for (int i = 0; i < tlist.size(); ++i) {
      cout << "," << std::setw(9) << tlist[i];
    }
    cout << "\n";
    for (int i = 0; i < plist.size(); ++i) {
      cout << std::setw(5) << plist[i];
      for (int j = 0; j < tlist.size(); ++j) {
        cout << "," << std::setw(9) << Yh2o[i][j];
      }
      cout << "\n";
    }
    cout << "\n";
  }
}



int main()
{
  // Test all functions in standalone mode.

  int status;

  test_pureFugacityCoefficient(CH4, "CH4");
  //test_pureFugacityCoefficient(H2O, "H2O");
  //test_pureFugacityCoefficient(CO2, "CO2");

  test_waterliquidReducedDensityAtSaturation();

  test_waterVapourReducedPressureAtSaturation();

  test_molFractionOfWaterInGas();
  
  // test_muCH4StandardLiquidByRT();
  // test_lambdaCH4_Na();
  // test_xiCH4_Na_Cl(); // Not really needed (the results are uniform)

  // test_molFractionOfMethaneInLiquid();

  // test_methaneSolubilityInLiquid();
  test_methaneSolubilityInLiquid2();

  // test_equilibriumMassFrac();
  // test_equilibriumMolFractions();

  return 0;
}

#endif // STANDALONE
