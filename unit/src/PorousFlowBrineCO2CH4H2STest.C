//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PorousFlowBrineCO2CH4H2STest.h"
#include "FluidPropertiesTestUtils.h"

/**
 * Test calculation of equilibrium mole fractions over entire temperature range
 */
TEST_F(PorousFlowBrineCO2CH4H2STest, equilibriumMoleFractions)
{
  // Test pure water (Xnacl = 0)
  // Low temperature regime
  Real p = 4.82e6;
  Real T = 310.95;
  Real Xnacl = 0.0;
  Real yco2 = 0.6;
  Real yh2s = 0.1;
  Real ych4 = 0.3;

  Real xco2, xch4, xh2s, yh2o;
  _fp->equilibriumMoleFractions(p, T, Xnacl, yco2, yh2s, ych4, xco2, xch4, xh2s, yh2o);

  ABS_TEST(xco2, 0.0093, 1.0e-3);
  ABS_TEST(xch4, 0.000276, 1.0e-3);
  ABS_TEST(xh2s, 0.00503, 1.0e-3);
  ABS_TEST(yh2o, 0.00191, 1.0e-3);
  // std::cout << xco2 << ", " << xch4 << ", " << xh2s << ", " << yh2o << std::endl;

  p = 16.93e6;
  T = 310.95;
  Xnacl = 0.0;

  yco2 = 0.6;
  yh2s = 0.1;
  ych4 = 0.3;

  _fp->equilibriumMoleFractions(p, T, Xnacl, yco2, yh2s, ych4, xco2, xch4, xh2s, yh2o);

  ABS_TEST(xco2, 0.0154, 1.0e-3);
  ABS_TEST(xch4, 0.00099, 1.0e-4);
  ABS_TEST(xh2s, 0.00608, 1.0e-3);
  ABS_TEST(yh2o, 0.00199, 1.0e-3);
  // std::cout << xco2 << ", " << xch4 << ", " << xh2s << ", " << yh2o << std::endl;
}

TEST_F(PorousFlowBrineCO2CH4H2STest, equilibriumMassFractions)
{
  // Test pure water (Xnacl = 0)
  // Low temperature regime
  Real p = 4.82e6;
  Real T = 310.95;
  Real Xnacl = 0.0;
  Real yco2 = 0.6;
  Real yh2s = 0.1;
  Real ych4 = 0.3;

  Real Xco2, Xch4, Xh2s, Yh2o;
  _fp->equilibriumMassFractions(p, T, Xnacl, yco2, yh2s, ych4, Xco2, Xch4, Xh2s, Yh2o);

  ABS_TEST(Xco2, 0.022064415, 1.0e-3);
  ABS_TEST(Xch4, 0.0002576449, 1.0e-3);
  ABS_TEST(Xh2s, 0.008173358, 1.0e-3);
  // ABS_TEST(Yh2o, 0.00191, 1.0e-3);
  // std::cout << Xco2 << ", " << Xch4 << ", " << Xh2s << std::endl;
}

TEST_F(PorousFlowBrineCO2CH4H2STest, gascompressibilty)
{
  Real p = 6.89476e6;
  Real T = 377.594;
  Real yco2 = 1.0;
  Real yh2s = 0.0;
  Real ych4 = 0.0;
  Real Zg;

  _fp->GasCompressibilityFactor(p, T, yco2, yh2s, ych4, Zg);

  ABS_TEST(Zg, 0.83, 1.0e-2);
}

TEST_F(PorousFlowBrineCO2CH4H2STest, gasProperties)
{

  const Real p = 5.0e6;
  const Real T = 350.0;
  const Real Xnacl = 0.1;

  FluidStatePhaseEnum phase_state;
  const unsigned int np = _fp->numPhases();
  const unsigned int nc = _fp->numComponents();
  const unsigned int nz = 3;
  std::vector<FluidStateProperties> fsp(np, FluidStateProperties(nc, nz));

  // Gas region
  Real Zco2 = 1.0;
  Real Zh2s = 0.00;
  Real Zch4 = 0.00;

  _fp->massFractions(p, T, Xnacl, Zco2, Zh2s, Zch4, phase_state, fsp);

  EXPECT_EQ(phase_state, FluidStatePhaseEnum::GAS);

  // Verify fluid density, viscosity and enthalpy
  _fp->gasProperties(p, T, fsp);
  Real gas_density = fsp[1].density;
  Real gas_viscosity = fsp[1].viscosity;
  Real gas_enthalpy = fsp[1].enthalpy;

  // Real density = _co2_fp->rho_from_p_T(p, T);
  Real viscosity = _co2_fp->mu_from_p_T(p, T);
  Real enthalpy = _co2_fp->h_from_p_T(p, T);

  ABS_TEST(gas_density, 91.2126, 1.0e-1);
  ABS_TEST(gas_viscosity, viscosity, 1.0e-8);
  // ABS_TEST(gas_enthalpy, enthalpy, 1.0e-8);

  // std::cout << gas_density <<  std::endl;

  // Verify derivatives
  Real ddensity_dp = fsp[1].ddensity_dp;
  Real ddensity_dT = fsp[1].ddensity_dT;
  Real ddensity_dZ = fsp[1].ddensity_dZ[0];
  Real dviscosity_dp = fsp[1].dviscosity_dp;
  Real dviscosity_dT = fsp[1].dviscosity_dT;
  Real dviscosity_dZ = fsp[1].dviscosity_dZ[0];
  Real denthalpy_dp = fsp[1].denthalpy_dp;
  Real denthalpy_dT = fsp[1].denthalpy_dT;
  Real denthalpy_dZ = fsp[1].denthalpy_dZ[0];

  // Derivatives wrt pressure
  const Real dp = 1.0e-2;
  _fp->gasProperties(p + dp, T, fsp);
  Real rho1 = fsp[1].density;
  Real mu1 = fsp[1].viscosity;
  Real h1 = fsp[1].enthalpy;

  _fp->gasProperties(p - dp, T, fsp);
  Real rho2 = fsp[1].density;
  Real mu2 = fsp[1].viscosity;
  Real h2 = fsp[1].enthalpy;

  REL_TEST(ddensity_dp, (rho1 - rho2) / (2.0 * dp), 1.0e-6);
  REL_TEST(dviscosity_dp, (mu1 - mu2) / (2.0 * dp), 1.0e-4);
  REL_TEST(denthalpy_dp, (h1 - h2) / (2.0 * dp), 1.0e-6);

  // Note: mass fraction changes with Z
  const Real dZ = 1.0e-8;
  _fp->massFractions(p, T, Xnacl, Zco2 + dZ, Zh2s + dZ, Zch4 + dZ, phase_state, fsp);
  _fp->gasProperties(p, T, fsp);
  rho1 = fsp[1].density;
  mu1 = fsp[1].viscosity;
  h1 = fsp[1].enthalpy;

  _fp->massFractions(p, T, Xnacl, Zco2 - dZ, Zh2s - dZ, Zch4 - dZ, phase_state, fsp);
  _fp->gasProperties(p, T, fsp);
  rho2 = fsp[1].density;
  mu2 = fsp[1].viscosity;
  h2 = fsp[1].enthalpy;

  // ABS_TEST(ddensity_dZ, (rho1 - rho2) / (2.0 * dZ), 1.0e-8);
  ABS_TEST(dviscosity_dZ, (mu1 - mu2) / (2.0 * dZ), 1.0e-8);
  ABS_TEST(denthalpy_dZ, (h1 - h2) / (2.0 * dZ), 1.0e-6);
}

TEST_F(PorousFlowBrineCO2CH4H2STest, liquidProperties)
{

  const Real p = 10.0e6;
  const Real T = 310.15;
  const Real Xnacl = 0.0;

  FluidStatePhaseEnum phase_state;
  const unsigned int np = _fp->numPhases();
  const unsigned int nc = _fp->numComponents();
  const unsigned int nz = 3;
  std::vector<FluidStateProperties> fsp(np, FluidStateProperties(nc, nz));

  // Liquid region
  Real Zco2 = 0.01;
  Real Zh2s = 0.00;
  Real Zch4 = 0.0009103;

  _fp->massFractions(p, T, Xnacl, Zco2, Zh2s, Zch4, phase_state, fsp);
  // std::cout << "Xco2 " << fsp[0].mass_fraction[1] << std::endl;

  EXPECT_EQ(phase_state, FluidStatePhaseEnum::LIQUID);

  // Verify fluid density, viscosity and enthalpy
  _fp->liquidProperties(p, T, Xnacl, fsp);
  Real liquid_density = fsp[0].density;
  Real liquid_viscosity = fsp[0].viscosity;
  Real liquid_enthalpy = fsp[0].enthalpy;

  // Real density = _co2_fp->rho_from_p_T(p, T);
  // Real viscosity = _co2_fp->mu_from_p_T(p, T);
  // Real enthalpy = _co2_fp->h_from_p_T(p, T);

  ABS_TEST(liquid_density, 998.72, 1.0e-1);
  // ABS_TEST(gas_viscosity, viscosity, 1.0e-8);
  // ABS_TEST(gas_enthalpy, enthalpy, 1.0e-8);

  // std::cout << "liquid density " << liquid_density <<  std::endl;

  // Verify derivatives

  Real ddensity_dp = fsp[0].ddensity_dp;
  Real ddensity_dT = fsp[0].ddensity_dT;
  Real ddensity_dZ = fsp[0].ddensity_dZ[0];
  Real dviscosity_dp = fsp[0].dviscosity_dp;
  Real dviscosity_dT = fsp[0].dviscosity_dT;
  Real dviscosity_dZ = fsp[0].dviscosity_dZ[0];
  Real denthalpy_dp = fsp[0].denthalpy_dp;
  Real denthalpy_dT = fsp[0].denthalpy_dT;
  Real denthalpy_dZ = fsp[0].denthalpy_dZ[0];

  // Derivatives wrt pressure

  const Real dp = 1.0e-1;
  _fp->gasProperties(p + dp, T, fsp);
  Real rho1 = fsp[0].density;
  Real mu1 = fsp[0].viscosity;
  Real h1 = fsp[0].enthalpy;

  _fp->gasProperties(p - dp, T, fsp);
  Real rho2 = fsp[0].density;
  Real mu2 = fsp[0].viscosity;
  Real h2 = fsp[0].enthalpy;

  REL_TEST(ddensity_dp, (rho1 - rho2) / (2.0 * dp), 1.0e-6);
  REL_TEST(dviscosity_dp, (mu1 - mu2) / (2.0 * dp), 1.0e-6);
  REL_TEST(denthalpy_dp, (h1 - h2) / (2.0 * dp), 1.0e-6);

  // Derivatives wrt Z
  const Real dZ = 1.0e-8;
  _fp->massFractions(p, T, Xnacl, Zco2 + dZ, Zh2s, Zch4, phase_state, fsp);
  _fp->liquidProperties(p, T, Xnacl, fsp);
  rho1 = fsp[0].density;
  mu1 = fsp[0].viscosity;
  h1 = fsp[0].enthalpy;

  _fp->massFractions(p, T, Xnacl, Zco2 - dZ, Zh2s, Zch4, phase_state, fsp);
  _fp->liquidProperties(p, T, Xnacl, fsp);
  rho2 = fsp[0].density;
  mu2 = fsp[0].viscosity;
  h2 = fsp[0].enthalpy;

  // REL_TEST(ddensity_dZ, (rho1 - rho2) / (2.0 * dZ), 1.0e-6);
  ABS_TEST(dviscosity_dZ, (mu1 - mu2) / (2.0 * dZ), 1.0e-6);
  REL_TEST(denthalpy_dZ, (h1 - h2) / (2.0 * dZ), 1.0e-6);
}

TEST_F(PorousFlowBrineCO2CH4H2STest, massFractions)
{
  const Real tol = 1.0e-8;
  Real p = 4.82e6;
  Real T = 310.95;
  Real Xnacl = 0.0;

  FluidStatePhaseEnum phase_state;
  const unsigned int np = _fp->numPhases();
  const unsigned int nc = _fp->numComponents();
  const unsigned int nz = 3;
  std::vector<FluidStateProperties> fsp(np, FluidStateProperties(nc, nz));

  // Single phase liquid region
  Real Zco2 = 0.002;
  Real Zh2s = 0.03;
  Real Zch4 = 0.000;

  _fp->massFractions(p, T, Xnacl, Zco2, Zh2s, Zch4, phase_state, fsp);
  EXPECT_EQ(phase_state, FluidStatePhaseEnum::LIQUID);

  // Verfify mass fraction values in liquid region
  Real Xh2o = fsp[0].mass_fraction[0];
  Real Yh2o = fsp[1].mass_fraction[0];
  Real Xco2 = fsp[0].mass_fraction[1];
  Real Yco2 = fsp[1].mass_fraction[1];
  Real Xch4 = fsp[0].mass_fraction[2];
  Real Ych4 = fsp[1].mass_fraction[2];
  Real Xh2s = fsp[0].mass_fraction[3];
  Real Yh2s = fsp[1].mass_fraction[3];
  Real Xnacl2 = fsp[0].mass_fraction[4];
  ABS_TEST(Xco2, Zco2, tol);
  ABS_TEST(Yco2, 0.0, tol);
  ABS_TEST(Xch4, Zch4, tol);
  ABS_TEST(Ych4, 0.0, tol);
  ABS_TEST(Xh2s, Zh2s, tol);
  ABS_TEST(Yh2s, 0.0, tol);
  ABS_TEST(Xh2o, 1.0 - Zco2 - Zch4 - Zh2s, tol);
  ABS_TEST(Yh2o, 0.0, tol);
  ABS_TEST(Xnacl2, Xnacl, tol);

  // Two phase region
  Zco2 = 0.1;
  Zh2s = 0.0;
  Zch4 = 0.0;

  _fp->massFractions(p, T, Xnacl, Zco2, Zh2s, Zch4, phase_state, fsp);
  EXPECT_EQ(phase_state, FluidStatePhaseEnum::TWOPHASE);

  // std::cout << "Xco2 " << fsp[0].mass_fraction[1] << std::endl;
}
