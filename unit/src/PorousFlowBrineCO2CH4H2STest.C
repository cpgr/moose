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
  _fp->equilibriumMoleFractions(p, T, Xnacl, yco2, ych4, yh2s, xco2, xch4, xh2s, yh2o);

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

  _fp->equilibriumMoleFractions(p, T, Xnacl, yco2, ych4, yh2s, xco2, xch4, xh2s, yh2o);

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

  const Real Mco2 = 0.0440098;
  const Real Mch4 = 0.01604;
  const Real Mh2s = 0.03408;

  const Real Yco2 = yco2 * Mco2 / (yco2 * Mco2 + yh2s * Mh2s + ych4 * Mch4);
  const Real Ych4 = ych4 * Mch4 / (yco2 * Mco2 + yh2s * Mh2s + ych4 * Mch4);
  const Real Yh2s = yh2s * Mh2s / (yco2 * Mco2 + yh2s * Mh2s + ych4 * Mch4);

  Real Xco2, Xch4, Xh2s, Yh2o;
  _fp->equilibriumMassFractions(p, T, Xnacl, Yco2, Ych4, Yh2s, Xco2, Xch4, Xh2s, Yh2o);

  ABS_TEST(Xco2, 0.022064415, 1.0e-3);
  ABS_TEST(Xch4, 0.0002576449, 1.0e-3);
  ABS_TEST(Xh2s, 0.008173358, 1.0e-3);
  ABS_TEST(Yh2o, 0.00191, 1.0e-3);
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
  Real Zco2 = 0.995;
  Real Zh2s = 0.001;
  Real Zch4 = 0.004;

  _fp->massFractions(p, T, Xnacl, Zco2, Zch4, Zh2s, phase_state, fsp);

  EXPECT_EQ(phase_state, FluidStatePhaseEnum::GAS);

  // Verify fluid density and viscosity
  _fp->gasProperties(p, T, fsp);
  Real gas_density = fsp[1].density;
  Real gas_viscosity = fsp[1].viscosity;

  Real viscosity = _co2_fp->mu_from_p_T(p, T);

  ABS_TEST(gas_density, 89.9580148861, 1.0e-8);
  ABS_TEST(gas_viscosity, viscosity, 1.0e-8);

  // Verify derivatives
  Real ddensity_dp = fsp[1].ddensity_dp;
  Real ddensity_dT = fsp[1].ddensity_dT;
  Real ddensity_dZco2 = fsp[1].ddensity_dZ[0];
  Real ddensity_dZch4 = fsp[1].ddensity_dZ[1];
  Real ddensity_dZh2s = fsp[1].ddensity_dZ[2];

  Real dviscosity_dp = fsp[1].dviscosity_dp;
  Real dviscosity_dT = fsp[1].dviscosity_dT;
  Real dviscosity_dZco2 = fsp[1].dviscosity_dZ[0];
  Real dviscosity_dZch4 = fsp[1].dviscosity_dZ[1];
  Real dviscosity_dZh2s = fsp[1].dviscosity_dZ[2];

  // Derivatives wrt pressure
  const Real dp = 1.0e-2;
  _fp->gasProperties(p + dp, T, fsp);
  Real rho1 = fsp[1].density;
  Real mu1 = fsp[1].viscosity;

  _fp->gasProperties(p - dp, T, fsp);
  Real rho2 = fsp[1].density;
  Real mu2 = fsp[1].viscosity;

  REL_TEST(ddensity_dp, (rho1 - rho2) / (2.0 * dp), 1.0e-6);
  REL_TEST(dviscosity_dp, (mu1 - mu2) / (2.0 * dp), 1.0e-4);

  // Derivatives wrt temperature
  const Real dT = 1.0e-6;
  _fp->gasProperties(p, T + dT, fsp);
  rho1 = fsp[1].density;
  mu1 = fsp[1].viscosity;

  _fp->gasProperties(p, T - dT, fsp);
  rho2 = fsp[1].density;
  mu2 = fsp[1].viscosity;

  REL_TEST(ddensity_dT, (rho1 - rho2) / (2.0 * dT), 1.0e-6);
  REL_TEST(dviscosity_dT, (mu1 - mu2) / (2.0 * dT), 1.0e-4);

  // Note: mass fraction changes with Z
  const Real dZ = 1.0e-8;
  _fp->massFractions(p, T, Xnacl, Zco2 + dZ, Zch4, Zh2s, phase_state, fsp);
  _fp->gasProperties(p, T, fsp);
  rho1 = fsp[1].density;
  mu1 = fsp[1].viscosity;

  _fp->massFractions(p, T, Xnacl, Zco2 - dZ, Zch4, Zh2s, phase_state, fsp);
  _fp->gasProperties(p, T, fsp);
  rho2 = fsp[1].density;
  mu2 = fsp[1].viscosity;

  REL_TEST(ddensity_dZco2, (rho1 - rho2) / (2.0 * dZ), 1.0e-6);
  REL_TEST(dviscosity_dZco2, (mu1 - mu2) / (2.0 * dZ), 1.0e-6);

  _fp->massFractions(p, T, Xnacl, Zco2, Zch4 + dZ, Zh2s, phase_state, fsp);
  _fp->gasProperties(p, T, fsp);
  rho1 = fsp[1].density;
  mu1 = fsp[1].viscosity;

  _fp->massFractions(p, T, Xnacl, Zco2, Zch4 - dZ, Zh2s, phase_state, fsp);
  _fp->gasProperties(p, T, fsp);
  rho2 = fsp[1].density;
  mu2 = fsp[1].viscosity;

  REL_TEST(ddensity_dZch4, (rho1 - rho2) / (2.0 * dZ), 1.0e-6);
  REL_TEST(dviscosity_dZch4, (mu1 - mu2) / (2.0 * dZ), 1.0e-6);

  _fp->massFractions(p, T, Xnacl, Zco2, Zch4, Zh2s + dZ, phase_state, fsp);
  _fp->gasProperties(p, T, fsp);
  rho1 = fsp[1].density;
  mu1 = fsp[1].viscosity;

  _fp->massFractions(p, T, Xnacl, Zco2, Zch4, Zh2s - dZ, phase_state, fsp);
  _fp->gasProperties(p, T, fsp);
  rho2 = fsp[1].density;
  mu2 = fsp[1].viscosity;

  REL_TEST(ddensity_dZh2s, (rho1 - rho2) / (2.0 * dZ), 1.0e-6);
  REL_TEST(dviscosity_dZh2s, (mu1 - mu2) / (2.0 * dZ), 1.0e-6);
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
  Real Zco2 = 0.04;
  Real Zh2s = 0.005;
  Real Zch4 = 0.0001;

  _fp->massFractions(p, T, Xnacl, Zco2, Zch4, Zh2s, phase_state, fsp);

  EXPECT_EQ(phase_state, FluidStatePhaseEnum::LIQUID);

  // Verify fluid density, viscosity and enthalpy
  _fp->liquidProperties(p, T, Xnacl, fsp);
  Real liquid_density = fsp[0].density;
  Real liquid_viscosity = fsp[0].viscosity;

  const Real brine_viscosity = _brine_fp->mu(p, T, Xnacl);

  ABS_TEST(liquid_density, 1005.85962535, 1.0e-6);
  ABS_TEST(liquid_viscosity, brine_viscosity, 1.0e-8);

  // Verify derivatives
  Real ddensity_dp = fsp[0].ddensity_dp;
  Real ddensity_dT = fsp[0].ddensity_dT;
  Real ddensity_dZco2 = fsp[0].ddensity_dZ[0];
  Real ddensity_dZch4 = fsp[0].ddensity_dZ[1];
  Real ddensity_dZh2s = fsp[0].ddensity_dZ[2];

  Real dviscosity_dp = fsp[0].dviscosity_dp;
  Real dviscosity_dT = fsp[0].dviscosity_dT;
  Real dviscosity_dZco2 = fsp[0].dviscosity_dZ[0];
  Real dviscosity_dZch4 = fsp[0].dviscosity_dZ[1];
  Real dviscosity_dZh2s = fsp[0].dviscosity_dZ[2];

  // Derivatives wrt pressure
  const Real dp = 1.0e-1;
  _fp->liquidProperties(p + dp, T, Xnacl, fsp);
  Real rho1 = fsp[0].density;
  Real mu1 = fsp[0].viscosity;

  _fp->liquidProperties(p - dp, T, Xnacl, fsp);
  Real rho2 = fsp[0].density;
  Real mu2 = fsp[0].viscosity;

  REL_TEST(ddensity_dp, (rho1 - rho2) / (2.0 * dp), 1.0e-4);
  REL_TEST(dviscosity_dp, (mu1 - mu2) / (2.0 * dp), 1.0e-3);

  // Derivatives wrt temperature
  const Real dT = 1.0e-4;
  _fp->liquidProperties(p, T + dT, Xnacl, fsp);
  rho1 = fsp[0].density;
  mu1 = fsp[0].viscosity;

  _fp->liquidProperties(p, T - dT, Xnacl, fsp);
  rho2 = fsp[0].density;
  mu2 = fsp[0].viscosity;

  REL_TEST(ddensity_dT, (rho1 - rho2) / (2.0 * dT), 1.0e-4);
  REL_TEST(dviscosity_dT, (mu1 - mu2) / (2.0 * dT), 1.0e-4);

  // Derivatives wrt Z's
  const Real dZ = 1.0e-8;
  _fp->massFractions(p, T, Xnacl, Zco2 + dZ, Zch4, Zh2s, phase_state, fsp);
  _fp->liquidProperties(p, T, Xnacl, fsp);
  rho1 = fsp[0].density;
  mu1 = fsp[0].viscosity;

  _fp->massFractions(p, T, Xnacl, Zco2 - dZ, Zch4, Zh2s, phase_state, fsp);
  _fp->liquidProperties(p, T, Xnacl, fsp);
  rho2 = fsp[0].density;
  mu2 = fsp[0].viscosity;

  REL_TEST(ddensity_dZco2, (rho1 - rho2) / (2.0 * dZ), 1.0e-6);
  ABS_TEST(dviscosity_dZco2, (mu1 - mu2) / (2.0 * dZ), 1.0e-6);

  _fp->massFractions(p, T, Xnacl, Zco2, Zch4 + dZ, Zh2s, phase_state, fsp);
  _fp->liquidProperties(p, T, Xnacl, fsp);
  rho1 = fsp[0].density;
  mu1 = fsp[0].viscosity;

  _fp->massFractions(p, T, Xnacl, Zco2, Zch4 - dZ, Zh2s, phase_state, fsp);
  _fp->liquidProperties(p, T, Xnacl, fsp);
  rho2 = fsp[0].density;
  mu2 = fsp[0].viscosity;

  REL_TEST(ddensity_dZch4, (rho1 - rho2) / (2.0 * dZ), 1.0e-6);
  ABS_TEST(dviscosity_dZch4, (mu1 - mu2) / (2.0 * dZ), 1.0e-6);

  _fp->massFractions(p, T, Xnacl, Zco2, Zch4, Zh2s + dZ, phase_state, fsp);
  _fp->liquidProperties(p, T, Xnacl, fsp);
  rho1 = fsp[0].density;
  mu1 = fsp[0].viscosity;

  _fp->massFractions(p, T, Xnacl, Zco2, Zch4, Zh2s - dZ, phase_state, fsp);
  _fp->liquidProperties(p, T, Xnacl, fsp);
  rho2 = fsp[0].density;
  mu2 = fsp[0].viscosity;

  REL_TEST(ddensity_dZh2s, (rho1 - rho2) / (2.0 * dZ), 1.0e-6);
  ABS_TEST(dviscosity_dZh2s, (mu1 - mu2) / (2.0 * dZ), 1.0e-6);
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

  _fp->massFractions(p, T, Xnacl, Zco2, Zch4, Zh2s, phase_state, fsp);
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

  // Verify derivatives in the liquid region
  Real dXco2_dp = fsp[0].dmass_fraction_dp[1];
  Real dXco2_dT = fsp[0].dmass_fraction_dT[1];
  Real dXco2_dX = fsp[0].dmass_fraction_dX[1];
  Real dXco2_dZco2 = fsp[0].dmass_fraction_dZ[1][0];
  Real dXco2_dZch4 = fsp[0].dmass_fraction_dZ[1][1];
  Real dXco2_dZh2s = fsp[0].dmass_fraction_dZ[1][2];
  Real dYco2_dp = fsp[1].dmass_fraction_dp[1];
  Real dYco2_dT = fsp[1].dmass_fraction_dT[1];
  Real dYco2_dX = fsp[1].dmass_fraction_dX[1];
  Real dYco2_dZco2 = fsp[1].dmass_fraction_dZ[1][0];
  Real dYco2_dZch4 = fsp[1].dmass_fraction_dZ[1][1];
  Real dYco2_dZh2s = fsp[1].dmass_fraction_dZ[1][2];

  Real dXch4_dp = fsp[0].dmass_fraction_dp[2];
  Real dXch4_dT = fsp[0].dmass_fraction_dT[2];
  Real dXch4_dX = fsp[0].dmass_fraction_dX[2];
  Real dXch4_dZco2 = fsp[0].dmass_fraction_dZ[2][0];
  Real dXch4_dZch4 = fsp[0].dmass_fraction_dZ[2][1];
  Real dXch4_dZh2s = fsp[0].dmass_fraction_dZ[2][2];
  Real dYch4_dp = fsp[1].dmass_fraction_dp[2];
  Real dYch4_dT = fsp[1].dmass_fraction_dT[2];
  Real dYch4_dX = fsp[1].dmass_fraction_dX[2];
  Real dYch4_dZco2 = fsp[1].dmass_fraction_dZ[2][0];
  Real dYch4_dZch4 = fsp[1].dmass_fraction_dZ[2][1];
  Real dYch4_dZh2s = fsp[1].dmass_fraction_dZ[2][2];

  Real dXh2s_dp = fsp[0].dmass_fraction_dp[3];
  Real dXh2s_dT = fsp[0].dmass_fraction_dT[3];
  Real dXh2s_dX = fsp[0].dmass_fraction_dX[3];
  Real dXh2s_dZco2 = fsp[0].dmass_fraction_dZ[3][0];
  Real dXh2s_dZch4 = fsp[0].dmass_fraction_dZ[3][1];
  Real dXh2s_dZh2s = fsp[0].dmass_fraction_dZ[3][2];
  Real dYh2s_dp = fsp[1].dmass_fraction_dp[3];
  Real dYh2s_dT = fsp[1].dmass_fraction_dT[3];
  Real dYh2s_dX = fsp[1].dmass_fraction_dX[3];
  Real dYh2s_dZco2 = fsp[1].dmass_fraction_dZ[3][0];
  Real dYh2s_dZch4 = fsp[1].dmass_fraction_dZ[3][1];
  Real dYh2s_dZh2s = fsp[1].dmass_fraction_dZ[3][2];

  Real dXnacl_dX = fsp[0].dmass_fraction_dX[4];

  ABS_TEST(dXco2_dp, 0.0, tol);
  ABS_TEST(dXco2_dT, 0.0, tol);
  ABS_TEST(dXco2_dX, 0.0, tol);
  ABS_TEST(dXco2_dZco2, 1.0, tol);
  ABS_TEST(dXco2_dZch4, 0.0, tol);
  ABS_TEST(dXco2_dZh2s, 0.0, tol);
  ABS_TEST(dYco2_dp, 0.0, tol);
  ABS_TEST(dYco2_dT, 0.0, tol);
  ABS_TEST(dYco2_dX, 0.0, tol);
  ABS_TEST(dYco2_dZco2, 0.0, tol);
  ABS_TEST(dYco2_dZch4, 0.0, tol);
  ABS_TEST(dYco2_dZh2s, 0.0, tol);

  ABS_TEST(dXch4_dp, 0.0, tol);
  ABS_TEST(dXch4_dT, 0.0, tol);
  ABS_TEST(dXch4_dX, 0.0, tol);
  ABS_TEST(dXch4_dZco2, 0.0, tol);
  ABS_TEST(dXch4_dZch4, 1.0, tol);
  ABS_TEST(dXch4_dZh2s, 0.0, tol);
  ABS_TEST(dYch4_dp, 0.0, tol);
  ABS_TEST(dYch4_dT, 0.0, tol);
  ABS_TEST(dYch4_dX, 0.0, tol);
  ABS_TEST(dYch4_dZco2, 0.0, tol);
  ABS_TEST(dYch4_dZch4, 0.0, tol);
  ABS_TEST(dYch4_dZh2s, 0.0, tol);

  ABS_TEST(dXh2s_dp, 0.0, tol);
  ABS_TEST(dXh2s_dT, 0.0, tol);
  ABS_TEST(dXh2s_dX, 0.0, tol);
  ABS_TEST(dXh2s_dZco2, 0.0, tol);
  ABS_TEST(dXh2s_dZch4, 0.0, tol);
  ABS_TEST(dXh2s_dZh2s, 1.0, tol);
  ABS_TEST(dYh2s_dp, 0.0, tol);
  ABS_TEST(dYh2s_dT, 0.0, tol);
  ABS_TEST(dYh2s_dX, 0.0, tol);
  ABS_TEST(dYh2s_dZco2, 0.0, tol);
  ABS_TEST(dYh2s_dZch4, 0.0, tol);
  ABS_TEST(dYh2s_dZh2s, 0.0, tol);

  ABS_TEST(dXnacl_dX, 1.0, tol);

  // Single phase gas region
  Zco2 = 0.995;
  Zh2s = 0.003;
  Zch4 = 0.002;

  _fp->massFractions(p, T, Xnacl, Zco2, Zch4, Zh2s, phase_state, fsp);
  EXPECT_EQ(phase_state, FluidStatePhaseEnum::GAS);

  // Verfify mass fraction values in gas region
  Xh2o = fsp[0].mass_fraction[0];
  Yh2o = fsp[1].mass_fraction[0];
  Xco2 = fsp[0].mass_fraction[1];
  Yco2 = fsp[1].mass_fraction[1];
  Xch4 = fsp[0].mass_fraction[2];
  Ych4 = fsp[1].mass_fraction[2];
  Xh2s = fsp[0].mass_fraction[3];
  Yh2s = fsp[1].mass_fraction[3];
  ABS_TEST(Xco2, 0.0, tol);
  ABS_TEST(Yco2, Zco2, tol);
  ABS_TEST(Xch4, 0.0, tol);
  ABS_TEST(Ych4, Zch4, tol);
  ABS_TEST(Xh2s, 0.0, tol);
  ABS_TEST(Yh2s, Zh2s, tol);
  ABS_TEST(Xh2o, 0.0, tol);
  ABS_TEST(Yh2o, 1.0 - Zco2 - Zch4 - Zh2s, tol);

  // Verify derivatives in the gas region
  dXco2_dp = fsp[0].dmass_fraction_dp[1];
  dXco2_dT = fsp[0].dmass_fraction_dT[1];
  dXco2_dX = fsp[0].dmass_fraction_dX[1];
  dXco2_dZco2 = fsp[0].dmass_fraction_dZ[1][0];
  dXco2_dZch4 = fsp[0].dmass_fraction_dZ[1][1];
  dXco2_dZh2s = fsp[0].dmass_fraction_dZ[1][2];
  dYco2_dp = fsp[1].dmass_fraction_dp[1];
  dYco2_dT = fsp[1].dmass_fraction_dT[1];
  dYco2_dX = fsp[1].dmass_fraction_dX[1];
  dYco2_dZco2 = fsp[1].dmass_fraction_dZ[1][0];
  dYco2_dZch4 = fsp[1].dmass_fraction_dZ[1][1];
  dYco2_dZh2s = fsp[1].dmass_fraction_dZ[1][2];

  dXch4_dp = fsp[0].dmass_fraction_dp[2];
  dXch4_dT = fsp[0].dmass_fraction_dT[2];
  dXch4_dX = fsp[0].dmass_fraction_dX[2];
  dXch4_dZco2 = fsp[0].dmass_fraction_dZ[2][0];
  dXch4_dZch4 = fsp[0].dmass_fraction_dZ[2][1];
  dXch4_dZh2s = fsp[0].dmass_fraction_dZ[2][2];
  dYch4_dp = fsp[1].dmass_fraction_dp[2];
  dYch4_dT = fsp[1].dmass_fraction_dT[2];
  dYch4_dX = fsp[1].dmass_fraction_dX[2];
  dYch4_dZco2 = fsp[1].dmass_fraction_dZ[2][0];
  dYch4_dZch4 = fsp[1].dmass_fraction_dZ[2][1];
  dYch4_dZh2s = fsp[1].dmass_fraction_dZ[2][2];

  dXh2s_dp = fsp[0].dmass_fraction_dp[3];
  dXh2s_dT = fsp[0].dmass_fraction_dT[3];
  dXh2s_dX = fsp[0].dmass_fraction_dX[3];
  dXh2s_dZco2 = fsp[0].dmass_fraction_dZ[3][0];
  dXh2s_dZch4 = fsp[0].dmass_fraction_dZ[3][1];
  dXh2s_dZh2s = fsp[0].dmass_fraction_dZ[3][2];
  dYh2s_dp = fsp[1].dmass_fraction_dp[3];
  dYh2s_dT = fsp[1].dmass_fraction_dT[3];
  dYh2s_dX = fsp[1].dmass_fraction_dX[3];
  dYh2s_dZco2 = fsp[1].dmass_fraction_dZ[3][0];
  dYh2s_dZch4 = fsp[1].dmass_fraction_dZ[3][1];
  dYh2s_dZh2s = fsp[1].dmass_fraction_dZ[3][2];

  ABS_TEST(dXco2_dp, 0.0, tol);
  ABS_TEST(dXco2_dT, 0.0, tol);
  ABS_TEST(dXco2_dX, 0.0, tol);
  ABS_TEST(dXco2_dZco2, 0.0, tol);
  ABS_TEST(dXco2_dZch4, 0.0, tol);
  ABS_TEST(dXco2_dZh2s, 0.0, tol);
  ABS_TEST(dYco2_dp, 0.0, tol);
  ABS_TEST(dYco2_dT, 0.0, tol);
  ABS_TEST(dYco2_dX, 0.0, tol);
  ABS_TEST(dYco2_dZco2, 1.0, tol);
  ABS_TEST(dYco2_dZch4, 0.0, tol);
  ABS_TEST(dYco2_dZh2s, 0.0, tol);

  ABS_TEST(dXch4_dp, 0.0, tol);
  ABS_TEST(dXch4_dT, 0.0, tol);
  ABS_TEST(dXch4_dX, 0.0, tol);
  ABS_TEST(dXch4_dZco2, 0.0, tol);
  ABS_TEST(dXch4_dZch4, 0.0, tol);
  ABS_TEST(dXch4_dZh2s, 0.0, tol);
  ABS_TEST(dYch4_dp, 0.0, tol);
  ABS_TEST(dYch4_dT, 0.0, tol);
  ABS_TEST(dYch4_dX, 0.0, tol);
  ABS_TEST(dYch4_dZco2, 0.0, tol);
  ABS_TEST(dYch4_dZch4, 1.0, tol);
  ABS_TEST(dYch4_dZh2s, 0.0, tol);

  ABS_TEST(dXh2s_dp, 0.0, tol);
  ABS_TEST(dXh2s_dT, 0.0, tol);
  ABS_TEST(dXh2s_dX, 0.0, tol);
  ABS_TEST(dXh2s_dZco2, 0.0, tol);
  ABS_TEST(dXh2s_dZch4, 0.0, tol);
  ABS_TEST(dXh2s_dZh2s, 0.0, tol);
  ABS_TEST(dYh2s_dp, 0.0, tol);
  ABS_TEST(dYh2s_dT, 0.0, tol);
  ABS_TEST(dYh2s_dX, 0.0, tol);
  ABS_TEST(dYh2s_dZco2, 0.0, tol);
  ABS_TEST(dYh2s_dZch4, 0.0, tol);
  ABS_TEST(dYh2s_dZh2s, 1.0, tol);

  // Two phase region
  Zco2 = 0.25;
  Zh2s = 0.2;
  Zch4 = 0.1;

  _fp->massFractions(p, T, Xnacl, Zco2, Zch4, Zh2s, phase_state, fsp);
  EXPECT_EQ(phase_state, FluidStatePhaseEnum::TWOPHASE);

  Xh2o = fsp[0].mass_fraction[0];
  Yh2o = fsp[1].mass_fraction[0];
  Xco2 = fsp[0].mass_fraction[1];
  Yco2 = fsp[1].mass_fraction[1];
  Xch4 = fsp[0].mass_fraction[2];
  Ych4 = fsp[1].mass_fraction[2];
  Xh2s = fsp[0].mass_fraction[3];
  Yh2s = fsp[1].mass_fraction[3];

  ABS_TEST(Xco2, 0.0118612, 1.0e-4);
  ABS_TEST(Yco2, 0.458951, 1.0e-4);
  ABS_TEST(Xch4, 0.00033786, 1.0e-4);
  ABS_TEST(Ych4, 0.187447, 1.0e-4);
  ABS_TEST(Xh2s, 0.0263567, 1.0e-4);
  ABS_TEST(Yh2s, 0.352359, 1.0e-4);
  ABS_TEST(Xh2o, 0.961444, 1.0e-4);
  ABS_TEST(Yh2o, 1.0 - Yco2 - Ych4 - Yh2s, tol);

  // Verify derivatives in the two-phase region
  dXco2_dp = fsp[0].dmass_fraction_dp[1];
  dXco2_dT = fsp[0].dmass_fraction_dT[1];
  dXco2_dX = fsp[0].dmass_fraction_dX[1];
  dXco2_dZco2 = fsp[0].dmass_fraction_dZ[1][0];
  dXco2_dZch4 = fsp[0].dmass_fraction_dZ[1][1];
  dXco2_dZh2s = fsp[0].dmass_fraction_dZ[1][2];
  dYco2_dp = fsp[1].dmass_fraction_dp[1];
  dYco2_dT = fsp[1].dmass_fraction_dT[1];
  dYco2_dX = fsp[1].dmass_fraction_dX[1];
  dYco2_dZco2 = fsp[1].dmass_fraction_dZ[1][0];
  dYco2_dZch4 = fsp[1].dmass_fraction_dZ[1][1];
  dYco2_dZh2s = fsp[1].dmass_fraction_dZ[1][2];

  dXch4_dp = fsp[0].dmass_fraction_dp[2];
  dXch4_dT = fsp[0].dmass_fraction_dT[2];
  dXch4_dX = fsp[0].dmass_fraction_dX[2];
  dXch4_dZco2 = fsp[0].dmass_fraction_dZ[2][0];
  dXch4_dZch4 = fsp[0].dmass_fraction_dZ[2][1];
  dXch4_dZh2s = fsp[0].dmass_fraction_dZ[2][2];
  dYch4_dp = fsp[1].dmass_fraction_dp[2];
  dYch4_dT = fsp[1].dmass_fraction_dT[2];
  dYch4_dX = fsp[1].dmass_fraction_dX[2];
  dYch4_dZco2 = fsp[1].dmass_fraction_dZ[2][0];
  dYch4_dZch4 = fsp[1].dmass_fraction_dZ[2][1];
  dYch4_dZh2s = fsp[1].dmass_fraction_dZ[2][2];

  dXh2s_dp = fsp[0].dmass_fraction_dp[3];
  dXh2s_dT = fsp[0].dmass_fraction_dT[3];
  dXh2s_dX = fsp[0].dmass_fraction_dX[3];
  dXh2s_dZco2 = fsp[0].dmass_fraction_dZ[3][0];
  dXh2s_dZch4 = fsp[0].dmass_fraction_dZ[3][1];
  dXh2s_dZh2s = fsp[0].dmass_fraction_dZ[3][2];
  dYh2s_dp = fsp[1].dmass_fraction_dp[3];
  dYh2s_dT = fsp[1].dmass_fraction_dT[3];
  dYh2s_dX = fsp[1].dmass_fraction_dX[3];
  dYh2s_dZco2 = fsp[1].dmass_fraction_dZ[3][0];
  dYh2s_dZch4 = fsp[1].dmass_fraction_dZ[3][1];
  dYh2s_dZh2s = fsp[1].dmass_fraction_dZ[3][2];

  // Derivatives wrt pressure
  const Real dp = 1.0e-2;
  _fp->massFractions(p + dp, T, Xnacl, Zco2, Zch4, Zh2s, phase_state, fsp);
  REL_TEST(dXco2_dp, (fsp[0].mass_fraction[1] - Xco2) / dp, tol);
  REL_TEST(dXch4_dp, (fsp[0].mass_fraction[2] - Xch4) / dp, tol);
  REL_TEST(dXh2s_dp, (fsp[0].mass_fraction[3] - Xh2s) / dp, tol);
  REL_TEST(dYco2_dp, (fsp[1].mass_fraction[1] - Yco2) / dp, tol);
  REL_TEST(dYch4_dp, (fsp[1].mass_fraction[2] - Ych4) / dp, tol);
  REL_TEST(dYh2s_dp, (fsp[1].mass_fraction[3] - Yh2s) / dp, tol);

  // Derivatives wrt temprature
  const Real dT = 1.0e-4;
  _fp->massFractions(p, T + dT, Xnacl, Zco2, Zch4, Zh2s, phase_state, fsp);
  REL_TEST(dXco2_dT, (fsp[0].mass_fraction[1] - Xco2) / dT, tol);
  REL_TEST(dXch4_dT, (fsp[0].mass_fraction[2] - Xch4) / dT, tol);
  REL_TEST(dXh2s_dT, (fsp[0].mass_fraction[3] - Xh2s) / dT, tol);
  REL_TEST(dYco2_dT, (fsp[1].mass_fraction[1] - Yco2) / dT, tol);
  REL_TEST(dYch4_dT, (fsp[1].mass_fraction[2] - Ych4) / dT, tol);
  REL_TEST(dYh2s_dT, (fsp[1].mass_fraction[3] - Yh2s) / dT, tol);

  // Derivatives wrt salt mass fraction
  const Real dX = 1.0e-8;
  _fp->massFractions(p, T, Xnacl + dX, Zco2, Zch4, Zh2s, phase_state, fsp);
  REL_TEST(dXco2_dX, (fsp[0].mass_fraction[1] - Xco2) / dX, tol);
  REL_TEST(dXch4_dX, (fsp[0].mass_fraction[2] - Xch4) / dX, tol);
  REL_TEST(dXh2s_dX, (fsp[0].mass_fraction[3] - Xh2s) / dX, tol);
  REL_TEST(dYco2_dX, (fsp[1].mass_fraction[1] - Yco2) / dX, tol);
  REL_TEST(dYch4_dX, (fsp[1].mass_fraction[2] - Ych4) / dX, tol);
  REL_TEST(dYh2s_dX, (fsp[1].mass_fraction[3] - Yh2s) / dX, tol);

  // Derivatives wrt Z's
  const Real dZ = 1.0e-8;
  _fp->massFractions(p, T, Xnacl, Zco2 + dZ, Zch4, Zh2s, phase_state, fsp);
  REL_TEST(dXco2_dZco2, (fsp[0].mass_fraction[1] - Xco2) / dZ, tol);
  REL_TEST(dXch4_dZco2, (fsp[0].mass_fraction[2] - Xch4) / dZ, tol);
  REL_TEST(dXh2s_dZco2, (fsp[0].mass_fraction[3] - Xh2s) / dZ, tol);
  REL_TEST(dYco2_dZco2, (fsp[1].mass_fraction[1] - Yco2) / dZ, tol);
  REL_TEST(dYch4_dZco2, (fsp[1].mass_fraction[2] - Ych4) / dZ, tol);
  REL_TEST(dYh2s_dZco2, (fsp[1].mass_fraction[3] - Yh2s) / dZ, tol);

  _fp->massFractions(p, T, Xnacl, Zco2, Zch4 + dZ, Zh2s, phase_state, fsp);
  REL_TEST(dXco2_dZch4, (fsp[0].mass_fraction[1] - Xco2) / dZ, tol);
  REL_TEST(dXch4_dZch4, (fsp[0].mass_fraction[2] - Xch4) / dZ, tol);
  REL_TEST(dXh2s_dZch4, (fsp[0].mass_fraction[3] - Xh2s) / dZ, tol);
  REL_TEST(dYco2_dZch4, (fsp[1].mass_fraction[1] - Yco2) / dZ, tol);
  REL_TEST(dYch4_dZch4, (fsp[1].mass_fraction[2] - Ych4) / dZ, tol);
  REL_TEST(dYh2s_dZch4, (fsp[1].mass_fraction[3] - Yh2s) / dZ, tol);

  _fp->massFractions(p, T, Xnacl, Zco2, Zch4, Zh2s + dZ, phase_state, fsp);
  REL_TEST(dXco2_dZh2s, (fsp[0].mass_fraction[1] - Xco2) / dZ, tol);
  REL_TEST(dXch4_dZh2s, (fsp[0].mass_fraction[2] - Xch4) / dZ, tol);
  REL_TEST(dXh2s_dZh2s, (fsp[0].mass_fraction[3] - Xh2s) / dZ, tol);
  REL_TEST(dYco2_dZh2s, (fsp[1].mass_fraction[1] - Yco2) / dZ, tol);
  REL_TEST(dYch4_dZh2s, (fsp[1].mass_fraction[2] - Ych4) / dZ, tol);
  REL_TEST(dYh2s_dZh2s, (fsp[1].mass_fraction[3] - Yh2s) / dZ, tol);
}

/*
 * Verify calculation of gas saturation and derivatives in the two-phase region
 */
TEST_F(PorousFlowBrineCO2CH4H2STest, twoPhaseProperties)
{
  const Real p = 1.0e6;
  const Real T = 350.0;
  const Real Xnacl = 0.1;

  FluidStatePhaseEnum phase_state;
  const unsigned int np = _fp->numPhases();
  const unsigned int nc = _fp->numComponents();
  const unsigned int nz = 3;
  std::vector<FluidStateProperties> fsp(np, FluidStateProperties(nc, nz));

  // Mass fractions that correspond to a gas saturation of 0.25
  Real Zco2 = 0.25;
  Real Zh2s = 0.2;
  Real Zch4 = 0.1;

  _fp->massFractions(p, T, Xnacl, Zco2, Zch4, Zh2s, phase_state, fsp);
  EXPECT_EQ(phase_state, FluidStatePhaseEnum::TWOPHASE);

  // Calculate Z that gives a saturation of 0.25
  Real gas_saturation = 0.25;
  Real liquid_pressure = p + _pc->capillaryPressure(1.0 - gas_saturation);
  // Calculate gas density and liquid density
  _fp->gasProperties(p, T, fsp);
  _fp->liquidProperties(liquid_pressure, T, Xnacl, fsp);

  // The mass fractions that correspond to a gas_saturation = 0.25
  Zco2 = (gas_saturation * fsp[1].density * fsp[1].mass_fraction[1] +
          (1.0 - gas_saturation) * fsp[0].density * fsp[0].mass_fraction[1]) /
         (gas_saturation * fsp[1].density + (1.0 - gas_saturation) * fsp[0].density);

  Zch4 = (gas_saturation * fsp[1].density * fsp[1].mass_fraction[2] +
          (1.0 - gas_saturation) * fsp[0].density * fsp[0].mass_fraction[2]) /
         (gas_saturation * fsp[1].density + (1.0 - gas_saturation) * fsp[0].density);

  Zh2s = (gas_saturation * fsp[1].density * fsp[1].mass_fraction[3] +
          (1.0 - gas_saturation) * fsp[0].density * fsp[0].mass_fraction[3]) /
         (gas_saturation * fsp[1].density + (1.0 - gas_saturation) * fsp[0].density);

  // Calculate the gas saturation and derivatives
  _fp->twoPhaseProperties(p, T, Xnacl, Zco2, Zch4, Zh2s, fsp);

  ABS_TEST(fsp[1].saturation, gas_saturation, 5.0e-2);

  // Test the derivatives
  const Real dp = 1.0e-1;
  gas_saturation = fsp[1].saturation;
  Real dgas_saturation_dp = fsp[1].dsaturation_dp;
  Real dgas_saturation_dT = fsp[1].dsaturation_dT;
  Real dgas_saturation_dX = fsp[1].dsaturation_dX;
  Real dgas_saturation_dZco2 = fsp[1].dsaturation_dZ[0];
  Real dgas_saturation_dZch4 = fsp[1].dsaturation_dZ[1];
  Real dgas_saturation_dZh2s = fsp[1].dsaturation_dZ[2];

  // Derivative wrt pressure
  _fp->massFractions(p + dp, T, Xnacl, Zco2, Zch4, Zh2s, phase_state, fsp);
  _fp->gasProperties(p + dp, T, fsp);
  _fp->twoPhaseProperties(p + dp, T, Xnacl, Zco2, Zch4, Zh2s, fsp);
  Real gsat1 = fsp[1].saturation;

  _fp->massFractions(p - dp, T, Xnacl, Zco2, Zch4, Zh2s, phase_state, fsp);
  _fp->gasProperties(p - dp, T, fsp);
  _fp->twoPhaseProperties(p - dp, T, Xnacl, Zco2, Zch4, Zh2s, fsp);
  Real gsat2 = fsp[1].saturation;

  REL_TEST(dgas_saturation_dp, (gsat1 - gsat2) / (2.0 * dp), 1.0e-6);

  // Derivative wrt temperature
  const Real dT = 1.0e-4;
  _fp->massFractions(p, T + dT, Xnacl, Zco2, Zch4, Zh2s, phase_state, fsp);
  _fp->gasProperties(p, T + dT, fsp);
  _fp->twoPhaseProperties(p, T + dT, Xnacl, Zco2, Zch4, Zh2s, fsp);
  gsat1 = fsp[1].saturation;

  _fp->massFractions(p, T - dT, Xnacl, Zco2, Zch4, Zh2s, phase_state, fsp);
  _fp->gasProperties(p, T - dT, fsp);
  _fp->twoPhaseProperties(p, T - dT, Xnacl, Zco2, Zch4, Zh2s, fsp);
  gsat2 = fsp[1].saturation;

  REL_TEST(dgas_saturation_dT, (gsat1 - gsat2) / (2.0 * dT), 1.0e-4);

  // Derivative wrt Xnacl
  const Real dX = 1.0e-8;
  _fp->massFractions(p, T, Xnacl + dX, Zco2, Zch4, Zh2s, phase_state, fsp);
  _fp->gasProperties(p, T, fsp);
  _fp->twoPhaseProperties(p, T, Xnacl + dX, Zco2, Zch4, Zh2s, fsp);
  gsat1 = fsp[1].saturation;

  _fp->massFractions(p, T, Xnacl - dX, Zco2, Zch4, Zh2s, phase_state, fsp);
  _fp->gasProperties(p, T, fsp);
  _fp->twoPhaseProperties(p, T, Xnacl - dX, Zco2, Zch4, Zh2s, fsp);
  gsat2 = fsp[1].saturation;

  REL_TEST(dgas_saturation_dX, (gsat1 - gsat2) / (2.0 * dX), 1.0e-4);

  // Derivative wrt Z
  const Real dZ = 1.0e-8;

  _fp->massFractions(p, T, Xnacl, Zco2 + dZ, Zch4, Zh2s, phase_state, fsp);
  _fp->gasProperties(p, T, fsp);
  _fp->twoPhaseProperties(p, T, Xnacl, Zco2 + dZ, Zch4, Zh2s, fsp);
  gsat1 = fsp[1].saturation;

  _fp->massFractions(p, T, Xnacl, Zco2 - dZ, Zch4, Zh2s, phase_state, fsp);
  _fp->twoPhaseProperties(p, T, Xnacl, Zco2 - dZ, Zch4, Zh2s, fsp);
  gsat2 = fsp[1].saturation;

  REL_TEST(dgas_saturation_dZco2, (gsat1 - gsat2) / (2.0 * dZ), 1.0e-4);

  _fp->massFractions(p, T, Xnacl, Zco2, Zch4 + dZ, Zh2s, phase_state, fsp);
  _fp->gasProperties(p, T, fsp);
  _fp->twoPhaseProperties(p, T, Xnacl, Zco2, Zch4 + dZ, Zh2s, fsp);
  gsat1 = fsp[1].saturation;

  _fp->massFractions(p, T, Xnacl, Zco2, Zch4 - dZ, Zh2s, phase_state, fsp);
  _fp->twoPhaseProperties(p, T, Xnacl, Zco2, Zch4 - dZ, Zh2s, fsp);
  gsat2 = fsp[1].saturation;

  REL_TEST(dgas_saturation_dZch4, (gsat1 - gsat2) / (2.0 * dZ), 1.0e-4);

  _fp->massFractions(p, T, Xnacl, Zco2, Zch4, Zh2s + dZ, phase_state, fsp);
  _fp->gasProperties(p, T, fsp);
  _fp->twoPhaseProperties(p, T, Xnacl, Zco2, Zch4, Zh2s + dZ, fsp);
  gsat1 = fsp[1].saturation;

  _fp->massFractions(p, T, Xnacl, Zco2, Zch4, Zh2s - dZ, phase_state, fsp);
  _fp->twoPhaseProperties(p, T, Xnacl, Zco2, Zch4, Zh2s - dZ, fsp);
  gsat2 = fsp[1].saturation;

  REL_TEST(dgas_saturation_dZh2s, (gsat1 - gsat2) / (2.0 * dZ), 1.0e-4);
}

/*
 * Verify calculation of total mass fraction given a gas saturation
 */
TEST_F(PorousFlowBrineCO2CH4H2STest, totalMassFraction)
{
  const Real p = 1.0e6;
  const Real T = 350.0;
  const Real Xnacl = 0.1;
  const Real s = 0.25;
  std::vector<Real> Yi{0.8, 0.1, 0.1};
  const unsigned qp = 0;

  const unsigned int np = _fp->numPhases();
  const unsigned int nc = _fp->numComponents();
  const unsigned int nz = 3;
  std::vector<FluidStateProperties> fsp(np, FluidStateProperties(nc, nz));

  // Calculate Z based on composition of gas and saturation
  const Real Zco2 = _fp->totalMassFraction(p, T, Xnacl, s, Yi, fsp, 0, qp);
  const Real Zch4 = _fp->totalMassFraction(p, T, Xnacl, s, Yi, fsp, 1, qp);
  const Real Zh2s = _fp->totalMassFraction(p, T, Xnacl, s, Yi, fsp, 2, qp);

  // Test that the saturation calculated in this fluid state using Z's is equal to s
  FluidStatePhaseEnum phase_state;
  _fp->clearFluidStateProperties(fsp);

  _fp->massFractions(p, T, Xnacl, Zco2, Zch4, Zh2s, phase_state, fsp);
  EXPECT_EQ(phase_state, FluidStatePhaseEnum::TWOPHASE);

  // Must call gas properties before calculation gas saturation
  _fp->gasProperties(p, T, fsp);

  const Real gas_saturation = _fp->gasSaturation(p, T, Xnacl, Zco2, Zch4, Zh2s, fsp);
  // The tolerance is quite loose as the calculated saturation isn't exact due to the
  // presence of water vopor in the gas phase that is not included in the input
  ABS_TEST(gas_saturation, s, 2.0e-2);
}
