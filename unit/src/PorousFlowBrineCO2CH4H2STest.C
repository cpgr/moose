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
  std::cout << xco2 << ", " << xch4 << ", " << xh2s << ", " << yh2o << std::endl;

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
  std::cout << xco2 << ", " << xch4 << ", " << xh2s << ", " << yh2o << std::endl;
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
  std::cout << Xco2 << ", " << Xch4 << ", " << Xh2s << std::endl;
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
  std::vector<FluidStateProperties> fsp(np, FluidStateProperties(nc));

  // Single phase liquid region
  Real Zco2 = 0.002;
  Real Zh2s = 0.001;
  Real Zch4 = 0.003;

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
  Zco2 = 0.2;
  Zh2s = 0.1;
  Zch4 = 0.3;

  _fp->massFractions(p, T, Xnacl, Zco2, Zh2s, Zch4, phase_state, fsp);
  EXPECT_EQ(phase_state, FluidStatePhaseEnum::TWOPHASE);
}
