//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PorousFlowWaterNCGADTest.h"
#include "FluidPropertiesTestUtils.h"

/**
 * Verify that the correct name is supplied
 */
TEST_F(PorousFlowWaterNCGADTest, name) { EXPECT_EQ("water-ncg", _fp->fluidStateName()); }

/*
 * Verify calculation of equilibrium mass fraction and derivatives
 */
TEST_F(PorousFlowWaterNCGADTest, equilibriumMassFraction)
{
  DualReal p = 1.0e6;
  p.derivatives()[0] = 1.0;

  DualReal T = 350.0;
  T.derivatives()[1] = 1.0;

  const Real dp = 1.0e-1;
  const Real dT = 1.0e-6;

  DualReal Xncg, Yh2o, Xncg1, Yh2o1, Xncg2, Yh2o2;
  _fp->equilibriumMassFractions(p, T, Xncg, Yh2o);
  _fp->equilibriumMassFractions(p - dp, T, Xncg1, Yh2o1);
  _fp->equilibriumMassFractions(p + dp, T, Xncg2, Yh2o2);

  Real dXncg_dp_fd = (Xncg2.value() - Xncg1.value()) / (2.0 * dp);
  Real dYh2o_dp_fd = (Yh2o2.value() - Yh2o1.value()) / (2.0 * dp);

  REL_TEST(Xncg.derivatives()[0], dXncg_dp_fd, 1.0e-8);
  REL_TEST(Yh2o.derivatives()[0], dYh2o_dp_fd, 1.0e-7);

  _fp->equilibriumMassFractions(p, T - dT, Xncg1, Yh2o1);
  _fp->equilibriumMassFractions(p, T + dT, Xncg2, Yh2o2);

  Real dXncg_dT_fd = (Xncg2.value() - Xncg1.value()) / (2.0 * dT);
  Real dYh2o_dT_fd = (Yh2o2.value() - Yh2o1.value()) / (2.0 * dT);

  REL_TEST(Xncg.derivatives()[1], dXncg_dT_fd, 1.0e-7);
  REL_TEST(Yh2o.derivatives()[1], dYh2o_dT_fd, 1.0e-7);
}

/*
 * Verify calculation of actual mass fraction and derivatives depending on value of
 * total mass fraction
 */
TEST_F(PorousFlowWaterNCGADTest, massFraction)
{
  DualReal p = 1.0e6;
  p.derivatives()[0] = 1.0;

  DualReal T = 350.0;
  T.derivatives()[1] = 1.0;

  FluidStatePhaseEnum phase_state;
  std::vector<ADFluidStateProperties> adfsp(2, ADFluidStateProperties(2));

  // Liquid region
  DualReal Z = 0.0001;
  Z.derivatives()[2] = 1.0;

  _fp->massFractions(p, T, Z, phase_state, adfsp);
  EXPECT_EQ(phase_state, FluidStatePhaseEnum::LIQUID);

  // Verfify mass fraction values
  DualReal Xncg = adfsp[0].mass_fraction[1];
  DualReal Yncg = adfsp[1].mass_fraction[1];
  DualReal Xh2o = adfsp[0].mass_fraction[0];
  DualReal Yh2o = adfsp[1].mass_fraction[0];
  ABS_TEST(Xncg.value(), Z.value(), 1.0e-12);
  ABS_TEST(Yncg.value(), 0.0, 1.0e-12);
  ABS_TEST(Xh2o.value(), 1.0 - Z.value(), 1.0e-12);
  ABS_TEST(Yh2o.value(), 0.0, 1.0e-12);

  // Verify derivatives
  Real dXncg_dp = Xncg.derivatives()[0];
  Real dXncg_dT = Xncg.derivatives()[1];
  Real dXncg_dZ = Xncg.derivatives()[2];
  Real dYncg_dp = Yncg.derivatives()[0];
  Real dYncg_dT = Yncg.derivatives()[1];
  Real dYncg_dZ = Yncg.derivatives()[2];
  ABS_TEST(dXncg_dp, 0.0, 1.0e-12);
  ABS_TEST(dXncg_dT, 0.0, 1.0e-12);
  ABS_TEST(dXncg_dZ, 1.0, 1.0e-12);
  ABS_TEST(dYncg_dp, 0.0, 1.0e-12);
  ABS_TEST(dYncg_dT, 0.0, 1.0e-12);
  ABS_TEST(dYncg_dZ, 0.0, 1.0e-12);

  // Gas region
  Z = 0.995;
  Z.derivatives()[2] = 1.0;

  _fp->massFractions(p, T, Z, phase_state, adfsp);
  EXPECT_EQ(phase_state, FluidStatePhaseEnum::GAS);

  // Verfify mass fraction values
  Xncg = adfsp[0].mass_fraction[1];
  Yncg = adfsp[1].mass_fraction[1];
  Xh2o = adfsp[0].mass_fraction[0];
  Yh2o = adfsp[1].mass_fraction[0];
  ABS_TEST(Xncg.value(), 0.0, 1.0e-12);
  ABS_TEST(Yncg.value(), Z.value(), 1.0e-12);
  ABS_TEST(Xh2o.value(), 0.0, 1.0e-12);
  ABS_TEST(Yh2o.value(), 1.0 - Z.value(), 1.0e-12);

  // Verify derivatives
  dXncg_dp = Xncg.derivatives()[0];
  dXncg_dT = Xncg.derivatives()[1];
  dXncg_dZ = Xncg.derivatives()[2];
  dYncg_dp = Yncg.derivatives()[0];
  dYncg_dT = Yncg.derivatives()[1];
  dYncg_dZ = Yncg.derivatives()[2];
  ABS_TEST(dXncg_dp, 0.0, 1.0e-12);
  ABS_TEST(dXncg_dT, 0.0, 1.0e-12);
  ABS_TEST(dXncg_dZ, 0.0, 1.0e-12);
  ABS_TEST(dYncg_dp, 0.0, 1.0e-12);
  ABS_TEST(dYncg_dT, 0.0, 1.0e-12);
  ABS_TEST(dYncg_dZ, 1.0, 1.0e-12);

  // Two phase region. In this region, the mass fractions and derivatives can
  //  be verified using the equilibrium mass fraction derivatives that have
  // been verified above
  Z = 0.45;
  Z.derivatives()[2] = 1.0;

  _fp->massFractions(p, T, Z, phase_state, adfsp);
  EXPECT_EQ(phase_state, FluidStatePhaseEnum::TWOPHASE);

  // Equilibrium mass fractions and derivatives
  DualReal Xncg_eq, Yh2o_eq;
  _fp->equilibriumMassFractions(p, T, Xncg_eq, Yh2o_eq);

  // Verfify mass fraction values
  Xncg = adfsp[0].mass_fraction[1];
  Yncg = adfsp[1].mass_fraction[1];
  Xh2o = adfsp[0].mass_fraction[0];
  Yh2o = adfsp[1].mass_fraction[0];
  ABS_TEST(Xncg, Xncg_eq, 1.0e-12);
  ABS_TEST(Yncg, 1.0 - Yh2o_eq, 1.0e-12);
  ABS_TEST(Xh2o, 1.0 - Xncg_eq, 1.0e-12);
  ABS_TEST(Yh2o, Yh2o_eq, 1.0e-12);

  // Verify derivatives wrt p and T
  dXncg_dp = Xncg.derivatives()[0];
  dXncg_dT = Xncg.derivatives()[1];
  dXncg_dZ = Xncg.derivatives()[2];
  dYncg_dp = Yncg.derivatives()[0];
  dYncg_dT = Yncg.derivatives()[1];
  dYncg_dZ = Yncg.derivatives()[2];
  ABS_TEST(dXncg_dp, Xncg_eq.derivatives()[0], 1.0e-8);
  ABS_TEST(dXncg_dT, Xncg_eq.derivatives()[1], 1.0e-8);
  ABS_TEST(dXncg_dZ, 0.0, 1.0e-8);
  ABS_TEST(dYncg_dp, -Yh2o_eq.derivatives()[0], 1.0e-8);
  ABS_TEST(dYncg_dT, -Yh2o_eq.derivatives()[1], 1.0e-8);
  ABS_TEST(dYncg_dZ, 0.0, 1.0e-8);

  // Use finite differences to verify derivative wrt Z is unaffected by Z
  const Real dZ = 1.0e-8;
  _fp->massFractions(p, T, Z + dZ, phase_state, adfsp);
  DualReal Xncg1 = adfsp[0].mass_fraction[1];
  DualReal Yncg1 = adfsp[1].mass_fraction[1];
  _fp->massFractions(p, T, Z - dZ, phase_state, adfsp);
  DualReal Xncg2 = adfsp[0].mass_fraction[1];
  DualReal Yncg2 = adfsp[1].mass_fraction[1];

  ABS_TEST(dXncg_dZ, (Xncg1 - Xncg2).value() / (2.0 * dZ), 1.0e-8);
  ABS_TEST(dYncg_dZ, (Yncg1 - Yncg2).value() / (2.0 * dZ), 1.0e-8);
}

/*
 * Verify calculation of gas density, viscosity, enthalpy and derivatives
 */
TEST_F(PorousFlowWaterNCGADTest, gasProperties)
{
  DualReal p = 1.0e6;
  p.derivatives()[0] = 1.0;

  DualReal T = 350.0;
  T.derivatives()[1] = 1.0;

  FluidStatePhaseEnum phase_state;
  std::vector<ADFluidStateProperties> adfsp(2, ADFluidStateProperties(2));

  // Gas region
  DualReal Z = 0.995;
  Z.derivatives()[2] = 1.0;

  _fp->massFractions(p, T, Z, phase_state, adfsp);
  EXPECT_EQ(phase_state, FluidStatePhaseEnum::GAS);

  // Verify fluid density, viscosity and enthalpy
  _fp->gasProperties(p, T, adfsp);
  DualReal gas_density = adfsp[1].density;
  DualReal gas_viscosity = adfsp[1].viscosity;
  DualReal gas_enthalpy = adfsp[1].enthalpy;

  DualReal density =
      _ncg_fp->rho_from_p_T(Z * p, T) + _water_fp->rho_from_p_T(_water_fp->vaporPressure(T), T);
  DualReal viscosity = Z * _ncg_fp->mu_from_p_T(Z * p, T) +
                       (1.0 - Z) * _water_fp->mu_from_p_T(_water_fp->vaporPressure(T), T);
  DualReal enthalpy = Z * _ncg_fp->h_from_p_T(Z * p, T) +
                      (1.0 - Z) * _water_fp->h_from_p_T(_water_fp->vaporPressure(T), T);

  ABS_TEST(gas_density, density, 1.0e-10);
  ABS_TEST(gas_viscosity, viscosity, 1.0e-10);
  ABS_TEST(gas_enthalpy, enthalpy, 1.0e-10);

  // Verify derivatives
  Real ddensity_dp = gas_density.derivatives()[0];
  Real ddensity_dT = gas_density.derivatives()[1];
  Real ddensity_dZ = gas_density.derivatives()[2];
  Real dviscosity_dp = gas_viscosity.derivatives()[0];
  Real dviscosity_dT = gas_viscosity.derivatives()[1];
  Real dviscosity_dZ = gas_viscosity.derivatives()[2];
  Real denthalpy_dp = gas_enthalpy.derivatives()[0];
  Real denthalpy_dT = gas_enthalpy.derivatives()[1];
  Real denthalpy_dZ = gas_enthalpy.derivatives()[2];

  const Real dp = 1.0e-1;
  _fp->gasProperties(p + dp, T, adfsp);
  DualReal rho1 = adfsp[1].density;
  DualReal mu1 = adfsp[1].viscosity;
  DualReal h1 = adfsp[1].enthalpy;

  _fp->gasProperties(p - dp, T, adfsp);
  DualReal rho2 = adfsp[1].density;
  DualReal mu2 = adfsp[1].viscosity;
  DualReal h2 = adfsp[1].enthalpy;

  REL_TEST(ddensity_dp, (rho1 - rho2).value() / (2.0 * dp), 1.0e-7);
  REL_TEST(dviscosity_dp, (mu1 - mu2).value() / (2.0 * dp), 1.0e-7);
  REL_TEST(denthalpy_dp, (h1 - h2).value() / (2.0 * dp), 5.0e-7);

  const Real dT = 1.0e-3;
  _fp->gasProperties(p, T + dT, adfsp);
  rho1 = adfsp[1].density;
  mu1 = adfsp[1].viscosity;
  h1 = adfsp[1].enthalpy;

  _fp->gasProperties(p, T - dT, adfsp);
  rho2 = adfsp[1].density;
  mu2 = adfsp[1].viscosity;
  h2 = adfsp[1].enthalpy;

  REL_TEST(ddensity_dT, (rho1 - rho2).value() / (2.0 * dT), 1.0e-7);
  REL_TEST(dviscosity_dT, (mu1 - mu2).value() / (2.0 * dT), 1.0e-8);
  REL_TEST(denthalpy_dT, (h1 - h2).value() / (2.0 * dT), 1.0e-8);

  // Note: mass fraction changes with Z
  const Real dZ = 1.0e-8;
  _fp->massFractions(p, T, Z + dZ, phase_state, adfsp);
  _fp->gasProperties(p, T, adfsp);
  rho1 = adfsp[1].density;
  mu1 = adfsp[1].viscosity;
  h1 = adfsp[1].enthalpy;

  _fp->massFractions(p, T, Z - dZ, phase_state, adfsp);
  _fp->gasProperties(p, T, adfsp);
  rho2 = adfsp[1].density;
  mu2 = adfsp[1].viscosity;
  h2 = adfsp[1].enthalpy;

  REL_TEST(ddensity_dZ, (rho1 - rho2).value() / (2.0 * dZ), 1.0e-8);
  REL_TEST(dviscosity_dZ, (mu1 - mu2).value() / (2.0 * dZ), 1.0e-8);
  REL_TEST(denthalpy_dZ, (h1 - h2).value() / (2.0 * dZ), 1.0e-8);

  // Check derivatives in the two phase region as well. Note that the mass fractions
  // vary with pressure and temperature in this region
  Z = 0.45;
  Z.derivatives()[2] = 1.0;

  _fp->massFractions(p, T, Z, phase_state, adfsp);
  EXPECT_EQ(phase_state, FluidStatePhaseEnum::TWOPHASE);

  _fp->gasProperties(p, T, adfsp);
  gas_density = adfsp[1].density;
  gas_viscosity = adfsp[1].viscosity;
  gas_enthalpy = adfsp[1].enthalpy;
  ddensity_dp = gas_density.derivatives()[0];
  ddensity_dT = gas_density.derivatives()[1];
  ddensity_dZ = gas_density.derivatives()[2];
  dviscosity_dp = gas_viscosity.derivatives()[0];
  dviscosity_dT = gas_viscosity.derivatives()[1];
  dviscosity_dZ = gas_viscosity.derivatives()[2];
  denthalpy_dp = gas_enthalpy.derivatives()[0];
  denthalpy_dT = gas_enthalpy.derivatives()[1];
  denthalpy_dZ = gas_enthalpy.derivatives()[2];

  _fp->massFractions(p + dp, T, Z, phase_state, adfsp);
  _fp->gasProperties(p + dp, T, adfsp);
  rho1 = adfsp[1].density;
  mu1 = adfsp[1].viscosity;
  h1 = adfsp[1].enthalpy;

  _fp->massFractions(p - dp, T, Z, phase_state, adfsp);
  _fp->gasProperties(p - dp, T, adfsp);
  rho2 = adfsp[1].density;
  mu2 = adfsp[1].viscosity;
  h2 = adfsp[1].enthalpy;

  REL_TEST(ddensity_dp, (rho1 - rho2).value() / (2.0 * dp), 1.0e-7);
  REL_TEST(dviscosity_dp, (mu1 - mu2).value() / (2.0 * dp), 1.0e-7);
  REL_TEST(denthalpy_dp, (h1 - h2).value() / (2.0 * dp), 1.0e-7);

  _fp->massFractions(p, T + dT, Z, phase_state, adfsp);
  _fp->gasProperties(p, T + dT, adfsp);
  rho1 = adfsp[1].density;
  mu1 = adfsp[1].viscosity;
  h1 = adfsp[1].enthalpy;

  _fp->massFractions(p, T - dT, Z, phase_state, adfsp);
  _fp->gasProperties(p, T - dT, adfsp);
  rho2 = adfsp[1].density;
  mu2 = adfsp[1].viscosity;
  h2 = adfsp[1].enthalpy;

  REL_TEST(ddensity_dT, (rho1 - rho2).value() / (2.0 * dT), 1.0e-7);
  REL_TEST(dviscosity_dT, (mu1 - mu2).value() / (2.0 * dT), 1.0e-8);
  REL_TEST(denthalpy_dT, (h1 - h2).value() / (2.0 * dT), 1.0e-8);

  _fp->massFractions(p, T, Z + dZ, phase_state, adfsp);
  _fp->gasProperties(p, T, adfsp);
  rho1 = adfsp[1].density;
  mu1 = adfsp[1].viscosity;
  h1 = adfsp[1].enthalpy;

  _fp->massFractions(p, T, Z - dZ, phase_state, adfsp);
  _fp->gasProperties(p, T, adfsp);
  rho2 = adfsp[1].density;
  mu2 = adfsp[1].viscosity;
  h2 = adfsp[1].enthalpy;

  ABS_TEST(ddensity_dZ, (rho1 - rho2).value() / (2.0 * dZ), 1.0e-8);
  ABS_TEST(dviscosity_dT, (mu1 - mu2).value() / (2.0 * dZ), 1.0e-7);
  ABS_TEST(denthalpy_dZ, (h1 - h2).value() / (2.0 * dZ), 1.0e-8);
}

/*
 * Verify calculation of liquid density and derivatives
 */
TEST_F(PorousFlowWaterNCGADTest, liquidDensity)
{
  DualReal p = 1.0e6;
  p.derivatives()[0] = 1.0;

  DualReal T = 350.0;
  T.derivatives()[1] = 1.0;

  const DualReal liquid_density = _fp->liquidDensity(p, T);
  std::cout << "liquid density " << liquid_density << std::endl;
}

/*
 * Verify calculation of liquid density, viscosity, enthalpy and derivatives. Note that as
 * density and viscosity don't depend on mass fraction, only the liquid region needs to be
 * tested (the calculations are identical in the two phase region). The enthalpy does depend
 * on mass fraction, so should be tested in the two phase region as well as the liquid region
 */
TEST_F(PorousFlowWaterNCGADTest, liquidProperties)
{
  DualReal p = 1.0e6;
  p.derivatives()[0] = 1.0;

  DualReal T = 350.0;
  T.derivatives()[1] = 1.0;

  FluidStatePhaseEnum phase_state;
  std::vector<ADFluidStateProperties> adfsp(2, ADFluidStateProperties(2));

  // Dummy qp value that isn't needed to test this
  const unsigned int qp = 0;

  // Verify fluid density and viscosity
  _fp->liquidProperties(p, T, qp, adfsp);
  DualReal liquid_density = adfsp[0].density;
  DualReal liquid_viscosity = adfsp[0].viscosity;
  DualReal liquid_enthalpy = adfsp[0].enthalpy;

  DualReal density = _water_fp->rho_from_p_T(p, T);
  DualReal viscosity = _water_fp->mu_from_p_T(p, T);
  DualReal enthalpy = _water_fp->h_from_p_T(p, T);

  ABS_TEST(liquid_density, density, 1.0e-12);
  ABS_TEST(liquid_viscosity, viscosity, 1.0e-12);
  ABS_TEST(liquid_enthalpy, enthalpy, 1.0e-12);

  // Verify derivatives
  Real ddensity_dp = liquid_density.derivatives()[0];
  Real ddensity_dT = liquid_density.derivatives()[1];
  Real ddensity_dZ = liquid_density.derivatives()[2];
  Real dviscosity_dp = liquid_viscosity.derivatives()[0];
  Real dviscosity_dT = liquid_viscosity.derivatives()[1];
  Real dviscosity_dZ = liquid_viscosity.derivatives()[2];
  Real denthalpy_dp = liquid_enthalpy.derivatives()[0];
  Real denthalpy_dT = liquid_enthalpy.derivatives()[1];
  Real denthalpy_dZ = liquid_enthalpy.derivatives()[2];

  const Real dp = 10.0;
  _fp->liquidProperties(p + dp, T, qp, adfsp);
  DualReal rho1 = adfsp[0].density;
  DualReal mu1 = adfsp[0].viscosity;
  DualReal h1 = adfsp[0].enthalpy;

  _fp->liquidProperties(p - dp, T, qp, adfsp);
  DualReal rho2 = adfsp[0].density;
  DualReal mu2 = adfsp[0].viscosity;
  DualReal h2 = adfsp[0].enthalpy;

  REL_TEST(ddensity_dp, (rho1 - rho2) / (2.0 * dp), 1.0e-8);
  REL_TEST(dviscosity_dp, (mu1 - mu2) / (2.0 * dp), 1.0e-8);
  REL_TEST(denthalpy_dp, (h1 - h2) / (2.0 * dp), 1.0e-8);

  const Real dT = 1.0e-4;
  _fp->liquidProperties(p, T + dT, qp, adfsp);
  rho1 = adfsp[0].density;
  mu1 = adfsp[0].viscosity;
  h1 = adfsp[0].enthalpy;

  _fp->liquidProperties(p, T - dT, qp, adfsp);
  rho2 = adfsp[0].density;
  mu2 = adfsp[0].viscosity;
  h2 = adfsp[0].enthalpy;

  REL_TEST(ddensity_dT, (rho1 - rho2) / (2.0 * dT), 1.0e-8);
  REL_TEST(dviscosity_dT, (mu1 - mu2) / (2.0 * dT), 1.0e-8);
  REL_TEST(denthalpy_dT, (h1 - h2) / (2.0 * dT), 1.0e-8);

  DualReal Z = 0.0001;
  Z.derivatives()[2] = 1.0;

  const Real dZ = 1.0e-8;

  _fp->massFractions(p, T, Z, phase_state, adfsp);
  _fp->liquidProperties(p, T, qp, adfsp);
  denthalpy_dZ = adfsp[0].enthalpy.derivatives()[2];

  _fp->massFractions(p, T, Z + dZ, phase_state, adfsp);
  _fp->liquidProperties(p, T, qp, adfsp);
  h1 = adfsp[0].enthalpy;

  _fp->massFractions(p, T, Z - dZ, phase_state, adfsp);
  _fp->liquidProperties(p, T, qp, adfsp);
  h2 = adfsp[0].enthalpy;

  REL_TEST(denthalpy_dZ, (h1 - h2) / (2.0 * dZ), 1.0e-8);

  // Density and viscosity don't depend on Z, so derivatives should be 0
  ABS_TEST(ddensity_dZ, 0.0, 1.0e-12);
  ABS_TEST(dviscosity_dZ, 0.0, 1.0e-12);

  // Check enthalpy calculations in the two phase region as well. Note that the mass fractions
  // vary with pressure and temperature in this region
  Z = 0.45;
  _fp->massFractions(p, T, Z, phase_state, adfsp);
  _fp->liquidProperties(p, T, qp, adfsp);
  denthalpy_dp = adfsp[0].enthalpy.derivatives()[0];
  denthalpy_dT = adfsp[0].enthalpy.derivatives()[1];
  denthalpy_dZ = adfsp[0].enthalpy.derivatives()[2];

  _fp->massFractions(p + dp, T, Z, phase_state, adfsp);
  _fp->liquidProperties(p + dp, T, qp, adfsp);
  h1 = adfsp[0].enthalpy;

  _fp->massFractions(p - dp, T, Z, phase_state, adfsp);
  _fp->liquidProperties(p - dp, T, qp, adfsp);
  h2 = adfsp[0].enthalpy;

  REL_TEST(denthalpy_dp, (h1 - h2) / (2.0 * dp), 1.0e-8);

  _fp->massFractions(p, T + dT, Z, phase_state, adfsp);
  _fp->liquidProperties(p, T + dT, qp, adfsp);
  h1 = adfsp[0].enthalpy;

  _fp->massFractions(p, T - dT, Z, phase_state, adfsp);
  _fp->liquidProperties(p, T - dT, qp, adfsp);
  h2 = adfsp[0].enthalpy;

  REL_TEST(denthalpy_dT, (h1 - h2) / (2.0 * dT), 1.0e-8);

  _fp->massFractions(p, T, Z + dZ, phase_state, adfsp);
  _fp->liquidProperties(p, T, qp, adfsp);
  h1 = adfsp[0].enthalpy;

  _fp->massFractions(p, T, Z - dZ, phase_state, adfsp);
  _fp->liquidProperties(p, T, qp, adfsp);
  h2 = adfsp[0].enthalpy;

  ABS_TEST(denthalpy_dZ, (h1 - h2) / (2.0 * dZ), 1.0e-8);
}

/*
 * Verify calculation of gas saturation and derivatives in the two-phase region
 */
TEST_F(PorousFlowWaterNCGADTest, saturation)
{
  DualReal p = 1.0e6;
  p.derivatives()[0] = 1.0;

  DualReal T = 350.0;
  T.derivatives()[1] = 1.0;

  FluidStatePhaseEnum phase_state;
  std::vector<ADFluidStateProperties> adfsp(2, ADFluidStateProperties(2));

  // Dummy qp value that isn't needed to test this
  const unsigned int qp = 0;

  // In the two-phase region, the mass fractions are the equilibrium values, so
  // a temporary value of Z can be used (as long as it corresponds to the two-phase
  // region)
  DualReal Z = 0.45;
  Z.derivatives()[2] = 1.0;

  _fp->massFractions(p, T, Z, phase_state, adfsp);
  EXPECT_EQ(phase_state, FluidStatePhaseEnum::TWOPHASE);

  // Calculate Z that gives a saturation of 0.25
  DualReal target_gas_saturation = 0.25;
  adfsp[1].saturation = target_gas_saturation;

  // Calculate gas density and liquid density
  _fp->gasProperties(p, T, adfsp);
  _fp->liquidProperties(p, T, qp, adfsp);

  // The mass fraction that corresponds to a gas_saturation = 0.25
  DualReal Zc =
      (target_gas_saturation * adfsp[1].density * adfsp[1].mass_fraction[1] +
       (1.0 - target_gas_saturation) * adfsp[0].density * adfsp[0].mass_fraction[1]) /
      (target_gas_saturation * adfsp[1].density + (1.0 - target_gas_saturation) * adfsp[0].density);

  // Calculate the gas saturation and derivatives
  DualReal gas_saturation = _fp->saturation(p, T, Zc, adfsp);

  REL_TEST(gas_saturation, target_gas_saturation, 1.0e-5);

  // Test the derivatives of gas saturation
  gas_saturation = _fp->saturation(p, T, Z, adfsp);
  const Real dp = 1.0e-1;

  Real dgas_saturation_dp = gas_saturation.derivatives()[0];
  Real dgas_saturation_dT = gas_saturation.derivatives()[1];
  Real dgas_saturation_dZ = gas_saturation.derivatives()[2];

  _fp->massFractions(p + dp, T, Z, phase_state, adfsp);
  _fp->gasProperties(p + dp, T, adfsp);
  DualReal gsat1 = _fp->saturation(p + dp, T, Z, adfsp);

  _fp->massFractions(p - dp, T, Z, phase_state, adfsp);
  _fp->gasProperties(p - dp, T, adfsp);
  DualReal gsat2 = _fp->saturation(p - dp, T, Z, adfsp);

  REL_TEST(dgas_saturation_dp, (gsat1 - gsat2).value() / (2.0 * dp), 1.0e-7);

  const Real dT = 1.0e-4;

  _fp->massFractions(p, T + dT, Z, phase_state, adfsp);
  _fp->gasProperties(p, T + dT, adfsp);
  gsat1 = _fp->saturation(p, T + dT, Z, adfsp);

  _fp->massFractions(p, T - dT, Z, phase_state, adfsp);
  _fp->gasProperties(p, T - dT, adfsp);
  gsat2 = _fp->saturation(p, T - dT, Z, adfsp);

  REL_TEST(dgas_saturation_dT, (gsat1 - gsat2).value() / (2.0 * dT), 1.0e-7);

  const Real dZ = 1.0e-8;

  _fp->massFractions(p, T, Z + dZ, phase_state, adfsp);
  _fp->gasProperties(p, T, adfsp);
  gsat1 = _fp->saturation(p, T, Z + dZ, adfsp);

  _fp->massFractions(p, T, Z - dZ, phase_state, adfsp);
  _fp->gasProperties(p, T, adfsp);
  gsat2 = _fp->saturation(p, T, Z - dZ, adfsp);

  REL_TEST(dgas_saturation_dZ, (gsat1 - gsat2).value() / (2.0 * dZ), 1.0e-7);
}

/*
 * Verify calculation of gas saturation and derivatives in the two-phase region
 */
TEST_F(PorousFlowWaterNCGADTest, twoPhase)
{
  DualReal p = 1.0e6;
  p.derivatives()[0] = 1.0;

  DualReal T = 350.0;
  T.derivatives()[1] = 1.0;

  FluidStatePhaseEnum phase_state;
  std::vector<ADFluidStateProperties> adfsp(2, ADFluidStateProperties(2));

  // Dummy qp value that isn't needed to test this
  const unsigned int qp = 0;

  DualReal Z = 0.45;
  Z.derivatives()[2] = 1.0;

  _fp->massFractions(p, T, Z, phase_state, adfsp);
  EXPECT_EQ(phase_state, FluidStatePhaseEnum::TWOPHASE);

  _fp->twoPhaseProperties(p, T, Z, qp, adfsp);
  DualReal liquid_density = adfsp[0].density;
  DualReal liquid_viscosity = adfsp[0].viscosity;
  DualReal liquid_enthalpy = adfsp[0].enthalpy;
  DualReal gas_saturation = adfsp[1].saturation;
  DualReal gas_density = adfsp[1].density;
  DualReal gas_viscosity = adfsp[1].viscosity;
  DualReal gas_enthalpy = adfsp[1].enthalpy;

  // As all the values are provided by the gasProperties(), liquidProperties() and
  // saturation(), we can simply compare the DualReal results are the same (really
  // only need to make sure that liquid properties are calculated correctly given the saturation)
  std::vector<ADFluidStateProperties> adfsp2(2, ADFluidStateProperties(2));

  _fp->massFractions(p, T, Z, phase_state, adfsp2);
  _fp->gasProperties(p, T, adfsp2);
  adfsp2[1].saturation = _fp->saturation(p, T, Z, adfsp2);
  _fp->liquidProperties(p, T, qp, adfsp2);

  ABS_TEST(gas_density, adfsp2[1].density, 1.0e-12);
  ABS_TEST(gas_viscosity, adfsp2[1].viscosity, 1.0e-12);
  ABS_TEST(gas_enthalpy, adfsp2[1].enthalpy, 1.0e-12);

  ABS_TEST(liquid_density, adfsp2[0].density, 1.0e-12);
  ABS_TEST(liquid_viscosity, adfsp2[0].viscosity, 1.0e-12);
  ABS_TEST(liquid_enthalpy, adfsp2[0].enthalpy, 1.0e-12);

  // const Real dp = 1.0e-1;
  // Real dgas_saturation_dp = gas_saturation.derivatives()[0];
  // Real dgas_saturation_dT = gas_saturation.derivatives()[1];
  // Real dgas_saturation_dZ = gas_saturation.derivatives()[2];
  //
  // _fp->massFractions(p + dp, T, Z, phase_state, adfsp);
  // _fp->twoPhaseProperties(p + dp, T, Z, adfsp);
  // DualReal gsat1 = adfsp[1].saturation;
  //
  // _fp->massFractions(p - dp, T, Z, phase_state, adfsp);
  // _fp->twoPhaseProperties(p - dp, T, Z, adfsp);
  // DualReal gsat2 = adfsp[1].saturation;
  //
  // REL_TEST(dgas_saturation_dp, (gsat1 - gsat2).value() / (2.0 * dp), 1.0e-6);
  // std::cout << "dgas_saturation_dp " << dgas_saturation_dp << ", fd "
  //           << (gsat1 - gsat2).value() / (2.0 * dp) << std::endl;
  // // std::cout << "d2s_dp2 " << fsp[1].d2saturation_dp2 << std::endl;
  // // std::cout << "fd " << (gsat1 - 2.0 * gas_saturation + gsat2) / (dp * dp) << std::endl;
  //
  // // Derivative wrt T
  // const Real dT = 1.0e-4;
  // _fp->massFractions(p, T + dT, Z, phase_state, adfsp);
  // _fp->twoPhaseProperties(p, T + dT, Z, adfsp);
  // gsat1 = adfsp[1].saturation;
  //
  // _fp->massFractions(p, T - dT, Z, phase_state, adfsp);
  // _fp->twoPhaseProperties(p, T - dT, Z, adfsp);
  // gsat2 = adfsp[1].saturation;
  //
  // REL_TEST(dgas_saturation_dT, (gsat1 - gsat2).value() / (2.0 * dT), 1.0e-6);
  // std::cout << "dgas_saturation_dT " << dgas_saturation_dT << ", fd "
  //           << (gsat1 - gsat2).value() / (2.0 * dT) << std::endl;
  //
  // // Derivative wrt Z
  // const Real dZ = 1.0e-8;
  //
  // _fp->massFractions(p, T, Z, phase_state, fsp);
  // _fp->twoPhaseProperties(p, T, Z + dZ, fsp);
  // gsat1 = fsp[1].saturation;
  //
  // _fp->twoPhaseProperties(p, T, Z - dZ, fsp);
  // gsat2 = fsp[1].saturation;
  //
  // REL_TEST(dgas_saturation_dZ, (gsat1 - gsat2) / (2.0 * dZ), 1.0e-6);
  // std::cout << "d2s_dZ2 " << fsp[1].d2saturation_dZ2 << std::endl;
  // std::cout << "fd " << (gsat1 - 2.0 * gas_saturation + gsat2) / (dZ * dZ) << std::endl;
}

/*
 * Verify calculation of total mass fraction given a gas saturation
 */
TEST_F(PorousFlowWaterNCGADTest, totalMassFraction)
{
  const Real p = 1.0e6;
  const Real T = 350.0;
  const Real s = 0.2;
  const Real Xnacl = 0.1;
  const unsigned qp = 0;

  Real Z = _fp->totalMassFraction(p, T, Xnacl, s, qp);

  // Test that the saturation calculated in this fluid state using Z is equal to s
  FluidStatePhaseEnum phase_state;
  std::vector<ADFluidStateProperties> adfsp(2, ADFluidStateProperties(2));

  DualReal pressure = p;
  DualReal temperature = T;
  DualReal z = Z;
  _fp->massFractions(pressure, temperature, z, phase_state, adfsp);
  EXPECT_EQ(phase_state, FluidStatePhaseEnum::TWOPHASE);

  _fp->gasProperties(pressure, temperature, adfsp);
  DualReal gas_saturation = _fp->saturation(pressure, temperature, z, adfsp);
  REL_TEST(gas_saturation, s, 1.0e-5);
}

/*
 * Verify calculation of enthalpy of dissolution. Note: the values calculated compare
 * well by eye to the values presented in Figure 4 of Battistelli et al, "A fluid property
 * module for the TOUGH2 simulator for saline brines with non-condensible gas"
 */
TEST_F(PorousFlowWaterNCGADTest, enthalpyOfDissolution)
{
  // T = 50C
  DualReal T = 323.15;
  T.derivatives()[1] = 1.0;

  // Enthalpy of dissolution of NCG in water
  DualReal hdis = _fp->enthalpyOfDissolution(T);
  REL_TEST(hdis.value(), -3.45731e5, 1.0e-3);

  // T = 350C
  T = 623.15;
  T.derivatives()[1] = 1.0;

  // Enthalpy of dissolution of NCG in water
  hdis = _fp->enthalpyOfDissolution(T);
  REL_TEST(hdis.value(), 1.23423e+06, 1.0e-3);

  // Test the derivative wrt temperature
  const Real dT = 1.0e-4;
  DualReal hdis2 = _fp->enthalpyOfDissolution(T + dT);
  DualReal hdis3 = _fp->enthalpyOfDissolution(T - dT);

  REL_TEST(hdis.derivatives()[1], (hdis2 - hdis3) / (2.0 * dT), 1.0e-5);
}
