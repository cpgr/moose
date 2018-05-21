//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PorousFlowBrineMethaneTest.h"
#include <vector>

/**
 * Verify that the correct name is supplied
 */
TEST_F(PorousFlowBrineMethaneTest, name) { EXPECT_EQ("brine-methane", _fp->fluidStateName()); }

/*
 * Verify calculation of the equilibrium constants and their derivatives wrt temperature
 */
// TEST_F(PorousFlowBrineMethaneTest, equilibriumConstantsMethane)
// {
//   const Real T = 350.0;
//   const Real dT = 1.0e-6;

//   Real K0Methane, dK0Methane_dT;
//   // _fp->equilibriumConstantH2O(T, K0H2O, dK0H2O_dT);
//   _fp->equilibriumConstantMethane(T, K0Methane, dK0Methane_dT);

//   // ABS_TEST("K0H2O", K0H2O, 0.412597711705, 1.0e-10);
//   ABS_TEST("K0Methane", K0Methane, 74.0435888596, 1.0e-10);

//   Real K0Methane_2, dK0Methane_2_dT;
//   // _fp->equilibriumConstantH2O(T + dT, K0H2O_2, dK0H2O_2_dT);
//   _fp->equilibriumConstantMethane(T + dT, K0Methane_2, dK0Methane_2_dT);

//   // Real dK0H2O_dT_fd = (K0H2O_2 - K0H2O) / dT;
//   Real dK0Methane_dT_fd = (K0Methane_2 - K0Methane) / dT;
//   // REL_TEST("dK0H2O_dT", dK0H2O_dT, dK0H2O_dT_fd, 1.0e-6);
//   REL_TEST("dK0Methane_dT", dK0Methane_dT, dK0Methane_dT_fd, 1.0e-6);
// }

/*
 * Verify calculation of the fugacity coefficients of methane 
 * and their derivatives wrt pressure and temperature.
 * Fugacity is checked against values in Table 4, at 4 corners of the table 
 * (combination of pressure/temperature min/max) and at somewhere in the middle.
 * Fugacity derivative is checked against the value at the centre of the table.
 * Pressure range: 1 ~ 10,000 bar
 * Temperature range: 0 ~ 1200 degree celsius
 * The number of digits behind decimal point varies between 3 and 4 in the table.
 * Reference:
 * An equation of state for the CH4-CO2-H2O system: I. Pure systems from 
 * 0 to 1000 degree C and 0 to 8000 bar, 
 * Z. Duan, N, Moller, JH. Weare,
 * Geochimica et Cosmochimica Acta Vol. 56, 1992, pp. 2605-2617.
 */
TEST_F(PorousFlowBrineMethaneTest, fugacityCoefficientMethane)
{
  Real P, t;
  Real phi, dphi_dp, dphi_dT;
  Real phi2, dphi_dp2, dphi_dT2;

  _fp->fugacityCoefficientMethane(barToPascal(1.0), celsiusToKelvin(0.0), phi, dphi_dp, dphi_dT);
  EXPECT_NEAR(0.9977, phi, 1.0e-4) << "Min pressure and temperature";

  // Note that the reference says 29584.790, which I suspect is wrong
  _fp->fugacityCoefficientMethane(barToPascal(1.0e4), celsiusToKelvin(0.0), phi, dphi_dp, dphi_dT);
  EXPECT_NEAR(29548.790, phi, 1.0e-3) << "Max pressure and min temperature";

  _fp->fugacityCoefficientMethane(barToPascal(1.0e4), celsiusToKelvin(1200.0), phi, dphi_dp, dphi_dT);
  EXPECT_NEAR(11.158, phi, 1.0e-3) << "Max pressure and temperature";

  _fp->fugacityCoefficientMethane(barToPascal(1.0), celsiusToKelvin(1200.0), phi, dphi_dp, dphi_dT);
  EXPECT_NEAR(1.0002, phi, 1.0e-4) << "Min pressure and max temperature";

  // Pressure and temperature are somewhere in the middle
  P = 5000.0;
  t = 600.0;
  _fp->fugacityCoefficientMethane(barToPascal(P), celsiusToKelvin(t), phi, dphi_dp, dphi_dT);
  EXPECT_NEAR(6.7646, phi, 1.0e-4) << "Pressure: " << P << " bar, temperature: " << t << " deg C";
  
  const Real dp = 1.0; // Pa
  _fp->fugacityCoefficientMethane(barToPascal(P) + dp, celsiusToKelvin(t), phi2, dphi_dp2, dphi_dT2);
  EXPECT_NEAR((phi2-phi)/dp, dphi_dp, 1.0e-4) << "Derivative wrt. pressure";

  const Real dT = 1.0; // K
  _fp->fugacityCoefficientMethane(barToPascal(P), celsiusToKelvin(t) + dT, phi2, dphi_dp2, dphi_dT2);
  EXPECT_NEAR((phi2-phi)/dT, dphi_dT, 1.0e-4) << "Derivative wrt. temperature";
}


/*
 * Verify calculation of the fugacity coefficients of water and their derivatives wrt
 * pressure and temperature.
 * Values are compared against precomputed values using the same function.
 * Pressure range: 1 ~ 2000 bar
 * Temperature range: 273 ~ 523 Kelvin
 */
TEST_F(PorousFlowBrineMethaneTest, fugacityCoefficientH2O)
{
  const Real tolerance = 1.0e-5; // absolute tolerance

  Real phi, dphi_dp, dphi_dT;
  Real phi2, dphi_dp2, dphi_dT2;

  _fp->fugacityCoefficientH2O(barToPascal(1.0), 273.15, phi, dphi_dp, dphi_dT);
  EXPECT_NEAR(0.98226511105978742, phi, tolerance) << "Min pressure and temperature";

  // Note that the reference says 29584.790, which I suspect is wrong
  _fp->fugacityCoefficientH2O(barToPascal(2000.0), 273.15, phi, dphi_dp, dphi_dT);
  EXPECT_NEAR(6.0379021049117689, phi, tolerance) << "Max pressure and min temperature";

  _fp->fugacityCoefficientH2O(barToPascal(2000.0), 523.15, phi, dphi_dp, dphi_dT);
  EXPECT_NEAR(0.82959242425207924, phi, tolerance) << "Max pressure and temperature";

  _fp->fugacityCoefficientH2O(barToPascal(1.0), 523.15, phi, dphi_dp, dphi_dT);
  EXPECT_NEAR(0.98494687092214472, phi, tolerance) << "Min pressure and max temperature";

  // Pressure and temperature are somewhere in the middle
  Real P = 1000.0;
  Real T = 398.15;
  _fp->fugacityCoefficientH2O(barToPascal(P), T, phi, dphi_dp, dphi_dT);
  EXPECT_NEAR(0.66770302582191587, phi, tolerance) << "Pressure: " << P << " bar, temperature: " << T << " K";
  
  const Real dp = 1.0; // Pa
  _fp->fugacityCoefficientH2O(barToPascal(P) + dp, T, phi2, dphi_dp2, dphi_dT2);
  EXPECT_NEAR((phi2-phi)/dp, dphi_dp, tolerance) << "Derivative wrt. pressure";

  const Real dT = 1.0; // K
  _fp->fugacityCoefficientH2O(barToPascal(P), T + dT, phi2, dphi_dp2, dphi_dT2);
  EXPECT_NEAR((phi2-phi)/dT, dphi_dT, 10.0 * tolerance) << "Derivative wrt. temperature";
}


/* 
 * Check the calculation of methane solubility in liquid by comparing the results with 
 * Table 4 ~ 8 from the reference (shown below).
 * Test for each table (which has a specific value of NaCl concentration) is done on 5 locations:
 * 4 corners (combinations of min/max pressure and temperature) and somewhere in the middle. 
 * Note that the values in the tables have 5 digits behind decimal point.
 * By default, the values are compared using strict tolerance, where the function is required to
 * reproduce the values in the table exactly. In many locations this is not realistic, 
 * so the tolerance is made larger.
 * Pressure range: 1 ~ 2000 bar
 * Temperature range: 273 ~ 523 Kelvin
 * Reference:
 * A thermodynamics model for calculating methane solubility, density and gas phase
 * composition of methane-bearing aqueous fluids from 273 to 523 K and from 1 to 2000 bar,
 * Z. Duan and S. Mao, 
 * Geochimica et Cosmochimica Acta 70, 2006, pp. 3369-3386.
 */
TEST_F(PorousFlowBrineMethaneTest, methaneSolubilityInLiquid)
{
  const Real mh2o = 0.01801528; // Molar mass of H2O in kg/mol
  const Real strictTolerance = 1.0e-5; // The same as values in the tables
  const int n = 25; // 5 tables X 5 locations

  Real P[n];     // Pressure in Bar
  Real T[n];     // Temperature in Kelvin
  Real bnacl[n]; // NaCl molality (mol/kg)
  Real bch4[n];  // CH4 solubility in liquid (mol/kg)
  Real tol[n];   // tolerance

  // Compare with Table 4: pure water
  P[0] = 1.0;
  T[0] = 273.15;
  bnacl[0] = 0.0;
  bch4[0] = 0.00247;
  tol[0] = strictTolerance;

  P[1] = 2000.0;
  T[1] = 333.15;
  bnacl[1] = 0.0;
  bch4[1] = 0.35486;
  tol[1] = strictTolerance;

  P[2] = 2000.0;
  T[2] = 573.15;
  bnacl[2] = 0.0;
  bch4[2] = 2.83679;
  tol[2] = 100.0 * strictTolerance;

  P[3] = 100.0;
  T[3] = 573.15;
  bnacl[3] = 0.0;
  bch4[3] = 0.02009;
  tol[3] = 100.0 * strictTolerance;

  P[4] = 1000.0;
  T[4] = 423.15;
  bnacl[4] = 0.0;
  bch4[4] = 0.43466;
  tol[4] = strictTolerance;

  // Compare with Table 5: NaCl of 1 mol/kg
  P[5] = 1.0;
  T[5] = 273.15;
  bnacl[5] = 1.0;
  bch4[5] = 0.00177;
  tol[5] = strictTolerance;

  P[6] = 2000.0;
  T[6] = 333.15;
  bnacl[6] = 1.0;
  bch4[6] = 0.27167;
  tol[6] = strictTolerance;

  P[7] = 2000.0;
  T[7] = 573.15;
  bnacl[7] = 1.0;
  bch4[7] = 2.36348;
  tol[7] = 100.0 * strictTolerance;

  P[8] = 100.0;
  T[8] = 573.15;
  bnacl[8] = 1.0;
  bch4[8] = 0.04693;
  tol[8] = 10.0 * strictTolerance;

  P[9] = 1000.0;
  T[9] = 423.15;
  bnacl[9] = 1.0;
  bch4[9] = 0.35113;
  tol[9] = strictTolerance;

  // Compare with Table 6: NaCl of 2 mol/kg
  P[10] = 1.0;
  T[10] = 273.15;
  bnacl[10] = 2.0;
  bch4[10] = 0.00127;
  tol[10] = strictTolerance;

  P[11] = 2000.0;
  T[11] = 333.15;
  bnacl[11] = 2.0;
  bch4[11] = 0.2092;
  tol[11] = 10.0 * strictTolerance;

  P[12] = 2000.0;
  T[12] = 573.15;
  bnacl[12] = 2.0;
  bch4[12] = 1.97693;
  tol[12] = 100.0 * strictTolerance;

  P[13] = 100.0;
  T[13] = 573.15;
  bnacl[13] = 2.0;
  bch4[13] = 0.05986;
  tol[13] = 100.0 * strictTolerance;

  P[14] = 1000.0;
  T[14] = 423.15;
  bnacl[14] = 2.0;
  bch4[14] = 0.28535;
  tol[14] = strictTolerance;

  // Compare with Table 7: NaCl of 4 mol/kg
  P[15] = 1.0;
  T[15] = 273.15;
  bnacl[15] = 4.0;
  bch4[15] = 0.00067;
  tol[15] = strictTolerance;

  P[16] = 2000.0;
  T[16] = 303.15;
  bnacl[16] = 4.0;
  bch4[16] = 0.09776;
  tol[16] = strictTolerance;

  P[17] = 2000.0;
  T[17] = 573.15;
  bnacl[17] = 4.0;
  bch4[17] = 1.39884;
  tol[17] = 1000.0 * strictTolerance;

  P[18] = 100.0;
  T[18] = 573.15;
  bnacl[18] = 4.0;
  bch4[18] = 0.06309;
  tol[18] = 1000.0 * strictTolerance;

  P[19] = 1000.0;
  T[19] = 423.15;
  bnacl[19] = 4.0;
  bch4[19] = 0.19183;
  tol[19] = 10.0 * strictTolerance;

  // Compare with Table 8: NaCl of 6 mol/kg
  P[20] = 1.0;
  T[20] = 273.15;
  bnacl[20] = 6.0;
  bch4[20] = 0.00036;
  tol[20] = strictTolerance;

  P[21] = 2000.0;
  T[21] = 303.15;
  bnacl[21] = 6.0;
  bch4[21] = 0.05366;
  tol[21] = strictTolerance;

  P[22] = 2000.0;
  T[22] = 573.15;
  bnacl[22] = 6.0;
  bch4[22] = 1.00244;
  tol[22] = 100.0 * strictTolerance;

  P[23] = 100.0;
  T[23] = 573.15;
  bnacl[23] = 6.0;
  bch4[23] = 0.05306;
  tol[23] = 10.0 * strictTolerance;

  P[24] = 1000.0;
  T[24] = 423.15;
  bnacl[24] = 6.0;
  bch4[24] = 0.13205;
  tol[24] = strictTolerance;
  
  for (int i = 0; i < n; ++i) {
    EXPECT_NEAR(bch4[i], _fp->methaneSolubilityInLiquid(barToPascal(P[i]), T[i], bnacl[i]), tol[i]) << 
      "Entry number " << i << ", P: " << P[i] << " bar, T: " << T[i] << " K, bnacl: " << bnacl[i] << " mol/kg\n";
  }
}


// /*
//  * Verify calculation of the activity coefficient and its derivatives wrt
//  * pressure and temperature
//  */
// TEST_F(PorousFlowBrineCO2Test, activityCoefficient)
// {
//   const Real p = 10.0e6;
//   const Real T = 350.0;
//   const Real xnacl = 0.1;
//   const Real dp = 1.0e-1;
//   const Real dT = 1.0e-6;

//   Real gamma, dgamma_dp, dgamma_dT;
//   _fp->activityCoefficient(p, T, xnacl, gamma, dgamma_dp, dgamma_dT);
//   ABS_TEST("gamma", gamma, 1.43, 1.0e-2);

//   Real gamma_2, dgamma_2_dp, dgamma_2_dT;
//   _fp->activityCoefficient(p + dp, T, xnacl, gamma_2, dgamma_2_dp, dgamma_2_dT);

//   Real dgamma_dp_fd = (gamma_2 - gamma) / dp;
//   REL_TEST("dgamma_dp", dgamma_dp, dgamma_dp_fd, 1.0e-6);

//   _fp->activityCoefficient(p, T + dT, xnacl, gamma_2, dgamma_2_dp, dgamma_2_dT);

//   Real dgamma_dT_fd = (gamma_2 - gamma) / dT;
//   REL_TEST("dgamma_dT", dgamma_dT, dgamma_dT_fd, 1.0e-6);
// }

/*
 * Verify calculation of the activitty coefficient and its derivatives wrt
 * pressure and temperature.
 * Values are compared against precomputed values using the same function.
 * Pressure range: 1 ~ 2000 bar
 * Temperature range: 273 ~ 523 Kelvin
 */
TEST_F(PorousFlowBrineMethaneTest, DISABLED_activityCoefficient)
{
  const Real tolerance = 1.0e-5; // absolute tolerance
  const Real xnacl = 0.1; 
  Real gamma, dgamma_dp, dgamma_dT;
  Real gamma2, dgamma_dp2, dgamma_dT2;

  _fp->activityCoefficient(barToPascal(1.0), 273.15, xnacl, gamma, dgamma_dp, dgamma_dT);
  EXPECT_NEAR(2.7739050439816912, gamma, tolerance) << "Min pressure and temperature";

  // Note that the reference says 29584.790, which I suspect is wrong
  _fp->activityCoefficient(barToPascal(2000.0), 273.15, xnacl, gamma, dgamma_dp, dgamma_dT);
  EXPECT_NEAR(3.4857366863558608, gamma, tolerance) << "Max pressure and min temperature";

  _fp->activityCoefficient(barToPascal(2000.0), 523.15, xnacl, gamma, dgamma_dp, dgamma_dT);
  EXPECT_NEAR(1.7527494176406562, gamma, tolerance) << "Max pressure and temperature";

  _fp->activityCoefficient(barToPascal(1.0), 523.15, xnacl, gamma, dgamma_dp, dgamma_dT);
  EXPECT_NEAR(1.8615175599120284, gamma, tolerance) << "Min pressure and max temperature";

  // Pressure and temperature are somewhere in the middle
  Real P = 1000.0; // bar
  Real T = 398.15; // Kelvin
  _fp->activityCoefficient(barToPascal(P), T, xnacl, gamma, dgamma_dp, dgamma_dT);
  EXPECT_NEAR(1.9302180048382156, gamma, tolerance) << "Pressure: " << P << " bar, temperature: " << T << " K";
  
  const Real dp = 1.0; // Pa
  _fp->activityCoefficient(barToPascal(P) + dp, T, xnacl, gamma2, dgamma_dp2, dgamma_dT2);
  EXPECT_NEAR((gamma2-gamma)/dp, dgamma_dp, tolerance) << "Derivative wrt. pressure";

  const Real dT = 1.0; // K
  _fp->activityCoefficient(barToPascal(P), T + dT, xnacl, gamma2, dgamma_dp2, dgamma_dT2);
  EXPECT_NEAR((gamma2-gamma)/dT, dgamma_dT, tolerance) << "Derivative wrt. temperature";
}


// /*
//  * Verify calculation of the partial density of CO2 and its derivative wrt temperature
//  */
// TEST_F(PorousFlowBrineCO2Test, partialDensity)
// {
//   const Real T = 473.15;
//   const Real dT = 1.0e-6;

//   Real partial_density, dpartial_density_dT;
//   _fp->partialDensityCO2(T, partial_density, dpartial_density_dT);
//   ABS_TEST("partialDensity", partial_density, 893.332, 1.0e-3);

//   Real partial_density_2, dpartial_density_2_dT;
//   _fp->partialDensityCO2(T + dT, partial_density_2, dpartial_density_2_dT);

//   Real dpartial_density_dT_fd = (partial_density_2 - partial_density) / dT;
//   REL_TEST("dpartial_density_dT", dpartial_density_dT, dpartial_density_dT_fd, 1.0e-6);
// }

// /*
//  * Verify calculation of equilibrium mass fraction and derivatives
//  */
// TEST_F(PorousFlowBrineCO2Test, equilibriumMassFraction)
// {
//   const Real p = 1.0e6;
//   const Real T = 350.0;
//   const Real xnacl = 0.1;
//   const Real dp = 1.0e-2;
//   const Real dT = 1.0e-6;

//   Real Xco2, dXco2_dp, dXco2_dT, Yh2o, dYh2o_dp, dYh2o_dT;
//   Real Xco21, dXco21_dp, dXco21_dT, Yh2o1, dYh2o1_dp, dYh2o1_dT;
//   Real Xco22, dXco22_dp, dXco22_dT, Yh2o2, dYh2o2_dp, dYh2o2_dT;
//   _fp->equilibriumMassFractions(p, T, xnacl, Xco2, dXco2_dp, dXco2_dT, Yh2o, dYh2o_dp, dYh2o_dT);
//   _fp->equilibriumMassFractions(
//       p - dp, T, xnacl, Xco21, dXco21_dp, dXco21_dT, Yh2o1, dYh2o1_dp, dYh2o1_dT);
//   _fp->equilibriumMassFractions(
//       p + dp, T, xnacl, Xco22, dXco22_dp, dXco22_dT, Yh2o2, dYh2o2_dp, dYh2o2_dT);

//   Real dXco2_dp_fd = (Xco22 - Xco21) / (2.0 * dp);
//   Real dYh2o_dp_fd = (Yh2o2 - Yh2o1) / (2.0 * dp);

//   REL_TEST("dXco2_dp", dXco2_dp, dXco2_dp_fd, 1.0e-6);
//   REL_TEST("dYh2o_dp", dYh2o_dp, dYh2o_dp_fd, 1.0e-6);

//   _fp->equilibriumMassFractions(
//       p, T - dT, xnacl, Xco21, dXco21_dp, dXco21_dT, Yh2o1, dYh2o1_dp, dYh2o1_dT);
//   _fp->equilibriumMassFractions(
//       p, T + dT, xnacl, Xco22, dXco22_dp, dXco22_dT, Yh2o2, dYh2o2_dp, dYh2o2_dT);

//   Real dXco2_dT_fd = (Xco22 - Xco21) / (2.0 * dT);
//   Real dYh2o_dT_fd = (Yh2o2 - Yh2o1) / (2.0 * dT);

//   REL_TEST("dXco2_dT", dXco2_dT, dXco2_dT_fd, 1.0e-6);
//   REL_TEST("dYh2o_dT", dYh2o_dT, dYh2o_dT_fd, 1.0e-6);
// }

// /*
//  * Verify calculation of actual mass fraction and derivatives depending on value of
//  * total mass fraction
//  */
// TEST_F(PorousFlowBrineCO2Test, MassFraction)
// {
//   const Real p = 1.0e6;
//   const Real T = 350.0;
//   const Real xnacl = 0.1;

//   FluidStatePhaseEnum phase_state;
//   std::vector<FluidStateProperties> fsp(2, FluidStateProperties(2));

//   // Liquid region
//   Real z = 0.0001;
//   _fp->massFractions(p, T, xnacl, z, phase_state, fsp);
//   EXPECT_EQ(phase_state, FluidStatePhaseEnum::LIQUID);

//   // Verfify mass fraction values
//   Real Xco2 = fsp[0].mass_fraction[1];
//   Real Yco2 = fsp[1].mass_fraction[1];
//   Real Xh2o = fsp[0].mass_fraction[0];
//   Real Yh2o = fsp[1].mass_fraction[0];
//   ABS_TEST("Xco2", Xco2, z, 1.0e-8);
//   ABS_TEST("Yco2", Yco2, 0.0, 1.0e-8);
//   ABS_TEST("Xh2o", Xh2o, 1.0 - z, 1.0e-8);
//   ABS_TEST("Yh2o", Yh2o, 0.0, 1.0e-8);

//   // Verify derivatives
//   Real dXco2_dp = fsp[0].dmass_fraction_dp[1];
//   Real dXco2_dT = fsp[0].dmass_fraction_dT[1];
//   Real dXco2_dz = fsp[0].dmass_fraction_dz[1];
//   Real dYco2_dp = fsp[1].dmass_fraction_dp[1];
//   Real dYco2_dT = fsp[1].dmass_fraction_dT[1];
//   Real dYco2_dz = fsp[1].dmass_fraction_dz[1];
//   ABS_TEST("dXco2_dp", dXco2_dp, 0.0, 1.0e-8);
//   ABS_TEST("dXco2_dT", dXco2_dT, 0.0, 1.0e-8);
//   ABS_TEST("dXco2_dz", dXco2_dz, 1.0, 1.0e-8);
//   ABS_TEST("dYco2_dp", dYco2_dp, 0.0, 1.0e-8);
//   ABS_TEST("dYco2_dT", dYco2_dT, 0.0, 1.0e-8);
//   ABS_TEST("dYco2_dz", dYco2_dz, 0.0, 1.0e-8);

//   // Gas region
//   z = 0.995;
//   _fp->massFractions(p, T, xnacl, z, phase_state, fsp);
//   EXPECT_EQ(phase_state, FluidStatePhaseEnum::GAS);

//   // Verfify mass fraction values
//   Xco2 = fsp[0].mass_fraction[1];
//   Yco2 = fsp[1].mass_fraction[1];
//   Xh2o = fsp[0].mass_fraction[0];
//   Yh2o = fsp[1].mass_fraction[0];
//   ABS_TEST("Xco2", Xco2, 0.0, 1.0e-8);
//   ABS_TEST("Yco2", Yco2, z, 1.0e-8);
//   ABS_TEST("Xh2o", Xh2o, 0.0, 1.0e-8);
//   ABS_TEST("Yh2o", Yh2o, 1.0 - z, 1.0e-8);

//   // Verify derivatives
//   dXco2_dp = fsp[0].dmass_fraction_dp[1];
//   dXco2_dT = fsp[0].dmass_fraction_dT[1];
//   dXco2_dz = fsp[0].dmass_fraction_dz[1];
//   dYco2_dp = fsp[1].dmass_fraction_dp[1];
//   dYco2_dT = fsp[1].dmass_fraction_dT[1];
//   dYco2_dz = fsp[1].dmass_fraction_dz[1];
//   ABS_TEST("dXco2_dp", dXco2_dp, 0.0, 1.0e-8);
//   ABS_TEST("dXco2_dT", dXco2_dT, 0.0, 1.0e-8);
//   ABS_TEST("dXco2_dz", dXco2_dz, 0.0, 1.0e-8);
//   ABS_TEST("dYco2_dp", dYco2_dp, 0.0, 1.0e-8);
//   ABS_TEST("dYco2_dT", dYco2_dT, 0.0, 1.0e-8);
//   ABS_TEST("dYco2_dz", dYco2_dz, 1.0, 1.0e-8);

//   // Two phase region. In this region, the mass fractions and derivatives can
//   //  be verified using the equilibrium mass fraction derivatives that have
//   // been verified above
//   z = 0.45;
//   _fp->massFractions(p, T, xnacl, z, phase_state, fsp);
//   EXPECT_EQ(phase_state, FluidStatePhaseEnum::TWOPHASE);

//   // Equilibrium mass fractions and derivatives
//   Real Xco2_eq, dXco2_dp_eq, dXco2_dT_eq, Yh2o_eq, dYh2o_dp_eq, dYh2o_dT_eq;
//   _fp->equilibriumMassFractions(
//       p, T, xnacl, Xco2_eq, dXco2_dp_eq, dXco2_dT_eq, Yh2o_eq, dYh2o_dp_eq, dYh2o_dT_eq);

//   // Verfify mass fraction values
//   Xco2 = fsp[0].mass_fraction[1];
//   Yco2 = fsp[1].mass_fraction[1];
//   Xh2o = fsp[0].mass_fraction[0];
//   Yh2o = fsp[1].mass_fraction[0];
//   ABS_TEST("Xco2", Xco2, Xco2_eq, 1.0e-8);
//   ABS_TEST("Yco2", Yco2, 1.0 - Yh2o_eq, 1.0e-8);
//   ABS_TEST("Xh2o", Xh2o, 1.0 - Xco2_eq, 1.0e-8);
//   ABS_TEST("Yh2o", Yh2o, Yh2o_eq, 1.0e-8);

//   // Verify derivatives wrt p and T
//   dXco2_dp = fsp[0].dmass_fraction_dp[1];
//   dXco2_dT = fsp[0].dmass_fraction_dT[1];
//   dXco2_dz = fsp[0].dmass_fraction_dz[1];
//   dYco2_dp = fsp[1].dmass_fraction_dp[1];
//   dYco2_dT = fsp[1].dmass_fraction_dT[1];
//   dYco2_dz = fsp[1].dmass_fraction_dz[1];
//   ABS_TEST("dXco2_dp", dXco2_dp, dXco2_dp_eq, 1.0e-8);
//   ABS_TEST("dXco2_dT", dXco2_dT, dXco2_dT_eq, 1.0e-8);
//   ABS_TEST("dXco2_dz", dXco2_dz, 0.0, 1.0e-8);
//   ABS_TEST("dYco2_dp", dYco2_dp, -dYh2o_dp_eq, 1.0e-8);
//   ABS_TEST("dYco2_dT", dYco2_dT, -dYh2o_dT_eq, 1.0e-8);
//   ABS_TEST("dYco2_dz", dYco2_dz, 0.0, 1.0e-8);

//   // Use finite differences to verify derivative wrt z is unaffected by z
//   const Real dz = 1.0e-8;
//   _fp->massFractions(p, T, xnacl, z + dz, phase_state, fsp);
//   Real Xco21 = fsp[0].mass_fraction[1];
//   Real Yco21 = fsp[1].mass_fraction[1];
//   _fp->massFractions(p, T, xnacl, z - dz, phase_state, fsp);
//   Real Xco22 = fsp[0].mass_fraction[1];
//   Real Yco22 = fsp[1].mass_fraction[1];

//   ABS_TEST("dXco2_dz", dXco2_dz, (Xco21 - Xco22) / (2.0 * dz), 1.0e-8);
//   ABS_TEST("dXYco2_dz", dYco2_dz, (Yco21 - Yco22) / (2.0 * dz), 1.0e-8);
// }

// /*
//  * Verify calculation of gas density and viscosity, and derivatives. Note that as
//  * these properties don't depend on mass fraction, only the gas region needs to be
//  * tested (the calculations are identical in the two phase region)
//  */
// TEST_F(PorousFlowBrineCO2Test, gasProperties)
// {
//   const Real p = 1.0e6;
//   const Real T = 350.0;
//   const Real xnacl = 0.1;

//   FluidStatePhaseEnum phase_state;
//   std::vector<FluidStateProperties> fsp(2, FluidStateProperties(2));

//   // Gas region
//   Real z = 0.995;
//   _fp->massFractions(p, T, xnacl, z, phase_state, fsp);
//   EXPECT_EQ(phase_state, FluidStatePhaseEnum::GAS);

//   // Verify fluid density and viscosity
//   _fp->gasProperties(p, T, fsp);
//   Real gas_density = fsp[1].density;
//   Real gas_viscosity = fsp[1].viscosity;
//   Real density = _co2_fp->rho(p, T);
//   Real viscosity = _co2_fp->mu(p, T);

//   ABS_TEST("gas density", gas_density, density, 1.0e-8);
//   ABS_TEST("gas viscosity", gas_viscosity, viscosity, 1.0e-8);

//   // Verify derivatives
//   Real ddensity_dp = fsp[1].ddensity_dp;
//   Real ddensity_dT = fsp[1].ddensity_dT;
//   Real ddensity_dz = fsp[1].ddensity_dz;
//   Real dviscosity_dp = fsp[1].dviscosity_dp;
//   Real dviscosity_dT = fsp[1].dviscosity_dT;
//   Real dviscosity_dz = fsp[1].dviscosity_dz;

//   const Real dp = 1.0e-1;
//   _fp->gasProperties(p + dp, T, fsp);
//   Real rho1 = fsp[1].density;
//   Real mu1 = fsp[1].viscosity;

//   _fp->gasProperties(p - dp, T, fsp);
//   Real rho2 = fsp[1].density;
//   Real mu2 = fsp[1].viscosity;

//   REL_TEST("ddensity_dp", ddensity_dp, (rho1 - rho2) / (2.0 * dp), 1.0e-6);
//   REL_TEST("dviscosity_dp", dviscosity_dp, (mu1 - mu2) / (2.0 * dp), 1.0e-6);

//   const Real dT = 1.0e-3;
//   _fp->gasProperties(p, T + dT, fsp);
//   rho1 = fsp[1].density;
//   mu1 = fsp[1].viscosity;
//   _fp->gasProperties(p, T - dT, fsp);
//   rho2 = fsp[1].density;
//   mu2 = fsp[1].viscosity;

//   REL_TEST("ddensity_dT", ddensity_dT, (rho1 - rho2) / (2.0 * dT), 1.0e-6);
//   REL_TEST("dviscosity_dT", dviscosity_dT, (mu1 - mu2) / (2.0 * dT), 1.0e-6);

//   // Note: mass fraction changes with z
//   const Real dz = 1.0e-8;
//   _fp->massFractions(p, T, xnacl, z + dz, phase_state, fsp);
//   _fp->gasProperties(p, T, fsp);
//   rho1 = fsp[1].density;
//   mu1 = fsp[1].viscosity;

//   _fp->massFractions(p, T, xnacl, z - dz, phase_state, fsp);
//   _fp->gasProperties(p, T, fsp);
//   rho2 = fsp[1].density;
//   mu2 = fsp[1].viscosity;

//   ABS_TEST("ddensity_dz", ddensity_dz, (rho1 - rho2) / (2.0 * dz), 1.0e-8);
//   ABS_TEST("dviscosity_dz", dviscosity_dz, (mu1 - mu2) / (2.0 * dz), 1.0e-8);
// }

// /*
//  * Verify calculation of liquid density and viscosity, and derivatives
//  */
// TEST_F(PorousFlowBrineCO2Test, liquidProperties)
// {
//   const Real p = 1.0e6;
//   const Real T = 350.0;
//   const Real xnacl = 0.1;

//   FluidStatePhaseEnum phase_state;
//   std::vector<FluidStateProperties> fsp(2, FluidStateProperties(2));

//   // Liquid region
//   Real z = 0.0001;
//   _fp->massFractions(p, T, xnacl, z, phase_state, fsp);
//   EXPECT_EQ(phase_state, FluidStatePhaseEnum::LIQUID);

//   // Verify fluid density and viscosity
//   _fp->liquidProperties(p, T, xnacl, fsp);

//   Real liquid_density = fsp[0].density;
//   Real liquid_viscosity = fsp[0].viscosity;

//   Real co2_partial_density, dco2_partial_density_dT;
//   _fp->partialDensityCO2(T, co2_partial_density, dco2_partial_density_dT);
//   Real brine_density = _brine_fp->rho(p, T, xnacl);

//   Real density = 1.0 / (z / co2_partial_density + (1.0 - z) / brine_density);

//   Real water_density = _water_fp->rho(p, T);
//   Real viscosity = _brine_fp->mu_from_rho_T(water_density, T, xnacl);

//   ABS_TEST("liquid density", liquid_density, density, 1.0e-12);
//   ABS_TEST("liquid viscosity", liquid_viscosity, viscosity, 1.0e-12);

//   // Verify fluid density and viscosity derivatives
//   Real ddensity_dp = fsp[0].ddensity_dp;
//   Real ddensity_dT = fsp[0].ddensity_dT;
//   Real ddensity_dz = fsp[0].ddensity_dz;
//   Real dviscosity_dp = fsp[0].dviscosity_dp;
//   Real dviscosity_dT = fsp[0].dviscosity_dT;
//   Real dviscosity_dz = fsp[0].dviscosity_dz;

//   const Real dp = 1.0;
//   _fp->liquidProperties(p + dp, T, xnacl, fsp);
//   Real rho1 = fsp[0].density;
//   Real mu1 = fsp[0].viscosity;

//   _fp->liquidProperties(p - dp, T, xnacl, fsp);
//   Real rho2 = fsp[0].density;
//   Real mu2 = fsp[0].viscosity;

//   REL_TEST("ddensity_dp", ddensity_dp, (rho1 - rho2) / (2.0 * dp), 1.0e-4);
//   REL_TEST("dviscosity_dp", dviscosity_dp, (mu1 - mu2) / (2.0 * dp), 1.0e-5);

//   const Real dT = 1.0e-4;
//   _fp->liquidProperties(p, T + dT, xnacl, fsp);
//   rho1 = fsp[0].density;
//   mu1 = fsp[0].viscosity;

//   _fp->liquidProperties(p, T - dT, xnacl, fsp);
//   rho2 = fsp[0].density;
//   mu2 = fsp[0].viscosity;

//   REL_TEST("ddensity_dT", ddensity_dT, (rho1 - rho2) / (2.0 * dT), 1.0e-6);
//   REL_TEST("dviscosity_dT", dviscosity_dT, (mu1 - mu2) / (2.0 * dT), 1.0e-6);

//   const Real dz = 1.0e-8;
//   _fp->massFractions(p, T, xnacl, z + dz, phase_state, fsp);
//   _fp->liquidProperties(p, T, xnacl, fsp);
//   rho1 = fsp[0].density;
//   mu1 = fsp[0].viscosity;

//   _fp->massFractions(p, T, xnacl, z - dz, phase_state, fsp);
//   _fp->liquidProperties(p, T, xnacl, fsp);
//   rho2 = fsp[0].density;
//   mu2 = fsp[0].viscosity;

//   REL_TEST("ddensity_dz", ddensity_dz, (rho1 - rho2) / (2.0 * dz), 1.0e-6);
//   ABS_TEST("dviscosity_dz", dviscosity_dz, (mu1 - mu2) / (2.0 * dz), 1.0e-6);

//   // Two-phase region
//   z = 0.045;
//   _fp->massFractions(p, T, xnacl, z, phase_state, fsp);
//   EXPECT_EQ(phase_state, FluidStatePhaseEnum::TWOPHASE);

//   // Verify fluid density and viscosity derivatives
//   _fp->liquidProperties(p, T, xnacl, fsp);

//   ddensity_dp = fsp[0].ddensity_dp;
//   ddensity_dT = fsp[0].ddensity_dT;
//   ddensity_dz = fsp[0].ddensity_dz;
//   dviscosity_dp = fsp[0].dviscosity_dp;
//   dviscosity_dT = fsp[0].dviscosity_dT;
//   dviscosity_dz = fsp[0].dviscosity_dz;

//   _fp->massFractions(p + dp, T, xnacl, z, phase_state, fsp);
//   _fp->liquidProperties(p + dp, T, xnacl, fsp);
//   rho1 = fsp[0].density;
//   mu1 = fsp[0].viscosity;

//   _fp->massFractions(p - dp, T, xnacl, z, phase_state, fsp);
//   _fp->liquidProperties(p - dp, T, xnacl, fsp);
//   rho2 = fsp[0].density;
//   mu2 = fsp[0].viscosity;

//   REL_TEST("ddensity_dp", ddensity_dp, (rho1 - rho2) / (2.0 * dp), 1.0e-4);
//   REL_TEST("dviscosity_dp", dviscosity_dp, (mu1 - mu2) / (2.0 * dp), 1.0e-5);

//   _fp->massFractions(p, T + dT, xnacl, z, phase_state, fsp);
//   _fp->liquidProperties(p, T + dT, xnacl, fsp);
//   rho1 = fsp[0].density;
//   mu1 = fsp[0].viscosity;

//   _fp->massFractions(p, T - dT, xnacl, z, phase_state, fsp);
//   _fp->liquidProperties(p, T - dT, xnacl, fsp);
//   rho2 = fsp[0].density;
//   mu2 = fsp[0].viscosity;

//   REL_TEST("ddensity_dT", ddensity_dT, (rho1 - rho2) / (2.0 * dT), 1.0e-6);
//   REL_TEST("dviscosity_dT", dviscosity_dT, (mu1 - mu2) / (2.0 * dT), 1.0e-6);

//   _fp->massFractions(p, T, xnacl, z + dz, phase_state, fsp);
//   _fp->liquidProperties(p, T, xnacl, fsp);
//   rho1 = fsp[0].density;
//   mu1 = fsp[0].viscosity;

//   _fp->massFractions(p, T, xnacl, z - dz, phase_state, fsp);
//   _fp->liquidProperties(p, T, xnacl, fsp);
//   rho2 = fsp[0].density;
//   mu2 = fsp[0].viscosity;

//   ABS_TEST("ddensity_dz", ddensity_dz, (rho1 - rho2) / (2.0 * dz), 1.0e-6);
//   ABS_TEST("dviscosity_dz", dviscosity_dz, (mu1 - mu2) / (2.0 * dz), 1.0e-6);
// }

// /*
//  * Verify calculation of gas saturation and derivatives in the two-phase region
//  */
// TEST_F(PorousFlowBrineCO2Test, saturationTwoPhase)
// {
//   const Real p = 1.0e6;
//   const Real T = 350.0;
//   const Real xnacl = 0.1;

//   FluidStatePhaseEnum phase_state;
//   std::vector<FluidStateProperties> fsp(2, FluidStateProperties(2));

//   // In the two-phase region, the mass fractions are the equilibrium values, so
//   // a temporary value of z can be used (as long as it corresponds to the two-phase
//   // region)
//   Real z = 0.45;
//   _fp->massFractions(p, T, xnacl, z, phase_state, fsp);
//   EXPECT_EQ(phase_state, FluidStatePhaseEnum::TWOPHASE);

//   // Calculate z that gives a saturation of 0.25
//   Real gas_saturation = 0.25;
//   Real liquid_pressure = p + _pc->capillaryPressure(1.0 - gas_saturation);
//   // Calculate gas density and liquid density
//   _fp->gasProperties(p, T, fsp);
//   _fp->liquidProperties(liquid_pressure, T, xnacl, fsp);

//   // The mass fraction that corresponds to a gas_saturation = 0.25
//   z = (gas_saturation * fsp[1].density * fsp[1].mass_fraction[1] +
//        (1.0 - gas_saturation) * fsp[0].density * fsp[0].mass_fraction[1]) /
//       (gas_saturation * fsp[1].density + (1.0 - gas_saturation) * fsp[0].density);

//   // Calculate the gas saturation and derivatives
//   _fp->saturationTwoPhase(p, T, xnacl, z, fsp);

//   ABS_TEST("gas saturation", fsp[1].saturation, gas_saturation, 1.0e-8);

//   // Test the derivatives
//   const Real dp = 1.0e-1;
//   gas_saturation = fsp[1].saturation;
//   Real dgas_saturation_dp = fsp[1].dsaturation_dp;
//   Real dgas_saturation_dT = fsp[1].dsaturation_dT;
//   Real dgas_saturation_dz = fsp[1].dsaturation_dz;

//   _fp->massFractions(p + dp, T, xnacl, z, phase_state, fsp);
//   _fp->gasProperties(p + dp, T, fsp);
//   _fp->saturationTwoPhase(p + dp, T, xnacl, z, fsp);
//   Real gsat1 = fsp[1].saturation;

//   _fp->massFractions(p - dp, T, xnacl, z, phase_state, fsp);
//   _fp->gasProperties(p - dp, T, fsp);
//   _fp->saturationTwoPhase(p - dp, T, xnacl, z, fsp);
//   Real gsat2 = fsp[1].saturation;

//   REL_TEST("dgas_saturation_dp", dgas_saturation_dp, (gsat1 - gsat2) / (2.0 * dp), 1.0e-6);

//   // Derivative wrt T
//   const Real dT = 1.0e-4;
//   _fp->massFractions(p, T + dT, xnacl, z, phase_state, fsp);
//   _fp->gasProperties(p, T + dT, fsp);
//   _fp->saturationTwoPhase(p, T + dT, xnacl, z, fsp);
//   gsat1 = fsp[1].saturation;

//   _fp->massFractions(p, T - dT, xnacl, z, phase_state, fsp);
//   _fp->gasProperties(p, T - dT, fsp);
//   _fp->saturationTwoPhase(p, T - dT, xnacl, z, fsp);
//   gsat2 = fsp[1].saturation;

//   REL_TEST("dgas_saturation_dT", dgas_saturation_dT, (gsat1 - gsat2) / (2.0 * dT), 1.0e-6);

//   // Derivative wrt z
//   const Real dz = 1.0e-8;

//   _fp->massFractions(p, T, xnacl, z, phase_state, fsp);
//   _fp->gasProperties(p, T, fsp);
//   _fp->saturationTwoPhase(p, T, xnacl, z + dz, fsp);
//   gsat1 = fsp[1].saturation;

//   _fp->saturationTwoPhase(p, T, xnacl, z - dz, fsp);
//   gsat2 = fsp[1].saturation;

//   REL_TEST("dgas_saturation_dz", dgas_saturation_dz, (gsat1 - gsat2) / (2.0 * dz), 1.0e-6);
// }

// /*
//  * Verify calculation of total mass fraction given a gas saturation
//  */
// TEST_F(PorousFlowBrineCO2Test, totalMassFraction)
// {
//   const Real p = 1.0e6;
//   const Real T = 350.0;
//   const Real xnacl = 0.1;
//   const Real s = 0.2;

//   Real z = _fp->totalMassFraction(p, T, xnacl, s);

//   // Test that the saturation calculated in this fluid state using z is equal to s
//   FluidStatePhaseEnum phase_state;
//   std::vector<FluidStateProperties> fsp(2, FluidStateProperties(2));

//   _fp->massFractions(p, T, xnacl, z, phase_state, fsp);
//   EXPECT_EQ(phase_state, FluidStatePhaseEnum::TWOPHASE);

//   _fp->gasProperties(p, T, fsp);
//   Real liquid_pressure = p + _pc->capillaryPressure(1.0 - s);
//   _fp->liquidProperties(liquid_pressure, T, xnacl, fsp);
//   _fp->saturationTwoPhase(p, T, xnacl, z, fsp);
//   ABS_TEST("gas saturation", fsp[1].saturation, s, 1.0e-8);
// }
