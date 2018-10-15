//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PorousFlowBrineCO2CH4H2S.h"
#include "BrineFluidProperties.h"
#include "SinglePhaseFluidPropertiesPT.h"
#include "MathUtils.h"
#include "Conversion.h"

registerMooseObject("PorousFlowApp", PorousFlowBrineCO2CH4H2S);

template <>
InputParameters
validParams<PorousFlowBrineCO2CH4H2S>()
{
  InputParameters params = validParams<PorousFlowFluidStateBase>();
  params.addRequiredParam<UserObjectName>("brine_fp", "The name of the user object for brine");
  params.addRequiredParam<UserObjectName>("co2_fp", "The name of the user object for CO2");
  params.addRequiredParam<UserObjectName>("ch4_fp", "The name of the user object for CH4");
  params.addClassDescription("Fluid state class for brine and CO2");
  return params;
}

PorousFlowBrineCO2CH4H2S::PorousFlowBrineCO2CH4H2S(const InputParameters & parameters)
  : PorousFlowFluidStateBase(parameters),
    _co2_component(1),
    _ch4_component(2),
    _h2s_component(3),
    _salt_component(4),
    _zco2_idx(0),
    _zch4_idx(1),
    _zh2s_idx(2),
    _brine_fp(getUserObject<BrineFluidProperties>("brine_fp")),
    _co2_fp(getUserObject<SinglePhaseFluidPropertiesPT>("co2_fp")),
    _ch4_fp(getUserObject<SinglePhaseFluidPropertiesPT>("ch4_fp")),
    _water_fp(_brine_fp.getComponent(BrineFluidProperties::WATER)),
    _Mh2o(_brine_fp.molarMassH2O()),
    _invMh2o(1.0 / _Mh2o),
    _Mco2(_co2_fp.molarMass()),
    _Mh2s(0.03408),
    _Mch4(0.01604),
    _Mnacl(_brine_fp.molarMassNaCl())
{
  // Set the number of phases and components, and their indexes
  _num_phases = 2;
  _num_components = 5;
  _num_zvars = 3;
  _gas_phase_number = 1 - _aqueous_phase_number;
}

std::string
PorousFlowBrineCO2CH4H2S::fluidStateName() const
{
  return "brine-co2-ch4-h2s";
}

void
PorousFlowBrineCO2CH4H2S::equilibriumMassFractions(Real pressure,
                                                   Real temperature,
                                                   Real Xnacl,
                                                   Real Yco2,
                                                   Real Ych4,
                                                   Real Yh2s,
                                                   Real & Xco2,
                                                   Real & Xch4,
                                                   Real & Xh2s,
                                                   Real & Yh2o) const
{
  // Calculate mole fraction of each component in the gas phase
  const Real Yh2otmp = 1.0 - Yco2 - Ych4 - Yh2s;
  Real denominator = Yco2 / _Mco2 + Yh2s / _Mh2s + Ych4 / _Mch4 + Yh2otmp / _Mh2o;
  const Real yco2 = Yco2 / _Mco2 / denominator;
  const Real yh2s = Yh2s / _Mh2s / denominator;
  const Real ych4 = Ych4 / _Mch4 / denominator;

  // Solve for the mole fractions in the liquid phase
  Real xco2, xch4, xh2s, yh2o;
  equilibriumMoleFractions(pressure, temperature, Xnacl, yco2, ych4, yh2s, xco2, xch4, xh2s, yh2o);

  denominator = xco2 * _Mco2 + xh2s * _Mh2s + xch4 * _Mch4 + (1.0 - xco2 - xh2s - xch4) * _Mh2o;
  Xco2 = xco2 * _Mco2 / denominator;
  Xh2s = xh2s * _Mh2s / denominator;
  Xch4 = xch4 * _Mch4 / denominator;

  Yh2o = yh2o * _Mh2o / (yco2 * _Mco2 + yh2s * _Mh2s + ych4 * _Mch4 + yh2o * _Mh2o);
}

void
PorousFlowBrineCO2CH4H2S::equilibriumMoleFractions(Real pressure,
                                                   Real temperature,
                                                   Real /*Xnacl*/,
                                                   Real yco2,
                                                   Real ych4,
                                                   Real yh2s,
                                                   Real & xco2,
                                                   Real & xch4,
                                                   Real & xh2s,
                                                   Real & yh2o) const
{
  // Equilibrium mass fraction of CO2,CH4 and H2S in liquid and H2O in gas phases (Zirrahi et. al
  // 2013)
  Real y_hs, y_co, y_ch;

  y_hs = yh2s;
  y_co = yco2;
  y_ch = ych4;

  const Real t = temperature;       // k
  const Real p = pressure * 1.0e-5; // bar

  // Intermolecular intraction parameters(ap_ij)
  const Real ap_ch_co = 3480693.2119225934;
  const Real ap_ch_hs = 4957015.776228365;
  const Real ap_co_hs = 1413558.6307468852;
  const Real ap_hs_co = -2728833.965482667;

  // Binary intraction parameters(Kij)
  const Real k_co_ch = 0.15196821446633596;
  const Real k_co_hs = -0.037849313973062265;
  const Real k_hs_ch = 1.5716638150922622;

  //  gas constant
  const Real r = 83.14;

  //  critical temperature and pressure
  const Real tc_ch = 190.6;
  const Real tc_hs = 373.53;
  const Real tc_co = 31.05 + 273.15;

  const Real pc_ch = 46.41;
  const Real pc_hs = 89.63;
  const Real pc_co = 73.825;

  // acentric factors
  const Real w_ch = 0.012;
  const Real w_co = 0.224;
  const Real w_hs = 0.094;

  // calculation of bi
  const Real b_ch = 0.0778 * r * tc_ch / pc_ch;
  const Real b_co = 0.0778 * r * tc_co / pc_co;
  const Real b_hs = 0.0778 * r * tc_hs / pc_hs;

  // calculation of beta
  const Real bet_ch = 0.37464 + 1.54226 * w_ch - 0.26992 * std::pow(w_ch, 2);
  const Real bet_co = 0.37464 + 1.54226 * w_co - 0.26992 * std::pow(w_co, 2);
  const Real bet_hs = 0.37464 + 1.54226 * w_hs - 0.26992 * std::pow(w_hs, 2);

  // calculation of alpha_t
  const Real at_ch = std::pow((1 + bet_ch * (1 - std::pow((t / tc_ch), 0.5))), 2);
  const Real at_co = std::pow((1 + bet_co * (1 - std::pow((t / tc_co), 0.5))), 2);
  const Real at_hs = std::pow((1 + bet_hs * (1 - std::pow((t / tc_hs), 0.5))), 2);

  // calculation of a_tc
  const Real atc_ch = 0.45724 * std::pow((r * tc_ch), 2) / pc_ch;
  const Real atc_co = 0.45724 * std::pow((r * tc_co), 2) / pc_co;
  const Real atc_hs = 0.45724 * std::pow((r * tc_hs), 2) / pc_hs;

  // calculation of ai
  const Real a_ch = atc_ch * at_ch;
  const Real a_co = atc_co * at_co;
  const Real a_hs = atc_hs * at_hs;

  // Initial gas composition
  Real y1 = y_ch;
  Real a1 = a_ch;
  Real y2 = y_hs;
  Real a2 = a_hs;
  Real y3 = y_co;
  Real a3 = a_co;

  //
  const Real k13 = k_co_ch;
  const Real k23 = k_co_hs;
  const Real k12 = k_hs_ch;

  //  non- random mixing rule(equation-A8)

  const Real a_mix1 = y1 * y1 * std::pow((a1 * a1), 0.5) +
                      (1 - k12) * 2 * y1 * y2 * std::pow((a1 * a2), 0.5) +
                      (1 - k13) * 2 * y1 * y3 * std::pow((a1 * a3), 0.5);
  const Real a_mix2 = y2 * y2 * std::pow((a2 * a2), 0.5) +
                      (1 - k23) * 2 * y2 * y3 * std::pow((a2 * a3), 0.5) +
                      y3 * y3 * std::pow((a3 * a3), 0.5);
  const Real a_mix_c = a_mix1 + a_mix2;

  //  non- random mixing rule(equation-A9)
  const Real a_mix_p =
      y1 * y2 * ap_ch_hs + y1 * y3 * ap_ch_co + y2 * y3 * ap_hs_co + y3 * y2 * ap_co_hs;

  //  non- random mixing rule(equation-A7)
  const Real a_mix = a_mix_c + a_mix_p;

  //  non- random mixing rule(equation-A10)
  const Real b_mix = y_co * b_co + y_ch * b_ch + y_hs * b_hs;

  //  non- random mixing rule(equation-A13-A14)
  const Real A = a_mix * p / (r * t * r * t);
  const Real B = b_mix * p / (r * t);

  // calculation of fugacity of liquid phase
  // calculation of z-factor
  // calculation of z-factor of gas phase

  const Real bet = b_mix * p / (r * t);
  const Real q = a_mix / (b_mix * r * t * std::pow(t, 0.5));

  const Real state = 1;
  Real er, zn, z, ff;
  if (state == 1)
  {

    {
      er = 1;
      z = 1;
      int kk = 1;
      while (er > 0.0001)
      {
        zn = 1 + bet -
             q * bet *
                 ((z - bet) / ((z + bet * (1 - std::pow(2, 0.5))) * (z + bet)) *
                  (z + bet * (1 + std::pow(2, 0.5))));
        er = abs(z - zn);
        z = zn;
        kk = kk + 1;
        if (kk > 500)
        {
          er = 0;
        }
      }
    }
  }
  else
  {
    // calculation of z-factor of liquid phase

    {
      er = 1;
      z = 0;
      // int kk = 1;
      while (er > 0.0001)
      {
        z = z + 0.000005;
        a1 = -(1 - B);
        a2 = (A - 2 * B - 3 * B * B);
        a3 = -(A * B - B * B - std::pow(B, 3));
        ff = std::pow(z, 3) + a1 * z * z + a2 * z + a3;
        er = abs(ff);
      }
    }
  }
  {
    zn = z;
  }

  // Gas phase volume

  const Real vm = zn * r * t / p;
  // fugacity coefficient calculation(equation A-11)

  Real fi, fi_ch, fi_hs, fi_co, ak, bk;
  // a_ik used  in A-11

  const Real a_co_ch = std::pow((a_co * a_ch), 0.5) * (1 - k_co_ch);
  const Real a_co_hs = std::pow((a_co * a_hs), 0.5) * (1 - k_co_hs);
  const Real a_hs_ch = std::pow((a_hs * a_ch), 0.5) * (1 - k_hs_ch);

  // fi of CH4
  Real am = a_mix;
  Real bm = b_mix;
  bk = b_ch;
  ak = y_co * (a_co_ch) + y_ch * a_ch + y_hs * a_hs_ch;

  a1 = ((bk / bm) * (zn - 1)) - log(p * (vm - bm) / (r * t));
  a2 = log((vm + (1 + std::pow(2, 0.5)) * bm) / (vm + (1 - std::pow(2, 0.5)) * bm));
  a3 = am / (2 * std::pow(2, 0.5) * bm * r * t);
  Real a4 = (2 * ak / am) - (bk / bm);

  fi = a1 - a2 * a3 * a4;
  fi_ch = exp(fi);

  // fi of CO2

  bk = b_co;
  ak = y_ch * (a_co_ch + ap_ch_co) + y_co * a_co + y_hs * (a_co_hs + ap_hs_co);

  a1 = ((bk / bm) * (zn - 1)) - log(p * (vm - bm) / (r * t));
  a2 = log((vm + (1 + std::pow(2, 0.5)) * bm) / (vm + (1 - std::pow(2, 0.5)) * bm));
  a3 = am / (2 * std::pow(2, 0.5) * bm * r * t);
  a4 = (2 * ak / am) - (bk / bm);

  fi = a1 - a2 * a3 * a4;
  fi_co = exp(fi);

  // fi of h2s

  bk = b_hs;
  ak = y_co * (a_co_hs + ap_co_hs) + y_hs * a_hs + y_ch * (a_hs_ch + ap_ch_hs);

  a1 = ((bk / bm) * (zn - 1)) - log(p * (vm - bm) / (r * t));
  a2 = log((vm + (1 + std::pow(2, 0.5)) * bm) / (vm + (1 - std::pow(2, 0.5)) * bm));
  a3 = am / (2 * std::pow(2, 0.5) * bm * r * t);
  a4 = (2 * ak / am) - (bk / bm);

  fi = a1 - a2 * a3 * a4;
  fi_hs = exp(fi);

  // Partial molar volume CH4 (equation 11)
  Real k, tt, v_ch, k_ch;

  v_ch = 3.541 + 1.23 * 0.001 * (t - 273.15);
  v_ch = exp(v_ch);

  // Henry's constant of CH4 (equation 12)
  tt = t * 0.01;
  k = 127.173804 - 155.575631 / tt - 65.2552591 * log(tt) + 6.16975729 * tt;
  k_ch = exp(k);
  k_ch = k_ch / 100000;

  // Henry's constant of gaseous and liquid CO2 (equation 7&8)
  Real k_co;
  Real v_co;
  if (state == 1)
  {
    v_co = 32.6; // Partial molar volume CO2
    k = 1.189 + 1.304 * std::pow(10, -2) * (t - 273) -
        (5.446 * std::pow(10, -5)) * std::pow((t - 273), 2);
    k_co = std::pow(10, k);
  }
  else
  {
    v_co = 32;
    k = 1.169 + 1.368 * std::pow(10, -2) * (t - 273) -
        (5.38 * std::pow(10, -5)) * std::pow((t - 273), 2);
    k_co = std::pow(10, k);
  }

  {
    k_co = 55.508 * k_co;
  }
  // Henry's constant of gaseous and liquid H2S (equation 9&10)

  const Real v_hs = 35.5; // Partial molar volume H2S

  Real k_hs;
  if (state == 1)
  {
    const Real bk1 = 640.1901754766855;
    const Real bk2 = 0.27643029942021946;
    const Real bk3 = 0.11099712781812077;
    const Real bk4 = 0.1646266738977218;
    const Real bk5 = 0.26532569367900316;
    k = bk1 + bk2 * (t) - (bk3 * std::pow(10, -3)) * std::pow(t, 2) -
        (bk4 * std::pow(10, 5)) / (t) - (bk5 * std::pow(10, 3)) * log10(t);
    k_hs = std::pow(10, k);
    k_hs = 55.508 * k_hs;
  }
  else
  {
    const Real br1 = 358.0300140642526;
    const Real br2 = -12604.137769112192;
    const Real br3 = 0.07327349335024513;
    const Real br4 = -59.26961904310672;

    k_hs = br1 + br2 / (t) + (br3) * (t) + br4 * log(t);
    k_hs = 10 * exp(k_hs);
    k_hs = 55.508 * k_hs;
  }

  // mole fraction of dissolved gases in the aqueous phase (equation 6)
  // activity coefficient is not known

  const Real bb_ch = y_ch * (fi_ch * p * exp(-v_ch * (p - 1) / (r * t))) / k_ch;
  const Real bb_co = y_co * (fi_co * p * exp(-v_co * (p - 1) / (r * t))) / k_co;
  const Real bb_hs = y_hs * (fi_hs * p * exp(-v_hs * (p - 1) / (r * t))) / k_hs;

  // mole fraction of water in gas phase (equation 13)
  am = a_mix;
  bm = b_mix;
  bk = 18.18;

  // Binary and intermolecular interaction Parameters

  Real ak_ch, ak_co, ak_hs, ak1, ak2, ak3;

  if (state == 1)
  {
    // Binary and intermolecular interaction H2O-Gas
    ak_ch = 2160707.911971258;
    ak_co = 1.6658915512132246E7 - 35955.73073289321 * t;
    ak_hs = 5903682.145090856;

    // Binary and intermolecular interaction H2O-Gas-Gas
    ak1 = -1863733.8738954344;
    ak2 = 1.9524270682006843E7 - 64566.01299897294 * t + p * 145139.97070876797;
    ak3 = -2702471.307888152;
  }
  else
  {
    // Binary and intermolecular interaction H2O-Gas
    ak_ch = 1.1012594358645048E7;
    ak_hs = 7478118.8440403305 - 4143.303720151969 * t;
    ak_co = 7722804.441254632 - 8958.797691006377 * t;

    // Binary and intermolecular interaction H2O-Gas-Gas
    ak1 = -1.2357484794511057E7;
    ak2 = -3.2996010087810606E8;
    ak3 = -1.456473568862047E7;
  }

  // intermolecular interaction parameter between water and Acid gas(equation 15)
  Real fi_h2o, yc, yn_ch, yn_co, yn_hs;
  {
    ak = y_ch * ak_ch + y_co * ak_co + y_hs * ak_hs + y_co * y_ch * ak1 + y_co * y_hs * (ak2) +
         y_hs * y_ch * (ak3);
  }
  // fugacity coefficient calculation(equation A-11)
  a1 = ((bk / bm) * (zn - 1)) - log(p * (vm - bm) / (r * t));
  a2 = log((vm + (1 + std::pow(2, 0.5)) * bm) / (vm + (1 - std::pow(2, 0.5)) * bm));
  a3 = am / (2 * std::pow(2, 0.5) * bm * r * t);
  a4 = (2 * ak / am) - (bk / bm);
  fi = a1 - a2 * a3 * a4;
  fi_h2o = exp(fi);

  // mole fraction of water in gas phase (equation 13)
  // equilibrium constant of water at reference pressure( equation 14)

  const Real vh = 18.18;
  k = -2.209 + 3.097 * std::pow(10, -2) * (t - 273) -
      (1.098 * std::pow(10, -4)) * std::pow((t - 273), 2) +
      (2.048 * std::pow(10, -7)) * std::pow((t - 273), 3);
  k = std::pow(10, k);
  Real aa_h2o = (k / (fi_h2o * p)) * exp(vh * (p - 1) / (r * t));

  // yc is water mole fraction in gas phase (equation 13)

  yc = (1 - bb_co - bb_ch - bb_hs) / ((1 / aa_h2o) - bb_co - bb_ch - bb_hs);

  // normolization process for gas composition

  yn_ch = y_ch / (y_ch + y_co + y_hs + yc);
  yn_co = y_co / (y_ch + y_co + y_hs + yc);
  yn_hs = y_hs / (y_ch + y_co + y_hs + yc);
  y_ch = yn_ch;
  y_co = yn_co;
  y_hs = yn_hs;

  Real x_ch, x_co, x_hs;
  // mole fraction of dissolved gases in the aqueous phase (equation 6)
  x_ch = y_ch * (fi_ch * p * exp(-v_ch * (p - 1) / (r * t))) / k_ch;
  x_co = y_co * (fi_co * p * exp(-v_co * (p - 1) / (r * t))) / k_co;
  x_hs = y_hs * (fi_hs * p * exp(-v_hs * (p - 1) / (r * t))) / k_hs;
  xch4 = x_ch;
  xh2s = x_hs;
  xco2 = x_co;
  yh2o = yc;
}

void
PorousFlowBrineCO2CH4H2S::thermophysicalProperties(Real pressure,
                                                   Real temperature,
                                                   Real Xnacl,
                                                   std::vector<const VariableValue *> & Z,
                                                   unsigned int qp,
                                                   std::vector<FluidStateProperties> & fsp) const
{
  // Extract the Z's for each component TODO: fix numbering
  const Real ZCO2 = (*Z[_zco2_idx])[qp];
  const Real ZCH4 = (*Z[_zch4_idx])[qp];
  const Real ZH2S = (*Z[_zh2s_idx])[qp];

  FluidStateProperties & liquid = fsp[_aqueous_phase_number];
  FluidStateProperties & gas = fsp[_gas_phase_number];

  // Check whether the input temperature is within the region of validity
  checkVariables(pressure, temperature);

  // Clear all of the FluidStateProperties data
  clearFluidStateProperties(fsp);

  FluidStatePhaseEnum phase_state;
  massFractions(pressure, temperature, Xnacl, ZCO2, ZCH4, ZH2S, phase_state, fsp);

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
      liquidProperties(liquid_pressure, temperature, Xnacl, fsp);

      break;
    }

    case FluidStatePhaseEnum::TWOPHASE:
    {
      // Calculate the gas properties
      gasProperties(pressure, temperature, fsp);

      // Calculate the saturation
      twoPhaseProperties(pressure, temperature, Xnacl, ZCO2, ZCH4, ZH2S, fsp);

      // Calculate the liquid properties
      Real liquid_pressure = pressure - _pc_uo.capillaryPressure(1.0 - gas.saturation);
      liquidProperties(liquid_pressure, temperature, Xnacl, fsp);

      break;
    }
  }

  // Liquid saturations can now be set
  liquid.saturation = 1.0 - gas.saturation;
  liquid.dsaturation_dp = -gas.dsaturation_dp;
  liquid.dsaturation_dT = -gas.dsaturation_dT;
  liquid.dsaturation_dX = -gas.dsaturation_dX;
  for (unsigned int z = 0; z < _num_zvars; ++z)
    liquid.dsaturation_dZ[z] = -gas.dsaturation_dZ[z];

  // Save pressures to FluidStateProperties object
  gas.pressure = pressure;
  liquid.pressure = pressure - _pc_uo.capillaryPressure(liquid.saturation);
}

void
PorousFlowBrineCO2CH4H2S::massFractions(Real pressure,
                                        Real temperature,
                                        Real Xnacl,
                                        Real ZCO2,
                                        Real ZCH4,
                                        Real ZH2S,
                                        Real & Xco2,
                                        Real & Yco2,
                                        Real & Xch4,
                                        Real & Ych4,
                                        Real & Xh2s,
                                        Real & Yh2s) const
{
  // Determine the vapor mass fraction from Rachford Rice equation
  // Note: assume that gas composition is just Zi initially with no H2O
  const Real Zco2m = ZCO2 / (ZCO2 + ZH2S + ZCH4);
  const Real Zh2sm = ZH2S / (ZCO2 + ZH2S + ZCH4);
  const Real Zch4m = ZCH4 / (ZCO2 + ZH2S + ZCH4);

  Real Yh2o;
  equilibriumMassFractions(
      pressure, temperature, Xnacl, Zco2m, Zch4m, Zh2sm, Xco2, Xch4, Xh2s, Yh2o);

  // Equilibrium constants (Yi/Xi)
  std::vector<Real> Zi{ZCO2, ZCH4, ZH2S};

  Real Ki_CO2, Ki_CH4, Ki_H2S;
  if (ZCO2 != 0)
    Ki_CO2 = Zco2m / Xco2;
  else
    Ki_CO2 = 0;

  if (ZCH4 != 0)
    Ki_CH4 = Zch4m / Xch4;
  else
    Ki_CH4 = 0;

  if (ZH2S != 0)
    Ki_H2S = Zh2sm / Xh2s;
  else
    Ki_H2S = 0;

  std::vector<Real> Ki{Ki_CO2, Ki_CH4, Ki_H2S, 0.0};

  // Iteratively solve for V as Ki is updated
  Real V0 = vaporMassFraction(Zi, Ki);

  Real V = 1.0;
  Real tol = 1.0e-8;
  unsigned int iter = 0;

  while (std::abs(V - V0) > tol)
  {
    V = V0;

    if (Zi[0] != 0)
      Xco2 = Zi[0] / ((Ki[0] - 1.0) * V + 1);
    else
      Xco2 = 0;

    if (Zi[1] != 0)
      Xch4 = Zi[1] / ((Ki[1] - 1.0) * V + 1);
    else
      Xch4 = 0;

    if (Zi[2] != 0)
      Xh2s = Zi[2] / ((Ki[2] - 1.0) * V + 1);
    else
      Xh2s = 0;

    Yco2 = Xco2 * Ki[0];
    Ych4 = Xch4 * Ki[1];
    Yh2s = Xh2s * Ki[2];

    equilibriumMassFractions(
        pressure, temperature, Xnacl, Yco2, Ych4, Yh2s, Xco2, Xch4, Xh2s, Yh2o);

    if (Yco2 != 0)
      Ki[0] = Yco2 / Xco2;
    else
      Ki[0] = 0;

    if (Ych4 != 0)
      Ki[1] = Ych4 / Xch4;
    else
      Ki[1] = 0;

    if (Yh2s != 0)
      Ki[2] = Yh2s / Xh2s;
    else
      Ki[2] = 0;

    Ki[3] = Yh2o / (1.0 - Xco2 - Xch4 - Xh2s);

    V0 = vaporMassFraction(Zi, Ki);

    iter++;

    if (iter > _nr_max_its)
      break;
  }

  // Update the equilibrium mass fractions
  if (Zi[0] != 0)
    Xco2 = Zi[0] / ((Ki[0] - 1.0) * V + 1);
  else
    Xco2 = 0;

  if (Zi[1] != 0)
    Xch4 = Zi[1] / ((Ki[1] - 1.0) * V + 1);
  else
    Xch4 = 0;

  if (Zi[2] != 0)
    Xh2s = Zi[2] / ((Ki[2] - 1.0) * V + 1);
  else
    Xh2s = 0;
}

void
PorousFlowBrineCO2CH4H2S::massFractions(Real pressure,
                                        Real temperature,
                                        Real Xnacl,
                                        Real ZCO2,
                                        Real ZCH4,
                                        Real ZH2S,
                                        FluidStatePhaseEnum & phase_state,
                                        std::vector<FluidStateProperties> & fsp) const
{
  FluidStateProperties & liquid = fsp[_aqueous_phase_number];
  FluidStateProperties & gas = fsp[_gas_phase_number];

  Real Xco2, Yco2, Xch4, Ych4, Xh2s, Yh2s;
  massFractions(pressure, temperature, Xnacl, ZCO2, ZCH4, ZH2S, Xco2, Yco2, Xch4, Ych4, Xh2s, Yh2s);

  // Determine which phases are present based on the value of z
  FluidStatePhaseEnum phase_state_co2, phase_state_ch4, phase_state_h2s;
  phaseState(ZCO2, Xco2, Yco2, phase_state_co2);
  phaseState(ZCH4, Xch4, Ych4, phase_state_ch4);
  phaseState(ZH2S, Xh2s, Yh2s, phase_state_h2s);

  // If all of these are liquid, then only a liquid phase is present.
  // If one of these is gas, then all the liquid component has been evaporated
  // and all phase states should be set to gas.
  // Otherwise, both fluid phases are present
  if (phase_state_co2 == phase_state_ch4 && phase_state_co2 == phase_state_h2s)
    phase_state = phase_state_co2;

  else if (phase_state_co2 == FluidStatePhaseEnum::GAS ||
           phase_state_ch4 == FluidStatePhaseEnum::GAS ||
           phase_state_h2s == FluidStatePhaseEnum::GAS)
  {
    phase_state_co2 = FluidStatePhaseEnum::GAS;
    phase_state_ch4 = FluidStatePhaseEnum::GAS;
    phase_state_h2s = FluidStatePhaseEnum::GAS;
    phase_state = FluidStatePhaseEnum::GAS;
  }

  else
    phase_state = FluidStatePhaseEnum::TWOPHASE;

  // The equilibrium mass fractions calculated above are only correct in the two phase
  // state. If only liquid or gas phases are present, the mass fractions are given by
  // the total mass fraction z
  Real Xh2o = 0.0, Yh2o = 0.0;
  Real dXco2_dp = 0.0, dYco2_dp = 0.0, dXco2_dT = 0.0, dYco2_dT = 0.0;
  Real dXco2_dX = 0.0, dYco2_dX = 0.0, dXco2_dZco2 = 0.0, dYco2_dZco2 = 0.0;
  Real dXco2_dZch4 = 0.0, dYco2_dZch4 = 0.0, dXco2_dZh2s = 0.0, dYco2_dZh2s = 0.0;
  Real dXch4_dp = 0.0, dYch4_dp = 0.0, dXch4_dT = 0.0, dYch4_dT = 0.0;
  Real dXch4_dX = 0.0, dYch4_dX = 0.0, dXch4_dZco2 = 0.0, dYch4_dZco2 = 0.0;
  Real dXch4_dZch4 = 0.0, dYch4_dZch4 = 0.0, dXch4_dZh2s = 0.0, dYch4_dZh2s = 0.0;
  Real dXh2s_dp = 0.0, dYh2s_dp = 0.0, dXh2s_dT = 0.0, dYh2s_dT = 0.0;
  Real dXh2s_dX = 0.0, dYh2s_dX = 0.0, dXh2s_dZco2 = 0.0, dYh2s_dZco2 = 0.0;
  Real dXh2s_dZch4 = 0.0, dYh2s_dZch4 = 0.0, dXh2s_dZh2s = 0.0, dYh2s_dZh2s = 0.0;

  // Check each component individually
  switch (phase_state_co2)
  {
    case FluidStatePhaseEnum::LIQUID:
    {
      Xco2 = ZCO2;
      Yco2 = 0.0;
      dXco2_dZco2 = 1.0;
      break;
    }

    case FluidStatePhaseEnum::GAS:
    {
      Xco2 = 0.0;
      Yco2 = ZCO2;
      dYco2_dZco2 = 1.0;
      break;
    }

    default:
      break;
  }

  switch (phase_state_ch4)
  {
    case FluidStatePhaseEnum::LIQUID:
    {
      Xch4 = ZCH4;
      Ych4 = 0.0;
      dXch4_dZch4 = 1.0;
      break;
    }

    case FluidStatePhaseEnum::GAS:
    {
      Xch4 = 0.0;
      Ych4 = ZCH4;
      dYch4_dZch4 = 1.0;
      break;
    }

    default:
      break;
  }

  switch (phase_state_h2s)
  {
    case FluidStatePhaseEnum::LIQUID:
    {
      Xh2s = ZH2S;
      Yh2s = 0.0;
      dXh2s_dZh2s = 1.0;
      break;
    }

    case FluidStatePhaseEnum::GAS:
    {
      Xh2s = 0.0;
      Yh2s = ZH2S;
      dYh2s_dZh2s = 1.0;
      break;
    }

    default:
      break;
  }

  // Set the H2O mass fractions depending on what the overall phase state is
  Xh2o = 1.0 - Xco2 - Xch4 - Xh2s;
  Yh2o = 1.0 - Yco2 - Ych4 - Yh2s;

  if (phase_state == FluidStatePhaseEnum::LIQUID)
    Yh2o = 0.0;

  else if (phase_state == FluidStatePhaseEnum::GAS)
    Xh2o = 0.0;
  else // phase_state == FluidStatePhaseEnum::TWOPHASE
  {
    // Finite difference derivatives
    // Derivatives wrt pressure
    Real dXco2, dXch4, dXh2s, dYco2, dYch4, dYh2s;
    const Real dp = 1.0e-2;
    massFractions(pressure + dp,
                  temperature,
                  Xnacl,
                  ZCO2,
                  ZCH4,
                  ZH2S,
                  dXco2,
                  dYco2,
                  dXch4,
                  dYch4,
                  dXh2s,
                  dYh2s);

    dXco2_dp = (dXco2 - Xco2) / dp;
    dXch4_dp = (dXch4 - Xch4) / dp;
    dXh2s_dp = (dXh2s - Xh2s) / dp;

    dYco2_dp = (dYco2 - Yco2) / dp;
    dYch4_dp = (dYch4 - Ych4) / dp;
    dYh2s_dp = (dYh2s - Yh2s) / dp;

    // Derivatives wrt temperature
    const Real dT = 1.0e-4;
    massFractions(pressure,
                  temperature + dT,
                  Xnacl,
                  ZCO2,
                  ZCH4,
                  ZH2S,
                  dXco2,
                  dYco2,
                  dXch4,
                  dYch4,
                  dXh2s,
                  dYh2s);

    dXco2_dT = (dXco2 - Xco2) / dT;
    dXch4_dT = (dXch4 - Xch4) / dT;
    dXh2s_dT = (dXh2s - Xh2s) / dT;

    dYco2_dT = (dYco2 - Yco2) / dT;
    dYch4_dT = (dYch4 - Ych4) / dT;
    dYh2s_dT = (dYh2s - Yh2s) / dT;

    // Derivatives wrt salt mass fraction
    const Real dX = 1.0e-8;
    massFractions(pressure,
                  temperature,
                  Xnacl + dX,
                  ZCO2,
                  ZCH4,
                  ZH2S,
                  dXco2,
                  dYco2,
                  dXch4,
                  dYch4,
                  dXh2s,
                  dYh2s);

    dXco2_dX = (dXco2 - Xco2) / dX;
    dXch4_dX = (dXch4 - Xch4) / dX;
    dXh2s_dX = (dXh2s - Xh2s) / dX;

    dYco2_dX = (dYco2 - Yco2) / dX;
    dYch4_dX = (dYch4 - Ych4) / dX;
    dYh2s_dX = (dYh2s - Yh2s) / dX;

    // Derivatives wrt Z's in the two phase region
    const Real dZ = 1.0e-8;
    massFractions(pressure,
                  temperature,
                  Xnacl,
                  ZCO2 + dZ,
                  ZCH4,
                  ZH2S,
                  dXco2,
                  dYco2,
                  dXch4,
                  dYch4,
                  dXh2s,
                  dYh2s);

    dXco2_dZco2 = (dXco2 - Xco2) / dZ;
    dXch4_dZco2 = (dXch4 - Xch4) / dZ;
    dXh2s_dZco2 = (dXh2s - Xh2s) / dZ;

    dYco2_dZco2 = (dYco2 - Yco2) / dZ;
    dYch4_dZco2 = (dYch4 - Ych4) / dZ;
    dYh2s_dZco2 = (dYh2s - Yh2s) / dZ;

    massFractions(pressure,
                  temperature,
                  Xnacl,
                  ZCO2,
                  ZCH4 + dZ,
                  ZH2S,
                  dXco2,
                  dYco2,
                  dXch4,
                  dYch4,
                  dXh2s,
                  dYh2s);

    dXco2_dZch4 = (dXco2 - Xco2) / dZ;
    dXch4_dZch4 = (dXch4 - Xch4) / dZ;
    dXh2s_dZch4 = (dXh2s - Xh2s) / dZ;

    dYco2_dZch4 = (dYco2 - Yco2) / dZ;
    dYch4_dZch4 = (dYch4 - Ych4) / dZ;
    dYh2s_dZch4 = (dYh2s - Yh2s) / dZ;

    massFractions(pressure,
                  temperature,
                  Xnacl,
                  ZCO2,
                  ZCH4,
                  ZH2S + dZ,
                  dXco2,
                  dYco2,
                  dXch4,
                  dYch4,
                  dXh2s,
                  dYh2s);

    dXco2_dZh2s = (dXco2 - Xco2) / dZ;
    dXch4_dZh2s = (dXch4 - Xch4) / dZ;
    dXh2s_dZh2s = (dXh2s - Xh2s) / dZ;

    dYco2_dZh2s = (dYco2 - Yco2) / dZ;
    dYch4_dZh2s = (dYch4 - Ych4) / dZ;
    dYh2s_dZh2s = (dYh2s - Yh2s) / dZ;
  }

  // Save the mass fractions in the FluidStateProperties object
  liquid.mass_fraction[_aqueous_fluid_component] = Xh2o;
  liquid.mass_fraction[_co2_component] = Xco2;
  liquid.mass_fraction[_ch4_component] = Xch4;
  liquid.mass_fraction[_h2s_component] = Xh2s;
  liquid.mass_fraction[_salt_component] = Xnacl;

  gas.mass_fraction[_aqueous_fluid_component] = Yh2o;
  gas.mass_fraction[_co2_component] = Yco2;
  gas.mass_fraction[_ch4_component] = Ych4;
  gas.mass_fraction[_h2s_component] = Yh2s;

  // Save the mole fractions in the FluidStateProperties object
  const Real ldenom = Xco2 / _Mco2 + Xh2s / _Mh2s + Xch4 / _Mch4 + Xh2o / _Mh2o;
  liquid.mole_fraction[_aqueous_fluid_component] = Xh2o / _Mh2o / ldenom;
  liquid.mole_fraction[_co2_component] = Xco2 / _Mco2 / ldenom;
  liquid.mole_fraction[_ch4_component] = Xch4 / _Mch4 / ldenom;
  liquid.mole_fraction[_h2s_component] = Xh2s / _Mh2s / ldenom;

  const Real gdenom = Yco2 / _Mco2 + Yh2s / _Mh2s + Ych4 / _Mch4 + Yh2o / _Mh2o;
  gas.mole_fraction[_aqueous_fluid_component] = Yh2o / _Mh2o / gdenom;
  gas.mole_fraction[_co2_component] = Yco2 / _Mco2 / gdenom;
  gas.mole_fraction[_ch4_component] = Ych4 / _Mch4 / gdenom;
  gas.mole_fraction[_h2s_component] = Yh2s / _Mh2s / gdenom;

  // Save the derivatives wrt PorousFlow variables
  // Liquid component
  liquid.dmass_fraction_dp[_aqueous_fluid_component] = -dXco2_dp - dXch4_dp - dXh2s_dp;
  liquid.dmass_fraction_dT[_aqueous_fluid_component] = -dXco2_dT - dXch4_dT - dXh2s_dT;
  liquid.dmass_fraction_dX[_aqueous_fluid_component] = -dXco2_dX - dXch4_dX - dXh2s_dX;
  liquid.dmass_fraction_dZ[_aqueous_fluid_component][_zco2_idx] =
      -dXco2_dZco2 - dXch4_dZco2 - dXh2s_dZco2;
  liquid.dmass_fraction_dZ[_aqueous_fluid_component][_zch4_idx] =
      -dXco2_dZch4 - dXch4_dZch4 - dXh2s_dZch4;
  liquid.dmass_fraction_dZ[_aqueous_fluid_component][_zh2s_idx] =
      -dXco2_dZh2s - dXch4_dZh2s - dXh2s_dZh2s;
  gas.dmass_fraction_dp[_aqueous_fluid_component] = -dYco2_dp - dYch4_dp - dYh2s_dp;
  gas.dmass_fraction_dT[_aqueous_fluid_component] = -dYco2_dT - dYch4_dT - dYh2s_dT;
  gas.dmass_fraction_dX[_aqueous_fluid_component] = -dYco2_dX - dYch4_dX - dYh2s_dX;
  gas.dmass_fraction_dZ[_aqueous_fluid_component][_zco2_idx] =
      -dYco2_dZco2 - dYch4_dZco2 - dYh2s_dZco2;
  gas.dmass_fraction_dZ[_aqueous_fluid_component][_zch4_idx] =
      -dYco2_dZch4 - dYch4_dZch4 - dYh2s_dZch4;
  gas.dmass_fraction_dZ[_aqueous_fluid_component][_zh2s_idx] =
      -dYco2_dZh2s - dYch4_dZh2s - dYh2s_dZh2s;
  // CO2 component
  liquid.dmass_fraction_dp[_co2_component] = dXco2_dp;
  liquid.dmass_fraction_dT[_co2_component] = dXco2_dT;
  liquid.dmass_fraction_dX[_co2_component] = dXco2_dX;
  liquid.dmass_fraction_dZ[_co2_component][_zco2_idx] = dXco2_dZco2;
  liquid.dmass_fraction_dZ[_co2_component][_zch4_idx] = dXco2_dZch4;
  liquid.dmass_fraction_dZ[_co2_component][_zh2s_idx] = dXco2_dZh2s;
  gas.dmass_fraction_dp[_co2_component] = dYco2_dp;
  gas.dmass_fraction_dT[_co2_component] = dYco2_dT;
  gas.dmass_fraction_dX[_co2_component] = dYco2_dX;
  gas.dmass_fraction_dZ[_co2_component][_zco2_idx] = dYco2_dZco2;
  gas.dmass_fraction_dZ[_co2_component][_zch4_idx] = dYco2_dZch4;
  gas.dmass_fraction_dZ[_co2_component][_zh2s_idx] = dYco2_dZh2s;

  // CH4 component
  liquid.dmass_fraction_dp[_ch4_component] = dXch4_dp;
  liquid.dmass_fraction_dT[_ch4_component] = dXch4_dT;
  liquid.dmass_fraction_dX[_ch4_component] = dXch4_dX;
  liquid.dmass_fraction_dZ[_ch4_component][_zco2_idx] = dXch4_dZco2;
  liquid.dmass_fraction_dZ[_ch4_component][_zch4_idx] = dXch4_dZch4;
  liquid.dmass_fraction_dZ[_ch4_component][_zh2s_idx] = dXch4_dZh2s;
  gas.dmass_fraction_dp[_ch4_component] = dYch4_dp;
  gas.dmass_fraction_dT[_ch4_component] = dYch4_dT;
  gas.dmass_fraction_dX[_ch4_component] = dYch4_dX;
  gas.dmass_fraction_dZ[_ch4_component][_zco2_idx] = dYch4_dZco2;
  gas.dmass_fraction_dZ[_ch4_component][_zch4_idx] = dYch4_dZch4;
  gas.dmass_fraction_dZ[_ch4_component][_zh2s_idx] = dYch4_dZh2s;

  // H2S component
  liquid.dmass_fraction_dp[_h2s_component] = dXh2s_dp;
  liquid.dmass_fraction_dT[_h2s_component] = dXh2s_dT;
  liquid.dmass_fraction_dX[_h2s_component] = dXh2s_dX;
  liquid.dmass_fraction_dZ[_h2s_component][_zco2_idx] = dXh2s_dZco2;
  liquid.dmass_fraction_dZ[_h2s_component][_zch4_idx] = dXh2s_dZch4;
  liquid.dmass_fraction_dZ[_h2s_component][_zh2s_idx] = dXh2s_dZh2s;
  gas.dmass_fraction_dp[_h2s_component] = dYh2s_dp;
  gas.dmass_fraction_dT[_h2s_component] = dYh2s_dT;
  gas.dmass_fraction_dX[_h2s_component] = dYh2s_dX;
  gas.dmass_fraction_dZ[_h2s_component][_zco2_idx] = dYh2s_dZco2;
  gas.dmass_fraction_dZ[_h2s_component][_zch4_idx] = dYh2s_dZch4;
  gas.dmass_fraction_dZ[_h2s_component][_zh2s_idx] = dYh2s_dZh2s;

  // Salt component
  liquid.dmass_fraction_dX[_salt_component] = 1.0;
}

Real
PorousFlowBrineCO2CH4H2S::gasDensity(
    Real pressure, Real temperature, Real Yco2, Real Ych4, Real Yh2s) const
{
  const Real gdenom =
      Yco2 / _Mco2 + Yh2s / _Mh2s + Ych4 / _Mch4 + (1.0 - Yco2 - Ych4 - Yh2s) / _Mh2o;
  const Real yco2 = Yco2 / _Mco2 / gdenom;
  const Real ych4 = Ych4 / _Mco2 / gdenom;
  const Real yh2s = Yh2s / _Mco2 / gdenom;

  Real Zg;
  Real Mw_mix = yco2 * _Mco2 + ych4 * _Mch4 + yh2s * _Mh2s; //(gr/mol)

  GasCompressibilityFactor(pressure, temperature, yco2, ych4, yh2s, Zg);

  return pressure * Mw_mix / (Zg * _R * temperature); // (gr/m3)
}

void
PorousFlowBrineCO2CH4H2S::gasProperties(Real pressure,
                                        Real temperature,
                                        std::vector<FluidStateProperties> & fsp) const
{
  FluidStateProperties & gas = fsp[_gas_phase_number];

  // Mole fraction of gas components in gas phase
  // const Real Yh2o = gas.mass_fraction[_aqueous_fluid_component];
  const Real Yco2 = gas.mass_fraction[_co2_component];
  const Real Ych4 = gas.mass_fraction[_ch4_component];
  const Real Yh2s = gas.mass_fraction[_h2s_component];

  const Real dYco2_dZco2 = gas.dmass_fraction_dZ[_co2_component][_zco2_idx];
  const Real dYco2_dZch4 = gas.dmass_fraction_dZ[_co2_component][_zch4_idx];
  const Real dYco2_dZh2s = gas.dmass_fraction_dZ[_co2_component][_zh2s_idx];

  const Real dYch4_dZco2 = gas.dmass_fraction_dZ[_ch4_component][_zco2_idx];
  const Real dYch4_dZch4 = gas.dmass_fraction_dZ[_ch4_component][_zch4_idx];
  const Real dYch4_dZh2s = gas.dmass_fraction_dZ[_ch4_component][_zh2s_idx];

  const Real dYh2s_dZco2 = gas.dmass_fraction_dZ[_h2s_component][_zco2_idx];
  const Real dYh2s_dZch4 = gas.dmass_fraction_dZ[_h2s_component][_zch4_idx];
  const Real dYh2s_dZh2s = gas.dmass_fraction_dZ[_h2s_component][_zh2s_idx];

  // Gas density and derivatives
  const Real rho = gasDensity(pressure, temperature, Yco2, Ych4, Yh2s);

  const Real dp = 1.0e-2;
  Real drho = gasDensity(pressure + dp, temperature, Yco2, Ych4, Yh2s);
  const Real drho_dp = (drho - rho) / dp;

  const Real dT = 1.0e-4;
  drho = gasDensity(pressure, temperature + dT, Yco2, Ych4, Yh2s);
  const Real drho_dT = (drho - rho) / dT;

  // Derivative wrt Z by the chain rule drho/dZ = drho/dY * dY/dZ
  const Real dY = 1.0e-8;
  drho = gasDensity(pressure, temperature, Yco2 + dY, Ych4, Yh2s);
  const Real drho_dYco2 = (drho - rho) / dY;
  drho = gasDensity(pressure, temperature, Yco2, Ych4 + dY, Yh2s);
  const Real drho_dYch4 = (drho - rho) / dY;
  drho = gasDensity(pressure, temperature, Yco2, Ych4, Yh2s + dY);
  const Real drho_dYh2s = (drho - rho) / dY;

  const Real drho_dZco2 =
      drho_dYco2 * dYco2_dZco2 + drho_dYch4 * dYch4_dZco2 + drho_dYh2s * dYh2s_dZco2;

  const Real drho_dZch4 =
      drho_dYco2 * dYco2_dZch4 + drho_dYch4 * dYch4_dZch4 + drho_dYh2s * dYh2s_dZch4;

  const Real drho_dZh2s =
      drho_dYco2 * dYco2_dZh2s + drho_dYch4 * dYch4_dZh2s + drho_dYh2s * dYh2s_dZh2s;

  // Save the values to the FluidStateProperties object
  gas.density = rho;
  gas.ddensity_dp = drho_dp;
  gas.ddensity_dT = drho_dT;
  gas.ddensity_dZ[_zco2_idx] = drho_dZco2;
  gas.ddensity_dZ[_zch4_idx] = drho_dZch4;
  gas.ddensity_dZ[_zh2s_idx] = drho_dZh2s;

  // Gas viscosity is approximated with pure CO2
  Real co2_viscosity, dco2_viscosity_dp, dco2_viscosity_dT;
  _co2_fp.mu_from_p_T(pressure, temperature, co2_viscosity, dco2_viscosity_dp, dco2_viscosity_dT);

  gas.viscosity = co2_viscosity;
  gas.dviscosity_dp = dco2_viscosity_dp;
  gas.dviscosity_dT = dco2_viscosity_dT;
}

Real
PorousFlowBrineCO2CH4H2S::liquidDensity(Real pressure,
                                        Real temperature,
                                        Real Xnacl,
                                        std::vector<FluidStateProperties> & fsp) const
{
  FluidStateProperties & liquid = fsp[_aqueous_phase_number];

  // The liquid density includes the density increase due to dissolved CO2
  Real brine_density, dbrine_density_dp, dbrine_density_dT, dbrine_density_dX;
  _brine_fp.rho_dpTx(pressure,
                     temperature,
                     Xnacl,
                     brine_density,
                     dbrine_density_dp,
                     dbrine_density_dT,
                     dbrine_density_dX);

  // Mass fraction of Gas components in liquid phase
  const Real Xco2 = liquid.mass_fraction[_co2_component];
  const Real Xch4 = liquid.mass_fraction[_ch4_component];
  const Real Xh2s = liquid.mass_fraction[_h2s_component];

  // The liquid density
  Real partial_density_co2, partial_density_ch4, partial_density_h2s, dpartial_density_dT_co2,
      dpartial_density_dT_ch4, dpartial_density_dT_h2s;

  partialDensityCO2CH4H2S(temperature,
                          partial_density_co2,
                          partial_density_ch4,
                          partial_density_h2s,
                          dpartial_density_dT_co2,
                          dpartial_density_dT_ch4,
                          dpartial_density_dT_h2s);

  const Real liquid_density =
      1.0 / (Xco2 / partial_density_co2 + Xch4 / partial_density_ch4 + Xh2s / partial_density_h2s +
             (1.0 - Xco2 - Xch4 - Xh2s - Xnacl) / brine_density);

  return liquid_density;
}

void
PorousFlowBrineCO2CH4H2S::liquidProperties(Real pressure,
                                           Real temperature,
                                           Real Xnacl,
                                           std::vector<FluidStateProperties> & fsp) const
{
  FluidStateProperties & liquid = fsp[_aqueous_phase_number];

  // The liquid density includes the density increase due to dissolved gas components
  Real brine_density, dbrine_density_dp, dbrine_density_dT, dbrine_density_dX;
  _brine_fp.rho_dpTx(pressure,
                     temperature,
                     Xnacl,
                     brine_density,
                     dbrine_density_dp,
                     dbrine_density_dT,
                     dbrine_density_dX);

  // Mass fraction of gas components in liquid phase
  const Real Xco2 = liquid.mass_fraction[_co2_component];
  const Real Xch4 = liquid.mass_fraction[_ch4_component];
  const Real Xh2s = liquid.mass_fraction[_h2s_component];

  // The liquid density
  Real partial_density_co2, partial_density_ch4, partial_density_h2s, dpartial_density_co2_dT,
      dpartial_density_ch4_dT, dpartial_density_h2s_dT;

  partialDensityCO2CH4H2S(temperature,
                          partial_density_co2,
                          partial_density_ch4,
                          partial_density_h2s,
                          dpartial_density_co2_dT,
                          dpartial_density_ch4_dT,
                          dpartial_density_h2s_dT);

  const Real liquid_density =
      1.0 / (Xco2 / partial_density_co2 + Xch4 / partial_density_ch4 + Xh2s / partial_density_h2s +
             (1.0 - Xco2 - Xch4 - Xh2s) / brine_density);

  // Derivatives are all of the same form
  auto drho =
      [Xco2,
       Xch4,
       Xh2s,
       partial_density_co2,
       partial_density_ch4,
       partial_density_h2s,
       brine_density](
          Real dXco2, Real dXch4, Real dXh2s, Real dco2pd, Real dch4pd, Real dh2spd, Real drhob) {
        // The general form of the derivative
        const Real denom =
            (Xco2 / partial_density_co2 + Xch4 / partial_density_ch4 + Xh2s / partial_density_h2s +
             (1.0 - Xco2 - Xch4 - Xh2s) / brine_density);
        Real deriv = -dXco2 / partial_density_co2 - dXch4 / partial_density_ch4 -
                     dXh2s / partial_density_h2s + (dXco2 + dXch4 + dXh2s) / brine_density;
        deriv += Xco2 * dco2pd / partial_density_co2 / partial_density_co2;
        deriv += Xch4 * dch4pd / partial_density_ch4 / partial_density_ch4;
        deriv += Xh2s * dh2spd / partial_density_h2s / partial_density_h2s;
        deriv += (1.0 - Xco2 - Xch4 - Xh2s) * drhob / Utility::pow<2>(brine_density);

        return deriv / denom / denom;
      };

  // Derivative wrt pressure
  const Real dXco2_dp = liquid.dmass_fraction_dp[_co2_component];
  const Real dXch4_dp = liquid.dmass_fraction_dp[_ch4_component];
  const Real dXh2s_dp = liquid.dmass_fraction_dp[_h2s_component];

  const Real dliquid_density_dp =
      drho(dXco2_dp, dXch4_dp, dXh2s_dp, 0.0, 0.0, 0.0, dbrine_density_dp);

  // Derivative wrt temperature
  const Real dXco2_dT = liquid.dmass_fraction_dT[_co2_component];
  const Real dXch4_dT = liquid.dmass_fraction_dT[_ch4_component];
  const Real dXh2s_dT = liquid.dmass_fraction_dT[_h2s_component];

  const Real dliquid_density_dT = drho(dXco2_dT,
                                       dXch4_dT,
                                       dXh2s_dT,
                                       dpartial_density_co2_dT,
                                       dpartial_density_ch4_dT,
                                       dpartial_density_h2s_dT,
                                       dbrine_density_dT);

  // Derivative wrt salt mass fraction
  const Real dXco2_dX = liquid.dmass_fraction_dX[_co2_component];
  const Real dXch4_dX = liquid.dmass_fraction_dX[_ch4_component];
  const Real dXh2s_dX = liquid.dmass_fraction_dX[_h2s_component];

  const Real dliquid_density_dX =
      drho(dXco2_dX, dXch4_dX, dXh2s_dX, 0.0, 0.0, 0.0, dbrine_density_dX);

  // Derivative wrt Zco2
  const Real dXco2_dZco2 = liquid.dmass_fraction_dZ[_co2_component][_zco2_idx];
  const Real dXch4_dZco2 = liquid.dmass_fraction_dZ[_ch4_component][_zco2_idx];
  const Real dXh2s_dZco2 = liquid.dmass_fraction_dZ[_h2s_component][_zco2_idx];

  const Real dliquid_density_dZco2 =
      drho(dXco2_dZco2, dXch4_dZco2, dXh2s_dZco2, 0.0, 0.0, 0.0, 0.0);

  // Derivative wrt Zch4
  const Real dXco2_dZch4 = liquid.dmass_fraction_dZ[_co2_component][_zch4_idx];
  const Real dXch4_dZch4 = liquid.dmass_fraction_dZ[_ch4_component][_zch4_idx];
  const Real dXh2s_dZch4 = liquid.dmass_fraction_dZ[_h2s_component][_zch4_idx];

  const Real dliquid_density_dZch4 =
      drho(dXco2_dZch4, dXch4_dZch4, dXh2s_dZch4, 0.0, 0.0, 0.0, 0.0);

  // Derivative wrt Zh2s
  const Real dXco2_dZh2s = liquid.dmass_fraction_dZ[_co2_component][_zh2s_idx];
  const Real dXch4_dZh2s = liquid.dmass_fraction_dZ[_ch4_component][_zh2s_idx];
  const Real dXh2s_dZh2s = liquid.dmass_fraction_dZ[_h2s_component][_zh2s_idx];

  const Real dliquid_density_dZh2s =
      drho(dXco2_dZh2s, dXch4_dZh2s, dXh2s_dZh2s, 0.0, 0.0, 0.0, 0.0);

  // Save the values to the FluidStateProperties object
  liquid.density = liquid_density;
  liquid.ddensity_dp = dliquid_density_dp;
  liquid.ddensity_dT = dliquid_density_dT;
  liquid.ddensity_dX = dliquid_density_dX;
  liquid.ddensity_dZ[_zco2_idx] = dliquid_density_dZco2;
  liquid.ddensity_dZ[_zch4_idx] = dliquid_density_dZch4;
  liquid.ddensity_dZ[_zh2s_idx] = dliquid_density_dZh2s;

  // Assume that liquid viscosity is just the brine viscosity
  Real liquid_viscosity, dliquid_viscosity_dp, dliquid_viscosity_dT, dliquid_viscosity_dX;
  _brine_fp.mu_dpTx(pressure,
                    temperature,
                    Xnacl,
                    liquid_viscosity,
                    dliquid_viscosity_dp,
                    dliquid_viscosity_dT,
                    dliquid_viscosity_dX);

  liquid.viscosity = liquid_viscosity;
  liquid.dviscosity_dp = dliquid_viscosity_dp;
  liquid.dviscosity_dT = dliquid_viscosity_dT;
  liquid.dviscosity_dX = dliquid_viscosity_dX;
}

Real
PorousFlowBrineCO2CH4H2S::gasSaturation(Real pressure,
                                        Real temperature,
                                        Real Xnacl,
                                        Real ZCO2,
                                        Real ZCH4,
                                        Real ZH2S,
                                        std::vector<FluidStateProperties> & fsp) const
{
  // Calculate mass fractions
  FluidStatePhaseEnum phase_state;
  massFractions(pressure, temperature, Xnacl, ZCO2, ZCH4, ZH2S, phase_state, fsp);

  // Calculate gas properties
  gasProperties(pressure, temperature, fsp);

  FluidStateProperties & liquid = fsp[_aqueous_phase_number];
  FluidStateProperties & gas = fsp[_gas_phase_number];

  // For a gas phase to exist, there must be at least one Z greater than zero,
  // and at least one component with equilibrium values of mutual solubility
  const Real Xco2 = liquid.mass_fraction[_co2_component];
  const Real Xch4 = liquid.mass_fraction[_ch4_component];
  const Real Xh2s = liquid.mass_fraction[_h2s_component];
  const Real Yco2 = gas.mass_fraction[_co2_component];
  const Real Ych4 = gas.mass_fraction[_ch4_component];
  const Real Yh2s = gas.mass_fraction[_h2s_component];

  // If Y * X > 0, then there is a gas phase
  Real K, V;
  if (Xco2 * Yco2 != 0.0)
  {
    K = Yco2 / Xco2;
    V = (ZCO2 / Xco2 - 1.0) / (K - 1.0);
  }
  else if (Xch4 * Ych4 != 0.0)
  {
    K = Ych4 / Xch4;
    V = (ZCH4 / Xch4 - 1.0) / (K - 1.0);
  }
  else // Xch4 * Ych4 != 0.0
  {
    K = Yh2s / Xh2s;
    V = (ZH2S / Xh2s - 1.0) / (K - 1.0);
  }

  // The liquid density
  const Real liquid_density = liquidDensity(pressure, temperature, Xnacl, fsp);

  // The gas saturation in the two phase case
  return V * liquid_density / (gas.density + V * (liquid_density - gas.density));
}

void
PorousFlowBrineCO2CH4H2S::twoPhaseProperties(Real pressure,
                                             Real temperature,
                                             Real Xnacl,
                                             Real ZCO2,
                                             Real ZCH4,
                                             Real ZH2S,
                                             std::vector<FluidStateProperties> & fsp) const
{
  FluidStateProperties & gas = fsp[_gas_phase_number];

  const Real gas_saturation = gasSaturation(pressure, temperature, Xnacl, ZCO2, ZCH4, ZH2S, fsp);

  // Finite difference derivatives
  const Real dp = 1.0e-2;
  Real ds = gasSaturation(pressure + dp, temperature, Xnacl, ZCO2, ZCH4, ZH2S, fsp);
  const Real ds_dp = (ds - gas_saturation) / dp;

  const Real dT = 1.0e-4;
  ds = gasSaturation(pressure, temperature + dT, Xnacl, ZCO2, ZCH4, ZH2S, fsp);
  const Real ds_dT = (ds - gas_saturation) / dT;

  const Real dX = 1.0e-6;
  ds = gasSaturation(pressure, temperature, Xnacl + dX, ZCO2, ZCH4, ZH2S, fsp);
  const Real ds_dX = (ds - gas_saturation) / dX;

  const Real dZ = 1.0e-8;
  ds = gasSaturation(pressure, temperature, Xnacl, ZCO2 + dZ, ZCH4, ZH2S, fsp);
  const Real ds_dZco2 = (ds - gas_saturation) / dZ;

  ds = gasSaturation(pressure, temperature, Xnacl, ZCO2, ZCH4 + dZ, ZH2S, fsp);
  const Real ds_dZch4 = (ds - gas_saturation) / dZ;

  ds = gasSaturation(pressure, temperature, Xnacl, ZCO2, ZCH4, ZH2S + dZ, fsp);
  const Real ds_dZh2s = (ds - gas_saturation) / dZ;

  // The gas saturation in the two phase case
  gas.saturation = gas_saturation;

  gas.dsaturation_dp = ds_dp;
  gas.dsaturation_dT = ds_dT;
  gas.dsaturation_dZ[_zco2_idx] = ds_dZco2;
  gas.dsaturation_dZ[_zch4_idx] = ds_dZch4;
  gas.dsaturation_dZ[_zh2s_idx] = ds_dZh2s;
  gas.dsaturation_dX = ds_dX;
}

void PorousFlowBrineCO2CH4H2S::checkVariables(Real /*pressure*/, Real /*temperature*/) const {}

void
PorousFlowBrineCO2CH4H2S::partialDensityCO2CH4H2S(Real temperature,
                                                  Real & partial_density_co2,
                                                  Real & partial_density_ch4,
                                                  Real & partial_density_h2s,
                                                  Real & dpartial_density_dT_co2,
                                                  Real & dpartial_density_dT_ch4,
                                                  Real & dpartial_density_dT_h2s) const
{
  // This correlation uses temperature in C
  const Real Tc = temperature - _T_c2k;
  // The parial molar volume CO2
  const Real V_co2 = 35.663 + Tc * (-0.059603 + Tc * 0.000630783);
  // 37.51 - 9.585e-2 * Tc + 8.74e-4 * Tc * Tc - 5.044e-7 * Tc * Tc * Tc; // Garcia 2001(cm3/mol)
  const Real dV_dT_co2 = -0.059603 + 2 * Tc * 0.000630783;
  //-9.585e-2 + 1.748e-3 * Tc - 1.5132e-6 * Tc * Tc;

  partial_density_co2 = 1.0e6 * _Mco2 / V_co2;
  dpartial_density_dT_co2 = -1.0e6 * _Mco2 * dV_dT_co2 / V_co2 / V_co2;

  // The parial molar volume CH4

  Real V_ch4 = 3.541 + 1.23 * 0.001 * Tc;
  V_ch4 = std::exp(V_ch4); // Rettich et. al 1981(cm3/mol)

  const Real dV_dT_ch4 = 1.23 * 0.001 * std::exp(3.541 + 1.23 * 0.001 * Tc);

  partial_density_ch4 = 1.0e6 * _Mch4 / V_ch4;
  dpartial_density_dT_ch4 = -1.0e6 * _Mch4 * dV_dT_ch4 / V_ch4 / V_ch4;

  // The parial molar volume H2S

  const Real V_h2s = 38.15; // Brelvi & O'Connell et. al 1972(cm3/mol)

  const Real dV_dT_h2s = 0;

  partial_density_h2s = 1.0e6 * _Mh2s / V_h2s;
  dpartial_density_dT_h2s = -1.0e6 * _Mh2s * dV_dT_h2s / V_h2s / V_h2s;
}

void
PorousFlowBrineCO2CH4H2S::GasCompressibilityFactor(
    Real pressure, Real temperature, Real yco2, Real ych4, Real yh2s, Real & Zg) const
{
  // Gas phase compressibility factor calculation

  Real y_hs, y_co, y_ch;

  y_hs = yh2s;
  y_co = yco2;
  y_ch = ych4;

  const Real t = temperature;       // k
  const Real p = pressure * 1.0e-5; // bar

  // Intermolecular intraction parameters(ap_ij)
  const Real ap_ch_co = 3480693.2119225934;
  const Real ap_ch_hs = 4957015.776228365;
  const Real ap_co_hs = 1413558.6307468852;
  const Real ap_hs_co = -2728833.965482667;

  // Binary intraction parameters(Kij)
  const Real k_co_ch = 0.15196821446633596;
  const Real k_co_hs = -0.037849313973062265;
  const Real k_hs_ch = 1.5716638150922622;

  //  gas constant
  const Real r = 83.14;

  //  critical temperature and pressure
  const Real tc_ch = 190.6;
  const Real tc_hs = 373.53;
  const Real tc_co = 31.05 + 273.15;

  const Real pc_ch = 46.41;
  const Real pc_hs = 89.63;
  const Real pc_co = 73.825;

  // acentric factors
  const Real w_ch = 0.012;
  const Real w_co = 0.224;
  const Real w_hs = 0.094;

  // calculation of bi
  const Real b_ch = 0.0778 * r * tc_ch / pc_ch;
  const Real b_co = 0.0778 * r * tc_co / pc_co;
  const Real b_hs = 0.0778 * r * tc_hs / pc_hs;

  // calculation of beta
  const Real bet_ch = 0.37464 + 1.54226 * w_ch - 0.26992 * std::pow(w_ch, 2);
  const Real bet_co = 0.37464 + 1.54226 * w_co - 0.26992 * std::pow(w_co, 2);
  const Real bet_hs = 0.37464 + 1.54226 * w_hs - 0.26992 * std::pow(w_hs, 2);

  // calculation of alpha_t
  const Real at_ch = std::pow((1 + bet_ch * (1 - std::pow((t / tc_ch), 0.5))), 2);
  const Real at_co = std::pow((1 + bet_co * (1 - std::pow((t / tc_co), 0.5))), 2);
  const Real at_hs = std::pow((1 + bet_hs * (1 - std::pow((t / tc_hs), 0.5))), 2);

  // calculation of a_tc
  const Real atc_ch = 0.45724 * std::pow((r * tc_ch), 2) / pc_ch;
  const Real atc_co = 0.45724 * std::pow((r * tc_co), 2) / pc_co;
  const Real atc_hs = 0.45724 * std::pow((r * tc_hs), 2) / pc_hs;

  // calculation of ai
  const Real a_ch = atc_ch * at_ch;
  const Real a_co = atc_co * at_co;
  const Real a_hs = atc_hs * at_hs;

  // Initial gas composition
  Real y1 = y_ch;
  Real a1 = a_ch;
  Real y2 = y_hs;
  Real a2 = a_hs;
  Real y3 = y_co;
  Real a3 = a_co;

  //
  const Real k13 = k_co_ch;
  const Real k23 = k_co_hs;
  const Real k12 = k_hs_ch;

  //  non- random mixing rule(equation-A8)

  const Real a_mix1 = y1 * y1 * std::pow((a1 * a1), 0.5) +
                      (1 - k12) * 2 * y1 * y2 * std::pow((a1 * a2), 0.5) +
                      (1 - k13) * 2 * y1 * y3 * std::pow((a1 * a3), 0.5);
  const Real a_mix2 = y2 * y2 * std::pow((a2 * a2), 0.5) +
                      (1 - k23) * 2 * y2 * y3 * std::pow((a2 * a3), 0.5) +
                      y3 * y3 * std::pow((a3 * a3), 0.5);
  const Real a_mix_c = a_mix1 + a_mix2;

  //  non- random mixing rule(equation-A9)
  const Real a_mix_p =
      y1 * y2 * ap_ch_hs + y1 * y3 * ap_ch_co + y2 * y3 * ap_hs_co + y3 * y2 * ap_co_hs;

  //  non- random mixing rule(equation-A7)
  const Real a_mix = a_mix_c + a_mix_p;

  //  non- random mixing rule(equation-A10)
  const Real b_mix = y_co * b_co + y_ch * b_ch + y_hs * b_hs;

  //  non- random mixing rule(equation-A13-A14)
  const Real A = a_mix * p / (r * t * r * t);
  const Real B = b_mix * p / (r * t);

  // calculation of fugacity of liquid phase
  // calculation of z-factor
  // calculation of z-factor of gas phase

  // const Real bet = b_mix * p / (r * t);
  // const Real q = a_mix / (b_mix * r * t * std::pow(t, 0.5));

  // const Real state = 1;
  Real er, z, ff;

  er = 1;
  z = 0;
  // int kk = 1;
  while (er > 0.0001)
  {
    z = z + 0.000005;
    a1 = -(1 - B);
    a2 = (A - 2 * B - 3 * B * B);
    a3 = -(A * B - B * B - std::pow(B, 3));
    ff = std::pow(z, 3) + a1 * z * z + a2 * z + a3;
    er = abs(ff);
  }

  Zg = z;
}

Real
PorousFlowBrineCO2CH4H2S::totalMassFraction(Real pressure,
                                            Real temperature,
                                            Real Xnacl,
                                            Real saturation,
                                            std::vector<Real> Yi,
                                            std::vector<FluidStateProperties> & fsp,
                                            unsigned int component,
                                            unsigned int qp) const
{
  // Saturation must be greater than zero
  if (saturation == 0.0)
    mooseError("Saturation must be greater than zero in totalMassFraction");

  FluidStateProperties & liquid = fsp[_aqueous_phase_number];

  // Mass fractions of gas components in gas phase
  const Real Yco2 = Yi[0];
  const Real Ych4 = Yi[1];
  const Real Yh2s = Yi[2];

  // Calculate the equilibrium solubilities of the gas components
  // in the liquid phase
  Real Xco2, Xch4, Xh2s, Yh2o;
  equilibriumMassFractions(pressure, temperature, Xnacl, Yco2, Ych4, Yh2s, Xco2, Xch4, Xh2s, Yh2o);

  // To calculate the gas and liquid densities, the FluidStateProperties class
  // must be populated with Xi's and yi's
  liquid.mass_fraction[_aqueous_fluid_component] = 1.0 - Xco2 - Xch4 - Xh2s;
  liquid.mass_fraction[_co2_component] = Xco2;
  liquid.mass_fraction[_ch4_component] = Xch4;
  liquid.mass_fraction[_h2s_component] = Xh2s;
  liquid.mass_fraction[_salt_component] = Xnacl;

  const Real gas_density = gasDensity(pressure, temperature, Yco2, Ych4, Yh2s);
  const Real liquid_pressure = pressure - _pc_uo.capillaryPressure(saturation, qp);
  const Real liquid_density = liquidDensity(liquid_pressure, temperature, Xnacl, fsp);

  // The mass fraction of vapor can then be calculated directly
  const Real V =
      saturation * gas_density / (saturation * gas_density + (1.0 - saturation) * liquid_density);

  Real Z;

  switch (component)
  {
    case 0:
      if (Yco2 > 0.0)
        Z = Xco2 * ((Yco2 / Xco2 - 1.0) * V + 1.0);
      else
        Z = 0.0;
      break;

    case 1:
      if (Ych4 > 0.0)
        Z = Xch4 * ((Ych4 / Xch4 - 1.0) * V + 1.0);
      else
        Z = 0.0;
      break;

    case 2:
      if (Yh2s > 0.0)
        Z = Xh2s * ((Yh2s / Xh2s - 1.0) * V + 1.0);
      else
        Z = 0.0;
      break;
  }

  return Z;
}

Real
PorousFlowBrineCO2CH4H2S::totalMassFraction(Real /*pressure*/,
                                            Real /*temperature*/,
                                            Real /*Xnacl*/,
                                            Real /*saturation*/,
                                            unsigned int /*qp*/) const
{
  mooseError("Call totalMassFraction(pressure, temperature, Xnacl, saturation, Yi, qp) instead");
}
