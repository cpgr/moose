//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef POROUSFLOWBRINECO2CH4H2S_H
#define POROUSFLOWBRINECO2CH4H2S_H

#include "PorousFlowFluidStateBase.h"

class BrineFluidProperties;
class SinglePhaseFluidPropertiesPT;
class PorousFlowBrineCO2CH4H2S;

template <>
InputParameters validParams<PorousFlowBrineCO2CH4H2S>();

/**
 * Specialized class for brine, CO2, CH4 and H2S including calculation of mutual
 * solubility of the all components.
 *
 * Notation convention
 * Throughout this class, both mole fractions and mass fractions will be used.
 * The following notation will be used:
 * yk: mole fraction of component k in the gas phase
 * xk: mole fraction of component k in the liquid phase
 * Yk: mass fraction of component k in the gas phase
 * Xk: mass fraction of component k in the liquid phase
 */
class PorousFlowBrineCO2CH4H2S : public PorousFlowFluidStateBase
{
public:
  PorousFlowBrineCO2CH4H2S(const InputParameters & parameters);

  virtual std::string fluidStateName() const override;

  virtual void thermophysicalProperties(Real pressure,
                                        Real temperature,
                                        Real Xnacl,
                                        std::vector<const VariableValue *> & Z,
                                        unsigned int qp,
                                        std::vector<FluidStateProperties> & fsp) const override;

  /**
   * Mole fractions of CO2, CH4 and H2S in brine and water vapor in gas phase at equilibrium.
   *
   *
   * @param pressure phase pressure (Pa)
   * @param temperature phase temperature (K)
   * @param Xnacl NaCl mass fraction (kg/kg)
   * @param[out] xco2 mole fraction of CO2 in liquid
   * @param[out] xch4 mole fraction of CO2 in liquid
   * @param[out] xh2s mole fraction of CO2 in liquid
   * @param[out] yh2o mass fraction of mole in gas
   */
  void equilibriumMoleFractions(Real pressure,
                                Real temperature,
                                Real Xnacl,
                                Real yco2,
                                Real ych4,
                                Real yh2s,
                                Real & xco2,
                                Real & xch4,
                                Real & xh2s,
                                Real & yh2o) const;

  /**
   * Mass fractions of CO2, CH4 and H2S in brine and water vapor in gas phase at equilibrium.
   *
   *
   * @param pressure phase pressure (Pa)
   * @param temperature phase temperature (K)
   * @param Xnacl NaCl mass fraction (kg/kg)
   * @param Yco2 mass fraction of CO2 in gas
   * @param Ych4 mass fraction of CO2 in gas
   * @param Yh2s mass fraction of CO2 in gas
   * @param[out] Xco2 mass fraction of CO2 in liquid
   * @param[out] Xch4 mass fraction of CO2 in liquid
   * @param[out] Xh2s mass fraction of CO2 in liquid
   * @param[out] Yh2o mass fraction of mole in gas
   */
  void equilibriumMassFractions(Real pressure,
                                Real temperature,
                                Real Xnacl,
                                Real Yco2,
                                Real Ych4,
                                Real Yh2s,
                                Real & Xco2,
                                Real & Xch4,
                                Real & Xh2s,
                                Real & Yh2o) const;

  /**
   * Mass fractions of CO2, CH4 and H2S in brine and water vapor in gas phase
   *
   * @param pressure phase pressure (Pa)
   * @param temperature phase temperature (K)
   * @param Xnacl NaCl mass fraction (kg/kg)
   * @param Z total mass fraction of CO2 component
   * @param[out] PhaseStateEnum current phase state
   * @param[out] FluidStateMassFractions data structure
   */
  void massFractions(Real pressure,
                     Real temperature,
                     Real Xnacl,
                     Real ZCO2,
                     Real ZCH4,
                     Real ZH2S,
                     FluidStatePhaseEnum & phase_state,
                     std::vector<FluidStateProperties> & fsp) const;

  void massFractions(Real pressure,
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
                     Real & Yh2s) const;

  /**
   * Gas density
   *
   * @param pressure gas pressure (Pa)
   * @param temperature temperature (K)
   * @param yco2 mole fraction of CO2 in gas phase
   * @param ych4 mole fraction of CH4 in gas phase
   * @param yh2s mole fraction of H2S in gas phase
   * @return gas density (kg/m^3)
   */
  Real gasDensity(Real pressure, Real temperature, Real yco2, Real ych4, Real yh2s) const;

  /**
   * Thermophysical properties of the gaseous state
   *
   * @param pressure gas pressure (Pa)
   * @param temperature temperature (K)
   * @param[out] FluidStateDensity data structure
   */
  void
  gasProperties(Real pressure, Real temperature, std::vector<FluidStateProperties> & fsp) const;

  /**
   * Liquid density
   *
   * @param pressure gas pressure (Pa)
   * @param temperature temperature (K)
   * @param Xnacl NaCl mass fraction (kg/kg)
   * @param FluidStateDensity data structure
   * @return liquid density (kg/m^3)
   */
  Real liquidDensity(Real pressure,
                     Real temperature,
                     Real Xnacl,
                     std::vector<FluidStateProperties> & fsp) const;

  /**
   * Thermophysical properties of the liquid state
   *
   * @param pressure liquid pressure (Pa)
   * @param temperature temperature (K)
   * @param Xnacl NaCl mass fraction (kg/kg)
   * @param[out] FluidStateDensity data structure
   */
  void liquidProperties(Real pressure,
                        Real temperature,
                        Real Xnacl,
                        std::vector<FluidStateProperties> & fsp) const;

  /**
   * Gas and liquid properties in the two-phase region
   *
   * @param pressure gas pressure (Pa)
   * @param temperature phase temperature (K)
   * @param Xnacl NaCl mass fraction (kg/kg)
   * @param Z total mass fraction of CO2 component
   * @param[out] FluidStateSaturation data structure
   */
  void twoPhaseProperties(Real pressure,
                          Real temperature,
                          Real Xnacl,
                          Real ZCO2,
                          Real ZCH4,
                          Real ZH2S,
                          std::vector<FluidStateProperties> & fsp) const;

  /**
   * Gas saturation in the two phase region
   *
   * @param pressure gas pressure (Pa)
   * @param temperature phase temperature (K)
   * @param Xnacl NaCl mass fraction (kg/kg)
   * @param FluidStateSaturation data structure
   * @return gas saturation
   */
  Real gasSaturation(Real pressure,
                     Real temperature,
                     Real Xnacl,
                     Real ZCO2,
                     Real ZCH4,
                     Real ZH2S,
                     std::vector<FluidStateProperties> & fsp) const;

  void GasCompressibilityFactor(
      Real pressure, Real temperature, Real yco2, Real yh2s, Real ych4, Real & Zg) const;

  void partialDensityCO2CH4H2S(Real temperature,
                               Real & partial_density_co2,
                               Real & partial_density_ch4,
                               Real & partial_density_h2s,
                               Real & dpartial_density_dT_co2,
                               Real & dpartial_density_dT_ch4,
                               Real & dpartial_density_dT_h2s) const;

  virtual Real totalMassFraction(
      Real pressure, Real temperature, Real Xnacl, Real saturation, unsigned int qp) const override;

  virtual Real totalMassFraction(Real pressure,
                                 Real temperature,
                                 Real Xnacl,
                                 Real saturation,
                                 std::vector<Real> Yi,
                                 std::vector<FluidStateProperties> & fsp,
                                 unsigned int component,
                                 unsigned int qp) const;

protected:
  /// Check the input variables
  void checkVariables(Real pressure, Real temperature) const;

  /// CO2 fluid component index
  const unsigned int _co2_component;
  /// CH4 fluid component index
  const unsigned int _ch4_component;
  /// H2S fluid component index
  const unsigned int _h2s_component;
  /// Salt component index
  const unsigned int _salt_component;
  /// Z index for CO2
  const unsigned int _zco2_idx;
  /// Z index for CH4
  const unsigned int _zch4_idx;
  /// Z index for H2S
  const unsigned int _zh2s_idx;
  /// Fluid properties UserObject for water
  const BrineFluidProperties & _brine_fp;
  /// Fluid properties UserObject for the CO2
  const SinglePhaseFluidPropertiesPT & _co2_fp;
  /// Fluid properties UserObject for the CH4
  const SinglePhaseFluidPropertiesPT & _ch4_fp;
  /// Fluid properties UserObject for H20
  const SinglePhaseFluidPropertiesPT & _water_fp;
  /// Molar mass of water (kg/mol)
  const Real _Mh2o;
  /// Inverse of molar mass of H2O (mol/kg)
  const Real _invMh2o;
  /// Molar mass of CO2 (kg/mol)
  const Real _Mco2;
  /// Molar mass of H2S (kg/mol)
  const Real _Mh2s;
  /// Molar mass of CH4 (kg/mol)
  const Real _Mch4;
  /// Molar mass of NaCL
  const Real _Mnacl;
};

#endif // POROUSFLOWBRINECO2CH4H2S_H
