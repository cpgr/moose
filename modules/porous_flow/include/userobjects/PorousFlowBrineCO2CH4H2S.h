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

  /**
   * Name of FluidState
   * @return brine-co2
   */
  virtual std::string fluidStateName() const;

  /**
   * Determines the complete thermophysical state of the system for a given set of
   * primary variables
   *
   * @param pressure gas phase pressure (Pa)
   * @param temperature fluid temperature (K)
   * @param Xnacl mass fraction of NaCl
   * @param Z total mass fraction of CO2 component
   * @param[out] fsp the FluidStateProperties struct containing all properties
   */
  void thermophysicalProperties(Real pressure,
                                Real temperature,
                                Real Xnacl,
                                Real ZCO2,
                                Real ZH2S,
                                Real ZCH4,
                                std::vector<FluidStateProperties> & fsp) const;

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
                                Real yh2s,
                                Real ych4,
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
   * @param[out] xco2 mole fraction of CO2 in liquid
   * @param[out] xch4 mole fraction of CO2 in liquid
   * @param[out] xh2s mole fraction of CO2 in liquid
   * @param[out] yh2o mass fraction of mole in gas
   */
  void equilibriumMassFractions(Real pressure,
                                Real temperature,
                                Real Xnacl,
                                Real yco2,
                                Real yh2s,
                                Real ych4,
                                Real & xco2,
                                Real & xch4,
                                Real & xh2s,
                                Real & yh2o) const;

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
                     Real ZH2S,
                     Real ZCH4,
                     FluidStatePhaseEnum & phase_state,
                     std::vector<FluidStateProperties> & fsp) const;

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
   * Gas and liquid saturations for the two-phase region
   *
   * @param pressure gas pressure (Pa)
   * @param temperature phase temperature (K)
   * @param Xnacl NaCl mass fraction (kg/kg)
   * @param Z total mass fraction of CO2 component
   * @param[out] FluidStateSaturation data structure
   */
  void saturationTwoPhase(Real pressure,
                          Real temperature,
                          Real Xnacl,
                          Real Z,
                          std::vector<FluidStateProperties> & fsp) const;

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
  /// Fluid properties UserObject for water
  const BrineFluidProperties & _brine_fp;
  /// Fluid properties UserObject for the CO2
  const SinglePhaseFluidPropertiesPT & _co2_fp;
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
