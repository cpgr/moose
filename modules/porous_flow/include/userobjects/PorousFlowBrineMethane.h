//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef POROUSFLOWBRINEMETHANE_H
#define POROUSFLOWBRINEMETHANE_H

#include <vector>
#include <math.h>
#include "PorousFlowFluidStateBase.h"

class BrineFluidProperties;
class SinglePhaseFluidPropertiesPT;
class PorousFlowBrineMethane;

template <>
InputParameters validParams<PorousFlowBrineMethane>();

/**
 * Specialized class for brine and Methane including calculation of mutual
 * solubility of the two fluids using the high-accuracy formulation of
 * Duan and Mao, A thermodynamic model for calculating methane solubility, density
 * and gas phase composition of methane-bearing aqueous fluids from 273 to 523 K
 * and from 1 to 2000 bar, Geochimica et Cosmochimica Acta, 70, 3369-3386 (2006)
 */
class PorousFlowBrineMethane : public PorousFlowFluidStateBase
{
public:
  PorousFlowBrineMethane(const InputParameters & parameters);

  /**
   * Name of FluidState
   * @return brine-methane
   */
  virtual std::string fluidStateName() const;

  /**
   * Convert NaCl mass fraction to mol fraction in brine solution
   * @param massFrac NaCl mass fraction (g/g)
   * @return NaCl mol fraction (mol/mol)
   */
  Real brineMassToMolFraction(const Real massFrac) const;

  /**
   * Convert NaCl mol fraction to mass fraction in brine solution
   * @param molFrac NaCl mol fraction (mol/mol)
   * @return NaCl mass fraction (g/g)
   */
  Real brineMolToMassFraction(const Real molFrac) const;

  /**
   * Convert NaCl mol fraction to molality in brine solution
   * @param molFrac NaCl mol fraction (mol/mol)
   * @return NaCl molality (mol/kg of H2O)
   */
  Real brineMolFractionToMolality(const Real molFrac) const;

  /**
   * Convert NaCl molality to mol fraction in brine solution
   * @param molality NaCl molality (mol/kg of H2O)
   * @return NaCl mol fraction (mol/mol)
   */
  Real brineMolalityToMolFraction(const Real molality) const;

  void thermophysicalProperties(Real pressure,
                                Real temperature,
                                Real Xnacl,
                                Real Z,
                                std::vector<FluidStateProperties> & fsp) const;
  /**
   * Mass fractions of methane and brine
   *
   * @param pressure phase pressure (Pa)
   * @param temperature phase temperature (C)
   * @param Xnacl NaCl mass fraction (kg/kg)
   * @param[out] Xch4l mass fraction of methane in liquid (kg/kg)
   * @param[out] dXch4l_dp derivative of mass fraction of methane in liquid wrt pressure
   * @param[out] dXch4l_dT derivative of mass fraction of methane in liqiud wrt temperature
   * @param[out] Xh2og mass fraction of H2O in gas (kg/kg)
   * @param[out] dXh2ogl_dp derivative of mass fraction of H2O in gas wrt pressure
   * @param[out] dXh2og_dT derivative of mass fraction of H2O in gas wrt temperature
   */
  void equilibriumMassFractions(Real pressure,
                                Real temperature,
                                Real Xnacl,
                                Real & Xch4l,
                                Real & dXch4l_dp,
                                Real & dXch4l_dT,
                                Real & Xh2og,
                                Real & dXh2og_dp,
                                Real & dXh2og_dT) const;

  /**
   * Mass fractions of methane and brine
   *
   * @param pressure phase pressure (Pa)
   * @param temperature phase temperature (C)
   * @param Xnacl NaCl mass fraction (kg/kg)
   * @param[out] Xch4l mass fraction of methane in liquid (kg/kg)
   * @param[out] Xh2og mass fraction of H2O in gas (kg/kg)
   */
  void equilibriumMassFractions(
      Real pressure, Real temperature, Real Xnacl, Real & Xch4l, Real & Xh2og) const;

  /**
   * Mass fractions of methane and H2O in both phases, as well as derivatives wrt
   * PorousFlow variables. Values depend on the phase state (liquid, gas or two phase)
   *
   * @param pressure phase pressure (Pa)
   * @param temperature phase temperature (C)
   * @param Xnacl NaCl mass fraction (kg/kg)
   * @param Z total mass fraction of methane component
   * @param[out] PhaseStateEnum current phase state
   * @param[out] FluidStateMassFractions data structure
   */
  void massFractions(Real pressure,
                     Real temperature,
                     Real Xnacl,
                     Real Z,
                     FluidStatePhaseEnum & phase_state,
                     std::vector<FluidStateProperties> & fsp) const;

  /**
   * Thermophysical properties of the gaseous state
   *
   * @param pressure gas pressure (Pa)
   * @param temperature temperature (C)
   * @param Xnacl NaCl mass fraction (kg/kg)
   * @param[out] FluidStateDensity data structure
   */
  void
  gasProperties(Real pressure, Real temperature, std::vector<FluidStateProperties> & fsp) const;

  /**
   * Thermophysical properties of the liquid state
   *
   * @param pressure liquid pressure (Pa)
   * @param temperature temperature (C)
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
   * @param temperature phase temperature (C)
   * @param Xnacl NaCl mass fraction (kg/kg)
   * @param Z total mass fraction of methane component
   * @param[out] FluidStateSaturation data structure
   */
  void saturationTwoPhase(Real pressure,
                          Real temperature,
                          Real Xnacl,
                          Real Z,
                          std::vector<FluidStateProperties> & fsp) const;

  /**
   * Fugacity coefficient for methane
   *
   * @param pressure gas pressure (Pa)
   * @param temperature temperature (K)
   * @param[out] fch4 fugacity coefficient for Methane
   * @param[out] dfch4_dp derivative of fugacity coefficient wrt pressure
   * @param[out] dfch4_dT derivative of fugacity coefficient wrt temperature
   */
  void fugacityCoefficientMethane(
      Real pressure, Real temperature, Real & fch4, Real & dfch4_dp, Real & dfch4_dT) const;

  /**
   * Fugacity coefficient for H2O
   *
   * @param pressure gas pressure (Pa)
   * @param temperature temperature (K)
   * @param[out] fh2o fugacity coefficient for H2O
   * @param[out] dfh2o_dp derivative of fugacity coefficient wrt pressure
   * @param[out] dfh2o_dT derivative of fugacity coefficient wrt temperature
   */
  void fugacityCoefficientH2O(
      Real pressure, Real temperature, Real & fh2o, Real & dfh2o_dp, Real & dfh2o_dT) const;

  Real molFractionOfWaterInGas(const Real pressure, const Real temperature, const Real Xh2o) const;

  /**
   * Activity coefficient for methane in brine
   *
   * @param pressure phase pressure (Pa)
   * @param temperature phase temperature (K)
   * @param bnacl molality of salt (mol/kg)
   */
  Real activityCoefficientMolality(Real pressure, Real temperature, Real bnacl) const;

  /**
   * Activity coefficient for methane in brine
   *
   * @param pressure phase pressure (Pa)
   * @param temperature phase temperature (K)
   * @param Xnacl mol fraction (g/g)
   */
  Real activityCoefficientMolFraction(Real pressure, Real temperature, Real Xnacl) const;

  /**
   * Activity coefficient for methane in brine, with molality as NaCl concentration.
   *
   * @param pressure phase pressure (Pa)
   * @param temperature phase temperature (K)
   * @param Xnacl salt mass fraction (g/g)
   * @param[out] gamma activity coefficient for Methane in brine (output)
   * @param[out] dgamma_dp derivative of activity coefficient wrt pressure
   * @param[out] dgamma_dT derivative of activity coefficient wrt temperature
   */
  void activityCoefficient(Real pressure,
                           Real temperature,
                           Real Xnacl,
                           Real & gamma,
                           Real & dgamma_dp,
                           Real & dgamma_dT) const;

  /**
   * Equilibrium constant for H2O
   *
   * @param temperature temperature (C)
   * @param[out] kh2o equilibrium constant for H2O
   * @param[out] dkh2o_dT derivative of equilibrium constant wrt temperature
   */
  void equilibriumConstantH2O(Real temperature, Real & kh2o, Real & dkh2o_dT) const;

  /**
   * Equilibrium constant for methane
   *
   * @param temperature temperature (C)
   * @param[out] kch4 equilibrium constant for Methane
   * @param[out] dkch4_dT derivative of equilibrium constant wrt temperature
   */
  void equilibriumConstantMethane(Real temperature, Real & kch4, Real & dkch4_dT) const;

  /**
   * Total mass fraction of methane summed over all phases in the two-phase state
   *
   * @param pressure gas pressure (Pa)
   * @param temperature temperature (K)
   * @param Xnacl NaCl mass fraction (kg/kg)
   * @param saturation gas saturation (-)
   * @return total mass fraction z (-)
   */
  Real totalMassFraction(Real pressure, Real temperature, Real Xnacl, Real saturation) const;

  /**
   * Compute methane solubility in liquid, based on Ref. 1, Eq. 8.
   *
   * @param pressure gas pressure (Pa)
   * @param temperature temperature (K)
   * @param Xnacl NaCl molality (mol/kg)
   * @return methane solubility in the liquid (mol/kg)
   */
  Real methaneSolubilityInLiquid(Real pressure, Real temperature, Real bnacl) const;

protected:
  /// Check the input variables
  void checkVariables(Real pressure, Real temperature) const;

  /// Fluid properties UserObject for water
  const BrineFluidProperties & _brine_fp;
  /// Fluid properties UserObject for methane
  const SinglePhaseFluidPropertiesPT & _ch4_fp;
  /// Fluid properties UserObject for H20
  const SinglePhaseFluidPropertiesPT & _water_fp;
  /// Molar mass of water (kg/mol)
  const Real _Mh2o;
  /// Inverse of molar mass of H2O (mol/kg)
  const Real _invMh2o;
  /// Molar mass of methane (kg/mol)
  const Real _Mch4;
  /// Molar mass of NaCL
  const Real _Mnacl;
  /// Molar gas constant in bar L /(K mol)
  const Real _Rbar;
};

#endif // POROUSFLOWBRINEMETHANE_H
