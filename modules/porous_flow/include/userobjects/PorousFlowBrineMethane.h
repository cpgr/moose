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

  void thermophysicalProperties(Real pressure,
                                Real temperature,
                                Real xnacl,
                                Real z,
                                std::vector<FluidStateProperties> & fsp) const;
  /**
   * Mass fractions of methane and brine
   *
   * @param pressure phase pressure (Pa)
   * @param temperature phase temperature (C)
   * @param xnacl NaCl mass fraction (kg/kg)
   * @param[out] xch4l mass fraction of methane in liquid (kg/kg)
   * @param[out] dxch4l_dp derivative of mass fraction of methane in liquid wrt pressure
   * @param[out] dxch4l_dT derivative of mass fraction of methane in liqiud wrt temperature
   * @param[out] xh2og mass fraction of H2O in gas (kg/kg)
   * @param[out] dh2ogl_dp derivative of mass fraction of H2O in gas wrt pressure
   * @param[out] dh2og_dT derivative of mass fraction of H2O in gas wrt temperature
   */
  void equilibriumMassFractions(Real pressure,
                                Real temperature,
                                Real xnacl,
                                Real & xch4l,
                                Real & dxch4l_dp,
                                Real & dxch4l_dT,
                                Real & xh2og,
                                Real & dxh2og_dp,
                                Real & dxh2og_dT) const;

  /**
   * Mass fractions of methane and H2O in both phases, as well as derivatives wrt
   * PorousFlow variables. Values depend on the phase state (liquid, gas or two phase)
   *
   * @param pressure phase pressure (Pa)
   * @param temperature phase temperature (C)
   * @param xnacl NaCl mass fraction (kg/kg)
   * @param z total mass fraction of methane component
   * @param[out] PhaseStateEnum current phase state
   * @param[out] FluidStateMassFractions data structure
   */
  void massFractions(Real pressure,
                     Real temperature,
                     Real xnacl,
                     Real z,
                     FluidStatePhaseEnum & phase_state,
                     std::vector<FluidStateProperties> & fsp) const;

  /**
   * Thermophysical properties of the gaseous state
   *
   * @param pressure gas pressure (Pa)
   * @param temperature temperature (C)
   * @param xnacl NaCl mass fraction (kg/kg)
   * @param[out] FluidStateDensity data structure
   */
  void
  gasProperties(Real pressure, Real temperature, std::vector<FluidStateProperties> & fsp) const;

  /**
   * Thermophysical properties of the liquid state
   *
   * @param pressure liquid pressure (Pa)
   * @param temperature temperature (C)
   * @param xnacl NaCl mass fraction (kg/kg)
   * @param[out] FluidStateDensity data structure
   */
  void liquidProperties(Real pressure,
                        Real temperature,
                        Real xnacl,
                        std::vector<FluidStateProperties> & fsp) const;

  /**
   * Gas and liquid saturations for the two-phase region
   *
   * @param pressure gas pressure (Pa)
   * @param temperature phase temperature (C)
   * @param xnacl NaCl mass fraction (kg/kg)
   * @param z total mass fraction of methane component
   * @param[out] FluidStateSaturation data structure
   */
  void saturationTwoPhase(Real pressure,
                          Real temperature,
                          Real xnacl,
                          Real z,
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

  /**
   * Activity coefficient for methane in brine
   *
   * @param pressure phase pressure (Pa)
   * @param temperature phase temperature (K)
   * @param xnacl salt mass fraction (kg/kg)
   * @param[out] gamma activity coefficient for Methane in brine (output)
   * @param[out] dgamma_dp derivative of activity coefficient wrt pressure
   * @param[out] dgamma_dT derivative of activity coefficient wrt temperature
   */
  void activityCoefficient(Real pressure,
                           Real temperature,
                           Real xnacl,
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
   * @param xnacl NaCl mass fraction (kg/kg)
   * @param saturation gas saturation (-)
   * @return total mass fraction z (-)
   */
  Real totalMassFraction(Real pressure, Real temperature, Real xnacl, Real saturation) const;

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
};

#endif // POROUSFLOWBRINEMETHANE_H
