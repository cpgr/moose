//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef POROUSFLOWWATERNCG_H
#define POROUSFLOWWATERNCG_H

#include "PorousFlowFluidStateBase.h"

class SinglePhaseFluidProperties;
class PorousFlowWaterNCG;

template <>
InputParameters validParams<PorousFlowWaterNCG>();

/**
 * Specialized class for water and a non-condensable gas (NCG)
 * Includes dissolution of gas in liquid water phase using Henry's law
 *
 * Notation convention
 * Throughout this class, both mole fractions and mass fractions will be used.
 * The following notation will be used:
 * yk: mole fraction of component k in the gas phase
 * xk: mole fraction of component k in the liquid phase
 * Yk: mass fraction of component k in the gas phase
 * Xk: mass fraction of component k in the liquid phase
 */
class PorousFlowWaterNCG : public PorousFlowFluidStateBase
{
public:
  PorousFlowWaterNCG(const InputParameters & parameters);

  virtual std::string fluidStateName() const override;

  void thermophysicalProperties(Real pressure,
                                Real temperature,
                                Real Xnacl,
                                Real Z,
                                unsigned int qp,
                                std::vector<FluidStateProperties> & fsp) const override;
  /**
   * Mass fractions of NCG in liquid phase and H2O in gas phase at thermodynamic
   * equilibrium. Calculated using Henry's law (for NCG component), and Raoult's
   * law (for water).
   *
   * @param pressure phase pressure (Pa)
   * @param temperature phase temperature (C)
   * @param[out] Xncg mass fraction of NCG in liquid (kg/kg)
   * @param[out] dXncg_dp derivative of mass fraction of NCG in liquid wrt pressure
   * @param[out] dXncg_dT derivative of mass fraction of NCG in liqiud wrt temperature
   * @param[out] Yh2o mass fraction of H2O in gas (kg/kg)
   * @param[out] dYh2o_dp derivative of mass fraction of NCG in gas wrt pressure
   * @param[out] dYh2o_dT derivative of mass fraction of NCG in gas wrt temperature
   */
  void equilibriumMassFractions(Real pressure,
                                Real temperature,
                                Real & Xncg,
                                Real & dXncg_dp,
                                Real & dXncg_dT,
                                Real & Yh2o,
                                Real & dYh2o_dp,
                                Real & dYh2o_dT) const;

  /**
   * Mass fractions of NCG and H2O in both phases, as well as derivatives wrt
   * PorousFlow variables. Values depend on the phase state (liquid, gas or two phase)
   *
   * @param pressure phase pressure (Pa)
   * @param temperature phase temperature (C)
   * @param Z total mass fraction of NCG component
   * @param[out] PhaseStateEnum current phase state
   * @param[out] FluidStateProperties data structure
   */
  void massFractions(Real pressure,
                     Real temperature,
                     Real Z,
                     FluidStatePhaseEnum & phase_state,
                     std::vector<FluidStateProperties> & fsp) const;

  /**
   * Gas properties - density, viscosity and enthalpy
   *
   * @param pressure gas pressure (Pa)
   * @param temperature temperature (C)
   * @param[out] FluidStateProperties data structure
   */
  void
  gasProperties(Real pressure, Real temperature, std::vector<FluidStateProperties> & fsp) const;

  /**
   * Liquid properties - density, viscosity and enthalpy
   * Note: The pressure here is the liquid pressure. In this class, enthalpy includes a
   * contribution due to the enthalpy of dissolution of the NCG into the liquid phase. As
   * a result, the derivatives can include a dependence on the capillary pressure, so this
   * method should be called after the saturation is calculated for the two phase case
   * ie: after calling saturationTwoPhase(). For the single phase liquid case, it is ok to
   * call this method by itself, as gas saturation is initialized to zero.
   *
   * @param pressure liquid pressure (Pa)
   * @param temperature temperature (C)
   * @param[out] FluidStateProperties data structure
   */
  void
  liquidProperties(Real pressure, Real temperature, std::vector<FluidStateProperties> & fsp) const;

  /**
   * Density of the gas state (assumes no effect due to evaporation of brine)
   *
   * @param pressure pressure (Pa)
   * @param temperature temperature (K)
   * @return gas density (kg/m^3)
   */
  Real gasDensity(Real pressure, Real temperature) const;

  /**
   * Density of the liquid state
   * Note: The pressure here is the liquid pressure. As a result, the derivatives can
   * include a dependence on saturation due to the capillary pressure, so this
   * method should be called after the saturation is calculated for the two phase case
   * ie: after calling saturation(). For the single phase liquid case, it is ok to
   * call this method by itself, as gas saturation is initialized to zero.
   *
   * @param pressure liquid pressure (Pa)
   * @param temperature temperature (K)
   * @return liquid density (kg/m^3)
   */
  Real liquidDensity(Real pressure, Real temperature) const;

  /**
   * Gas saturation in the two-phase region
   *
   * @param pressure gas pressure (Pa)
   * @param temperature phase temperature (K)
   * @param Z total mass fraction of CO2 component
   * @param Xncg NCG mass fraction in liquid (kg/kg)
   * @param Yncg mass fraction of NCG in gas (kg/kg)
   * @return gas saturation (-)
   */
  Real saturation(Real pressure, Real temperature, Real Z, Real Xncg, Real Yncg) const;

  /**
   * Gas and liquid properties in the two-phase region
   *
   * @param pressure gas pressure (Pa)
   * @param temperature phase temperature (C)
   * @param Z total mass fraction of NCG component
   * @param[out] FluidStateProperties data structure
   */
  void twoPhaseProperties(Real pressure,
                          Real temperature,
                          Real Z,
                          std::vector<FluidStateProperties> & fsp) const;

  /**
   * Enthalpy of dissolution of NCG in water calculated using Henry's constant
   * From Himmelblau, Partial molal heats and entropies of solution for gases dissolved
   * in water from the freezing to the near critical point, J. Phys. Chem. 63 (1959)
   *
   * @param temperature fluid temperature (K)
   * @param[out] hdis enthalpy of dissolution (J/kg)
   * @param[out] dhdis_dT derivative of enthalpy of dissolution wrt temperature
   */
  void enthalpyOfDissolution(Real temperature, Real & hdis, Real & dhdis_dT) const;

  virtual Real totalMassFraction(
      Real pressure, Real temperature, Real Xnacl, Real saturation, unsigned int qp) const override;

protected:
  /**
   * Convert mole fraction to mass fraction
   * @param xmol mole fraction
   * @return mass fraction
   */
  Real moleFractionToMassFraction(Real xmol) const;

  /**
   * Check that the temperature is between the triple and critical values
   * @param temperature fluid temperature (K)
   */
  void checkVariables(Real temperature) const;

  /// Fluid properties UserObject for water
  const SinglePhaseFluidProperties & _water_fp;
  /// Fluid properties UserObject for the NCG
  const SinglePhaseFluidProperties & _ncg_fp;
  /// Molar mass of water (kg/mol)
  const Real _Mh2o;
  /// Molar mass of non-condensable gas (kg/mol)
  const Real _Mncg;
  /// Triple point temperature of water (K)
  const Real _water_triple_temperature;
  /// Critical temperature of water (K)
  const Real _water_critical_temperature;
};

#endif // POROUSFLOWWATERNCG_H
