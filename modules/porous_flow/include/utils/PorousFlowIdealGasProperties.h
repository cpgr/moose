/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#ifndef POROUSFLOWIDEALGASPROPERTIES_H
#define POROUSFLOWIDEALGASPROPERTIES_H

#include "MooseTypes.h"
#include "MooseError.h"

namespace PorousFlowIdealGasProperties
{
  /**
   * Fluid name
   * @return fluid Name
   */
  std::string fluidName();

  /**
   * Ideal gas density as a function of  pressure, temperature and molar mass.
   *
   * @param pressure gas pressure (Pa)
   * @param temperature gas temperature (C)
   * @param molar mass gas molar mass (kg/mol)
   * @return density (kg/m^3)
   */
  Real density(Real pressure, Real temperature, Real molar_mass);

  /**
   * Derivative of the density of an ideal gas as a function of
   * pressure.
   *
   * @param temperature gas temperature (C)
   * @param molar mass gas molar mass (kg/mol)
   * @return derivative of CO2 density (kg/m^3) with respect to pressure
   */
  Real dDensity_dP(Real temperature, Real molar_mass);

  /**
   * Derivative of the density of an ideal gas as a function of
   * temperature.
   *
   * @param pressure gas pressure (Pa)
   * @param temperature gas temperature (C)
   * @param molar mass gas molar mass (kg/mol)
   * @return derivative of CO2 density (kg/m^3) with respect to pressure
   */
  Real dDensity_dT(Real pressure, Real temperature, Real molar_mass);

 /// Conversion of temperature from Celcius to Kelvin
 const Real _t_c2k = 273.15;
 /// Gas constant
 const Real _R = 8.3144621;
}

#endif // POROUSFLOWIDEALGASPROPERTIES_H
