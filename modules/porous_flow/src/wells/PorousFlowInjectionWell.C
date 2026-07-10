//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PorousFlowInjectionWell.h"

registerMooseObject("PorousFlowApp", PorousFlowInjectionWell);

InputParameters
PorousFlowInjectionWell::validParams()
{
  InputParameters params = PorousFlowWellBase::validParams();
  params.addRequiredParam<Point>("point", "The x,y,z coordinates of the injection point");
  params.addRequiredParam<FunctionName>(
      "rate",
      "The TOTAL injection rate in kg/s, positive for injection.  This is split between the fluid "
      "components according to 'injected_mass_fractions'.");
  params.addRequiredParam<std::vector<Real>>(
      "injected_mass_fractions",
      "The mass fraction of each fluid component in the injected stream, in fluid-component order. "
      " Must have one entry per fluid component and sum to 1.");
  params.addParam<VariableName>("temperature",
                                "The temperature variable.  Supply this, along with 'fp' and "
                                "'injection_temperature', for a nonisothermal simulation.");
  params.addParam<UserObjectName>(
      "fp", "The fluid-properties UserObject used to compute the enthalpy of the injected fluid");
  params.addParam<FunctionName>("injection_temperature",
                                "The temperature of the injected fluid, in K");

  // Declared so that a multi-completion injector gets an explanatory error rather than an
  // "unused parameter" one.  Completions along a trajectory need the Peaceman well of phase 1b.
  params.addParam<std::string>("point_file",
                               "File of completion coordinates.  Not supported for fixed-rate "
                               "injectors: see the error message.");

  params.addClassDescription(
      "A rate-controlled injection well, injecting a fluid of fixed composition at a point.  The "
      "rate is prescribed, not computed from a bottom-hole pressure, so the well has no "
      "injectivity limit.");
  return params;
}

PorousFlowInjectionWell::PorousFlowInjectionWell(const InputParameters & parameters)
  : PorousFlowWellBase(parameters)
{
}
