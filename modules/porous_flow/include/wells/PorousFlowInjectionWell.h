//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "PorousFlowWellBase.h"

/**
 * Surrogate object for a rate-controlled injection well.  PorousFlowWellAction expands this into
 * one PorousFlowPointSourceFromPostprocessor per fluid component (plus, when nonisothermal, a
 * PorousFlowPointEnthalpySourceFromPostprocessor) and the associated rate reporters.
 */
class PorousFlowInjectionWell : public PorousFlowWellBase
{
public:
  static InputParameters validParams();

  PorousFlowInjectionWell(const InputParameters & parameters);
};
