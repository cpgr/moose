//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "MooseObject.h"

/**
 * Surrogate object for a well in the [Wells] block.
 *
 * Like PorousFlowSinkBC, this object specifies the parameters that the user must provide, but
 * performs no computation itself: PorousFlowWellAction reads these parameters and builds the
 * concrete DiracKernels, UserObjects and Postprocessors that model the well.
 */
class PorousFlowWellBase : public MooseObject
{
public:
  static InputParameters validParams();

  /// Parameters shared by every well type (injectors and producers alike)
  static InputParameters validParamsCommon();

  /**
   * Parameters that describe a bottom-hole-pressure constraint.  These are declared so that an
   * input using them receives an explanatory error from PorousFlowWellAction rather than having
   * them silently ignored: they require the implicit scalar-BHP well of phase 1b.
   */
  static InputParameters validParamsNotYetImplemented();

  PorousFlowWellBase(const InputParameters & parameters);
};
