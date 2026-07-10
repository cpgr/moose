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
 * Surrogate object for a production well operating on a deliverability curve.  PorousFlowWellAction
 * expands this into one PorousFlowPolyLineSink per fluid component, each with its own
 * PorousFlowSumQuantity, and the associated rate and cumulative reporters.
 */
class PorousFlowProductionWell : public PorousFlowWellBase
{
public:
  static InputParameters validParams();

  PorousFlowProductionWell(const InputParameters & parameters);
};
