//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "PorousFlowLineGeometry.h"

/**
 * Line source (or sink) that adds (removes) fluid at a total mass flux rate specified by a
 * postprocessor, divided between the completions in proportion to their length
 */
class PorousFlowLineSourceFromPostprocessor : public PorousFlowLineGeometry
{
public:
  static InputParameters validParams();

  PorousFlowLineSourceFromPostprocessor(const InputParameters & parameters);

  virtual Real computeQpResidual() override;

protected:
  /// The total mass flux over the whole line (kg/s)
  const PostprocessorValue & _mass_flux;
};
