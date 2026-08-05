//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PorousFlowLineSourceFromPostprocessor.h"

registerMooseObject("PorousFlowApp", PorousFlowLineSourceFromPostprocessor);

InputParameters
PorousFlowLineSourceFromPostprocessor::validParams()
{
  InputParameters params = PorousFlowLineGeometry::validParams();
  params.addRequiredParam<PostprocessorName>(
      "mass_flux",
      "The postprocessor name holding the TOTAL mass flux over the whole line in kg/s (positive is "
      "flux in, negative is flux out).  This total is divided between the completions, so it is "
      "the same whether the line is resolved by one Dirac point or by many.");
  params.addClassDescription(
      "Line source (or sink) that adds (or removes) fluid at a total mass flux rate specified by a "
      "postprocessor.  The total is divided between the completions in proportion to the length of "
      "line each represents, so the rate is prescribed exactly and does not depend on the "
      "conditions in the reservoir.  The weight column of the point file is not used by this "
      "object: the flux is distributed by length alone.");
  return params;
}

PorousFlowLineSourceFromPostprocessor::PorousFlowLineSourceFromPostprocessor(
    const InputParameters & parameters)
  : PorousFlowLineGeometry(parameters), _mass_flux(getPostprocessorValue("mass_flux"))
{
}

Real
PorousFlowLineSourceFromPostprocessor::computeQpResidual()
{
  // Negative sign to make a positive mass_flux in the input file a source
  return -_test[_i][_qp] * _mass_flux * completionFraction(currentPointCachedID());
}
