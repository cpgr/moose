//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PorousFlowProductionWell.h"

registerMooseObject("PorousFlowApp", PorousFlowProductionWell);

InputParameters
PorousFlowProductionWell::validParams()
{
  InputParameters params = PorousFlowWellBase::validParams();

  // Geometry, passed straight through to PorousFlowLineGeometry
  params.addParam<std::string>("point_file",
                               "File containing the completion coordinates and their weights, one "
                               "'r x y z' per line, as for PorousFlowPolyLineSink");
  params.addParam<std::vector<Real>>("line_base",
                                     "Completion base point, in the same 'r x y z' format as a "
                                     "single-line point_file");
  params.addRangeCheckedParam<Real>(
      "line_length",
      0.0,
      "line_length>=0",
      "Length of the completion.  Only used if there is a single completion point.");
  params.addParam<RealVectorValue>(
      "line_direction",
      RealVectorValue(0.0, 0.0, 1.0),
      "Direction of the completion.  Only used if there is a single completion point.");

  // Deliverability curve, passed straight through to PorousFlowPolyLineSink
  params.addRequiredParam<std::vector<Real>>(
      "p_or_t_vals",
      "Tuple of pressure (or temperature) values.  Must be monotonically increasing.");
  params.addRequiredParam<std::vector<Real>>(
      "fluxes",
      "Tuple of flux values defining the deliverability curve, positive for production.  These are "
      "multiplied by the completion length, and by any of the 'use_*' weightings requested below, "
      "to give a rate in kg/s.");

  MooseEnum p_or_t_choice("pressure=0 temperature=1", "pressure");
  params.addParam<MooseEnum>(
      "function_of",
      p_or_t_choice,
      "Whether the deliverability curve is a function of pressure or of temperature");
  params.addParam<unsigned int>(
      "fluid_phase", 0, "The fluid phase whose pressure drives flow into this well");
  params.addParam<bool>(
      "use_relative_permeability", false, "Multiply the flux by the fluid relative permeability");
  params.addParam<bool>("use_mobility", false, "Multiply the flux by the fluid mobility");

  params.addClassDescription(
      "A production well whose rate follows a prescribed deliverability curve: the flux is a "
      "piecewise-linear function of the local porepressure (or temperature).  The produced stream "
      "carries the resident composition of the reservoir.");
  return params;
}

PorousFlowProductionWell::PorousFlowProductionWell(const InputParameters & parameters)
  : PorousFlowWellBase(parameters)
{
}
