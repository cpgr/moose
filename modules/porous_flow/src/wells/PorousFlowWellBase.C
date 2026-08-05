//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PorousFlowWellBase.h"

InputParameters
PorousFlowWellBase::validParamsCommon()
{
  InputParameters params = emptyInputParameters();
  params.addRequiredParam<UserObjectName>(
      "PorousFlowDictator", "The UserObject that holds the list of PorousFlow variable names");
  params.addRequiredParam<VariableName>(
      "porepressure",
      "The variable that holds the mass balance of the final fluid component.  This is the same "
      "porepressure variable that the PorousFlow Physics/Action uses.");
  params.addParam<std::vector<VariableName>>(
      "mass_fraction_vars",
      {},
      "List of variables holding the mass balances of the fluid components.  The i^th variable "
      "holds the equation for the i^th fluid component, and 'porepressure' holds the equation for "
      "the final fluid component, exactly as in the PorousFlow Physics/Action.");
  // Completion geometry, passed straight through to PorousFlowLineGeometry.  Every well is a
  // trajectory through the formation, whatever drives the flux at its completions, so these belong
  // to all well types alike
  params.addParam<std::string>("point_file",
                               "File containing the completion coordinates and their weights, one "
                               "'r x y z' per line, as for PorousFlowPolyLineSink");
  params.addParam<std::vector<Real>>("line_base",
                                     "Completion base point, in the same 'r x y z' format as a "
                                     "single-line point_file.  One completion is then generated "
                                     "for each element the completion interval passes through.");
  params.addRangeCheckedParam<Real>(
      "line_length",
      0.0,
      "line_length>=0",
      "Length of the completion interval.  Only used if there is a single completion point.");
  params.addParam<RealVectorValue>(
      "line_direction",
      RealVectorValue(0.0, 0.0, 1.0),
      "Direction of the completion interval.  Only used if there is a single completion point.");
  params.addParam<Real>("start_time", "The time at which this well starts operating");
  params.addParam<Real>("end_time", "The time at which this well stops operating");
  params.addParam<std::vector<Real>>(
      "sync_times",
      {},
      "Times at which the schedule changes.  These are added to the simulation's sync times so "
      "that a rate change never lands part-way through a timestep.");
  params.addParam<FunctionName>(
      "shut_if",
      "A Function that shuts this well in whenever it evaluates to a nonzero (true) value.  Use a "
      "ParsedFunction to test a Postprocessor, whose name may be given directly in 'symbol_values' "
      "(e.g. expression = 'wc > 0.95', symbol_names = 'wc', symbol_values = 'my_water_cut_pp').  "
      "The Postprocessor is read from the end of the previous timestep, so the well shuts one "
      "timestep after the condition is first met.");
  return params;
}

InputParameters
PorousFlowWellBase::validParamsNotYetImplemented()
{
  InputParameters params = emptyInputParameters();
  params.addParam<Real>("bhp", "Bottom-hole pressure target");
  params.addParam<Real>("bhp_max", "Maximum bottom-hole pressure (an injector constraint)");
  params.addParam<Real>("bhp_min", "Minimum bottom-hole pressure (a producer constraint)");
  params.addParam<std::string>("group", "The well group this well belongs to");
  params.addParam<std::string>("reopen_if", "Condition upon which this well is reopened");
  return params;
}

InputParameters
PorousFlowWellBase::validParams()
{
  InputParameters params = MooseObject::validParams();
  params += PorousFlowWellBase::validParamsCommon();
  params += PorousFlowWellBase::validParamsNotYetImplemented();
  // Must match the system name given to registerMooseObjectTask in PorousFlowApp.C
  params.registerBase("PorousFlowWellBase");
  return params;
}

PorousFlowWellBase::PorousFlowWellBase(const InputParameters & parameters) : MooseObject(parameters)
{
}
