//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PorousFlowConstrainedPeacemanBorehole.h"
#include "Function.h"

registerMooseObject("PorousFlowApp", PorousFlowConstrainedPeacemanBorehole);

InputParameters
PorousFlowConstrainedPeacemanBorehole::validParams()
{
  InputParameters params = PorousFlowPeacemanBorehole::validParams();
  // The bottom-hole pressure is solved for, so neither the function that would prescribe it nor
  // the gate that would set the direction of flow applies here
  params.suppressParameter<FunctionName>("bottom_p_or_t");
  params.suppressParameter<FunctionName>("character");
  // The temperature-dependent wellbore profile is integrated up from bottom_p_or_t, which this
  // borehole does not have, so only the constant 'unit_weight' is available here
  params.suppressParameter<UserObjectName>("unit_weight_fp");
  params.renameCoupledVar("scalar_variable",
                          "bottom_hole_pressure",
                          "The scalar variable holding the pressure at the bottom of this "
                          "borehole.  Its value is solved for, not prescribed: it is whatever "
                          "delivers 'rate'.");
  params.addRequiredParam<FunctionName>(
      "rate",
      "The rate this well is asked to deliver, summed over all of its completions, positive for "
      "production and negative for injection.  The units are those of the equation the borehole "
      "contributes to, so kg/s for a fluid mass equation and J/s for a heat equation.");
  params.addParam<FunctionName>(
      "bhp_max",
      "The bottom-hole pressure this well may not exceed.  Where delivering 'rate' would need more "
      "than this, the well holds this pressure instead and delivers whatever rate follows, which "
      "is "
      "the usual treatment of a well that cannot meet its target without over-pressuring the "
      "formation.  Without this parameter the rate is delivered whatever pressure it takes.");
  params.addParam<unsigned int>(
      "max_constraint_swaps",
      5,
      "The number of times in one timestep the well may switch between delivering 'rate' and "
      "holding 'bhp_max'.  A well whose operating point sits on the limit can otherwise switch "
      "back "
      "and forth without ever converging; once this many switches have been made the well keeps "
      "whichever constraint it holds until the next timestep.");
  params.addClassDescription(
      "A Peaceman borehole whose bottom-hole pressure is a scalar unknown, constrained so that the "
      "flux summed over the completions equals a requested rate.  The rate is therefore delivered "
      "exactly at every timestep, while the flux at each completion still follows the local "
      "Peaceman expression, so the allocation along the completions responds to permeability and "
      "to how far each layer has pressured up.");
  return params;
}

PorousFlowConstrainedPeacemanBorehole::PorousFlowConstrainedPeacemanBorehole(
    const InputParameters & parameters)
  : PorousFlowPeacemanBorehole(parameters),
    _rate(getFunction("rate")),
    _bhp_max(isParamValid("bhp_max") ? &getFunction("bhp_max") : nullptr),
    _max_swaps(getParam<unsigned int>("max_constraint_swaps")),
    _clamped(false),
    _n_swaps(0),
    _rate_mode_residual(0.0),
    _pbh_last(0.0)
{
  if (!_use_scalar)
    paramError("bottom_hole_pressure",
               "this borehole solves for its bottom-hole pressure, so it needs a scalar variable "
               "to hold it.  Supply 'bottom_hole_pressure', or use a PorousFlowPeacemanBorehole "
               "instead if the bottom-hole pressure is known in advance.");

  if (_k_order != 1)
    paramError("bottom_hole_pressure",
               "a borehole has one bottom-hole pressure, so this must be a scalar variable "
               "carrying a single degree of freedom, but the one supplied carries ",
               _k_order);
}

void
PorousFlowConstrainedPeacemanBorehole::timestepSetup()
{
  PorousFlowPeacemanBorehole::timestepSetup();
  _n_swaps = 0;

  // The constraint is decided here, once per timestep, so that the equation being solved never
  // changes part way through a solve.  Deciding it between Newton iterations instead was measured
  // to cost a failed timestep wherever the limit engages, while never engaging the limit within the
  // step that first breaches it, so it bought nothing for that cost
  if (_bhp_max)
  {
    _communicator.sum(_rate_mode_residual);
    updateConstraint();
  }
}

void
PorousFlowConstrainedPeacemanBorehole::residualSetup()
{
  PorousFlowPeacemanBorehole::residualSetup();

  // The decision is made in timestepSetup, so this only resets the accumulator: the value left
  // behind by the last evaluation of a timestep is the one that decides the next
  _rate_mode_residual = 0.0;
}

void
PorousFlowConstrainedPeacemanBorehole::updateConstraint()
{
  // Once the well has switched often enough in one timestep, stop: an operating point sitting on
  // the limit can otherwise chatter between the two constraints and never converge
  if (_n_swaps >= _max_swaps)
    return;

  const Real bhp_max = _bhp_max->value(_t, _bottom_point);

  if (!_clamped && _pbh_last > bhp_max)
  {
    // Delivering the rate needs a pressure the well is not allowed to reach
    _clamped = true;
    ++_n_swaps;
  }
  else if (_clamped && _rate_mode_residual < 0.0)
  {
    // The flux falls as the bottom-hole pressure rises, so a negative rate-constraint residual at
    // the limit means the rate would be met below it, and the limit no longer binds
    _clamped = false;
    ++_n_swaps;
  }
}

Real
PorousFlowConstrainedPeacemanBorehole::boreholePressure() const
{
  return _kappa[0] + _unit_weight * (_q_point[_qp] - _bottom_point);
}

Real
PorousFlowConstrainedPeacemanBorehole::conductance(unsigned current_dirac_ptid) const
{
  return wellConstantSum(current_dirac_ptid) * _test[_i][_qp];
}

Real
PorousFlowConstrainedPeacemanBorehole::computeQpBaseOutflow(unsigned current_dirac_ptid) const
{
  return conductance(current_dirac_ptid) * (ptqp() - boreholePressure());
}

void
PorousFlowConstrainedPeacemanBorehole::computeQpBaseOutflowJacobian(unsigned jvar,
                                                                    unsigned current_dirac_ptid,
                                                                    Real & outflow,
                                                                    Real & outflowp) const
{
  outflow = 0.0;
  outflowp = 0.0;

  if (_dictator.notPorousFlowVariable(jvar))
    return;
  const unsigned pvar = _dictator.porousFlowVariableNum(jvar);

  const Real c = conductance(current_dirac_ptid);
  outflow = c * (ptqp() - boreholePressure());
  outflowp = c * dptqp(pvar) * _phi[_j][_qp];
}

Real
PorousFlowConstrainedPeacemanBorehole::computeQpScalarResidual()
{
  // The flux of the whole completion is the sum over the element's nodes, since the residual of
  // each node carries only its own share via the test function
  const unsigned current_dirac_ptid = currentPointCachedID();
  const Real fraction = completionFraction(current_dirac_ptid);

  const auto i_save = _i;
  Real flux = 0.0;
  for (_i = 0; _i < _test.size(); ++_i)
    flux += weightedOutflow();
  _i = i_save;

  // The rate is a property of the well, not of any one completion, so it must enter the constraint
  // exactly once however the completions are spread over elements and ranks.  Subtracting each
  // completion's share of it sums to the rate by construction, needs no communication, and does
  // not depend on how many ranks the well is divided between
  const Real rate_mode = flux - _rate.value(_t, _q_point[_qp]) * fraction;

  // Recorded even while the limit is being held, because it is what says when the limit stops
  // binding.  residualSetup() sums it across ranks and acts on it between Newton iterations
  _rate_mode_residual += rate_mode;
  _pbh_last = _kappa[0];

  if (_clamped)
    // The limit has taken over, so the well no longer delivers the rate and its bottom-hole
    // pressure is pinned.  Sharing this out between the completions the same way the rate is
    // shared keeps it entering the constraint exactly once
    return (_kappa[0] - _bhp_max->value(_t, _bottom_point)) * fraction;

  return rate_mode;
}

Real
PorousFlowConstrainedPeacemanBorehole::computeQpJacobianScalar()
{
  // outflow = weighting * conductance * (porepressure - boreholePressure()), and neither the
  // weighting nor the conductance depends on the borehole pressure
  return -weightingFactor() * conductance(currentPointCachedID());
}

Real
PorousFlowConstrainedPeacemanBorehole::computeQpScalarJacobian()
{
  if (_clamped)
    // d(P_bh - bhp_max)/d(P_bh), shared between the completions like the constraint itself
    return completionFraction(currentPointCachedID());

  const auto i_save = _i;
  Real jacobian = 0.0;
  for (_i = 0; _i < _test.size(); ++_i)
    jacobian += computeQpJacobianScalar();
  _i = i_save;

  return jacobian;
}

Real
PorousFlowConstrainedPeacemanBorehole::computeQpScalarOffDiagJacobian(const unsigned int jvar_num)
{
  if (_clamped)
    // A pinned bottom-hole pressure does not depend on the state of the reservoir
    return 0.0;

  // The constraint is the flux summed over the element's nodes, so its derivative is the derivative
  // of each node's residual summed the same way.  Reusing jac() keeps the constraint and the flux
  // on one code path: the product rule through the use_* weightings is applied in exactly one place
  const auto i_save = _i;
  Real jacobian = 0.0;
  for (_i = 0; _i < _test.size(); ++_i)
    jacobian += jac(jvar_num);
  _i = i_save;

  return jacobian;
}
