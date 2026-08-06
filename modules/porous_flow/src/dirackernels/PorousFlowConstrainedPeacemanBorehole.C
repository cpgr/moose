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
  : PorousFlowPeacemanBorehole(parameters), _rate(getFunction("rate"))
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
  const auto i_save = _i;
  Real flux = 0.0;
  for (_i = 0; _i < _test.size(); ++_i)
    flux += weightedOutflow();
  _i = i_save;

  // The rate is a property of the well, not of any one completion, so it must enter the constraint
  // exactly once however the completions are spread over elements and ranks.  Subtracting each
  // completion's share of it sums to the rate by construction, needs no communication, and does
  // not depend on how many ranks the well is divided between
  return flux - _rate.value(_t, _q_point[_qp]) * completionFraction(current_dirac_ptid);
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
