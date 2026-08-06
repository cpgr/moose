//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ConstrainedPointSource.h"

registerMooseObject("PorousFlowTestApp", ConstrainedPointSource);

InputParameters
ConstrainedPointSource::validParams()
{
  InputParameters params = PorousFlowDiracKernelScalarBase::validParams();
  params.addRequiredParam<Real>("flux_coefficient",
                                "The strength of the sink per unit (u - scalar_variable)");
  params.addRequiredParam<Real>("rate", "The rate that the constraint equation demands");
  params.addRequiredParam<Point>("point", "The location of the point sink");
  params.addClassDescription("A point sink whose strength is set by a scalar variable, together "
                             "with the constraint that fixes that strength to a prescribed rate");
  return params;
}

ConstrainedPointSource::ConstrainedPointSource(const InputParameters & parameters)
  : PorousFlowDiracKernelScalarBase(parameters),
    _flux_coefficient(getParam<Real>("flux_coefficient")),
    _rate(getParam<Real>("rate")),
    _point(getParam<Point>("point"))
{
  if (_use_scalar && _k_order != 1)
    mooseError("ConstrainedPointSource requires a scalar_variable carrying a single degree of "
               "freedom, but the supplied one carries ",
               _k_order);
}

void
ConstrainedPointSource::addPoints()
{
  addPoint(_point, 0);
}

Real
ConstrainedPointSource::computeQpResidual()
{
  return _flux_coefficient * (_u[_qp] - _kappa[0]) * _test[_i][_qp];
}

Real
ConstrainedPointSource::computeQpJacobian()
{
  return _flux_coefficient * _phi[_j][_qp] * _test[_i][_qp];
}

Real
ConstrainedPointSource::computeQpJacobianScalar()
{
  return -_flux_coefficient * _test[_i][_qp];
}

Real
ConstrainedPointSource::computeQpScalarResidual()
{
  return _flux_coefficient * (_u[_qp] - _kappa[0]) - _rate;
}

Real
ConstrainedPointSource::computeQpScalarJacobian()
{
  return -_flux_coefficient;
}

Real
ConstrainedPointSource::computeQpScalarOffDiagJacobian(const unsigned int jvar_num)
{
  if (jvar_num != _var.number())
    return 0.0;

  return _flux_coefficient * _phi[_j][_qp];
}
