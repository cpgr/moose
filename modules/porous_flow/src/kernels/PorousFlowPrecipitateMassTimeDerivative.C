//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PorousFlowPrecipitateMassTimeDerivative.h"

#include "MooseVariable.h"

registerMooseObject("PorousFlowApp", PorousFlowPrecipitateMassTimeDerivative);

InputParameters
PorousFlowPrecipitateMassTimeDerivative::validParams()
{
  InputParameters params = TimeKernel::validParams();
  params.addRequiredParam<UserObjectName>(
      "PorousFlowDictator", "The UserObject that holds the list of PorousFlow variable names.");
  params.addRangeCheckedParam<Real>(
      "halite_density",
      2165.0,
      "halite_density > 0",
      "Density of solid halite (kg/m^3). Must match the value used by the "
      "PorousFlowHaliteVolumeFraction material.");
  params.set<bool>("use_displaced_mesh") = false;
  params.suppressParameter<bool>("use_displaced_mesh");
  params.addClassDescription("Derivative of the precipitated (solid) halite mass with respect to "
                             "time, applied to the salt-component equation. Mass lumping to the "
                             "nodes is used.");
  return params;
}

PorousFlowPrecipitateMassTimeDerivative::PorousFlowPrecipitateMassTimeDerivative(
    const InputParameters & parameters)
  : TimeKernel(parameters),
    _dictator(getUserObject<PorousFlowDictator>("PorousFlowDictator")),
    _var_is_porflow_var(_dictator.isPorousFlowVariable(_var.number())),
    _halite_density(getParam<Real>("halite_density")),
    _halite_volume_fraction(getMaterialProperty<Real>("PorousFlow_halite_volume_fraction_nodal")),
    _halite_volume_fraction_old(
        getMaterialPropertyOld<Real>("PorousFlow_halite_volume_fraction_nodal")),
    _dhalite_volume_fraction_dvar(
        getMaterialProperty<std::vector<Real>>("dPorousFlow_halite_volume_fraction_nodal_dvar"))
{
}

Real
PorousFlowPrecipitateMassTimeDerivative::computeQpResidual()
{
  return _test[_i][_qp] * _halite_density *
         (_halite_volume_fraction[_i] - _halite_volume_fraction_old[_i]) / _dt;
}

Real
PorousFlowPrecipitateMassTimeDerivative::computeQpJacobian()
{
  if (!_var_is_porflow_var)
    return 0.0;
  return computeQpJac(_dictator.porousFlowVariableNum(_var.number()));
}

Real
PorousFlowPrecipitateMassTimeDerivative::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (_dictator.notPorousFlowVariable(jvar))
    return 0.0;
  return computeQpJac(_dictator.porousFlowVariableNum(jvar));
}

Real
PorousFlowPrecipitateMassTimeDerivative::computeQpJac(unsigned int pvar) const
{
  // The halite volume fraction is lumped to the nodes and depends only on the (nodal) PorousFlow
  // variables, so the only non-zero contribution is for _i == _j.
  if (_i != _j)
    return 0.0;
  return _test[_i][_qp] * _halite_density * _dhalite_volume_fraction_dvar[_i][pvar] / _dt;
}
