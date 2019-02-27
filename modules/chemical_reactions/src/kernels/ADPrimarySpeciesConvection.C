//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADPrimarySpeciesConvection.h"

registerADMooseObject("ChemicalReactionsApp", ADPrimarySpeciesConvection);

defineADValidParams(
    ADPrimarySpeciesConvection, ADKernel, params.addRequiredCoupledVar("p", "Pressure");
    RealVectorValue g(0, 0, 0);
    params.addParam<RealVectorValue>("gravity", g, "Gravity vector (default is (0, 0, 0))");
    params.addClassDescription("Convection of primary species"););

template <ComputeStage compute_stage>
ADPrimarySpeciesConvection<compute_stage>::ADPrimarySpeciesConvection(
    const InputParameters & parameters)
  : ADKernel<compute_stage>(parameters),
    _cond(adGetMaterialProperty<Real>("conductivity")),
    _gravity(adGetParam<RealVectorValue>("gravity")),
    _density(adGetMaterialProperty<Real>("density")),
    _grad_p(adCoupledGradient("p"))
{
}

template <ComputeStage compute_stage>
ADResidual
ADPrimarySpeciesConvection<compute_stage>::computeQpResidual()
{
  const auto darcy_vel = -_cond[_qp] * (_grad_p[_qp] - _density[_qp] * _gravity);

  return _test[_i][_qp] * (darcy_vel * _grad_u[_qp]);
}
