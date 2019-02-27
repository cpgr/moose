//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADDarcyFluxPressure.h"

registerADMooseObject("ChemicalReactionsApp", ADDarcyFluxPressure);

defineADValidParams(ADDarcyFluxPressure, ADKernel, params.addRequiredCoupledVar("p", "Pressure");
                    RealVectorValue g(0, 0, 0);
                    params.addParam<RealVectorValue>("gravity",
                                                     g,
                                                     "Gravity vector (default is (0, 0, 0))");
                    params.addClassDescription("Darcy flux"););

template <ComputeStage compute_stage>
ADDarcyFluxPressure<compute_stage>::ADDarcyFluxPressure(const InputParameters & parameters)
  : ADKernel<compute_stage>(parameters),
    _cond(adGetMaterialProperty<Real>("conductivity")),
    _gravity(adGetParam<RealVectorValue>("gravity")),
    _density(adGetMaterialProperty<Real>("density"))
{
}

template <ComputeStage compute_stage>
ADResidual
ADDarcyFluxPressure<compute_stage>::computeQpResidual()
{
  return _grad_test[_i][_qp] * _cond[_qp] * (_grad_u[_qp] - _density[_qp] * _gravity);
}
