//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PorousFlowFluidStateBrineCO2.h"
#include "PorousFlowCapillaryPressure.h"
#include "PorousFlowBrineCO2.h"

registerMooseObject("PorousFlowApp", PorousFlowFluidStateBrineCO2);

template <>
InputParameters
validParams<PorousFlowFluidStateBrineCO2>()
{
  InputParameters params = validParams<PorousFlowFluidState>();
  params.addClassDescription("Fluid state class for brine and CO2");
  return params;
}

PorousFlowFluidStateBrineCO2::PorousFlowFluidStateBrineCO2(const InputParameters & parameters)
  : PorousFlowFluidState(parameters)
{
}

void
PorousFlowFluidStateBrineCO2::thermophysicalProperties()
{
  // The FluidProperties objects use temperature in K
  Real Tk = _temperature[_qp] + _T_c2k;

  _fs.thermophysicalProperties(_gas_porepressure[_qp], Tk, _Xnacl[_qp], _Z, _qp, _fsp);
}
