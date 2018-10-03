//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PorousFlowFluidStateWaterNCG.h"
#include "PorousFlowCapillaryPressure.h"
#include "PorousFlowWaterNCG.h"

registerMooseObject("PorousFlowApp", PorousFlowFluidStateWaterNCG);

template <>
InputParameters
validParams<PorousFlowFluidStateWaterNCG>()
{
  InputParameters params = validParams<PorousFlowFluidState>();
  params.addClassDescription("Fluid state class for water and non-condensable gas");
  return params;
}

PorousFlowFluidStateWaterNCG::PorousFlowFluidStateWaterNCG(const InputParameters & parameters)
  : PorousFlowFluidState(parameters)
{
}

void
PorousFlowFluidStateWaterNCG::thermophysicalProperties()
{
  // The FluidProperties objects use temperature in K
  Real Tk = _temperature[_qp] + _T_c2k;

  _fs.thermophysicalProperties(_gas_porepressure[_qp], Tk, _Xnacl[_qp], _Z, _qp, _fsp);
}
