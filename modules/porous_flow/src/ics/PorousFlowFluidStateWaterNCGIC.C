//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PorousFlowFluidStateWaterNCGIC.h"
#include "PorousFlowWaterNCG.h"

registerMooseObject("PorousFlowApp", PorousFlowFluidStateWaterNCGIC);

template <>
InputParameters
validParams<PorousFlowFluidStateWaterNCGIC>()
{
  InputParameters params = validParams<PorousFlowFluidStateIC>();
  params.addRequiredParam<UserObjectName>("fluid_state", "Name of the FluidState UserObject");
  params.addClassDescription(
      "An initial condition to calculate z from saturation for water and non-condensable gas");
  return params;
}

PorousFlowFluidStateWaterNCGIC::PorousFlowFluidStateWaterNCGIC(const InputParameters & parameters)
  : PorousFlowFluidStateIC(parameters)
{
}

Real
PorousFlowFluidStateWaterNCGIC::value(const Point & /*p*/)
{
  // The water-ncg fluid state needs temperature in K
  Real Tk = _temperature[_qp] + _T_c2k;

  // The total mass fraction corresponding to the input saturation
  Real Z = _fs.totalMassFraction(_gas_porepressure[_qp], Tk, _xnacl[_qp], _saturation[_qp], _qp);

  return Z;
}
