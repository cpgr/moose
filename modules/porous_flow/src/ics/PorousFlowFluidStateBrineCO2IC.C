//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PorousFlowFluidStateBrineCO2IC.h"
#include "PorousFlowBrineCO2.h"

registerMooseObject("PorousFlowApp", PorousFlowFluidStateBrineCO2IC);

template <>
InputParameters
validParams<PorousFlowFluidStateBrineCO2IC>()
{
  InputParameters params = validParams<PorousFlowFluidStateIC>();
  params.addRequiredParam<UserObjectName>("fluid_state", "Name of the FluidState UserObject");
  params.addClassDescription(
      "An initial condition to calculate z from saturation for brine and CO2");
  return params;
}

PorousFlowFluidStateBrineCO2IC::PorousFlowFluidStateBrineCO2IC(const InputParameters & parameters)
  : PorousFlowFluidStateIC(parameters)
{
}

Real
PorousFlowFluidStateBrineCO2IC::value(const Point & /*p*/)
{
  // The brine-co2 fluid state needs temperature in K
  Real Tk = _temperature[_qp] + _T_c2k;

  // The total mass fraction corresponding to the input saturation
  Real Z = _fs.totalMassFraction(_gas_porepressure[_qp], Tk, _xnacl[_qp], _saturation[_qp], _qp);

  return Z;
}
