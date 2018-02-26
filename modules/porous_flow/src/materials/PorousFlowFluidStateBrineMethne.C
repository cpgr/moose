//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PorousFlowFluidStateBrineMethane.h"
#include "PorousFlowCapillaryPressure.h"
#include "PorousFlowBrineMethane.h"

template <>
InputParameters
validParams<PorousFlowFluidStateBrineMethane>()
{
  InputParameters params = validParams<PorousFlowFluidStateFlashBase>();
  params.addCoupledVar("xnacl", 0, "The salt mass fraction in the brine (kg/kg)");
  params.addClassDescription("Fluid state class for brine and methane");
  return params;
}

PorousFlowFluidStateBrineMethane::PorousFlowFluidStateBrineMethane(
    const InputParameters & parameters)
  : PorousFlowFluidStateFlashBase(parameters),
    _xnacl(_nodal_material ? coupledNodalValue("xnacl") : coupledValue("xnacl")),
    _fs_uo(getUserObject<PorousFlowBrineMethane>("fluid_state"))
{
  // Check that a valid Brine-Methane FluidState has been supplied in fluid_state
  if (_fs_uo.fluidStateName() != "brine-methane")
    mooseError("Only a valid Brine-Methane FluidState can be used in ", _name);
}

void
PorousFlowFluidStateBrineMethane::thermophysicalProperties()
{
  // The FluidProperty objects use temperature in K
  Real Tk = _temperature[_qp] + _T_c2k;

  _fs_uo.thermophysicalProperties(_gas_porepressure[_qp], Tk, _xnacl[_qp], (*_z[0])[_qp], _fsp);
}
