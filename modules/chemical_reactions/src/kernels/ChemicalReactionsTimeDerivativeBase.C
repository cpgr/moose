//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ChemicalReactionsTimeDerivativeBase.h"

template <>
InputParameters
validParams<ChemicalReactionsTimeDerivativeBase>()
{
  InputParameters params = validParams<TimeDerivative>();
  params.addParam<Real>(
      "weight", 1.0, "The weight of the equilibrium species in total concentration");
  MooseEnum concentration_enum("molarity molality default", "default");
  params.addParam<MooseEnum>("concentration_units",
                             concentration_enum,
                             "Units for concentration. Default is moles/volume");
  params.addClassDescription("Base class for chemical reactions time derivative kernels");
  return params;
}

ChemicalReactionsTimeDerivativeBase::ChemicalReactionsTimeDerivativeBase(
    const InputParameters & parameters)
  : TimeDerivative(parameters),
    _weight(getParam<Real>("weight")),
    _porosity(getMaterialProperty<Real>("porosity")),
    _porosity_old(getMaterialPropertyOld<Real>("porosity")),
    _u_old(valueOld()),
    _concentration_enum(getParam<MooseEnum>("concentration_units").getEnum<ConcentrationEnum>())
{
  // If concentration is in units of molality (moles / kg H2O), then get density
  if (_concentration_enum == ConcentrationEnum::MOLALITY)
    _density = &getMaterialProperty<Real>("density");
}
