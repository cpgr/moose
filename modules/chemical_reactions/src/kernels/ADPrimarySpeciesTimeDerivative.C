//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADPrimarySpeciesTimeDerivative.h"

registerADMooseObject("ChemicalReactionsApp", ADPrimarySpeciesTimeDerivative);

defineADValidParams(
    ADPrimarySpeciesTimeDerivative,
    ADTimeDerivative,
    params.addClassDescription("Derivative of primary species concentration wrt time"););

template <ComputeStage compute_stage>
ADPrimarySpeciesTimeDerivative<compute_stage>::ADPrimarySpeciesTimeDerivative(
    const InputParameters & parameters)
  : ADTimeDerivative<compute_stage>(parameters), _porosity(adGetMaterialProperty<Real>("porosity"))
{
}

template <ComputeStage compute_stage>
ADResidual
ADPrimarySpeciesTimeDerivative<compute_stage>::precomputeQpResidual()
{
  return _porosity[_qp] * ADTimeDerivative<compute_stage>::computeQpResidual();
}
