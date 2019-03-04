//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADPrimarySpeciesDiffusion.h"

registerADMooseObject("ChemicalReactionsApp", ADPrimarySpeciesDiffusion);

defineADValidParams(ADPrimarySpeciesDiffusion,
                    ADDiffusion,
                    params.addClassDescription("Diffusion of primary species"););

template <ComputeStage compute_stage>
ADPrimarySpeciesDiffusion<compute_stage>::ADPrimarySpeciesDiffusion(
    const InputParameters & parameters)
  : ADDiffusion<compute_stage>(parameters), _diffusivity(adGetMaterialProperty<Real>("diffusivity"))
{
}

template <ComputeStage compute_stage>
ADVectorResidual
ADPrimarySpeciesDiffusion<compute_stage>::precomputeQpResidual()
{
  return _diffusivity[_qp] * ADDiffusion<compute_stage>::computeQpResidual();
}
