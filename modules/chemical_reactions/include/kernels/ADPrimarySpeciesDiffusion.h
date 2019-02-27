//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#ifndef ADPRIMARYSPECIESDIFFUSION_H
#define ADPRIMARYSPECIESDIFFUSION_H

#include "ADDiffusion.h"

template <ComputeStage compute_stage>
class ADPrimarySpeciesDiffusion;

declareADValidParams(ADPrimarySpeciesDiffusion);

template <ComputeStage compute_stage>
class ADPrimarySpeciesDiffusion : public ADDiffusion<compute_stage>
{
public:
  ADPrimarySpeciesDiffusion(const InputParameters & parameters);

protected:
  virtual ADResidual computeQpResidual() override;

  /// Diffusion coefficient
  const MaterialProperty<Real> & _diffusivity;

  usingKernelMembers;
};

#endif // ADPRIMARYSPECIESDIFFUSION_H
