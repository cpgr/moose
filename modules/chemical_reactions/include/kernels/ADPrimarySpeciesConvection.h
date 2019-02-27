//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef ADPRIMARYSPECIESCONVECTION_H
#define ADPRIMARYSPECIESCONVECTION_H

#include "ADKernel.h"

template <ComputeStage compute_stage>
class ADPrimarySpeciesConvection;

declareADValidParams(ADPrimarySpeciesConvection);

template <ComputeStage compute_stage>
class ADPrimarySpeciesConvection : public ADKernel<compute_stage>
{
public:
  ADPrimarySpeciesConvection(const InputParameters & parameters);

protected:
  virtual ADResidual computeQpResidual() override;

  /// Hydraulic conductivity
  const MaterialProperty<Real> & _cond;
  /// Gravity
  const RealVectorValue _gravity;
  /// Fluid density
  const MaterialProperty<Real> & _density;
  /// Pressure gradient
  const ADVariableGradient & _grad_p;

  usingKernelMembers;
};

#endif // ADPRIMARYSPECIESCONVECTION_H
