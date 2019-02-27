//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#ifndef ADPRIMARYSPECIESTIMEDERIVATIVE_H
#define ADPRIMARYSPECIESTIMEDERIVATIVE_H

#include "ADTimeDerivative.h"

template <ComputeStage compute_stage>
class ADPrimarySpeciesTimeDerivative;

declareADValidParams(ADPrimarySpeciesTimeDerivative);

template <ComputeStage compute_stage>
class ADPrimarySpeciesTimeDerivative : public ADTimeDerivative<compute_stage>
{
public:
  ADPrimarySpeciesTimeDerivative(const InputParameters & parameters);

protected:
  virtual ADResidual computeQpResidual() override;

  /// Porosity
  const MaterialProperty<Real> & _porosity;

  usingTimeKernelMembers;
};

#endif // ADPRIMARYSPECIESTIMEDERIVATIVE_H
