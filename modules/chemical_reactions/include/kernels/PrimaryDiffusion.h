//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef PRIMARYDIFFUSION_H
#define PRIMARYDIFFUSION_H

#include "Diffusion.h"

class PrimaryDiffusion;

template <>
InputParameters validParams<PrimaryDiffusion>();

class PrimaryDiffusion : public Diffusion
{
public:
  PrimaryDiffusion(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

  /// Diffusion coefficient
  const MaterialProperty<Real> & _diffusivity;
};

#endif // PRIMARYDIFFUSION_H
