//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef PRIMARYCONVECTION_H
#define PRIMARYCONVECTION_H

#include "Kernel.h"
#include "DerivativeMaterialInterface.h"

class PrimaryConvection;

template <>
InputParameters validParams<PrimaryConvection>();

class PrimaryConvection : public DerivativeMaterialInterface<Kernel>
{
public:
  PrimaryConvection(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  /// Hydraulic conductivity
  const MaterialProperty<Real> & _cond;
  /// Gravity
  const RealVectorValue _gravity;
  /// Fluid density
  const MaterialProperty<Real> & _density;
  /// Pressure gradient
  const VariableGradient & _grad_p;
  /// Pressure variable number
  const unsigned int _pvar;
};

#endif // PRIMARYCONVECTION_H
