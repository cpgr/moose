//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef PRIMARYTIMEDERIVATIVE_H
#define PRIMARYTIMEDERIVATIVE_H

#include "TimeDerivative.h"

class PrimaryTimeDerivative;

template <>
InputParameters validParams<PrimaryTimeDerivative>();

class PrimaryTimeDerivative : public TimeDerivative
{
public:
  PrimaryTimeDerivative(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

  /// Porosity
  const MaterialProperty<Real> & _porosity;
};

#endif // PRIMARYTIMEDERIVATIVE_H
