//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef DARCYFLUXPRESSURE_H
#define DARCYFLUXPRESSURE_H

#include "Kernel.h"
#include "DerivativeMaterialInterface.h"

class DarcyFluxPressure;

template <>
InputParameters validParams<DarcyFluxPressure>();

class DarcyFluxPressure : public DerivativeMaterialInterface<Kernel>
{
public:
  DarcyFluxPressure(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

  /// Hydraulic conductivity
  const MaterialProperty<Real> & _cond;
  /// Gravity
  const RealVectorValue _gravity;
  /// Fluid density
  const MaterialProperty<Real> & _density;
};

#endif // DARCYFLUXPRESSURE_H
