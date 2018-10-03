//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef POROUSFLOWFLUIDSTATE_H
#define POROUSFLOWFLUIDSTATE_H

#include "PorousFlowFluidStateFlashBase.h"

class PorousFlowFluidState;

template <>
InputParameters validParams<PorousFlowFluidState>();

class PorousFlowFluidState : public PorousFlowFluidStateFlashBase
{
public:
  PorousFlowFluidState(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;
  virtual void thermophysicalProperties() override;

  /// Salt mass fraction (kg/kg)
  const VariableValue & _Xnacl;
  /// Gradient of salt mass fraction (only defined at the qps)
  const VariableGradient & _grad_Xnacl_qp;
  /// Salt mass fraction variable number
  const unsigned int _Xnacl_varnum;
  /// Salt mass fraction PorousFlow variable number
  const unsigned int _Xvar;
  /// Salt component index
  const unsigned int _salt_component;
};

#endif // POROUSFLOWFLUIDSTATE_H
