//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef POROUSFLOWFLUIDSTATEBRINEMETHANE_H
#define POROUSFLOWFLUIDSTATEBRINEMETHANE_H

#include "PorousFlowFluidStateFlashBase.h"

class PorousFlowBrineMethane;
class PorousFlowFluidStateBrineMethane;

template <>
InputParameters validParams<PorousFlowFluidStateBrineMethane>();

/**
 * Fluid state class for brine and Methane. from
 * Duan and Mao, A thermodynamic model for calculating methane solubility, density
 * and gas phase composition of methane-bearing aqueous fluids from 273
 * to 523 K and from 1 to 2000 bar, Geochimica et Cosmochimica Acta, 70, 3369-3386 (2006)
 */
class PorousFlowFluidStateBrineMethane : public PorousFlowFluidStateFlashBase
{
public:
  PorousFlowFluidStateBrineMethane(const InputParameters & parameters);

protected:
  virtual void thermophysicalProperties() override;

  /// Salt mass fraction (kg/kg)
  const VariableValue & _xnacl;
  /// FluidState UserObject
  const PorousFlowBrineMethane & _fs_uo;
};

#endif // POROUSFLOWFLUIDSTATEBRINEMETHANE_H
