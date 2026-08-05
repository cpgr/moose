//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "PorousFlowLineGeometry.h"

class SinglePhaseFluidProperties;

/**
 * Line source that adds the heat energy carried by fluid injected at a total mass flux rate
 * specified by a postprocessor, divided between the completions in proportion to their length
 */
class PorousFlowLineEnthalpySourceFromPostprocessor : public PorousFlowLineGeometry
{
public:
  static InputParameters validParams();

  PorousFlowLineEnthalpySourceFromPostprocessor(const InputParameters & parameters);

  virtual Real computeQpResidual() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

protected:
  /// The total mass flux of injected fluid over the whole line (kg/s)
  const PostprocessorValue & _mass_flux;

  /// Pressure at which the enthalpy of the injected fluid is evaluated
  const VariableValue & _pressure;

  /// Temperature of the injected fluid (K)
  const PostprocessorValue & _T_in;

  /// Fluid properties of the injected fluid
  const SinglePhaseFluidProperties & _fp;

  /// Variable number of the pressure variable
  const unsigned int _p_var_num;
};
