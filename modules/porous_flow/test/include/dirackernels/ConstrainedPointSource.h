//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "PorousFlowDiracKernelScalarBase.h"

/**
 * A point sink of strength flux_coefficient*(u - kappa), where kappa is a scalar variable, along
 * with the constraint that the strength equals a prescribed rate.  Solving for kappa is then the
 * question "what value of kappa delivers the requested rate", which is the structure of a
 * bottom-hole-pressure well stripped of all its physics.
 *
 * This exists to exercise PorousFlowDiracKernelScalarBase: all four of the blocks it assembles are
 * non-zero here, and kappa's only equation comes from a DiracKernel.
 */
class ConstrainedPointSource : public PorousFlowDiracKernelScalarBase
{
public:
  static InputParameters validParams();

  ConstrainedPointSource(const InputParameters & parameters);

  virtual void addPoints() override;

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpJacobianScalar() override;
  virtual Real computeQpScalarResidual() override;
  virtual Real computeQpScalarJacobian() override;
  virtual Real computeQpScalarOffDiagJacobian(unsigned int jvar_num) override;

  /// The strength of the sink per unit (u - kappa)
  const Real _flux_coefficient;

  /// The rate the constraint equation demands
  const Real _rate;

  /// Where the point sink sits
  const Point _point;
};
