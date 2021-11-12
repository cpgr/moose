//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "PorousFlowCapillaryPressure.h"

/**
 * Brooks-Corey effective saturation, capillary pressure and relative
 * permeability functions.
 *
 * From Brooks, R. H. and A. T. Corey (1966), Properties of porous media affecting
 * fluid flow, J. Irrig. Drain. Div., 6, 61
 */
class PorousFlowCapillaryPressureBC : public PorousFlowCapillaryPressure
{
public:
  static InputParameters validParams();

  PorousFlowCapillaryPressureBC(const InputParameters & parameters);

  virtual Real capillaryPressureCurve(Real saturation, unsigned int qp = 0) const override;
  virtual Real dCapillaryPressureCurve(Real saturation, unsigned int qp = 0) const override;
  virtual Real d2CapillaryPressureCurve(Real saturation, unsigned int qp = 0) const override;

  virtual Real effectiveSaturation(Real pc, unsigned int qp = 0) const override;
  virtual Real dEffectiveSaturation(Real pc, unsigned int qp = 0) const override;
  virtual Real d2EffectiveSaturation(Real pc, unsigned int qp = 0) const override;

protected:
  /// Brooks-Corey exponent lambda
  const Real _lambda;
  /// Threshold entry pressure
  const VariableValue & _pe;

private:
  /// Check that provided entry pressure is greater than or equal to zero
  void checkEntryPressure(unsigned int qp = 0) const;
};
