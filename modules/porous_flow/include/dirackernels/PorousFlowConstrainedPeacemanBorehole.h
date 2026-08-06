//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "PorousFlowPeacemanBorehole.h"

/**
 * A Peaceman borehole whose bottom-hole pressure is a scalar unknown rather than a prescribed
 * function, together with the constraint that fixes it: the flux summed over all the completions
 * equals a requested rate.  The well therefore delivers the rate exactly at every timestep, and the
 * pressure that achieves it is solved for implicitly alongside the reservoir equations.
 *
 * The flux at each completion is the same Peaceman expression the parent computes, so the layer
 * allocation along the completions still follows the local W*lambda*dP: high permeability layers
 * take more, layers that have pressured up take less, and the split shifts as the simulation runs.
 * A source of prescribed rate cannot reproduce that.
 *
 * Unlike the parent this carries no directional gate.  With the bottom-hole pressure an unknown,
 * the gate would make the flux non-differentiable exactly at the operating point the constraint
 * settles on, and it would suppress crossflow between completions, which is physical.  The
 * 'character' parameter is suppressed accordingly, and the direction of flow at a completion
 * follows the sign of the pressure difference.
 */
class PorousFlowConstrainedPeacemanBorehole : public PorousFlowPeacemanBorehole
{
public:
  static InputParameters validParams();

  PorousFlowConstrainedPeacemanBorehole(const InputParameters & parameters);

  virtual void timestepSetup() override;
  virtual void residualSetup() override;

protected:
  /**
   * Choose between delivering the rate and holding the pressure limit, using what the previous
   * residual evaluation recorded.  This is an active set decision taken between Newton iterations,
   * never within one, so that the residual and the Jacobian always describe the same constraint.
   */
  void updateConstraint();

  /// The pressure in the borehole at the current Dirac point, from the scalar unknown
  Real boreholePressure() const;

  /**
   * The flux out of the porespace at the current Dirac point per unit pressure difference and per
   * unit weighting, for the test function of node _i.  The outflow is this multiplied by
   * (porepressure - boreholePressure()), and nothing in it depends on the borehole pressure, which
   * is what makes the derivatives with respect to that unknown so simple.
   * @param current_dirac_ptid The ID assigned to the point by addPoints()
   */
  Real conductance(unsigned current_dirac_ptid) const;

  virtual Real computeQpBaseOutflow(unsigned current_dirac_ptid) const override;
  virtual void computeQpBaseOutflowJacobian(unsigned jvar,
                                            unsigned current_dirac_ptid,
                                            Real & outflow,
                                            Real & outflowp) const override;

  virtual Real computeQpScalarResidual() override;
  virtual Real computeQpScalarJacobian() override;
  virtual Real computeQpScalarOffDiagJacobian(unsigned int jvar_num) override;
  virtual Real computeQpJacobianScalar() override;

  /// The rate this well is asked to deliver, positive for production, or nullptr if it owns no
  /// constraint and so contributes flux alone
  const Function * const _rate;

  /// The bottom-hole pressure this well may not exceed, or nullptr if it is unlimited
  const Function * const _bhp_max;

  /// How many times per timestep the well may switch between its two constraints
  const unsigned int _max_swaps;

  /// Whether the well is currently holding _bhp_max rather than delivering _rate
  bool _clamped;

  /// Switches made so far this timestep
  unsigned int _n_swaps;

  /**
   * The residual the rate constraint would have, accumulated over the completions during the
   * previous residual evaluation and summed across ranks before use.  It is recorded even while
   * the well is holding its pressure limit, because it is what says when the limit stops binding:
   * the flux falls with rising bottom-hole pressure, so a negative value means the rate constraint
   * would settle below the limit.
   */
  Real _rate_mode_residual;

  /// The bottom-hole pressure seen during the previous residual evaluation
  Real _pbh_last;
};
