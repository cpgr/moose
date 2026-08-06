//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "DiracKernel.h"

class MooseVariableScalar;

/**
 * A DiracKernel that, as well as contributing to its own field variable, assembles the row of a
 * coupled scalar variable.  The scalar is "scalar_variable" in the input file and "kappa" in the
 * source, following KernelScalarBase, of which this is the point-source analogue.
 *
 * ComputeDiracKernel threading only ever calls computeResidual(), computeJacobian() and
 * computeOffDiagJacobian() on a DiracKernel.  Its Jacobian loop walks Assembly::couplingEntries(),
 * which holds pairs of *field* variables, so computeOffDiagJacobianScalar() is never reached for a
 * DiracKernel.  All four blocks (field-field, field-kappa, kappa-field, kappa-kappa) are therefore
 * assembled from inside the three calls that are made, and no framework loop needs to change.
 *
 * Without a "scalar_variable" this class behaves exactly as a plain DiracKernel, so it is safe as
 * the base of existing point sources.
 */
class PorousFlowDiracKernelScalarBase : public DiracKernel
{
public:
  static InputParameters validParams();

  PorousFlowDiracKernelScalarBase(const InputParameters & parameters);

  virtual void computeResidual() override;
  virtual void computeJacobian() override;
  virtual void computeOffDiagJacobian(unsigned int jvar_num) override;

protected:
  /// Contribution to the kappa row at the current Dirac point, for scalar component _h
  virtual Real computeQpScalarResidual();

  /// d(kappa row _h)/d(kappa _l) at the current Dirac point
  virtual Real computeQpScalarJacobian() { return 0.0; }

  /**
   * d(kappa row _h)/d(jvar _j) at the current Dirac point
   * @param jvar_num The number of the field variable being differentiated with respect to
   */
  virtual Real computeQpScalarOffDiagJacobian(unsigned int /*jvar_num*/) { return 0.0; }

  /// d(field row _i)/d(kappa _l) at the current Dirac point
  virtual Real computeQpJacobianScalar() { return 0.0; }

  /**
   * Run `body` at every Dirac point of this kernel that lies in the current element, having set
   * _qp, _current_point and `point_value`.  The three DiracKernel assembly routines each carry
   * their own copy of this loop; the four blocks assembled here would carry four more.
   * @param body Called once per active point, with the weight of that point
   */
  template <typename Body>
  void forEachPoint(Body && body);

  /// Assemble the kappa row of the residual
  void computeScalarResidual();

  /// Assemble d(kappa row)/d(kappa)
  void computeScalarJacobian();

  /// Assemble d(kappa row)/d(jvar)
  void computeScalarOffDiagJacobian(unsigned int jvar_num);

  /// Assemble d(field row)/d(kappa)
  void computeFieldScalarJacobian();

  /// Whether a scalar variable was supplied
  const bool _use_scalar;

  /// Whether this instance assembles the kappa row.  Only one object may own a shared constraint
  const bool _compute_scalar_residuals;

  /// Whether this instance assembles its field rows
  const bool _compute_field_residuals;

  /// (Pointer to) the scalar variable this kernel assembles into, or nullptr if there is none
  const MooseVariableScalar * const _kappa_var_ptr;

  /// The number of the scalar variable
  const unsigned int _kappa_var;

  /// Order of the scalar variable, ie how many dofs it carries
  const unsigned int _k_order;

  /// Current value of the scalar variable
  const VariableValue & _kappa;

  ///@{ Used to iterate over the components of the scalar variable
  unsigned int _h;
  unsigned int _l;
  ///@}
};

template <typename Body>
void
PorousFlowDiracKernelScalarBase::forEachPoint(Body && body)
{
  const std::vector<Real> * const point_values =
      _drop_duplicate_points ? nullptr
                             : &_local_dirac_kernel_info.getPoints()[_current_elem].second;
  unsigned int local_qp = 0;
  Real point_value = 1.0;

  for (_qp = 0; _qp < _qrule->n_points(); _qp++)
  {
    _current_point = _physical_point[_qp];
    if (isActiveAtPoint(_current_elem, _current_point))
    {
      if (!_drop_duplicate_points)
        point_value = (*point_values)[local_qp++];

      body(point_value);
    }
  }
}
