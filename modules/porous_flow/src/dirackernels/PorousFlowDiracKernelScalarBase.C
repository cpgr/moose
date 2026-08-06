//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PorousFlowDiracKernelScalarBase.h"

#include "Assembly.h"
#include "FEProblemBase.h"
#include "MooseVariableScalar.h"
#include "SystemBase.h"

#include "libmesh/quadrature.h"

InputParameters
PorousFlowDiracKernelScalarBase::validParams()
{
  InputParameters params = DiracKernel::validParams();
  // Derived classes are expected to rename this to something meaningful, as KernelScalarBase's
  // consumers do
  params.addCoupledVar("scalar_variable", "Primary coupled scalar variable");
  params.addParam<bool>("compute_scalar_residuals", true, "Whether to compute scalar residuals");
  params.addParam<bool>(
      "compute_field_residuals", true, "Whether to compute residuals for the field variable.");
  return params;
}

PorousFlowDiracKernelScalarBase::PorousFlowDiracKernelScalarBase(const InputParameters & parameters)
  : DiracKernel(parameters),
    _use_scalar(isParamValid("scalar_variable")),
    _compute_scalar_residuals(_use_scalar && getParam<bool>("compute_scalar_residuals")),
    _compute_field_residuals(getParam<bool>("compute_field_residuals")),
    _kappa_var_ptr(_use_scalar ? getScalarVar("scalar_variable", 0) : nullptr),
    _kappa_var(_use_scalar ? _kappa_var_ptr->number() : 0),
    _k_order(_use_scalar ? _kappa_var_ptr->order() : 0),
    _kappa(_use_scalar ? (_is_implicit ? _kappa_var_ptr->sln() : _kappa_var_ptr->slnOld()) : _zero)
{
}

Real
PorousFlowDiracKernelScalarBase::computeQpScalarResidual()
{
  mooseError("A scalar_variable has been set and compute_scalar_residuals=true, but the "
             "computeQpScalarResidual method was not overridden.  Accidental call of base class?");
  return 0.0;
}

void
PorousFlowDiracKernelScalarBase::computeResidual()
{
  if (_compute_field_residuals)
    DiracKernel::computeResidual();

  if (_compute_scalar_residuals)
    computeScalarResidual();
}

void
PorousFlowDiracKernelScalarBase::computeScalarResidual()
{
  std::vector<Real> scalar_residuals(_k_order, 0.0);

  forEachPoint(
      [&](const Real point_value)
      {
        for (_h = 0; _h < _k_order; ++_h)
          scalar_residuals[_h] += point_value * computeQpScalarResidual();
      });

  addResiduals(
      _assembly, scalar_residuals, _kappa_var_ptr->dofIndices(), _kappa_var_ptr->scalingFactor());

  // The Dirac Jacobian loop flushes cached Jacobians once it has finished, but the residual loop
  // has no such step: it only calls Assembly::addResidual, which visits field variables alone and
  // leaves the cache untouched.  So the scalar row has to be flushed here or it never reaches the
  // global vector.  Doing it per element does not double count, because
  // Assembly::addCachedResidualDirectly clears the cache as it adds.
  _fe_problem.addCachedResidual(_tid);
}

void
PorousFlowDiracKernelScalarBase::computeJacobian()
{
  if (_compute_field_residuals)
    DiracKernel::computeJacobian();

  if (!_use_scalar)
    return;

  if (_compute_field_residuals)
    computeFieldScalarJacobian();

  if (_compute_scalar_residuals)
    computeScalarJacobian();
}

void
PorousFlowDiracKernelScalarBase::computeOffDiagJacobian(const unsigned int jvar_num)
{
  // ComputeDiracThread reaches this kernel's own variable as one of the coupling pairs, so the
  // blocks that must be assembled exactly once per element are handled there
  if (jvar_num == _var.number())
    computeJacobian();
  else if (_compute_field_residuals)
    DiracKernel::computeOffDiagJacobian(jvar_num);

  if (_compute_scalar_residuals)
    computeScalarOffDiagJacobian(jvar_num);
}

void
PorousFlowDiracKernelScalarBase::computeScalarJacobian()
{
  DenseMatrix<Real> ke(_k_order, _k_order);

  forEachPoint(
      [&](const Real point_value)
      {
        for (_h = 0; _h < _k_order; ++_h)
          for (_l = 0; _l < _k_order; ++_l)
            ke(_h, _l) += point_value * computeQpScalarJacobian();
      });

  addJacobian(_assembly,
              ke,
              _kappa_var_ptr->dofIndices(),
              _kappa_var_ptr->dofIndices(),
              _kappa_var_ptr->scalingFactor());
}

void
PorousFlowDiracKernelScalarBase::computeScalarOffDiagJacobian(const unsigned int jvar_num)
{
  const auto & jvar = getVariable(jvar_num);
  DenseMatrix<Real> ke(_k_order, jvar.phiSize());

  forEachPoint(
      [&](const Real point_value)
      {
        for (_h = 0; _h < _k_order; ++_h)
          for (_j = 0; _j < jvar.phiSize(); ++_j)
            ke(_h, _j) += point_value * computeQpScalarOffDiagJacobian(jvar_num);
      });

  addJacobian(_assembly,
              ke,
              _kappa_var_ptr->dofIndices(),
              jvar.dofIndices(),
              _kappa_var_ptr->scalingFactor());
}

void
PorousFlowDiracKernelScalarBase::computeFieldScalarJacobian()
{
  DenseMatrix<Real> ke(_test.size(), _k_order);

  forEachPoint(
      [&](const Real point_value)
      {
        for (_i = 0; _i < _test.size(); ++_i)
          for (_l = 0; _l < _k_order; ++_l)
            ke(_i, _l) += point_value * computeQpJacobianScalar();
      });

  addJacobian(_assembly, ke, _var.dofIndices(), _kappa_var_ptr->dofIndices(), _var.scalingFactor());
}
