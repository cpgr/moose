/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/


#ifndef POROUSFLOWDICTATOR_H
#define POROUSFLOWDICTATOR_H

#include "GeneralUserObject.h"
#include "Coupleable.h"
#include "ZeroInterface.h"

class PorousFlowDictator;

template<>
InputParameters validParams<PorousFlowDictator>();

/**
 * This holds maps between the nonlinear variables
 * used in a PorousFlow simulation and the
 * variable number used internally by MOOSE, as
 * well as the number of fluid phases and
 * the number of fluid components.
 *
 * All PorousFlow Materials and Kernels calculate
 * and use derivatives with respect to all the variables
 * mentioned in this Object, at least in principal
 * (in practice they may be lazy and not compute
 * all derivatives).
 */
class PorousFlowDictator :
  public GeneralUserObject,
  public Coupleable,
  public ZeroInterface
{
 public:
  PorousFlowDictator(const InputParameters & parameters);

  void initialize();
  void execute();
  void finalize();

  /**
   * The number of PorousFlow variables.  Materials
   * and Kernels will calculate and use derivatives
   * with respect to these variables in the Jacobian
   */
  unsigned int num_variables() const;

  /// the number of fluid phases
  unsigned int num_phases() const;

  /// the number of fluid components
  unsigned int num_components() const;

  /**
   * the PorousFlow variable number
   * @param moose_var_num the MOOSE variable number
   * eg if porflow_vars = 'pwater pgas', and the variables in
   * the simulation are 'energy pwater pgas shape'
   * then porflow_var_num(2) = 1
   */
  unsigned int porflow_var_num(unsigned int moose_var_num) const;

  /**
   * returns true if moose_var_num is not a porflow var
   * @param moose_var_num the MOOSE variable number
   * eg if porflow_vars = 'pwater pgas', and the variables in
   * the simulation are 'energy pwater pgas shape'
   * then not_porflow_var(0) = true, not_porflow_var(1) = false
   */
  bool not_porflow_var(unsigned int moose_var_num) const;

  /**
   * a vector of pointers to VariableValues
   * @param porflow_var_num the porflow variable number
   * eg if porflow_vars = 'pwater pgas', then
   * (*porflow_vals(1))[qp] = pgas evaluated at quadpoint qp
   * Also porflow_vals(i) = &coupledValue
   */
  const VariableValue * porflow_vals(unsigned int porflow_var_num) const;

  /**
   * a vector of pointers to old VariableValues
   * @param porflow_var_num the porflow variable number
   * eg if porflow_vars = 'pwater pgas', then
   * (*porflow_vals_old(1))[qp] = old pgas evaluated at quadpoint qp
   * Also porflow_vals_old(i) = &coupledValueOld
   */
  const VariableValue * porflow_vals_old(unsigned int porflow_var_num) const;

  /**
   * a vector of pointers to grad(Variable)
   * @param porflow_var_num the porflow variable number
   * eg if porflow_vars = 'pwater pgas', then
   * (*grad_var(1))[qp] = grad(pgas) evaluated at quadpoint qp
   * Also grad_var(i) = &coupledGradient
   */
  const VariableGradient * grad_var(unsigned int porflow_var_num) const;

  /**
   * The nodal variable values for the given porflow_var_num
   * To extract the value of porflow variable "pvar", at
   * node i, use (*PorousFlowDictator.nodal_var(pvar))[i]
   * @param porflow_var_num the porflow variable number
   */
  const VariableValue * nodal_var(unsigned int porflow_var_num) const;

  /**
   * The old nodal variable values for the given porflow_var_num
   * @param porflow_var_num the porflow variable number
   */
  const VariableValue * nodal_var_old(unsigned int porflow_var_num) const;

  /**
   * Dummy pressure variable name for use in derivatives using the
   * DerivativeMaterialInterface
   */
  const VariableName pressureVariableNameDummy() const;

  /**
   * Dummy saturation variable name for use in derivatives using the
   * DerivativeMaterialInterface
   */
  const VariableName saturationVariableNameDummy() const;

  /**
   * Dummy temperature variable name for use in derivatives using the
   * DerivativeMaterialInterface
   */
  const VariableName temperatureVariableNameDummy() const;

  /**
   * Dummy mass fraction variable name for use in derivatives using the
   * DerivativeMaterialInterface
   */
  const VariableName massFractionVariableNameDummy() const;

 protected:

  /// number of porousflow variables
  unsigned int _num_variables;

  /// number of fluid phases
  unsigned int _num_phases;

  /// number of fluid components
  unsigned int _num_components;

  /// _moose_var_num[i] = the moose variable number corresponding to porflow variable i
  std::vector<unsigned int> _moose_var_num;

  /// _pf_var_num[i] = the porflow variable corresponding to moose variable i
  std::vector<unsigned int> _pf_var_num;

  /// moose_var_value[i] = values of porflow variable i
  std::vector<const VariableValue *> _moose_var_value; // this is a vector of pointers to VariableValues

  /// moose_var_value_old[i] = old values of porflow variable i
  std::vector<const VariableValue *> _moose_var_value_old;

  /// moose_var_value[i] = values of porflow variable i
  std::vector<const VariableValue *> _moose_nodal_var_value; // this is a vector of pointers to VariableValues

  /// moose_var_value_old[i] = old values of porflow variable i
  std::vector<const VariableValue *> _moose_nodal_var_value_old;

  /// moose_grad_var[i] = gradient values of porflow variable i
  std::vector<const VariableGradient *> _moose_grad_var;
};

#endif // POROUSFLOWDICTATOR_H
