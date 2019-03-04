//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADEquilibriumSpeciesTimeDerivative.h"

registerADMooseObject("ChemicalReactionsApp", ADEquilibriumSpeciesTimeDerivative);

defineADValidParams(
    ADEquilibriumSpeciesTimeDerivative,
    ADTimeDerivative,
    params.addParam<Real>("weight",
                          1.0,
                          "The weight of the equilibrium species in total concentration");
    params.addCoupledVar("log_k", 0.0, "The equilibrium constant of this equilibrium species");
    params.addParam<Real>(
        "sto_u",
        1.0,
        "The stoichiometric coefficient of the primary species this kernel operates on");
    params.addParam<std::vector<Real>>(
        "sto_v", "The stoichiometric coefficients of coupled primary species");
    params.addCoupledVar("v", "Coupled primary species constituting the equilibrium species");
    params.addClassDescription("Derivative of equilibrium species concentration wrt time"););

template <ComputeStage compute_stage>
ADEquilibriumSpeciesTimeDerivative<compute_stage>::ADEquilibriumSpeciesTimeDerivative(
    const InputParameters & parameters)
  : ADTimeDerivative<compute_stage>(parameters),
    _weight(adGetParam<Real>("weight")),
    _log_k(adCoupledValue("log_k")),
    _sto_u(adGetParam<Real>("sto_u")),
    _sto_v(adGetParam<std::vector<Real>>("sto_v")),
    // _gamma_u(adCoupledValue("gamma_u")),
    // _gamma_u_old(adCoupledValueOld("gamma_u")),
    // _gamma_eq(adCoupledValue("gamma_eq")),
    // _gamma_eq_old(adCoupledValueOld("gamma_eq")),
    _porosity(adGetMaterialProperty<Real>("porosity")),
    _porosity_old(adGetMaterialPropertyOld<Real>("porosity"))
{
  // const unsigned int n = coupledComponents("v");
  //
  // // Check that the correct number of coupled values have been provided
  // if (_sto_v.size() != n)
  //   mooseError("The number of stoichiometric coefficients in sto_v is not equal to the number of
  //   "
  //              "coupled species in ",
  //              _name);
  //
  // if (isCoupled("gamma_v"))
  //   if (coupledComponents("gamma_v") != n)
  //     mooseError("The number of activity coefficients in gamma_v is not equal to the number of "
  //                "coupled species in ",
  //                _name);
  //
  // _vars.resize(n);
  // _v_vals.resize(n);
  // _v_vals_old.resize(n);
  // _gamma_v.resize(n);
  // _gamma_v_old.resize(n);
  //
  // for (unsigned int i = 0; i < _vars.size(); ++i)
  // {
  //   _vars[i] = coupled("v", i);
  //   _v_vals[i] = &coupledValue("v", i);
  //   _v_vals_old[i] = &coupledValueOld("v", i);
  //   // If gamma_v has been supplied, use those values, but if not, use the default value
  //   _gamma_v[i] = (isCoupled("gamma_v") ? &coupledValue("gamma_v", i) :
  //   &coupledValue("gamma_v")); _gamma_v_old[i] =
  //       (isCoupled("gamma_v") ? &coupledValueOld("gamma_v", i) : &coupledValue("gamma_v"));
  // }
}

template <ComputeStage compute_stage>
ADResidual
ADEquilibriumSpeciesTimeDerivative<compute_stage>::precomputeQpResidual()
{
  // mooseAssert(_gamma_eq[_qp] > 0.0, "Activity coefficient must be greater than zero");
  //
  // // Contribution due to primary species that this kernel acts on
  // Real val_new =
  //     std::pow(10.0, _log_k[_qp]) * std::pow(_gamma_u[_qp] * _u[_qp], _sto_u) / _gamma_eq[_qp];
  // Real val_old = std::pow(10.0, _log_k[_qp]) * std::pow(_gamma_u_old[_qp] * _u_old[_qp], _sto_u)
  // /
  //                _gamma_eq_old[_qp];
  //
  // // Contribution due to coupled primary species
  // for (unsigned int i = 0; i < _vars.size(); ++i)
  // {
  //   val_new *= std::pow((*_gamma_v[i])[_qp] * (*_v_vals[i])[_qp], _sto_v[i]);
  //   val_old *= std::pow((*_gamma_v_old[i])[_qp] * (*_v_vals_old[i])[_qp], _sto_v[i]);
  // }
  //
  // return _porosity[_qp] * _weight * _test[_i][_qp] * (val_new - val_old) / _dt;
  return 0;
}
