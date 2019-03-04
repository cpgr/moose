//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADSecondarySpeciesDiffusion.h"

registerADMooseObject("ChemicalReactionsApp", ADSecondarySpeciesDiffusion);

defineADValidParams(
    ADSecondarySpeciesDiffusion,
    ADKernelGrad,
    params.addParam<Real>("weight",
                          1.0,
                          "The weight of the equilibrium species in total concentration");
    params.addCoupledVar("log_k", 0.0, "The equilibrium constant of this equilibrium species");
    params.addParam<Real>(
        "sto",
        1.0,
        "The stoichiometric coefficient of the primary species this kernel operates on");
    params.addParam<std::vector<Real>>(
        "sto_coupled", "The stoichiometric coefficients of coupled primary species");
    params.addRequiredCoupledVar(
        "coupled_primary", "Coupled primary species also present in this equilibrium species");
    params.addRequiredParam<AuxVariableName>("secondary_species",
                                             "The name of this secondary (equilibrium) species");
    params.addClassDescription("Diffusion of secondary (equilibrium) species"););

template <ComputeStage compute_stage>
ADSecondarySpeciesDiffusion<compute_stage>::ADSecondarySpeciesDiffusion(
    const InputParameters & parameters)
  : ADKernelGrad<compute_stage>(parameters),
    _n_primary(coupledComponents("coupled_primary")),
    _secondary_name(adGetParam<AuxVariableName>("secondary_species")),
    _diffusivity(adGetMaterialProperty<Real>("diffusivity")),
    _weight(adGetParam<Real>("weight")),
    _log_k(adCoupledValue("log_k")),
    _sto(adGetParam<Real>("sto")),
    _sto_coupled(adGetParam<std::vector<Real>>("sto_coupled")),
    _gamma_primary(adGetADMaterialProperty<Real>("activity_coefficient_" + _var.name())),
    _gamma_secondary(adGetADMaterialProperty<Real>("activity_coefficient_" + _secondary_name))
{
  // Check that the correct number of stoichiometric coefficients have been provided
  if (_sto_coupled.size() != _n_primary)
    paramError("sto_coupled",
               "The number of stoichiometric coefficients is not equal to the number of coupled "
               "primary species");

  _primary_species.resize(_n_primary);
  _grad_primary_species.resize(_n_primary);
  _gamma_coupled.resize(_n_primary);

  for (unsigned int i = 0; i < _n_primary; ++i)
  {
    _primary_species[i] = &adCoupledValue("coupled_primary", i);
    _grad_primary_species[i] = &adCoupledGradient("coupled_primary", i);
    const std::string psname = getVar("coupled_primary", i)->name();
    _gamma_coupled[i] = &adGetADMaterialProperty<Real>("activity_coefficient_" + psname);
  }
}

template <ComputeStage compute_stage>
ADVectorResidual
ADSecondarySpeciesDiffusion<compute_stage>::precomputeQpResidual()
{
  auto qpresult = _sto * _gamma_primary[_qp] * std::pow(_gamma_primary[_qp] * _u[_qp], _sto - 1.0) *
                  _grad_u[_qp];

  for (unsigned int i = 0; i < _n_primary; ++i)
    qpresult *= std::pow((*_gamma_coupled[i])[_qp] * (*_primary_species[i])[_qp], _sto_coupled[i]);

  const auto gamma_u = std::pow(_gamma_primary[_qp] * _u[_qp], _sto);

  for (unsigned int i = 0; i < _n_primary; ++i)
  {
    auto qptmp =
        gamma_u * _sto_coupled[i] * (*_gamma_coupled[i])[_qp] *
        std::pow((*_gamma_coupled[i])[_qp] * (*_primary_species[i])[_qp], _sto_coupled[i] - 1.0) *
        (*_grad_primary_species[i])[_qp];

    for (unsigned int j = 0; j < _n_primary; ++j)
      if (j != i)
        qptmp *= std::pow((*_gamma_coupled[i])[_qp] * (*_primary_species[j])[_qp], _sto_coupled[j]);

    qpresult += qptmp;
  }

  mooseAssert(_gamma_secondary[_qp] > 0.0, "Activity coefficient must be greater than zero");

  return _weight * std::pow(10.0, _log_k[_qp]) * _diffusivity[_qp] * qpresult /
         _gamma_secondary[_qp];
}
