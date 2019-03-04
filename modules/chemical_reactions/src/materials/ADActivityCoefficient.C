//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADActivityCoefficient.h"

registerADMooseObject("ChemicalReactionsApp", ADActivityCoefficient);

defineADValidParams(
    ADActivityCoefficient,
    ADMaterial,
    params.addRequiredCoupledVar("primary_species", "Primary species present in the model");
    params.addRequiredCoupledVar("secondary_species", "Equilibrium species present in the model");
    params.addClassDescription("Activity coefficients of primary and secondary species"););

template <ComputeStage compute_stage>
ADActivityCoefficient<compute_stage>::ADActivityCoefficient(const InputParameters & parameters)
  : ADMaterial<compute_stage>(parameters),
    _n_primary(coupledComponents("primary_species")),
    _n_secondary(coupledComponents("secondary_species"))
{
  // Primary species
  _primary_species.resize(_n_primary);
  _gamma_primary.resize(_n_primary);

  for (unsigned int i = 0; i < _n_primary; ++i)
  {
    _primary_species[i] = &adCoupledValue("primary_species", i);
    const std::string psname = getVar("primary_species", i)->name();
    _gamma_primary[i] = &adDeclareADProperty<Real>("activity_coefficient_" + psname);
  }

  // Secondary species
  _secondary_species.resize(_n_secondary);
  _gamma_secondary.resize(_n_secondary);

  for (unsigned int i = 0; i < _n_secondary; ++i)
  {
    _secondary_species[i] = &coupledValue("secondary_species", i);
    const std::string ssname = getVar("secondary_species", i)->name();
    _gamma_secondary[i] = &adDeclareADProperty<Real>("activity_coefficient_" + ssname);
  }
}

template <ComputeStage compute_stage>
void
ADActivityCoefficient<compute_stage>::computeQpProperties()
{
  // Default activity coefficient is unity for all primary and secondary species
  for (unsigned int i = 0; i < _n_primary; ++i)
    (*_gamma_primary[i])[_qp] = 1.0;

  for (unsigned int i = 0; i < _n_secondary; ++i)
    (*_gamma_secondary[i])[_qp] = 1.0;
}
