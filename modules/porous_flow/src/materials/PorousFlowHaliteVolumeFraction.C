//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PorousFlowHaliteVolumeFraction.h"

registerMooseObject("PorousFlowApp", PorousFlowHaliteVolumeFraction);

InputParameters
PorousFlowHaliteVolumeFraction::validParams()
{
  InputParameters params = PorousFlowMaterialVectorBase::validParams();
  params.addRangeCheckedParam<Real>(
      "halite_density",
      2165.0,
      "halite_density > 0",
      "Density of solid halite (kg/m^3) used to convert the precipitated salt mass to a volume "
      "fraction");
  params.addClassDescription(
      "Computes the solid halite volume fraction (m^3 halite / m^3 porous medium) from the "
      "local-equilibrium precipitated-salt mass fraction reported by a salt-precipitating fluid "
      "state");
  return params;
}

PorousFlowHaliteVolumeFraction::PorousFlowHaliteVolumeFraction(const InputParameters & parameters)
  : PorousFlowMaterialVectorBase(parameters),
    _halite_density(getParam<Real>("halite_density")),
    _precipitated_salt(_nodal_material
                           ? getMaterialProperty<Real>("PorousFlow_precipitated_salt_nodal")
                           : getMaterialProperty<Real>("PorousFlow_precipitated_salt_qp")),
    _saturation(_nodal_material
                    ? getMaterialProperty<std::vector<Real>>("PorousFlow_saturation_nodal")
                    : getMaterialProperty<std::vector<Real>>("PorousFlow_saturation_qp")),
    _fluid_density(
        _nodal_material
            ? getMaterialProperty<std::vector<Real>>("PorousFlow_fluid_phase_density_nodal")
            : getMaterialProperty<std::vector<Real>>("PorousFlow_fluid_phase_density_qp")),
    _porosity_old(_nodal_material ? getMaterialPropertyOld<Real>("PorousFlow_porosity_nodal")
                                  : getMaterialPropertyOld<Real>("PorousFlow_porosity_qp")),
    _halite_volume_fraction(_nodal_material
                                ? declareProperty<Real>("PorousFlow_halite_volume_fraction_nodal")
                                : declareProperty<Real>("PorousFlow_halite_volume_fraction_qp"))
{
}

void
PorousFlowHaliteVolumeFraction::initQpStatefulProperties()
{
  // Simulations start with no precipitated halite (the porous medium begins undersaturated). The
  // old porosity is unavailable at t = 0, so the algebraic conversion below cannot be applied here.
  _halite_volume_fraction[_qp] = 0.0;
}

void
PorousFlowHaliteVolumeFraction::computeQpProperties()
{
  // Fluid mass per unit pore volume, Sum_ph(S_ph rho_ph)
  Real fluid_mass = 0.0;
  for (const auto ph : make_range(_num_phases))
    fluid_mass += _saturation[_qp][ph] * _fluid_density[_qp][ph];

  _halite_volume_fraction[_qp] =
      _precipitated_salt[_qp] * _porosity_old[_qp] * fluid_mass / _halite_density;
}
