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
  params.addPrivateParam<std::string>("pf_material_type", "halite");
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
    _dprecipitated_salt_dvar(
        _nodal_material
            ? getMaterialProperty<std::vector<Real>>("dPorousFlow_precipitated_salt_nodal_dvar")
            : getMaterialProperty<std::vector<Real>>("dPorousFlow_precipitated_salt_qp_dvar")),
    _saturation(_nodal_material
                    ? getMaterialProperty<std::vector<Real>>("PorousFlow_saturation_nodal")
                    : getMaterialProperty<std::vector<Real>>("PorousFlow_saturation_qp")),
    _dsaturation_dvar(_nodal_material ? getMaterialProperty<std::vector<std::vector<Real>>>(
                                            "dPorousFlow_saturation_nodal_dvar")
                                      : getMaterialProperty<std::vector<std::vector<Real>>>(
                                            "dPorousFlow_saturation_qp_dvar")),
    _fluid_density(
        _nodal_material
            ? getMaterialProperty<std::vector<Real>>("PorousFlow_fluid_phase_density_nodal")
            : getMaterialProperty<std::vector<Real>>("PorousFlow_fluid_phase_density_qp")),
    _dfluid_density_dvar(_nodal_material ? getMaterialProperty<std::vector<std::vector<Real>>>(
                                               "dPorousFlow_fluid_phase_density_nodal_dvar")
                                         : getMaterialProperty<std::vector<std::vector<Real>>>(
                                               "dPorousFlow_fluid_phase_density_qp_dvar")),
    _porosity_old(_nodal_material ? getMaterialPropertyOld<Real>("PorousFlow_porosity_nodal")
                                  : getMaterialPropertyOld<Real>("PorousFlow_porosity_qp")),
    _halite_volume_fraction(_nodal_material
                                ? declareProperty<Real>("PorousFlow_halite_volume_fraction_nodal")
                                : declareProperty<Real>("PorousFlow_halite_volume_fraction_qp")),
    _dhalite_volume_fraction_dvar(
        _nodal_material
            ? declareProperty<std::vector<Real>>("dPorousFlow_halite_volume_fraction_nodal_dvar")
            : declareProperty<std::vector<Real>>("dPorousFlow_halite_volume_fraction_qp_dvar"))
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

  const Real coeff = _porosity_old[_qp] / _halite_density;
  _halite_volume_fraction[_qp] = _precipitated_salt[_qp] * coeff * fluid_mass;

  // Derivatives wrt the PorousFlow variables (phi_old carries no current-step derivative)
  _dhalite_volume_fraction_dvar[_qp].assign(_num_var, 0.0);
  for (const auto v : make_range(_num_var))
  {
    Real dfluid_mass = 0.0;
    for (const auto ph : make_range(_num_phases))
      dfluid_mass += _dsaturation_dvar[_qp][ph][v] * _fluid_density[_qp][ph] +
                     _saturation[_qp][ph] * _dfluid_density_dvar[_qp][ph][v];

    _dhalite_volume_fraction_dvar[_qp][v] = coeff * (_dprecipitated_salt_dvar[_qp][v] * fluid_mass +
                                                     _precipitated_salt[_qp] * dfluid_mass);
  }
}
