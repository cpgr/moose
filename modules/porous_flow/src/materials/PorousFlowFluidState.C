//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PorousFlowFluidState.h"
#include "PorousFlowCapillaryPressure.h"

registerMooseObject("PorousFlowApp", PorousFlowFluidState);

template <>
InputParameters
validParams<PorousFlowFluidState>()
{
  InputParameters params = validParams<PorousFlowVariableBase>();
  params.addRequiredCoupledVar("gas_porepressure",
                               "Variable that is the porepressure of the gas phase");
  params.addRequiredCoupledVar("z", "Total mass fraction of component i summed over all phases");
  params.addCoupledVar("xnacl", 0, "The salt mass fraction in the brine (kg/kg)");
  MooseEnum unit_choice("Kelvin=0 Celsius=1", "Kelvin");
  params.addParam<MooseEnum>(
      "temperature_unit", unit_choice, "The unit of the temperature variable");
  params.addRequiredParam<UserObjectName>("capillary_pressure",
                                          "Name of the UserObject defining the capillary pressure");
  params.addRequiredParam<UserObjectName>("fluid_state", "Name of the FluidState UserObject");
  params.addPrivateParam<std::string>("pf_material_type", "fluid_state");
  params.addClassDescription("Class for fluid state calculations using persistent primary "
                             "variables and a vapor-liquid flash");
  return params;
}

PorousFlowFluidState::PorousFlowFluidState(const InputParameters & parameters)
  : PorousFlowVariableBase(parameters),

    _gas_porepressure(_nodal_material ? coupledDofValues("gas_porepressure")
                                      : coupledValue("gas_porepressure")),
    _gas_gradp_qp(coupledGradient("gas_porepressure")),
    _gas_porepressure_varnum(coupled("gas_porepressure")),
    _pvar(_dictator.isPorousFlowVariable(_gas_porepressure_varnum)
              ? _dictator.porousFlowVariableNum(_gas_porepressure_varnum)
              : 0),

    _num_Z_vars(coupledComponents("z")),
    _Xnacl(_nodal_material ? coupledDofValues("xnacl") : coupledValue("xnacl")),
    _grad_Xnacl_qp(coupledGradient("xnacl")),
    _Xnacl_varnum(coupled("xnacl")),
    _Xvar(_dictator.isPorousFlowVariable(_Xnacl_varnum)
              ? _dictator.porousFlowVariableNum(_Xnacl_varnum)
              : 0),

    _fs(getUserObject<PorousFlowFluidStateMultiComponentBase>("fluid_state")),
    _aqueous_phase_number(_fs.aqueousPhaseIndex()),
    _gas_phase_number(_fs.gasPhaseIndex()),
    _aqueous_fluid_component(_fs.aqueousComponentIndex()),
    _gas_fluid_component(_fs.gasComponentIndex()),
    _salt_component(_fs.saltComponentIndex()),

    _temperature(_nodal_material ? getMaterialProperty<Real>("PorousFlow_temperature_nodal")
                                 : getMaterialProperty<Real>("PorousFlow_temperature_qp")),
    _gradT_qp(getMaterialProperty<RealGradient>("PorousFlow_grad_temperature_qp")),
    _dtemperature_dvar(
        _nodal_material
            ? getMaterialProperty<std::vector<Real>>("dPorousFlow_temperature_nodal_dvar")
            : getMaterialProperty<std::vector<Real>>("dPorousFlow_temperature_qp_dvar")),

    _mass_frac(_nodal_material
                   ? declareProperty<std::vector<std::vector<Real>>>("PorousFlow_mass_frac_nodal")
                   : declareProperty<std::vector<std::vector<Real>>>("PorousFlow_mass_frac_qp")),
    _grad_mass_frac_qp(_nodal_material ? nullptr
                                       : &declareProperty<std::vector<std::vector<RealGradient>>>(
                                             "PorousFlow_grad_mass_frac_qp")),
    _dmass_frac_dvar(_nodal_material ? declareProperty<std::vector<std::vector<std::vector<Real>>>>(
                                           "dPorousFlow_mass_frac_nodal_dvar")
                                     : declareProperty<std::vector<std::vector<std::vector<Real>>>>(
                                           "dPorousFlow_mass_frac_qp_dvar")),
    _saturation_old(_nodal_material
                        ? getMaterialPropertyOld<std::vector<Real>>("PorousFlow_saturation_nodal")
                        : getMaterialPropertyOld<std::vector<Real>>("PorousFlow_saturation_qp")),

    _fluid_density(_nodal_material
                       ? declareProperty<std::vector<Real>>("PorousFlow_fluid_phase_density_nodal")
                       : declareProperty<std::vector<Real>>("PorousFlow_fluid_phase_density_qp")),
    _dfluid_density_dvar(_nodal_material ? declareProperty<std::vector<std::vector<Real>>>(
                                               "dPorousFlow_fluid_phase_density_nodal_dvar")
                                         : declareProperty<std::vector<std::vector<Real>>>(
                                               "dPorousFlow_fluid_phase_density_qp_dvar")),
    _fluid_viscosity(_nodal_material
                         ? declareProperty<std::vector<Real>>("PorousFlow_viscosity_nodal")
                         : declareProperty<std::vector<Real>>("PorousFlow_viscosity_qp")),
    _dfluid_viscosity_dvar(
        _nodal_material
            ? declareProperty<std::vector<std::vector<Real>>>("dPorousFlow_viscosity_nodal_dvar")
            : declareProperty<std::vector<std::vector<Real>>>("dPorousFlow_viscosity_qp_dvar")),

    _fluid_enthalpy(
        _nodal_material
            ? declareProperty<std::vector<Real>>("PorousFlow_fluid_phase_enthalpy_nodal")
            : declareProperty<std::vector<Real>>("PorousFlow_fluid_phase_enthalpy_qp")),
    _dfluid_enthalpy_dvar(_nodal_material ? declareProperty<std::vector<std::vector<Real>>>(
                                                "dPorousFlow_fluid_phase_enthalpy_nodal_dvar")
                                          : declareProperty<std::vector<std::vector<Real>>>(
                                                "dPorousFlow_fluid_phase_enthalpy_qp_dvar")),

    _T_c2k(getParam<MooseEnum>("temperature_unit") == 0 ? 0.0 : 273.15),
    _is_initqp(false),
    _pc(getUserObject<PorousFlowCapillaryPressure>("capillary_pressure")),
    _pidx(_fs.getPressureIndex()),
    _Tidx(_fs.getTemperatureIndex()),
    _Zidx(_num_Z_vars),
    _Xidx(_fs.getXIndex())
{
  // Check that the number of phases in the fluidstate class is also provided in the Dictator
  if (_fs.numPhases() != _num_phases)
    mooseError(name(),
               ": only ",
               _fs.numPhases(),
               " phases are allowed. Please check the number of phases entered in the dictator is "
               "correct");

  // Check that the number of Z variables is correct for the chosen fluidstate
  if (_fs.numZVars() != _num_Z_vars)
    mooseError(name(), _fs.numZVars(), " variables must be supplied");

  // Z variable indexes
  for (unsigned int z = 0; z < _num_Z_vars; ++z)
    _Zidx[z] = _fs.getZIndex() + z;

  // Store all total mass fractions and associated variable numbers
  _Z.resize(_num_Z_vars);
  _Zqp.resize(_num_Z_vars);
  _gradZ_qp.resize(_num_Z_vars);
  _Z_varnum.resize(_num_Z_vars);
  _Zvar.resize(_num_Z_vars);

  for (unsigned int i = 0; i < _num_Z_vars; ++i)
  {
    _Z[i] = (_nodal_material ? &coupledDofValues("z", i) : &coupledValue("z", i));
    _gradZ_qp[i] = &coupledGradient("z", i);
    _Z_varnum[i] = coupled("z", i);
    _Zvar[i] = (_dictator.isPorousFlowVariable(_Z_varnum[i])
                    ? _dictator.porousFlowVariableNum(_Z_varnum[i])
                    : 0);
  }

  // Set the size of the FluidStateProperties vector
  _fsp.resize(_num_phases, FluidStateProperties(_num_components));
}

void
PorousFlowFluidState::thermophysicalProperties()
{
  // The FluidProperty objects use temperature in K
  const Real Tk = _temperature[_qp] + _T_c2k;

  // Z values at current qp
  for (unsigned int z = 0; z < _num_Z_vars; ++z)
    _Zqp[z] = (*_Z[z])[_qp];

  _fs.thermophysicalProperties(_gas_porepressure[_qp], Tk, _Xnacl[_qp], _Zqp, _qp, _fsp);
}

void
PorousFlowFluidState::initQpStatefulProperties()
{
  _is_initqp = true;
  // Set the size of pressure and saturation vectors
  PorousFlowVariableBase::initQpStatefulProperties();

  // Set the size of all other vectors
  setMaterialVectorSize();

  thermophysicalProperties();

  // Set the initial values of the properties at the nodes.
  // Note: not required for qp materials as no old values at the qps are requested
  if (_nodal_material)
  {
    thermophysicalProperties();

    for (unsigned int ph = 0; ph < _num_phases; ++ph)
    {
      _saturation[_qp][ph] = _fsp[ph].saturation.value();
      _porepressure[_qp][ph] = _fsp[ph].pressure.value();
      _fluid_density[_qp][ph] = _fsp[ph].density.value();
      _fluid_viscosity[_qp][ph] = _fsp[ph].viscosity.value();
      _fluid_enthalpy[_qp][ph] = _fsp[ph].enthalpy.value();

      for (unsigned int comp = 0; comp < _num_components; ++comp)
        _mass_frac[_qp][ph][comp] = _fsp[ph].mass_fraction[comp].value();
    }
  }
}

void
PorousFlowFluidState::computeQpProperties()
{
  _is_initqp = false;
  // Prepare the derivative vectors
  PorousFlowVariableBase::computeQpProperties();

  // Set the size of all other vectors
  setMaterialVectorSize();

  // Calculate all required thermophysical properties
  thermophysicalProperties();

  for (unsigned int ph = 0; ph < _num_phases; ++ph)
  {
    _saturation[_qp][ph] = _fsp[ph].saturation.value();
    _porepressure[_qp][ph] = _fsp[ph].pressure.value();
    _fluid_density[_qp][ph] = _fsp[ph].density.value();
    _fluid_viscosity[_qp][ph] = _fsp[ph].viscosity.value();
    _fluid_enthalpy[_qp][ph] = _fsp[ph].enthalpy.value();

    for (unsigned int comp = 0; comp < _num_components; ++comp)
      _mass_frac[_qp][ph][comp] = _fsp[ph].mass_fraction[comp].value();
  }

  // Derivative of saturation wrt variables
  for (unsigned int ph = 0; ph < _num_phases; ++ph)
  {
    _dsaturation_dvar[_qp][ph][_pvar] = _fsp[ph].saturation.derivatives()[_pidx];

    for (unsigned int z = 0; z < _num_Z_vars; ++z)
      _dsaturation_dvar[_qp][ph][_Zvar[z]] = _fsp[ph].saturation.derivatives()[_Zidx[z]];
  }
  // Derivative of capillary pressure
  const Real dpc = _pc.dCapillaryPressure(_fsp[_aqueous_phase_number].saturation.value());

  // Derivative of porepressure wrt variables
  if (_dictator.isPorousFlowVariable(_gas_porepressure_varnum))
  {
    for (unsigned int ph = 0; ph < _num_phases; ++ph)
    {
      _dporepressure_dvar[_qp][ph][_pvar] = 1.0;
      if (!_nodal_material)
        (*_dgradp_qp_dgradv)[_qp][ph][_pvar] = 1.0;
    }

    if (!_nodal_material)
    {
      (*_dgradp_qp_dgradv)[_qp][_aqueous_phase_number][_pvar] +=
          -dpc * _dsaturation_dvar[_qp][_aqueous_phase_number][_pvar];

      for (unsigned int z = 0; z < _num_Z_vars; ++z)
        (*_dgradp_qp_dgradv)[_qp][_aqueous_phase_number][_Zvar[z]] =
            -dpc * _dsaturation_dvar[_qp][_aqueous_phase_number][_Zvar[z]];
    }

    // The aqueous phase porepressure is also a function of liquid saturation,
    // which depends on both gas porepressure and Z
    _dporepressure_dvar[_qp][_aqueous_phase_number][_pvar] +=
        -dpc * _dsaturation_dvar[_qp][_aqueous_phase_number][_pvar];

    for (unsigned int z = 0; z < _num_Z_vars; ++z)
      _dporepressure_dvar[_qp][_aqueous_phase_number][_Zvar[z]] =
          -dpc * _dsaturation_dvar[_qp][_aqueous_phase_number][_Zvar[z]];
  }

  // Calculate derivatives of material properties wrt primary variables
  // Derivative of Z wrt variables
  std::vector<Real> dZ_dvar;
  dZ_dvar.assign(_num_pf_vars, 0.0);

  for (unsigned int z = 0; z < _num_Z_vars; ++z)
    if (_dictator.isPorousFlowVariable(_Z_varnum[z]))
      dZ_dvar[_Zvar[z]] = 1.0;

  // Derivatives of properties wrt primary variables
  for (unsigned int v = 0; v < _num_pf_vars; ++v)
    for (unsigned int ph = 0; ph < _num_phases; ++ph)
    {
      // Derivative of density in each phase
      _dfluid_density_dvar[_qp][ph][v] =
          _fsp[ph].density.derivatives()[_pidx] * _dporepressure_dvar[_qp][ph][v];
      _dfluid_density_dvar[_qp][ph][v] +=
          _fsp[ph].density.derivatives()[_Tidx] * _dtemperature_dvar[_qp][v];

      // Derivative of viscosity in each phase
      _dfluid_viscosity_dvar[_qp][ph][v] =
          _fsp[ph].viscosity.derivatives()[_pidx] * _dporepressure_dvar[_qp][ph][v];
      _dfluid_viscosity_dvar[_qp][ph][v] +=
          _fsp[ph].viscosity.derivatives()[_Tidx] * _dtemperature_dvar[_qp][v];

      // Derivative of enthalpy in each phase
      _dfluid_enthalpy_dvar[_qp][ph][v] =
          _fsp[ph].enthalpy.derivatives()[_pidx] * _dporepressure_dvar[_qp][ph][v];
      _dfluid_enthalpy_dvar[_qp][ph][v] +=
          _fsp[ph].enthalpy.derivatives()[_Tidx] * _dtemperature_dvar[_qp][v];
    }

  // Add derivatives wrt Z's
  for (unsigned int ph = 0; ph < _num_phases; ++ph)
    for (unsigned int z = 0; z < _num_Z_vars; ++z)
    {
      _dfluid_density_dvar[_qp][ph][_Zvar[z]] +=
          _fsp[ph].density.derivatives()[_Zidx[z]] * dZ_dvar[_Zvar[z]];
      _dfluid_viscosity_dvar[_qp][ph][_Zvar[z]] +=
          _fsp[ph].viscosity.derivatives()[_Zidx[z]] * dZ_dvar[_Zvar[z]];
      _dfluid_enthalpy_dvar[_qp][ph][_Zvar[z]] +=
          _fsp[ph].enthalpy.derivatives()[_Zidx[z]] * dZ_dvar[_Zvar[z]];
    }

  // The derivative of the mass fractions for each fluid component in each phase.
  // Note: these are all calculated in terms of gas pressuse, so there is no
  // capillary pressure effect, and hence no need to multiply by _dporepressure_dvar
  for (unsigned int ph = 0; ph < _num_phases; ++ph)
    for (unsigned int comp = 0; comp < _num_components; ++comp)
    {
      _dmass_frac_dvar[_qp][ph][comp][_pvar] = _fsp[ph].mass_fraction[comp].derivatives()[_pidx];

      for (unsigned int z = 0; z < _num_Z_vars; ++z)
        _dmass_frac_dvar[_qp][ph][comp][_Zvar[z]] =
            _fsp[ph].mass_fraction[comp].derivatives()[_Zidx[z]] * dZ_dvar[_Zvar[z]];
    }

  // If the material properties are being evaluated at the qps, calculate the
  // gradients as well. Note: only nodal properties are evaluated in
  // initQpStatefulProperties(), so no need to check _is_initqp flag for qp
  // properties
  if (!_nodal_material)
  {
    // Second derivative of capillary pressure
    Real d2pc = _pc.d2CapillaryPressure(_fsp[_aqueous_phase_number].saturation.value(), _qp);

    (*_grads_qp)[_qp][_aqueous_phase_number] =
        _dsaturation_dvar[_qp][_aqueous_phase_number][_pvar] * _gas_gradp_qp[_qp];

    for (unsigned int z = 0; z < _num_Z_vars; ++z)
      (*_grads_qp)[_qp][_aqueous_phase_number] +=
          _dsaturation_dvar[_qp][_aqueous_phase_number][_Zvar[z]] * (*_gradZ_qp[z])[_qp];

    (*_grads_qp)[_qp][_gas_phase_number] = -(*_grads_qp)[_qp][_aqueous_phase_number];

    (*_gradp_qp)[_qp][_gas_phase_number] = _gas_gradp_qp[_qp];
    (*_gradp_qp)[_qp][_aqueous_phase_number] =
        _gas_gradp_qp[_qp] - dpc * (*_grads_qp)[_qp][_aqueous_phase_number];

    (*_dgradp_qp_dv)[_qp][_aqueous_phase_number][_pvar] =
        -d2pc * (*_grads_qp)[_qp][_aqueous_phase_number] *
        _dsaturation_dvar[_qp][_aqueous_phase_number][_pvar];

    for (unsigned int z = 0; z < _num_Z_vars; ++z)
      (*_dgradp_qp_dv)[_qp][_aqueous_phase_number][_Zvar[z]] =
          -d2pc * (*_grads_qp)[_qp][_aqueous_phase_number] *
          _dsaturation_dvar[_qp][_aqueous_phase_number][_Zvar[z]];

    (*_grad_mass_frac_qp)[_qp][_aqueous_phase_number][_aqueous_fluid_component] =
        _fsp[_aqueous_phase_number].mass_fraction[_aqueous_fluid_component].derivatives()[_pidx] *
        _gas_gradp_qp[_qp];

    for (unsigned int z = 0; z < _num_Z_vars; ++z)
      (*_grad_mass_frac_qp)[_qp][_aqueous_phase_number][_aqueous_fluid_component] +=
          _fsp[_aqueous_phase_number]
              .mass_fraction[_aqueous_fluid_component]
              .derivatives()[_Zidx[z]] *
          (*_gradZ_qp[z])[_qp];

    (*_grad_mass_frac_qp)[_qp][_aqueous_phase_number][_gas_fluid_component] =
        -(*_grad_mass_frac_qp)[_qp][_aqueous_phase_number][_aqueous_fluid_component];

    (*_grad_mass_frac_qp)[_qp][_gas_phase_number][_aqueous_fluid_component] =
        _fsp[_gas_phase_number].mass_fraction[_aqueous_fluid_component].derivatives()[_pidx] *
        _gas_gradp_qp[_qp];

    for (unsigned int z = 0; z < _num_Z_vars; ++z)
      (*_grad_mass_frac_qp)[_qp][_gas_phase_number][_aqueous_fluid_component] +=
          _fsp[_gas_phase_number].mass_fraction[_aqueous_fluid_component].derivatives()[_Zidx[z]] *
          (*_gradZ_qp[z])[_qp];

    (*_grad_mass_frac_qp)[_qp][_gas_phase_number][_gas_fluid_component] =
        -(*_grad_mass_frac_qp)[_qp][_gas_phase_number][_aqueous_fluid_component];
  }

  // If Xnacl is a PorousFlow variable, add contribution to material properties
  // due to the Xnacl component that are not included in
  // PorousFlowFluidStateFlashBase::computeQpProperties();
  if (_dictator.isPorousFlowVariable(_Xnacl_varnum))
  {
    // Derivative of saturation wrt variables
    for (unsigned int ph = 0; ph < _num_phases; ++ph)
      _dsaturation_dvar[_qp][ph][_Xvar] = _fsp[ph].saturation.derivatives()[_Xidx];

    // Derivative of porepressure wrt variables
    if (_dictator.isPorousFlowVariable(_gas_porepressure_varnum))
    {
      if (!_nodal_material)
        (*_dgradp_qp_dgradv)[_qp][_aqueous_phase_number][_Xvar] =
            -dpc * _dsaturation_dvar[_qp][_aqueous_phase_number][_Xvar];

      // The aqueous phase porepressure is also a function of liquid saturation,
      // which depends on Xnacl
      _dporepressure_dvar[_qp][_aqueous_phase_number][_Xvar] =
          -dpc * _dsaturation_dvar[_qp][_aqueous_phase_number][_Xvar];
    }

    // Calculate derivatives of material properties wrt primary variables
    // Derivative of Xnacl wrt variables
    std::vector<Real> dX_dvar;
    dX_dvar.assign(_num_pf_vars, 0.0);

    if (_dictator.isPorousFlowVariable(_Xvar))
      dX_dvar[_Xvar] = 1.0;

    // Derivatives of properties wrt primary variables
    for (unsigned int v = 0; v < _num_pf_vars; ++v)
      for (unsigned int ph = 0; ph < _num_phases; ++ph)
      {
        // Derivative of density in each phase
        _dfluid_density_dvar[_qp][ph][v] += _fsp[ph].density.derivatives()[_Xidx] * dX_dvar[v];

        // Derivative of viscosity in each phase
        _dfluid_viscosity_dvar[_qp][ph][v] += _fsp[ph].viscosity.derivatives()[_Xidx] * dX_dvar[v];
      }

    // The derivative of the mass fractions for each fluid component in each phase.
    // Note: these are all calculated in terms of gas pressuse, so there is no
    // capillary pressure effect, and hence no need to multiply by _dporepressure_dvar
    for (unsigned int ph = 0; ph < _num_phases; ++ph)
      for (unsigned int comp = 0; comp < _num_components; ++comp)
        _dmass_frac_dvar[_qp][ph][comp][_Xvar] = _fsp[ph].mass_fraction[comp].derivatives()[_Xidx];

    // If the material properties are being evaluated at the qps, add the contribution
    // to the gradients as well. Note: only nodal properties are evaluated in
    // initQpStatefulProperties(), so no need to check _is_initqp flag for qp
    // properties
    if (!_nodal_material)
    {
      // Second derivative of capillary pressure
      Real d2pc = _pc.d2CapillaryPressure(_fsp[_aqueous_phase_number].saturation.value(), _qp);

      (*_grads_qp)[_qp][_aqueous_phase_number] +=
          _dsaturation_dvar[_qp][_aqueous_phase_number][_Xvar] * _grad_Xnacl_qp[_qp];

      (*_grads_qp)[_qp][_gas_phase_number] -=
          _dsaturation_dvar[_qp][_aqueous_phase_number][_Xvar] * _grad_Xnacl_qp[_qp];

      (*_gradp_qp)[_qp][_aqueous_phase_number] -=
          dpc * _dsaturation_dvar[_qp][_aqueous_phase_number][_Xvar] * _grad_Xnacl_qp[_qp];

      (*_dgradp_qp_dv)[_qp][_aqueous_phase_number][_Xvar] =
          -d2pc * (*_grads_qp)[_qp][_aqueous_phase_number] *
          _dsaturation_dvar[_qp][_aqueous_phase_number][_Xvar];

      (*_grad_mass_frac_qp)[_qp][_aqueous_phase_number][_salt_component] = _grad_Xnacl_qp[_qp];
      (*_grad_mass_frac_qp)[_qp][_aqueous_phase_number][_aqueous_fluid_component] +=
          _fsp[_aqueous_phase_number].mass_fraction[_aqueous_fluid_component].derivatives()[_Xidx] *
          _grad_Xnacl_qp[_qp];
      (*_grad_mass_frac_qp)[_qp][_aqueous_phase_number][_gas_fluid_component] -=
          _fsp[_aqueous_phase_number].mass_fraction[_aqueous_fluid_component].derivatives()[_Xidx] *
          _grad_Xnacl_qp[_qp];
      (*_grad_mass_frac_qp)[_qp][_gas_phase_number][_aqueous_fluid_component] +=
          _fsp[_gas_phase_number].mass_fraction[_aqueous_fluid_component].derivatives()[_Xidx] *
          _grad_Xnacl_qp[_qp];
      (*_grad_mass_frac_qp)[_qp][_gas_phase_number][_gas_fluid_component] -=
          _fsp[_gas_phase_number].mass_fraction[_aqueous_fluid_component].derivatives()[_Xidx] *
          _grad_Xnacl_qp[_qp];
    }
  }
}

void
PorousFlowFluidState::setMaterialVectorSize() const
{
  _fluid_density[_qp].assign(_num_phases, 0.0);
  _fluid_viscosity[_qp].assign(_num_phases, 0.0);
  _fluid_enthalpy[_qp].assign(_num_phases, 0.0);
  _mass_frac[_qp].resize(_num_phases);

  for (unsigned int ph = 0; ph < _num_phases; ++ph)
    _mass_frac[_qp][ph].resize(_num_components);

  // Derivatives and gradients are not required in initQpStatefulProperties
  if (!_is_initqp)
  {
    _dfluid_density_dvar[_qp].resize(_num_phases);
    _dfluid_viscosity_dvar[_qp].resize(_num_phases);
    _dfluid_enthalpy_dvar[_qp].resize(_num_phases);
    _dmass_frac_dvar[_qp].resize(_num_phases);

    if (!_nodal_material)
      (*_grad_mass_frac_qp)[_qp].resize(_num_phases);

    for (unsigned int ph = 0; ph < _num_phases; ++ph)
    {
      _dfluid_density_dvar[_qp][ph].assign(_num_pf_vars, 0.0);
      _dfluid_viscosity_dvar[_qp][ph].assign(_num_pf_vars, 0.0);
      _dfluid_enthalpy_dvar[_qp][ph].assign(_num_pf_vars, 0.0);
      _dmass_frac_dvar[_qp][ph].resize(_num_components);

      for (unsigned int comp = 0; comp < _num_components; ++comp)
        _dmass_frac_dvar[_qp][ph][comp].assign(_num_pf_vars, 0.0);

      if (!_nodal_material)
        (*_grad_mass_frac_qp)[_qp][ph].assign(_num_components, RealGradient());
    }
  }
}
