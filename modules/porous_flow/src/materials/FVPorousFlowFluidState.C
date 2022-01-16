//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FVPorousFlowFluidState.h"
#include "PorousFlowDictator.h"

registerMooseObject("PorousFlowApp", FVPorousFlowFluidState);

InputParameters
FVPorousFlowFluidState::validParams()
{
  InputParameters params = Material::validParams();
  params.addRequiredParam<UserObjectName>(
      "PorousFlowDictator", "The UserObject that holds the list of PorousFlow variable names");
  params.addRequiredCoupledVar("gas_porepressure",
                               "Variable that is the porepressure of the gas phase");
  params.addRequiredCoupledVar("z", "Total mass fraction of component i summed over all phases");
  params.addCoupledVar(
      "temperature", 293, "The fluid temperature (C or K, depending on temperature_unit)");
  params.addCoupledVar("xnacl", 0, "The salt mass fraction in the brine (kg/kg)");
  MooseEnum unit_choice("Kelvin=0 Celsius=1", "Kelvin");
  params.addParam<MooseEnum>(
      "temperature_unit", unit_choice, "The unit of the temperature variable");
  params.addRequiredParam<UserObjectName>("fluid_state", "Name of the FluidState UserObject");
  params.addClassDescription("Class for fluid state calculations using persistent primary "
                             "variables and a vapor-liquid flash");
  return params;
}

FVPorousFlowFluidState::FVPorousFlowFluidState(const InputParameters & parameters)
  : Material(parameters),

    _dictator(getUserObject<PorousFlowDictator>("PorousFlowDictator")),
    _num_phases(_dictator.numPhases()),
    _num_components(_dictator.numComponents()),

    _gas_porepressure(adCoupledValue("gas_porepressure")),
    _gas_porepressure_varnum(coupled("gas_porepressure")),
    _pvar(_dictator.isPorousFlowVariable(_gas_porepressure_varnum)
              ? _dictator.porousFlowVariableNum(_gas_porepressure_varnum)
              : 0),

    _num_Z_vars(coupledComponents("z")),
    _Xnacl(adCoupledValue("xnacl")),
    _Xnacl_varnum(coupled("xnacl")),
    _Xvar(_dictator.isPorousFlowVariable(_Xnacl_varnum)
              ? _dictator.porousFlowVariableNum(_Xnacl_varnum)
              : 0),

    _temperature(adCoupledValue("temperature")),
    _temperature_varnum(coupled("temperature")),
    _Tvar(_dictator.isPorousFlowVariable(_temperature_varnum)
              ? _dictator.porousFlowVariableNum(_temperature_varnum)
              : 0),

    _fs(getUserObject<PorousFlowBrineCO2>("fluid_state")),
    _aqueous_phase_number(_fs.aqueousPhaseIndex()),
    _gas_phase_number(_fs.gasPhaseIndex()),
    _aqueous_fluid_component(_fs.aqueousComponentIndex()),
    _gas_fluid_component(_fs.gasComponentIndex()),
    _salt_component(_fs.saltComponentIndex()),

    _pressure(declareADProperty<std::vector<Real>>("pressure")),
    _saturation(declareADProperty<std::vector<Real>>("saturation")),
    _mass_frac(declareADProperty<std::vector<std::vector<Real>>>("mass_fractions")),
    _fluid_density(declareADProperty<std::vector<Real>>("density")),
    _fluid_viscosity(declareADProperty<std::vector<Real>>("viscosity")),
    _fluid_enthalpy(declareADProperty<std::vector<Real>>("enthalpy")),
    _fluid_internal_energy(declareADProperty<std::vector<Real>>("internal_energy")),

    _T_c2k(getParam<MooseEnum>("temperature_unit") == 0 ? 0.0 : 273.15),
    _is_initqp(false),
    _pidx(_fs.getPressureIndex()),
    _Tidx(_fs.getTemperatureIndex()),
    _Zidx(_fs.getZIndex()),
    _Xidx(_fs.getXIndex())
{
  // Check that the number of phases in the fluidstate class is also provided in the Dictator
  if (_fs.numPhases() != _num_phases)
    mooseError(name(),
               ": only ",
               _fs.numPhases(),
               " phases are allowed. Please check the number of phases entered in the dictator is "
               "correct");

  // Store all total mass fractions and associated variable numbers
  _Z.resize(_num_Z_vars);
  _Z_varnum.resize(_num_Z_vars);
  _Zvar.resize(_num_Z_vars);

  for (unsigned int i = 0; i < _num_Z_vars; ++i)
  {
    _Z[i] = &adCoupledValue("z", i);
    _Z_varnum[i] = coupled("z", i);
    _Zvar[i] = (_dictator.isPorousFlowVariable(_Z_varnum[i])
                    ? _dictator.porousFlowVariableNum(_Z_varnum[i])
                    : 0);
  }

  // Set the size of the FluidStateProperties vector
  _fsp.resize(_num_phases, FluidStateProperties(_num_components));
}

void
FVPorousFlowFluidState::thermophysicalProperties()
{
  // The FluidProperty objects use temperature in K
  ADReal Tk = _temperature[_qp] + _T_c2k;

  _fs.thermophysicalProperties(_gas_porepressure[_qp], Tk, _Xnacl[_qp], (*_Z[0])[_qp], _qp, _fsp);
}

void
FVPorousFlowFluidState::initQpStatefulProperties()
{
  // Set the size of all other vectors
  setMaterialVectorSize();

  // Calculate all required thermophysical properties
  thermophysicalProperties();

  // Set the initial values of the properties
  for (unsigned int ph = 0; ph < _num_phases; ++ph)
  {
    _saturation[_qp][ph] = _fsp[ph].saturation;
    _pressure[_qp][ph] = _fsp[ph].pressure;
    _fluid_density[_qp][ph] = _fsp[ph].density;
    _fluid_viscosity[_qp][ph] = _fsp[ph].viscosity;
    _fluid_enthalpy[_qp][ph] = _fsp[ph].enthalpy;
    _fluid_internal_energy[_qp][ph] = _fsp[ph].internal_energy;

    for (unsigned int comp = 0; comp < _num_components; ++comp)
      _mass_frac[_qp][ph][comp] = _fsp[ph].mass_fraction[comp];
  }
}

void
FVPorousFlowFluidState::computeQpProperties()
{
  // Set the size of all other vectors
  setMaterialVectorSize();

  // Calculate all required thermophysical properties
  thermophysicalProperties();

  for (unsigned int ph = 0; ph < _num_phases; ++ph)
  {
    _saturation[_qp][ph] = _fsp[ph].saturation;
    _pressure[_qp][ph] = _fsp[ph].pressure;
    _fluid_density[_qp][ph] = _fsp[ph].density;
    _fluid_viscosity[_qp][ph] = _fsp[ph].viscosity;
    _fluid_enthalpy[_qp][ph] = _fsp[ph].enthalpy;
    _fluid_internal_energy[_qp][ph] = _fsp[ph].internal_energy;

    for (unsigned int comp = 0; comp < _num_components; ++comp)
      _mass_frac[_qp][ph][comp] = _fsp[ph].mass_fraction[comp];
  }
}

void
FVPorousFlowFluidState::setMaterialVectorSize() const
{
  _pressure[_qp].assign(_num_phases, 0.0);
  _saturation[_qp].assign(_num_phases, 0.0);
  _fluid_density[_qp].assign(_num_phases, 0.0);
  _fluid_viscosity[_qp].assign(_num_phases, 0.0);
  _fluid_enthalpy[_qp].assign(_num_phases, 0.0);
  _fluid_internal_energy[_qp].assign(_num_phases, 0.0);
  _mass_frac[_qp].resize(_num_phases);

  for (unsigned int ph = 0; ph < _num_phases; ++ph)
    _mass_frac[_qp][ph].resize(_num_components);
}
