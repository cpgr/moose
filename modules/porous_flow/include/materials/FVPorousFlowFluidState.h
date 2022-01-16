/*****************************************************************/
/*    FINCH - FINite volume Capillary Heterogeneity modelling    */
/*                                                               */
/*           All contents are licensed under MIT/BSD             */
/*              See LICENSE for full restrictions                */
/*****************************************************************/

#pragma once

#include "Material.h"
#include "PorousFlowBrineCO2.h"

class PorousFlowDictator;

/**
 * Fluid state class using a persistent set of primary variables for
 * the mutliphase, multicomponent case.
 *
 * Primary variables are: gas pressure, total mass fraction
 * of a component summed over all phases (and optionally temperature in a
 * non-isothermal case).
 *
 * The total mass fraction of component i summed over all phases, Z_i,
 * is defined as (for two phases)
 *
 * Z_i = (S_g rho_g Y_i + S_l rho_l X_i) / (S_g rho_g + S_l rho_l)
 *
 * where S is saturation, rho is density, and the subscripts correspond to gas
 * and liquid phases, respectively, and Y_i and X_i are the mass fractions of
 * the ith component in the gas and liquid phase, respectively.
 *
 * Depending on the phase conditions, the primary variable Z_i can represent either
 * a mass fraction (when only a single phase is present), or a saturation when
 * two phases are present, and hence it is a persistent variable.
 *
 * The PorousFlow kernels expect saturation and mass fractions (as well as pressure
 * and temperature), so these must be calculated from Z_i once the state of the
 * system is determined.
 *
 * A compositional flash calculation using the Rachford-Rice equation is solved
 * to determine vapor fraction (gas saturation), and subsequently the composition
 * of each phase.
 */
class FVPorousFlowFluidState : public Material
{
public:
  static InputParameters validParams();

  FVPorousFlowFluidState(const InputParameters & parameters);

protected:
  virtual void initialSetup() override;
  virtual void initQpStatefulProperties() override;
  virtual void computeQpProperties() override;

  /// Size material property vectors and initialise with zeros
  void setMaterialVectorSize() const;

  /**
   * Calculates all required thermophysical properties and derivatives for each phase
   * and fluid component
   */
  virtual void thermophysicalProperties();

  const PorousFlowDictator & _dictator;
  const unsigned int _num_phases;
  const unsigned int _num_components;

  /// Porepressure
  const ADVariableValue & _gas_porepressure;
  /// Moose variable number of the gas porepressure
  const unsigned int _gas_porepressure_varnum;
  /// PorousFlow variable number of the gas porepressure
  const unsigned int _pvar;
  /// Total mass fraction(s) of the gas component(s) summed over all phases
  std::vector<const ADVariableValue *> _Z;
  /// Moose variable number of Z
  std::vector<unsigned int> _Z_varnum;
  /// PorousFlow variable number of Z
  std::vector<unsigned int> _Zvar;
  /// Number of coupled total mass fractions. Should be _num_phases - 1
  const unsigned int _num_Z_vars;
  /// Salt mass fraction (kg/kg)
  const ADVariableValue & _Xnacl;
  /// Salt mass fraction variable number
  const unsigned int _Xnacl_varnum;
  /// Salt mass fraction PorousFlow variable number
  const unsigned int _Xvar;
  /// Temperature variable
  const ADVariableValue & _temperature;
  /// Moose variable number of the temperature
  const unsigned int _temperature_varnum;
  /// PorousFlow variable number of the temperature
  const unsigned int _Tvar;
  /// FluidState UserObject
  const PorousFlowBrineCO2 & _fs;
  /// Phase number of the aqueous phase
  const unsigned int _aqueous_phase_number;
  /// Phase number of the gas phase
  const unsigned int _gas_phase_number;
  /// Fluid component number of the aqueous component
  const unsigned int _aqueous_fluid_component;
  /// Fluid component number of the gas phase
  const unsigned int _gas_fluid_component;
  /// Salt component index
  const unsigned int _salt_component;

  ADMaterialProperty<std::vector<Real>> & _pressure;
  ADMaterialProperty<std::vector<Real>> & _saturation;
  /// Mass fraction matrix
  ADMaterialProperty<std::vector<std::vector<Real>>> & _mass_frac;
  /// Fluid density of each phase
  ADMaterialProperty<std::vector<Real>> & _fluid_density;
  /// Viscosity of each phase
  ADMaterialProperty<std::vector<Real>> & _fluid_viscosity;
  /// Enthalpy of each phase
  ADMaterialProperty<std::vector<Real>> & _fluid_enthalpy;
  /// Internal energy of each phase
  ADMaterialProperty<std::vector<Real>> & _fluid_internal_energy;

  /// Conversion from degrees Celsius to degrees Kelvin
  const Real _T_c2k;
  /// Flag to indicate whether to calculate stateful properties
  bool _is_initqp;
  /// FluidStateProperties data structure
  std::vector<FluidStateProperties> _fsp;
  /// Index of derivative wrt pressure
  const unsigned int _pidx;
  /// Index of derivative wrt temperature
  const unsigned int _Tidx;
  /// Index of derivative wrt total mass fraction Z
  const unsigned int _Zidx;
  /// Index of derivative wrt salt mass fraction X
  const unsigned int _Xidx;
};
