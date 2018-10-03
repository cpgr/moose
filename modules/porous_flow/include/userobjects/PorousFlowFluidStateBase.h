//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef POROUSFLOWFLUIDSTATEBASE_H
#define POROUSFLOWFLUIDSTATEBASE_H

#include "PorousFlowFluidStateFlash.h"
#include "PorousFlowCapillaryPressure.h"

class PorousFlowFluidStateBase;

/// Data structure to pass calculated thermophysical properties
struct FluidStateProperties
{
  FluidStateProperties(){};
  FluidStateProperties(unsigned int nc, unsigned int nz = 1)
    : pressure(0.0),
      saturation(0.0),
      density(0.0),
      viscosity(1.0), // to guard against division by zero
      enthalpy(0.0),
      mass_fraction(nc, 0.0),
      dsaturation_dp(0.0),
      dsaturation_dT(0.0),
      // dsaturation_dZ(0.0),
      dsaturation_dZ(nz, 0.0),
      dsaturation_dX(0.0),
      ddensity_dp(0.0),
      ddensity_dT(0.0),
      // ddensity_dZ(0.0),
      ddensity_dZ(nz, 0.0),
      ddensity_dX(0.0),
      dviscosity_dp(0.0),
      dviscosity_dT(0.0),
      // dviscosity_dZ(0.0),
      dviscosity_dZ(nz, 0.0),
      dviscosity_dX(0.0),
      denthalpy_dp(0.0),
      denthalpy_dT(0.0),
      // denthalpy_dZ(0.0),
      denthalpy_dZ(nz, 0.0),
      denthalpy_dX(0.0),
      dmass_fraction_dp(nc, 0.0),
      dmass_fraction_dT(nc, 0.0),
      // dmass_fraction_dZ(nc, 0.0),
      dmass_fraction_dZ(nc, std::vector<Real>(nz, 0.0)),
      dmass_fraction_dX(nc, 0.0){};

  Real pressure;
  Real saturation;
  Real density;
  Real viscosity;
  Real enthalpy;
  std::vector<Real> mass_fraction;
  Real dsaturation_dp;
  Real dsaturation_dT;
  // Real dsaturation_dZ;
  std::vector<Real> dsaturation_dZ;
  Real dsaturation_dX;
  Real ddensity_dp;
  Real ddensity_dT;
  // Real ddensity_dZ;
  std::vector<Real> ddensity_dZ;
  Real ddensity_dX;
  Real dviscosity_dp;
  Real dviscosity_dT;
  // Real dviscosity_dZ;
  std::vector<Real> dviscosity_dZ;
  Real dviscosity_dX;
  Real denthalpy_dp;
  Real denthalpy_dT;
  // Real denthalpy_dZ;
  std::vector<Real> denthalpy_dZ;
  Real denthalpy_dX;
  std::vector<Real> dmass_fraction_dp;
  std::vector<Real> dmass_fraction_dT;
  // std::vector<Real> dmass_fraction_dZ;
  std::vector<std::vector<Real>> dmass_fraction_dZ;
  std::vector<Real> dmass_fraction_dX;
};

template <>
InputParameters validParams<PorousFlowFluidStateBase>();

/**
 * Base class for fluid states for miscible multiphase flow in porous media.
 */
class PorousFlowFluidStateBase : public PorousFlowFluidStateFlash
{
public:
  PorousFlowFluidStateBase(const InputParameters & parameters);

  /**
   * Name of FluidState
   * @return name of fluid state class
   */
  virtual std::string fluidStateName() const = 0;

  /**
   * Determines the complete thermophysical state of the system for a given set of
   * primary variables
   *
   * @param pressure gas phase pressure (Pa)
   * @param temperature fluid temperature (K)
   * @param Xnacl mass fraction of NaCl
   * @param Z total mass fraction gas component(s)
   * @param qp quadpoint
   * @param[out] fsp the FluidStateProperties struct containing all properties
   */
  virtual void thermophysicalProperties(Real pressure,
                                        Real temperature,
                                        Real Xnacl,
                                        std::vector<const VariableValue *> & Z,
                                        unsigned int qp,
                                        std::vector<FluidStateProperties> & fsp) const = 0;
  /**
   * Total mass fraction of gas component summed over all phases in the two-phase state
   *
   * @param pressure gas pressure (Pa)
   * @param temperature temperature (K)
   * @param Xnacl NaCl mass fraction (kg/kg)
   * @param saturation gas saturation (-)
   * @param qp quadpoint, which may be used by the capillary-pressure userobject
   * @return total mass fraction Z (-)
   */
  virtual Real totalMassFraction(
      Real pressure, Real temperature, Real Xnacl, Real saturation, unsigned int qp) const = 0;

  /**
   * The maximum number of phases in this model
   * @return number of phases
   */
  unsigned int numPhases() const { return _num_phases; };

  /**
   * The maximum number of components in this model
   * @return number of components
   */
  unsigned int numComponents() const { return _num_components; };

  /**
   * The index of the aqueous phase
   * @return aqueous phase number
   */
  unsigned int aqueousPhaseIndex() const { return _aqueous_phase_number; };

  /**
   * The index of the gas phase
   * @return gas phase number
   */
  unsigned int gasPhaseIndex() const { return _gas_phase_number; };

  /**
   * The index of the aqueous fluid component
   * @return aqueous fluid component number
   */
  unsigned int aqueousComponentIndex() const { return _aqueous_fluid_component; };

  /**
   * The index of the gas fluid component
   * @return gas fluid component number
   */
  unsigned int gasComponentIndex() const { return _gas_fluid_component; };

  /**
   * The index of the salt component
   * @return salt component number
   */
  unsigned int saltComponentIndex() const { return _salt_component; };

  /**
   * Clears the contents of the FluidStateProperties data structure
   * @param[out] fsp FluidStateProperties data structure with all data initialized to 0
   */
  void clearFluidStateProperties(std::vector<FluidStateProperties> & fsp) const;

protected:
  /// Number of phases
  unsigned int _num_phases;
  /// Total number of components
  unsigned int _num_components;
  /// Number of nonlinear variables Zi
  unsigned int _num_zvars;
  /// Phase number of the aqueous phase
  const unsigned int _aqueous_phase_number;
  /// Phase number of the gas phase
  unsigned int _gas_phase_number;
  /// Fluid component number of the aqueous component
  const unsigned int _aqueous_fluid_component;
  /// Fluid component number of the gas phase
  unsigned int _gas_fluid_component;
  /// Salt component index
  unsigned int _salt_component;
  /// Universal gas constant (J/mol/K)
  const Real _R;
  /// Conversion from C to K
  const Real _T_c2k;
  /// Maximum number of iterations for the Newton-Raphson iterations
  const Real _nr_max_its;
  /// Tolerance for Newton-Raphson iterations
  const Real _nr_tol;
  /// Capillary pressure UserObject
  const PorousFlowCapillaryPressure & _pc_uo;
};

#endif // POROUSFLOWFLUIDSTATEBASE_H
