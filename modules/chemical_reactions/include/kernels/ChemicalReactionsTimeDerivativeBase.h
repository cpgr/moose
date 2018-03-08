//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef CHEMICALREACTIONSTIMEDERIVATIVEBASE_H
#define CHEMICALREACTIONSTIMEDERIVATIVEBASE_H

#include "TimeDerivative.h"

class ChemicalReactionsTimeDerivativeBase;

template <>
InputParameters validParams<ChemicalReactionsTimeDerivativeBase>();

/**
 * Base class for time derivative kernels in the chemical reactions module
 */
class ChemicalReactionsTimeDerivativeBase : public TimeDerivative
{
public:
  ChemicalReactionsTimeDerivativeBase(const InputParameters & parameters);

protected:
  /// Weight of the equilibrium species in the total primary species
  const Real _weight;
  /// Porosity
  const MaterialProperty<Real> & _porosity;
  /// Old value of porosity
  const MaterialProperty<Real> & _porosity_old;
  /// Old value of the primary species concentration.
  const VariableValue & _u_old;
  /// Solvent density (assumed to be water)
  const MaterialProperty<Real> * _density;
  /// Enum of concentration units
  const enum class ConcentrationEnum { MOLARITY, MOLALITY, DEFAULT } _concentration_enum;
};

#endif // CHEMICALREACTIONSTIMEDERIVATIVEBASE_H
