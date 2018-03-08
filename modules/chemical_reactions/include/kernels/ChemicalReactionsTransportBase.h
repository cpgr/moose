//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef CHEMICALREACTIONSTRANSPORTBASE_H
#define CHEMICALREACTIONSTRANSPORTBASE_H

#include "Kernel.h"

class ChemicalReactionsTransportBase;

template <>
InputParameters validParams<ChemicalReactionsTransportBase>();

/**
 * Base class for transport kernels in the chemical reactions module
 */
class ChemicalReactionsTransportBase : public Kernel
{
public:
  ChemicalReactionsTransportBase(const InputParameters & parameters);

protected:
  /// Weight of the equilibrium species in the total primary species
  const Real _weight;
  /// Solvent density (assumed to be water)
  const MaterialProperty<Real> * _density;
  /// Enum of concentration units
  const enum class ConcentrationEnum { MOLARITY, MOLALITY, DEFAULT } _concentration_enum;
};

#endif // CHEMICALREACTIONSTRANSPORTBASE_H
