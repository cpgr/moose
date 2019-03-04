//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef EQUILIBRIUMREACTIONS_H
#define EQUILIBRIUMREACTIONS_H

#include "MooseTypes.h"
#include "metaphysicl/numberarray.h"
#include "metaphysicl/dualnumber.h"

namespace EquilibriumReactions
{
/**
 * Equilibrium species concentration for aqueous equilibrium reaction
 * @param c_j[in] vector of primary species concentrations
 * @param gamma_j[in] vector of activity coefficients for primary species
 * @param sto_j stoichiometric coefficients for primary species
 * @param gamma_i activity coefficient of equilibrium species
 * @param K_i (logarithem of) equilibrium constant
 * @return concentration of aqueous equilibrium species
 */
Real equilibriumSpecies(std::vector<Real> & c_j,
                        std::vector<Real> & gamma_j,
                        std::vector<Real> & sto_j,
                        Real gamma_i,
                        Real log_Ki);

DualReal equilibriumSpecies(std::vector<DualReal> & c_j,
                            std::vector<DualReal> & gamma_j,
                            std::vector<DualReal> & sto_j,
                            DualReal gamma_i,
                            DualReal log_Ki);

Real equilibriumSpecies(std::vector<const VariableValue *> & c_j,
                        std::vector<const VariableValue *> & gamma_j,
                        const std::vector<Real> & sto_j,
                        const VariableValue & gamma_i,
                        const VariableValue & log_Ki,
                        unsigned int qp);

template <typename T>
T equilibriumSpeciesTempl(
    std::vector<T> & c_j, std::vector<T> & gamma_j, std::vector<T> & sto_j, T gamma_i, T log_Ki);
}

#endif // EQUILIBRIUMREACTIONS_H
