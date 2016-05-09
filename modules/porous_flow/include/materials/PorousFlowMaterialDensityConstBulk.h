/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#ifndef PORFLOWMATERIALDENSITYCONSTBULK_H
#define PORFLOWMATERIALDENSITYCONSTBULK_H

#include "PorousFlowMaterialFluidPropertiesBase.h"
#include "PorousFlowConstantBulkProperties.h"

class PorousFlowMaterialDensityConstBulk;

template<>
InputParameters validParams<PorousFlowMaterialDensityConstBulk>();

/**
 * Material designed to calculate fluid density
 * from porepressure, assuming constant bulk modulus
 * for the fluid.
 */
class PorousFlowMaterialDensityConstBulk : public PorousFlowMaterialFluidPropertiesBase
{
public:
  PorousFlowMaterialDensityConstBulk(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties();
  virtual void computeQpProperties();

  /// density at zero porepressure
  const Real _density_p0;
  /// constant bulk modulus
  const Real _bulk;
  /// Fluid phase density at the nodes
  MaterialProperty<Real> & _density_nodal;
  /// Old fluid phase density at the nodes
  MaterialProperty<Real> & _density_nodal_old;
  /// Derivative of fluid density wrt phase pore pressure at the nodes
  MaterialProperty<Real> & _ddensity_nodal_dp;
  /// Fluid phase density at the qps
  MaterialProperty<Real> & _density_qp;
  /// Derivative of fluid density wrt phase pore pressure at the qps
  MaterialProperty<Real> & _ddensity_qp_dp;
};

#endif //PORFLOWMATERIALDENSITYCONSTBULK_H
