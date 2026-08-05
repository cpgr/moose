//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PorousFlowLineEnthalpySourceFromPostprocessor.h"
#include "SinglePhaseFluidProperties.h"

registerMooseObject("PorousFlowApp", PorousFlowLineEnthalpySourceFromPostprocessor);

InputParameters
PorousFlowLineEnthalpySourceFromPostprocessor::validParams()
{
  InputParameters params = PorousFlowLineGeometry::validParams();
  params.addRequiredParam<PostprocessorName>(
      "mass_flux",
      "The postprocessor name holding the TOTAL mass flux of injected fluid over the whole line in "
      "kg/s (please ensure this is positive so that this object acts like a source)");
  params.addRequiredParam<UserObjectName>(
      "fp",
      "The name of the user object used to calculate the fluid properties of the injected fluid");
  params.addRequiredCoupledVar(
      "pressure", "Pressure used to calculate the injected fluid enthalpy (measured in Pa)");
  params.addRequiredParam<PostprocessorName>(
      "T_in", "The postprocessor name holding the temperature of injected fluid (measured in K)");
  params.addClassDescription(
      "Line source that adds heat energy corresponding to injection of a fluid at a total mass "
      "flux rate (specified by a postprocessor) and given temperature (specified by a "
      "postprocessor).  The total is divided between the completions in proportion to the length "
      "of line each represents, and the enthalpy of each is evaluated at the pressure there, so "
      "that a completion interval spanning a pressure gradient is treated correctly.");
  return params;
}

PorousFlowLineEnthalpySourceFromPostprocessor::PorousFlowLineEnthalpySourceFromPostprocessor(
    const InputParameters & parameters)
  : PorousFlowLineGeometry(parameters),
    _mass_flux(getPostprocessorValue("mass_flux")),
    _pressure(coupledValue("pressure")),
    _T_in(getPostprocessorValue("T_in")),
    _fp(getUserObject<SinglePhaseFluidProperties>("fp")),
    _p_var_num(coupled("pressure"))
{
}

Real
PorousFlowLineEnthalpySourceFromPostprocessor::computeQpResidual()
{
  // Negative sign to make a positive mass_flux in the input file a source
  const Real h = _fp.h_from_p_T(_pressure[_qp], _T_in);
  return -_test[_i][_qp] * _mass_flux * completionFraction(currentPointCachedID()) * h;
}

Real
PorousFlowLineEnthalpySourceFromPostprocessor::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _p_var_num)
  {
    Real h, dh_dp, dh_dT;
    _fp.h_from_p_T(_pressure[_qp], _T_in, h, dh_dp, dh_dT);
    return -_test[_i][_qp] * _phi[_j][_qp] * _mass_flux *
           completionFraction(currentPointCachedID()) * dh_dp;
  }
  else
    return 0.0;
}
