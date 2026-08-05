# PorousFlowLineEnthalpySourceFromPostprocessor

`PorousFlowLineEnthalpySourceFromPostprocessor` adds the heat energy carried by fluid that is
injected along a line at a total mass flux rate computed by a postprocessor, at a temperature also
computed by a postprocessor.

It is the line analogue of
[PorousFlowPointEnthalpySourceFromPostprocessor](PorousFlowPointEnthalpySourceFromPostprocessor.md),
and should be used together with
[PorousFlowLineSourceFromPostprocessor](PorousFlowLineSourceFromPostprocessor.md) sharing the same
`mass_flux` postprocessor and the same line, so that the heat added corresponds to the fluid added.

The fluid property object used by this Dirac kernel should be the same one used in the computational
domain where the object is located.

Parameter [!param](/DiracKernels/PorousFlowLineEnthalpySourceFromPostprocessor/pressure), along with
[!param](/DiracKernels/PorousFlowLineEnthalpySourceFromPostprocessor/T_in), is used to calculate the
enthalpy of the injected fluid.  Most frequently it is the PorousFlow porepressure variable, which
models fluid injected at a specified rate and temperature.  Alternately it may be the pressure of
the injected fluid as fixed by an external agent such as a pump, stored in an AuxVariable.

The enthalpy is evaluated separately at each completion, at the pressure prevailing there.  A
completion interval that spans an appreciable pressure gradient is therefore treated correctly,
which a single point source at one end of the interval could not do.  Unlike the injected mass, the
energy added consequently does depend on how the line is resolved.

For instance:

!listing modules/porous_flow/test/tests/dirackernels/line_source_thermal.i block=DiracKernels

Note that the `execute_on` parameter of the postprocessors is set to `timestep_begin` so that the
correct values are used within the timestep.

!listing modules/porous_flow/test/tests/dirackernels/line_source_thermal.i block=Postprocessors

!syntax parameters /DiracKernels/PorousFlowLineEnthalpySourceFromPostprocessor

!syntax inputs /DiracKernels/PorousFlowLineEnthalpySourceFromPostprocessor

!syntax children /DiracKernels/PorousFlowLineEnthalpySourceFromPostprocessor
