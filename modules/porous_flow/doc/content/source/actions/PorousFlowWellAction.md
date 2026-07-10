# PorousFlowWellAction

!syntax description /Wells/PorousFlowWellAction

This action expands each block of the [Wells](syntax/Wells/index.md) syntax into the concrete
objects that model the well.  Users do not add this action directly; it is created by the `[Wells]`
syntax.

The action reads only the [PorousFlowDictator](PorousFlowDictator.md), so it works with inputs whose
equations came from a PorousFlow Physics block and with hand-built inputs alike.  Depending on the
well type it builds:

- for a [PorousFlowInjectionWell](PorousFlowInjectionWell.md): one point source per fluid component,
  a rate postprocessor per component and, when nonisothermal, an enthalpy source;
- for a [PorousFlowProductionWell](PorousFlowProductionWell.md): one polyline sink per fluid
  component, each with its own [PorousFlowSumQuantity](PorousFlowSumQuantity.md), and the produced
  mass, rate and cumulative postprocessors.

Any parameter that requires the implicit scalar-BHP well raises an error rather than being silently
ignored.

!syntax parameters /Wells/PorousFlowWellAction
