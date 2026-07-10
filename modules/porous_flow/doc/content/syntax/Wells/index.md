# Wells

The `[Wells]` block is a convenience layer for specifying wells in a PorousFlow simulation.  Each
well is written once, and the action expands it into the DiracKernels, UserObjects and
Postprocessors that model it, including one kernel per fluid component and a standard set of rate
and cumulative-mass reporters.

`[Wells]` is a top-level block rather than part of a `[Physics]` block: the Physics blocks declare
the PDE terms, whereas wells are operational *forcing* on those equations, closer in spirit to
`[BCs]`.  Keeping it top-level means it works identically whether the equations were built by a
PorousFlow Physics block or written out by hand, since everything the expansion needs comes from the
[PorousFlowDictator](PorousFlowDictator.md).

A well is a point or line source, so `[Wells]` suits a well in a 2D or 3D reservoir model.  Areal or
radial injection across a boundary is a boundary-condition problem: use
[PorousFlowSink](PorousFlowSink.md) for that instead.

## Fluid components

The PorousFlowDictator records how many fluid components there are, but not which variable carries
the mass balance of each one.  `[Wells]` therefore takes `porepressure` and `mass_fraction_vars`,
with exactly the meaning they have in the PorousFlow Physics: the $i$th entry of
`mass_fraction_vars` holds the equation for the $i$th fluid component, and `porepressure` holds the
equation for the final one.

## Shutting a well in: gate the flux, never the object

A well is always shut in by driving its flux to zero, never by disabling its objects with a
`Control` such as `TimePeriod`. Enabling and disabling is the one thing that cannot express
"produce zero":

- a disabled `Postprocessor` is simply not executed, so it holds the last value it computed.  An
  injector gated this way goes on reporting the rate it would have had were it open, and its
  cumulative goes on integrating that rate;
- a disabled [PorousFlowPolyLineSink](PorousFlowPolyLineSink.md) never re-runs the `addPoints()`
  that zeroes its [PorousFlowSumQuantity](PorousFlowSumQuantity.md), so
  [PorousFlowPlotQuantity](PorousFlowPlotQuantity.md) reports the mass produced during the last open
  step, forever.

In both cases the physics is right and the reporting is wrong.  Instead, an injector's rate function
is multiplied by a gate that is 1 when open and 0 when shut, and a producer's polyline sinks are
given that same gate as their `multiplying_var`.  Everything downstream then reads zero honestly.

## Example

!listing modules/porous_flow/test/tests/wells/injection_2comp_action.i block=Wells

The equivalent hand-built input, which this reproduces exactly, is:

!listing modules/porous_flow/test/tests/wells/injection_2comp_hand.i block=DiracKernels

!syntax list /Wells objects=True actions=False subsystems=False

!syntax list /Wells objects=False actions=False subsystems=True
