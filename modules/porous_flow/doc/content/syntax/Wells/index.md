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

These two parameters are a map from fluid component to variable, and the order in which the
variables are listed *is* the map. They do not describe what the variables are.

In a compositional simulation the distinction does not show, because the variables carrying the
component equations are mass fractions, so `mass_fraction_vars` both maps the components and
describes its entries. Three components solved for pressure and two total mass fractions:

```
porepressure = pp
mass_fraction_vars = 'z_h2 z_ch4'
```

which puts component 0 on `z_h2`, component 1 on `z_ch4` and component 2 on `pp`.

An immiscible simulation separates the two meanings. Solving for pressure, gas saturation and the
mass fraction of hydrogen in the gas gives

```
porepressure = pp
mass_fraction_vars = 'sg y_h2'
```

which puts component 0 on `sg`, component 1 on `y_h2` and component 2 on `pp`. Only one of the two
entries is a mass fraction and the first is a saturation, so here the parameter's name describes
its entries wrongly while still mapping them correctly. What matters is the order, and it must
agree with the order the kernels were written in. Nothing checks this: the map is supplied by the
input rather than derived from the equations, so listing the variables in the wrong order silently
attaches each well's flux to the wrong component equation. Write the list by reading off which
variable each `PorousFlowMassTimeDerivative` acts on, in `fluid_component` order.

!alert note
Both parameters exist only because the PorousFlowDictator does not record which variable carries
each component's equation. Were it to do so, `[Wells]` would ask the dictator instead of the user,
the order would no longer be something to get wrong, and the naming difficulty above would not
arise. That is the intended direction, and these parameters should be expected to disappear.

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

## Wells that respect a pressure limit

A well given `bhp_max` no longer delivers its rate unconditionally.  Its bottom-hole pressure becomes
an unknown of the simulation, solved for alongside the reservoir equations, and the well delivers the
requested rate whenever it can do so beneath the limit and holds the limit when it cannot.  The rate
then declines as the reservoir pressures up, which a well of prescribed rate cannot represent.

Nothing else in the block changes: adding the one parameter is what selects the mechanism, and the
expansion is again the action's business.  It builds one borehole owning the rate constraint, one
per fluid component carrying the injected composition, and a single scalar variable holding the
pressure that all of them share, reported under `<well>_bhp`.

Such a well is modelled by a Peaceman borehole, whose well constant depends on the radius of the
hole, so its completions must be given by `point_file` or `line_base`, whose first entry per
completion is that radius.  A bare `point` carries no radius and is an error here.

The pressure at a completion is the bottom-hole pressure plus the weight of the fluid standing
between them, so a well completed over a tall interval does not see one pressure throughout.  That
gradient is `wellbore_density` times `gravity`, the latter normally set once in `GlobalParams` as the
Physics take it.  The density is the one in the hole, the injected stream rather than anything in the
formation, so it is the user's to give; left at its default of zero the wellbore carries no gradient,
which is exact for a well completed at a single point or along a horizontal interval.

## Example

!listing modules/porous_flow/test/tests/wells/injection_2comp_action.i block=Wells

The equivalent hand-built input, which this reproduces exactly, is:

!listing modules/porous_flow/test/tests/wells/injection_2comp_hand.i block=DiracKernels

!syntax list /Wells objects=True actions=False subsystems=False

!syntax list /Wells objects=False actions=False subsystems=True
