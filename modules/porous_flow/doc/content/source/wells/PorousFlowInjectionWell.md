# PorousFlowInjectionWell

!syntax description /Wells/PorousFlowInjectionWell

A rate-controlled injection well, specified inside the [Wells](syntax/Wells/index.md) block.  The
well injects a fluid of fixed composition at a single point, at a total mass rate given by the
`rate` function.

The [PorousFlowWellAction](/Wells/index.md) expands this object into:

- one [PorousFlowPointSourceFromPostprocessor](PorousFlowPointSourceFromPostprocessor.md) per fluid
  component, each carrying $q_i = w_i q_\mathrm{total}$, where $w_i$ is the $i$th entry of
  `injected_mass_fractions`.  Components with a zero mass fraction are skipped;
- when `temperature` is supplied, a
  [PorousFlowPointEnthalpySourceFromPostprocessor](PorousFlowPointEnthalpySourceFromPostprocessor.md)
  carrying the *total* mass rate at the injection temperature;
- the rate postprocessors `<well>_q_<i>` and the cumulative-mass postprocessors `<well>_cum_<i>`.

Because the rate is prescribed rather than computed from a bottom-hole pressure, this well has no
injectivity limit and cannot fail to bootstrap: it is safe from a cold start, even when injecting
gas into a fully brine-saturated reservoir.  The cost is that the well has no feedback, so it will
pressurise the near-well cell without limit if the prescribed rate exceeds what the formation can
accept.  Supplying `bhp_max` to obtain that feedback is not yet possible and raises an error.

## Sign convention

`rate` is positive in the well's natural sense: a positive rate injects.  Every reported rate and
cumulative is likewise positive during injection.

## Shutting the well in

`start_time` and `end_time` shut the well in outside the operating window, and are added to the
simulation's sync times so that the solve lands exactly on them.  The well is shut by zeroing its
rate, so its rate and cumulative reporters correctly read zero while it is shut.

`shut_if` shuts the well in whenever a supplied `Function` evaluates to a nonzero (true) value.
Because a `ParsedFunction` resolves a `Postprocessor` named in its `symbol_values`, the condition
can depend on the state of the simulation rather than only on time:

```
[too_wet_fn]
  type = ParsedFunction
  expression = 'wc > 0.95'
  symbol_names = 'wc'
  symbol_values = 'my_water_cut_pp'
[]
```

The postprocessor is read from the end of the previous timestep, so the well shuts one timestep
after the condition is first met.  Unlike `start_time`/`end_time`, the crossing time is not known in
advance and so cannot be a sync time.  There is no hysteresis: a well whose trigger quantity hovers
on the threshold will chatter from step to step.

Rates are sampled at the *end* of each timestep, so a well opening at `start_time` first flows
during the step that ends at `start_time + dt`.  A `PiecewiseConstant` rate schedule should
therefore use `direction = left`.

!syntax parameters /Wells/PorousFlowInjectionWell

!syntax inputs /Wells/PorousFlowInjectionWell

!syntax children /Wells/PorousFlowInjectionWell
