# PorousFlowProductionWell

!syntax description /Wells/PorousFlowProductionWell

A production well operating on a deliverability curve, specified inside the
[Wells](syntax/Wells/index.md) block.  The flux out of the reservoir is a piecewise-linear function
of the local porepressure (or temperature), given by the `p_or_t_vals`/`fluxes` pairs.

The [PorousFlowWellAction](/Wells/index.md) expands this object into, for each fluid component:

- a [PorousFlowSumQuantity](PorousFlowSumQuantity.md) user object (one per kernel: the line sink
  accumulates into it, so it cannot be shared);
- a [PorousFlowPolyLineSink](PorousFlowPolyLineSink.md) with `mass_fraction_component` set to that
  component, so the produced stream carries the resident composition of the reservoir;
- the produced-mass postprocessor `<well>_mass_<i>` (kg in the step), the rate `<well>_q_<i>`, and
  the cumulative `<well>_cum_<i>`.

The deliverability curve is a crude but honest model of well productivity.  It is deliberately
preferred over a fixed-rate producer, which has no feedback at all and will pull the near-well cell
to unphysical pressures once the prescribed rate outruns what the formation can deliver.

## Sign convention

`fluxes` and every reported rate and cumulative are positive in the well's natural sense: a positive
value produces fluid out of the reservoir.

## Shutting the well in

`start_time`/`end_time` and `shut_if` shut the well in, exactly as for an injection well.  A
producer's rate comes from its deliverability curve rather than from a rate function, so there is no
rate to zero: the action instead builds an AuxVariable carrying the gate and passes it to each
polyline sink as `multiplying_var`, which multiplies both the flux and its Jacobian.  The sink
therefore accumulates zero into its `PorousFlowSumQuantity` while the well is shut, and the produced
mass, rate and cumulative all correctly read zero.

## Limitations

Bottom-hole pressure control (`bhp`, `bhp_min`) is not yet supported and raises an error.

!syntax parameters /Wells/PorousFlowProductionWell

!syntax inputs /Wells/PorousFlowProductionWell

!syntax children /Wells/PorousFlowProductionWell
