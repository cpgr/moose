# PorousFlowConstrainedPeacemanBorehole

A `PorousFlowConstrainedPeacemanBorehole` is a [PorousFlowPeacemanBorehole](PorousFlowPeacemanBorehole.md)
whose bottom-hole pressure is not prescribed but solved for. The pressure is held in a scalar
variable, named by [!param](/DiracKernels/PorousFlowConstrainedPeacemanBorehole/bottom_hole_pressure),
and the borehole adds the constraint that fixes it: the flux summed over all of the completions
equals [!param](/DiracKernels/PorousFlowConstrainedPeacemanBorehole/rate). The well therefore
delivers the requested rate exactly at every timestep, and the pressure that achieves it is found
implicitly alongside the reservoir equations rather than chased by an outer controller.

The flux at each completion is the same Peaceman expression the parent computes, so the allocation
along the completions still follows the local $W \lambda \Delta P$: high permeability layers take
more, layers that have pressured up take less, and the split shifts as the simulation runs. A source
of prescribed rate imposes a fixed allocation instead, which is wrong in any heterogeneous or
partially pressurised interval.

Please see [sinks](sinks.md) for an extended discussion of line sinks, and
[PorousFlowDiracKernelScalarBase](PorousFlowDiracKernelScalarBase.md) for how a `DiracKernel`
assembles into a scalar variable at all.

## Differences from the parent

Neither [!param](/DiracKernels/PorousFlowPeacemanBorehole/bottom_p_or_t) nor
[!param](/DiracKernels/PorousFlowPeacemanBorehole/character) applies here, and both are suppressed.
The `character` parameter of the parent gates the flow, so that a production well does nothing where
the borehole pressure exceeds the porepressure and an injection well does nothing where it does not.
With the bottom-hole pressure an unknown, that gate would make the flux non-differentiable exactly
at the operating point the constraint settles on, and it would suppress crossflow between
completions, which is physical. The direction of flow at a completion instead follows the sign of
the pressure difference.

!alert note
The scalar variable's only equation is the one this object assembles, and the kernel coverage check
does not look at `DiracKernels`, so it will report the scalar as uncovered. Add a
[NullScalarKernel](NullScalarKernel.md) on it, as described in
[PorousFlowDiracKernelScalarBase](PorousFlowDiracKernelScalarBase.md).

!syntax parameters /DiracKernels/PorousFlowConstrainedPeacemanBorehole

!syntax inputs /DiracKernels/PorousFlowConstrainedPeacemanBorehole

!syntax children /DiracKernels/PorousFlowConstrainedPeacemanBorehole
