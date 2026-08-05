# PorousFlowLineSourceFromPostprocessor

`PorousFlowLineSourceFromPostprocessor` adds (or removes) fluid at a total mass flux rate computed
by a postprocessor, spread along a line rather than deposited at a single point.  The line is
specified exactly as it is for the other line objects described in [sinks](sinks.md), by
[!param](/DiracKernels/PorousFlowLineSourceFromPostprocessor/point_file) or by
[!param](/DiracKernels/PorousFlowLineSourceFromPostprocessor/line_base).

The rate is *prescribed*: it is whatever the postprocessor says, and does not depend on the pressure
or the mobility in the reservoir.  Use this object to model injection over a completion interval at
a known rate.  When the rate should instead respond to the reservoir, so that the well throttles
itself or honours a bottom-hole pressure, use
[PorousFlowPeacemanBorehole](PorousFlowPeacemanBorehole.md) or
[PorousFlowPolyLineSink](PorousFlowPolyLineSink.md).

The postprocessor gives the total over the whole line.  Each completion receives the share of that
total corresponding to the length of line it represents: half the distance back to the previous
point plus half the distance forward to the next.  The total injected is therefore the same however
finely the line is resolved, which matters because in
[!param](/DiracKernels/PorousFlowLineSourceFromPostprocessor/line_base) mode one point is generated
for each element the line passes through, and those points are regenerated when the mesh changes
under adaptivity or displacement.

!alert note
The weight column of the point file is not used by this object.  For
[PorousFlowPeacemanBorehole](PorousFlowPeacemanBorehole.md) that column is the borehole radius, and
for [PorousFlowPolyLineSink](PorousFlowPolyLineSink.md) it multiplies the flux, but here the total
is divided by length alone.  Similarly
[!param](/DiracKernels/PorousFlowLineSourceFromPostprocessor/line_length) is not needed when the
line has a single completion, because a prescribed rate does not scale with the length of the
completion.

For instance:

!listing modules/porous_flow/test/tests/dirackernels/line_source.i block=DiracKernels

Note that the `execute_on` parameter of the postprocessor is set to `timestep_begin` so that the
correct value is used within the timestep.

!listing modules/porous_flow/test/tests/dirackernels/line_source.i block=Postprocessors

!syntax parameters /DiracKernels/PorousFlowLineSourceFromPostprocessor

!syntax inputs /DiracKernels/PorousFlowLineSourceFromPostprocessor

!syntax children /DiracKernels/PorousFlowLineSourceFromPostprocessor
