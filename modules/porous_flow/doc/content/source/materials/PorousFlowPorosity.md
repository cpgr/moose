# PorousFlowPorosity

!syntax description /Materials/PorousFlowPorosity

This Material computes porosity (at the nodes or quadpoints, depending on the `at_nodes` flag):
\begin{equation}
\label{eq:poro_evolve}
\phi + M = \alpha_{B} + (\phi_{0} + M_{\mathrm{ref}} - \alpha_{B})\times \exp \left( \frac{\alpha_{B}'
  - 1}{K}(P_{f} - P_{f}^{\mathrm{ref}}) - \epsilon^{\mathrm{total}}_{ii} + \alpha_{T}(T - T^{\mathrm{ref}}) \right) \ ,
\end{equation}
A full description is provided in the [porosity documentation](/porous_flow/porosity.md)

Flags provided to `PorousFlowPorosity` control its evolution.

- If `mechanical = true` then the porosity will depend on $\epsilon^{\mathrm{total}}_{ii}$.
  Otherwise that term in [eq:poro_evolve] is ignored.

- If `fluid = true` then the porosity will depend on $(P_{f} - P_{f}^{\mathrm{ref}})$.  Otherwise
  that term in [eq:poro_evolve] is ignored.

- If `thermal = true` then the porosity will depend on $(T - T^{\mathrm{ref}})$.  Otherwise that term
  in [eq:poro_evolve] is ignored.

- If `chemical = true` then porosity will depend on $M$.  Otherwise that term in
  [eq:poro_evolve] is ignored.

- If `chemical_equilibrium = true` then porosity is additionally reduced by an equilibrium-mineral
  volume fraction (m$^3$ mineral / m$^3$ porous medium) supplied as a dedicated scalar material
  property, named by `equilibrium_mineral` (default `halite_volume_fraction`, as produced by
  `PorousFlowHaliteVolumeFraction` for salt-precipitating fluid states).  The volume fraction $c$
  enters as an extra contribution $w (c - c_{\mathrm{ref}})$ to
  $M$, with weight $w$ = `equilibrium_weight` and reference $c_{\mathrm{ref}}$ =
  `equilibrium_reference`.  This is independent of the kinetic `chemical` option (the two may be
  combined) and uses the OLD value of the volume fraction to break the porosity/precipitate cyclic
  dependency, so it contributes no current-step Jacobian term.  Combine with `porosity_min` (below)
  to keep porosity positive once the mineral fills the pore space.

!alert note title=Lower bound on porosity (`porosity_min`)
The optional parameter `porosity_min` places a hard lower bound on the computed
porosity: if [eq:poro_evolve] yields a value below `porosity_min`, the porosity is
set to `porosity_min` instead (its derivatives are scaled by `zero_modifier` rather
than zeroed, to aid Newton convergence).  By default no bound is imposed.  This is
primarily useful for `chemical = true` porosity: the exponential form only
guarantees a positive porosity for a positive decay (pore-pressure/strain/thermal)
term, so mineral precipitation - which enters through $M$, not the decay - can
otherwise drive the porosity negative once the precipitated mineral fills the pore
space.  Setting `porosity_min` to a small positive value prevents this.

!alert note
On the parameter `solid_bulk`: Selecting a functor with dependency on a model
variable $u$ (such as pore pressure) for the parameter `solid_bulk` $K$ in a way
that $\text{d}K / \text{d}u \neq 0$ without taking further precautions may result in an
inexact Jacobian matrix, which can lead to convergence problems. These problems
should not be triggered by a constant solid bulk or a pure dependency on the
model's time and/or spatial coordinates.

!syntax parameters /Materials/PorousFlowPorosity

!syntax inputs /Materials/PorousFlowPorosity

!syntax children /Materials/PorousFlowPorosity
