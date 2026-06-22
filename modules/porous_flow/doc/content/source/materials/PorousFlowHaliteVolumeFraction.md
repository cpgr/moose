# PorousFlowHaliteVolumeFraction

!syntax description /Materials/PorousFlowHaliteVolumeFraction

This Material computes the solid halite volume fraction $c_{\mathrm{halite}}$ (m$^3$ halite /
m$^3$ porous medium) from the local-equilibrium precipitated-salt mass fraction reported by a
salt-precipitating fluid state, namely [PorousFlowBrineCO2](PorousFlowBrineCO2.md) run with
`precipitate_salt = true`.  The fluid state reports `precipitated_salt` $= m_h / (m_g + m_l)$, the
mass of solid halite per unit fluid mass.  Multiplying by the fluid mass per unit medium volume,
$\phi \sum_{\mathrm{ph}} S_{\mathrm{ph}} \rho_{\mathrm{ph}}$, and dividing by the halite density
gives the volume fraction:
\begin{equation}
\label{eq:halite_vol_frac}
c_{\mathrm{halite}} = \mathrm{precipitated\_salt} \times \phi \,
  \frac{\sum_{\mathrm{ph}} S_{\mathrm{ph}} \rho_{\mathrm{ph}}}{\rho_{\mathrm{halite}}} \ ,
\end{equation}
where $S_{\mathrm{ph}}$ and $\rho_{\mathrm{ph}}$ are the saturation and density of phase
$\mathrm{ph}$, $\phi$ is the porosity, and $\rho_{\mathrm{halite}}$ is the solid halite density
supplied through the `halite_density` parameter (default 2165 kg/m$^3$).

The computation uses the *old* value of the porosity $\phi_{\mathrm{old}}$ in
[eq:halite_vol_frac].  This breaks the cyclic dependence between porosity and halite volume
fraction, exactly as [PorousFlowAqueousPreDisMineral](PorousFlowAqueousPreDisMineral.md) does for
kinetic minerals (see [Porosity](/porous_flow/porosity.md) for more details).  Because the old
porosity is fixed data carrying no current-step derivative, it contributes no current-step Jacobian
term.  At $t = 0$, where the old porosity is unavailable, the current porosity is used instead (the
two coincide) so that a simulation may start already oversaturated, with solid halite present,
without losing the excess salt at the first step.

The Material declares the scalar property `PorousFlow_halite_volume_fraction` (at the nodes or
quadpoints depending on the `at_nodes` flag).  This is the volume fraction consumed by
[PorousFlowPorosity](PorousFlowPorosity.md) when `chemical_equilibrium = true`, supplied through its
`equilibrium_mineral` parameter, to reduce porosity as halite fills the pore space.

A templated AD version, `ADPorousFlowHaliteVolumeFraction`, is available for the (AD-only)
finite-volume models.  The AD path propagates the derivatives automatically through the generic
precipitated-salt, saturation and density properties, while the non-AD version hand-codes them.

!syntax parameters /Materials/PorousFlowHaliteVolumeFraction

!syntax inputs /Materials/PorousFlowHaliteVolumeFraction

!syntax children /Materials/PorousFlowHaliteVolumeFraction
