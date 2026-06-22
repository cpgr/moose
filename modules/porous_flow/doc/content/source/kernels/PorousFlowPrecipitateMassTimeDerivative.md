# PorousFlowPrecipitateMassTimeDerivative

!syntax description /Kernels/PorousFlowPrecipitateMassTimeDerivative

This `Kernel` implements the weak form of the time derivative of the solid (precipitated) halite
mass per unit volume,
\begin{equation*}
  \frac{\partial}{\partial t}\left(\rho_{\mathrm{halite}}\, c_{\mathrm{halite}}\right) \ ,
\end{equation*}
where $c_{\mathrm{halite}}$ is the solid halite volume fraction (m$^3$ halite / m$^3$ porous
medium) computed by [PorousFlowHaliteVolumeFraction](PorousFlowHaliteVolumeFraction.md) and
$\rho_{\mathrm{halite}}$ is the halite density set by the `halite_density` parameter (default
2165 kg m$^{-3}$, which must match the value used by the material).

It is applied to the salt-component residual so that, together with the
[PorousFlowMassTimeDerivative](PorousFlowMassTimeDerivative.md) and
[PorousFlowAdvectiveFlux](PorousFlowAdvectiveFlux.md) terms for the dissolved salt, the conserved
salt is exactly repartitioned between the dissolved and solid phases as halite precipitates or
dissolves.  Unlike the kinetic [PorousFlowPreDis](PorousFlowPreDis.md) kernel, this term is not
premultiplied by the aqueous-phase saturation, so it does not vanish as a cell dries out - which
is what allows the local-equilibrium [PorousFlowBrineCO2](PorousFlowBrineCO2.md) path to remain
mass conservative through complete dry-out.

The halite mass is lumped to the nodes (as in
[PorousFlowMassTimeDerivative](PorousFlowMassTimeDerivative.md)), so the only non-zero Jacobian
contribution is the diagonal nodal term.  An end-to-end example using this kernel is given in
the [CO2 dry-out / halite precipitation example](halite_precipitation.md).

!syntax parameters /Kernels/PorousFlowPrecipitateMassTimeDerivative

!syntax inputs /Kernels/PorousFlowPrecipitateMassTimeDerivative

!syntax children /Kernels/PorousFlowPrecipitateMassTimeDerivative
