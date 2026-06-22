# CO2 injection with halite precipitation (dry-out)

When supercritical CO$_2$ is injected into a saline aquifer it evaporates water from the
brine.  The dissolved salt left behind concentrates the remaining brine until it reaches the
halite solubility limit, at which point solid halite precipitates, filling pore space and
choking the permeability near the wellbore.  This loss of injectivity ("dry-out" or "salting
out") is a well-known operational concern for geological CO$_2$ storage.

This example models the process with the [brine-co2](brineco2.md) equation of state run in its
salt-precipitating mode.  It extends the [1D radial intercomparison problem](1Dradial.md) so
that the salt is a conserved nonlinear variable that is allowed to drop out of solution as
solid halite once the brine saturates.

The input file is `modules/porous_flow/examples/halite_precipitation/1Dradial_dryout.i`.  See
[the comparison below](#comparison) for why the native local-equilibrium path is used here in
preference to the kinetic precipitation-dissolution chemistry.

## Method

Salt is fluid component 2 (`xnacl`), carried as a conserved total-salt variable with its own
mass balance.  The [PorousFlowBrineCO2](PorousFlowBrineCO2.md) fluid state is activated with
`precipitate_salt = true`, so that its flash imposes *local equilibrium*: the dissolved
salinity is clamped at the solubility $X_\mathrm{eq}(T)$ and any excess salt is reported as a
precipitated-salt mass fraction.  Precipitation is therefore instantaneous and exact at the
solubility limit, with no kinetic-rate, surface-area or solubility-constant parameters to tune.

!listing modules/porous_flow/examples/halite_precipitation/1Dradial_dryout.i block=UserObjects/fs

The solid halite is handled by three pieces:

- [PorousFlowHaliteVolumeFraction](PorousFlowHaliteVolumeFraction.md) converts the
  precipitated-salt mass fraction into a solid halite volume fraction $c_\mathrm{halite}$
  (m$^3$ halite / m$^3$ porous medium).
- [PorousFlowPrecipitateMassTimeDerivative](PorousFlowPrecipitateMassTimeDerivative.md) carries
  $\partial(\rho_\mathrm{halite}\, c_\mathrm{halite})/\partial t$ on the salt equation, removing
  the precipitated salt from the brine balance.  Because the same conserved salt either remains
  dissolved or becomes solid halite, the scheme is exactly mass conservative.
- [PorousFlowPorosity](PorousFlowPorosity.md) with `chemical_equilibrium = true` subtracts
  $c_\mathrm{halite}$ from the reference porosity, and
  [PorousFlowPermeabilityKozenyCarman](PorousFlowPermeabilityKozenyCarman.md) chokes the
  permeability as the porosity drops.

!listing modules/porous_flow/examples/halite_precipitation/1Dradial_dryout.i start=[Materials] end=[BCs]

### Reservoir properties

The reservoir properties match the [1D radial intercomparison problem](1Dradial.md) and are
summarized in [tab:halite_res].

!table id=tab:halite_res caption=Reservoir properties
| Property | Value |
| - | - |
| Pressure | 12 MPa |
| Temperature | 45 $^{\circ}$C |
| Permeability | $10^{-13}$ m$^2$ (100 md) |
| Porosity | 0.12 |
| Initial NaCl mass fraction | 0.15 |
| Halite solubility $X_\mathrm{eq}$(45 $^{\circ}$C) | 0.2672 |
| Halite density | 2165 kg m$^{-3}$ |

CO$_2$ is injected at the wellbore (inner radius 0.1 m) using a self-limiting flux that ramps
to zero as the well pressure approaches a bottom-hole limit, modelling injectivity loss as the
near-well rock clogs.  The outer boundary holds the reservoir pressure.  The model is
isothermal; the solubility limit $X_\mathrm{eq}(45\,^{\circ}\mathrm{C}) = 0.2672$ is evaluated
internally by the flash, so the local-equilibrium path is fully non-isothermal-capable with no
additional wiring.

## Results

[fig:halite_saturation] shows the gas saturation profile at several times.  A fully
gas-saturated region ($S_g = 1$) emerges near the well and grows radially outward as the brine
dries out.  [fig:halite_volfrac] shows the corresponding solid halite volume fraction: halite
forms a bank that sits exactly behind the dry-out front, dropping to zero precisely where $S_g$
comes off unity.

!media porous_flow/halite_dryout_saturation.png
       id=fig:halite_saturation
       style=width:60%;margin-left:10px;
       caption=Gas saturation profile at various times; the fully dried ($S_g = 1$) region advances radially outward from the well.

!media porous_flow/halite_dryout_volfrac.png
       id=fig:halite_volfrac
       style=width:60%;margin-left:10px;
       caption=Solid halite volume fraction at various times; the bank tracks the dry-out front and advances with it.

The bank is *distributed* across the dry-out zone rather than concentrated in a single cell,
advancing from one cell at 0.5 days to roughly 0.7 m by 10 days, as summarized in
[tab:halite_front].  This is the in-place equilibrium precipitation behaviour reported for
TOUGH2/ECO2N radial dry-out simulations, in which solid salt is distributed across the dried
zone rather than concentrated at the well.

!table id=tab:halite_front caption=Halite bank extent versus time
| Time | Dried extent ($S_g = 1$) | Peak $c_\mathrm{halite}$ |
| - | - | - |
| 0.5 days | $\sim$0.14 m | 0.0023 |
| 1 day | $\sim$0.22 m | 0.0050 |
| 2 days | $\sim$0.30 m | 0.0050 |
| 5 days | $\sim$0.46 m | 0.0050 |
| 10 days | $\sim$0.60 m | 0.0050 |

The halite plateau is essentially capped at $c_\mathrm{halite} \approx 0.005$: once a cell
reaches $S_g = 1$ there is no more brine to evaporate, so no further salt accumulates there and
the bank advances without thickening.  At the initial salinity of 0.15 this corresponds to a
mild clogging - porosity dips from 0.12 to about 0.115 and permeability from $10^{-13}$ to
$9.1\times 10^{-14}$ m$^2$.  A higher initial salinity (closer to $X_\mathrm{eq}$) produces a
correspondingly thicker bank and stronger clogging.

Because no salt crosses the boundary and the conserved salt is merely repartitioned between
dissolved and solid, the total salt mass is conserved to solver tolerance throughout the run -
the `total_salt_kg` postprocessor stays constant.

## Comparison with the kinetic chemistry stack id=comparison

The same physics can be approached with the native aqueous precipitation-dissolution (PreDis)
chemistry: [PorousFlowAqueousPreDisChemistry](PorousFlowAqueousPreDisChemistry.md),
[PorousFlowAqueousPreDisMineral](PorousFlowAqueousPreDisMineral.md) and the
[PorousFlowPreDis](PorousFlowPreDis.md) kernel, driven off `xnacl` against a solubility
"constant".

The kinetic stack premultiplies its reaction rate by the aqueous-phase saturation, so the
precipitation rate vanishes as a cell dries out.  This makes it unable to cap the dissolved
salinity through full dry-out: the brine cannot shed its excess salt fast enough, and keeping
the simulation in the valid salinity range requires a variational-inequality bound on `xnacl`.
Once that bound becomes active it supplies salt that the reaction keeps converting to halite,
so the scheme spuriously creates salt and the halite collects in the single near-well cell.
The local-equilibrium path used in `1Dradial_dryout.i` avoids this entirely - the flash caps the
dissolved salinity at $X_\mathrm{eq}$ exactly, needs no bound, and produces the distributed,
mass-conservative bank shown above.

!bibtex bibliography
