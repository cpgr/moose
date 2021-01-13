---
title: 'The geochemistry module of MOOSE'
tags:
  - C++
  - geochemistry
  - reactive transport
  - THMC
authors:
  - name: Andy Wilkins
    orcid: 0000-0001-6472-9560
    affiliation: 1
  - name: Christopher P. Green
    orcid: 0000-0002-7315-6597
    affiliation: 1
  - name: Logan Harbour
    orcid: TODO
    affiliation: 2
  - name: Robert Podgorney
    orcid: 0000-0003-4354-4967 TODO
    affiliation: 2
affiliations:
  - name: CSIRO (Commonwealth Scientific and Industrial Research Organisation)
    index: 1
  - name: Idaho National Laboratory
    index: 2
date: 14 January 2021
bibliography: paper.bib
---

# Summary

Geochemical models are chemical models of aqueous solutions within the Earth's subsurface.  Such models are used to understand water chemistry and its interaction with the subsurface in, for example, geothermal wells and reservoirs, contaminant flow through aquifers, ore creation and mining techniques, petroleum and gas exploitation, and rock diagenesis.  The models are frequently extremely complicated, with thousands of interacting chemical species, and often require computer software to find the chemical concentrations.  Reactive-transport simulations involve predicting spatio-temporal chemical-species flows and temperature flows in a spatial domain (such as an aquifer) in addition to solving a chemical model at each point, at each time, in the domain.  Such models are used to explore spatio-temporal dependence of concentration, precipitation, dissolution, sorption, etc, and require significant computational resources.  Models that also consider the deformation of the subsurface in response to the chemistry and flows are even more complex.

# Existing software

A number of commercial and free (open or closed source) packages can solve reactive-transport models.  Some of these are: [CrunchFlow](https://www.netl.doe.gov/sites/default/files/netl-file/CrunchFlow-Manual.pdf) [@beisman2015]; [EQ3/6](https://www-gs.llnl.gov/energy-homeland-security/geochemistry) [@osti_1231666]; the [Geochemist's Workbench](https://www.gwb.com/); [MIN3P](https://www.min3p.com/) [@su2020; @maher2019]; [OpenGeoSys](https://www.opengeosys.org/) [@kolditz2012; @bilke2019]; [PHREEQC](https://www.usgs.gov/software/phreeqc-version-3/) [@pankhurst1995; @pankhurst1999] and [PHAST](https://www.usgs.gov/software/phast-a-computer-program-simulating-groundwater-flow-solute-transport-and-multicomponent) [@pankhurst2010] and [HP1 and HP2](https://www.pc-progress.com/en/Default.aspx?h3d2-hp2); [PHT3D](http://www.pht3d.org/); [Reaktoro](https://reaktoro.org/) [@leal]; and [TOUGHREACT](https://tough.lbl.gov/software/toughreact/) [@xu2004].  Many of these have been developed over decades by teams of researchers, so have sophisticated reaction-modelling capabilities, along with well-established GUIs to assist in model creation and the analysis of results.  Many of them have been used in hundreds of studies by thousands of researchers.

# Statement of need

Researchers are using numerical models to answer increasingly complicated questions.  The questions often involve coupling of different physics, perhaps acting at different scales.  For instance, understanding the impact of subsurface CO$_{2}$ injection (in carbon sequestration schemes) can involve assessing the chemical changes resulting from temperature changes and the presence of CO$_{2}$, and how these impact the porosity and permeability of the subsurface, and how the subsurface deforms and fractures due to the injection.  Another example is understanding ore deposition, which can involve assessing how fluid and heat flows interact with geochemistry within shear bands within the Earth's crust.  Similar remarks hold for hydraulic fracturing, aquifer thermal energy storage, in-situ leaching, and petroleum/gas extraction using enhanced recovery techniques.

In the past, most researchers have answered such questions using separate models, or loosely-coupled models solved using different software packages.  The geochemistry MOOSE module allows researchers to perform stand-alone geochemical modelling, but also draw upon the power of other MOOSE modules to build complicated coupled models and utilise a single code to solve it.

# Overview and comparisons with other software

The geochemistry module is built upon the open-source, massively parallel, fully implicit multiphysics simulation framework MOOSE (Multiphysics Object-Oriented Simulation Environment) [@permann2020moose].  MOOSE is an open-source library from Idaho National Laboratory that provides a high-level interface to the libMesh finite element library [@libmesh] and PETSc nonlinear solvers [@petsc-web-page; @petsc-user-ref; @petsc-efficient].  MOOSE and the geochemistry module follow [strict quality controls](https://mooseframework.org/sqa/index.html).  The geochemistry module's [test suite](https://mooseframework.inl.gov/application_development/test_system.html) contains over 350 tests, ranging from simple unit tests to fully-fledged benchmarks against other codes.

As outlined in this article, the geochemistry module of MOOSE can solve models involving aqueous geochemistry, including aqueous equilibrium, oxidation-reduction, sorption and surface complexation, dissolution-precipitation and gas buffering.  One aspect that makes the geochemistry module different to the codes mentioned above is the ease of coupling additional physics to the geochemical functionality.  In particular, when used in conjuction with the MOOSE's PorousFlow module [@Wilkins2020], sophisticated reactive-transport simulations may be performed, including multi-phase and unsaturated fluid flows, high-precision equations of state for fluids, dynamic porosity and permeability distributions, and sophisticated numerical stabilization.  Geomechanics may also be incorporated by using MOOSE's TensorMechanics module, to explore the interplay between geomechanical stresses and strains, and fluids, heat and reactions.  Fracturing may also be incorporated by using MOOSE's XFEM module [@jiang2020ceramic; @zhang2018modified].

Moreover, harnessing the power of MOOSE means that geochemistry simulations efficiently use multiprocessors and automatically utilize threading/OpenMP, without developers or users having to worry about the software technicalities (MOOSE scales well to over 100,000 cores).  From the user's perspective, a variety of integration schemes and unstructured meshes may be used.  Various time-stepping schemes may be employed.  A variety of I/O formats are available, such as Exodus and VTK files as well as CSV plaintext.  Users may utilize the "MultiApp" system that allows simulations to be embedded inside other simulations.  A generic GUI helps users create, run and analyse models.  Finally, the MOOSE APIs are specifically designed to enable the rapid and simple development of new capability.

Having said this, the aforementioned codes offer advantages too.  Some offer specialist features not available in the geochemistry module.  Their prepostprocessing and postprocessing capabilities are far superior to the generic GUI offered by MOOSE.  Being focussed on reactive-transport alone, their learning curve is much easier.

The geochemistry module's functionality is a subset of that described in the authoratative, pedagogical textbook by @bethke_2007.  It is not possible to describe all geochemical-modelling concepts in this short article, so the reader is referred to @bethke_2007 for further information.

For convenience, the source code for ``PorousFlow`` is bundled within the MOOSE framework.  Detailed documentation found at https://mooseframework.org/modules/geochemistry/index.html.

# Reaction functionality

This section describes the type of geochemical reactions that the geochemistry module can solve.  Both equilibrium and kinetic reactions may be used in models solved by the geochemistry module.

## Database

The geochemical database defines the possible scope of geochemical models.  The geochemistry module assumes the database is in a specific JSON format.  The default database used in all geochemistry-module examples and tests is JSON version of the January 2019 [LLNL thermodynamic database](https://www.gwb.com/thermo.php).    A python script is provided to convert the popular [Geochemist's Workbench](https://www.gwb.com/) and EQ3/6 [@osti_1231666] databases to the JSON format.

## Types of species and reactions

The user may select from the database a subset of chemical components and reactions that are relevant for their particular model.  These are:

- aqueous species in equilibrium
- minerals in equilibrium with the aqueous solution that may precipitate or dissolve
- dissolved gases and gases in equilibrium with an external buffer such as the atmosphere
- species in alternative oxidation states (redox couples)
- sorbed species in equilibrium with the aqueous solution
- minerals, redox species and sorbed species whose dynamics are governed by kinetic rates

## Redox species

Redox species may be "coupled" which means they are in equilibrium with their redox couple, so are treated like any other secondary species.  Alternatively they may be in disequilibrium ("uncoupled") and treated on an equal footing to their primary redox couple (eg, they are provided with an initial condition).  Alternatively, they may be in disequilibrium and governed by a kinetic law.

## Sorption

Sorption may be modelled using a Langmuir approach or a surface-complexation approach.  In the Langmuir approach, a particular mineral's surface contains a number of sorbing sites, and all sorbing species compete to adsorb onto this single type of site.  The more sophisticated approach to surface complexation accounts for the electrical state of the porous-skeleton surface, which can vary sharply with pH, ionic strength and solution composition.  A surface potential for each type of sorbing site is introduced, and appears in the mass-action for each sorbed species.  The surface potential, $\Psi$ \[units: J.C$^{-1}$ = V\] obeys
\begin{equation}
\label{eq:psi}
\frac{A_{\mathrm{sf}}}{Fn_{w}}\sqrt{8RT\epsilon\epsilon_{0}\rho_{w} I}\sinh \left(\frac{\Psi F}{2RT}\right) = \sum_{q}z_{q}m_{q} \ .
\end{equation}
Here:

- $A_{\mathrm{sf}}$ \[units: m$^{2}$\] is the surface area of precipitated mineral upon which the surface site resides.  The user may specify this, or a specific surface area in \[m$^{2}$/g(of mineral)\].
- $n_{w}$ \[units: kg\] is the mass of solvent water
- $F = 96485\ldots \,$C.mol$^{-1}$ is the Faraday constant
- $R = 8.314\ldots\,$J.K$^{-1}$.mol$^{-1}$ is the gas constant
- $T$ \[units: K\] is temperature
- $\epsilon$ \[units: dimensionless\] is the dielectric constant: $\epsilon = 78.5$ at 25$^{\circ}$C for water.
- $\epsilon_{0} = 8.854\times 10^{-12}\,$F.m$^{-1}$ is the permittivity of free space.
- $\rho_{w}=1000\,$kg.m$^{-3}$ is the density of water.
- $I$ \[units: mol.kg$^{-1}$\] is the ionic strength of the aqueous solution
- $q$ represents all sorbed species that sorb to the sorbing site
- $z_{q}$ \[units: dimensionless\] is the charge of the sorbing species
- $m_{q}$ \[units: mol.kg$^{-1}$\] is the molality of the sorbing species, which is a function of $\Psi$ through the mass-action equation

## Kinetics

Minerals, redox species and sorbed species may be governed by kinetic laws:
\begin{equation}
\frac{\mathrm{d}n}{\mathrm{d}t} = - \sum_{\beta}r_{\beta} \ .
\end{equation}
Here:

- $n$ \[units: mol\] is the mole number of the kinetic species
- $t$ \[units: s\] is time
- $r_{\beta}$ \[units: mol.s$^{-1}$\] is a kinetic rate of dissolution

Each term in the sum, $r_{\beta}$, is of the form
\begin{equation}
\label{eqn.indiv.rate}
r = kA[M] \left( \prod_{\alpha}m_{\alpha}^{P_{\alpha}} \right) \left|1 - \left(Q/K\right)^{\theta}\right|^{\eta} \exp\left( \frac{E_{a}}{R} \left(\frac{1}{T_{0}} - \frac{1}{T}\right)\right) {\mathrm{sgn}}(1 - (Q/K)) \ .
\end{equation}
Each $r_{\beta}$ may have different parameters $k$, $A$, $P_{\alpha}$, etc, so acid-neutral-alkali promotion as listed in the correlations prepared by @palandri may be used in geochemistry models.  In this equation:

- $r$ \[units: mol.s$^{-1}$\] is the rate of the kinetic reaction.  If it is positive then the kinetic species' mass will be decreasing (eg, dissolution).  If it is negative then the kinetic species' mass will be increasing (eg, precipitation).
- $k$ is the intrinsic rate constant.  The product $kA[M]$ has units mol.s$^{-1}$
- $[M]$ \[units: g\] is the mass of the kinetic species.  It is optional.
- $A$ is either the surface area \[units: m$^{2}$\] for the kinetic species (if $r$ does not contain $M$), or the specific surface area \[units: m$^{2}$.g$^{-1}$\] (if $r$ contains $M$)
- $\alpha$ is a label denoting a promoting species
- $m_{\alpha}$ is either: mass of solvent water if the promoting species is H$_{2}$O; fugacity of a gas if the promoting species is a gas; activity if the promoting species is either H$^{+}$ or OH$^{-}$; mobility, otherwise.
- $P_{\alpha}$ is a dimensionless power
- $Q$ is the activity product of the kinetic species
- $K$ is the reaction's equilibrium constant defined through the database file
- $\theta$ and $\eta$ are dimensionless exponents
- $E_{a}$ \[units: J.mol$^{-1}$\] is the activation energy
- $R = 8.314\ldots\,$J.K$^{-1}$.mol$^{-1}$ is the gas constant
- $T_{0}$ \[units: K\] is a reference temperature
- $T$ \[units: K\] is the temperature
- ${\mathrm{sgn}}(1-(Q/K))$ is -1 if $Q>K$ (kinetic species mass will increase with time), 1 if $Q<K$ (kinetic species mass will decrease with time).

## Activity and fugacity models

Only the Debye-Huckel B-dot model along with the related formulae for neutral species and water are coded into the geochemistry module.  The virial Pitzer/HMW models are not included.  The activity of all mineral species is assumed to be unity.  The Spycher-Reed [@spycher1988] fugacity formula [@toughreact; @prausnitz] is used in the geochemistry module.

## Constraints

A constraint must be supplied by the user for each member of the equilibrium basis in order that the geochemical system has a unique mathematical solution.  The following types of constraints are possible in the geochemistry module:

- For water: the mass of solvent water, or the total bulk mole number, or the activity
- For aqueous basis species: the free molality, the total bulk number, or the activity (so that pH may be controlled, for instance)
- For precipitated minerals: the free (precipitated) mole number, or the total bulk mole number
- Fugacity for gases with fixed fugacity
- For sorbing sites: the free molality of unoccupied sites, or the total bulk mole number.

In addition, all kinetic species, $A_{\bar{k}}$, must be provided with an initial mole number.

## Reaction paths

The following reaction paths are available in the geochemistry module:

- Adding reactants at user-defined rates.  These reactants may be basis species, secondary species, minerals, etc.
- Controlling temperature with a user-defined function.
- Changing temperature and altering the composition by adding cooler/hotter reactants at user-defined rates.
- Removing one or more species activity constraints or gas fugacity constraints at user-supplied times.
- Controlling species activity (such as the pH) or gas fugacity with user-supplied functions.
- Discarding masses of any minerals present in the equilibrium solution (called a "dump" by @bethke_2007).
- Removing mineral masses at the end of each time-step (called "flow-through" by @bethke_2007)
- Adding pure H$_{2}$O and removing an equal mass of water components and solutes it contains (called "flush" by @bethke_2007)

Combinations of these may be used.  For instance, changing temperature while controlling species activity.

## Mathematical solution strategy

At each time-step, the geochemistry module calculates: the mass of solvent water; the molality of all aqueous species and dissolved gases; the mineral activity products and saturation indices; the mole number of precipitated minerals; the molality of unoccupied surface sites; the surface potentials; the molality of sorbed species; the mole number of kinetic species; and all activities and activity coefficients.  To achieve this, a fully-implicit, under-relaxed, iterative Newton-Raphson procedure is used.

During this procedure, care is taken to avoid overflows or underflows with numbers such as $10^{\pm 1000}$.  Multiple sub-steps may be performed in one time-step to ensure good convergence.  Minerals are allowed to precipitate or dissolve during the procedure.  Charge neutrality is enforced.

# Reactive transport

One of the notable features of the MOOSE framework is the ability to couple different physics together.  Different types of couplings are available: tight (full) coupling; operator-splitting; etc.  To perform reactive-transport simulations using the geochemistry module, we recommend employing MOOSE's PorousFlow module [@Wilkins2020] in an operator-split.  This means that one time-step involves a two-step process.

1. The PorousFlow module transports the solutes in the aqueous phase, as well as the temperature distribution.  This provides the temperature and extra chemical source terms at each point in the spatial domain.
2. With these sources, the geochemistry module solves for the mole numbers of each species at each point in the spatial domain, using the Newton-Raphson method mentioned in the previous section.

It is easy to use multiple sub-steps within (1) and/or (2) if desired.  The operator-split approach is also used by other reactive-transport solvers such the [Geochemist's Workbench](https://www.gwb.com/).  The PorousFlow module is a sophisticated multi-component, multi-phase transport code, and employing it means:

- pressure and temperature tightly couple with fluid flows
- densities, viscosities, etc, may depend on solute concentrations, temperature and pressure
- porosity and permeability can change with precipitation and dissolution
- multiphase flows can be used
- coupling with sophisticated geomechanics (including plasticity [@adhikary2016robust; @wilkins2020method], fracture and large strains) is straightforward
- sophisticated numerical stabilization is available

Nevertheless, rudimentary transport capability has been included into the geochemistry module, which models Darcy flow of the species' mobile concentrations, with hydrodynamic dispersion.  It is assumed that precipitated minerals, sorbed species and sorption sites are immobile.  The temperature, porosity, Darcy flux vector and dispersion tensor are fixed by the user.  All of these may vary spatially and temporally, but the geochemistry module provides them with no dynamics (in contrast with the PorousFlow module, where temperature evolution is governed by an advection-diffusion equation, for instance).  An operator-splitting method may be used to provide the mathematical solution, as described above.

# Computation aspects

Spatially-varying geochemistry simulations use a large amount of memory since they store information about a complete geochemical system at each finite-element node.  On the other hand, ignoring transport, they are almost embarrasingly parallel, so may be solved efficiently using a large number of processors.  Even reactive-transport models scale well with multiple processors, since multi-component porous flow is essentially multiple coupled diffusion equations.  \autoref{fig:scaling_eg} shows an example of scaling in which a model of fixed size is run on multiple processors.  Further experiments involving memory usage, solver choices and cpu-scaling can be found in the [online documentation](https://mooseframework.org/modules/geochemistry/theory/compute_efficiencies.html).

![CPU time required to solve a reactive-transport simulation using the PorousFlow and geochemistry modules.\label{fig:scaling_eg}](joss_paper_scaling.png){ width=60% }

# Example: cooling with feldspars

One of the geochemistry [tests and examples](https://mooseframework.org/modules/geochemistry/tests_and_examples/index.html) involves slowly cooling an aqueous solution from 300$^{\circ}$C to 25$^{\circ}$C.  The aqueous solution is in equilibrium contact with the feldspars albite, maximum microcline, muscovite and quartz.  This example is documented in section 14.1 of @bethke_2007.  \autoref{fig:feldspar_eg} shows the comparison between MOOSE's geochemistry module and the Geochemist's Workbench software.

![Precipitated volumes as a function of temperature when an aqueous solution in contact with feldspars is cooled.\label{fig:feldspar_eg}](joss_paper_feldspar.png){ width=60% }

# Example: reactive-transport in an aquifer thermal energy storage scenario

It has been proposed to use the Weber-Tensleep formation in the USA to store renewably-generated, high-temperature water, for later use in electricity generation.  During the storage phase, formation water is produced from the reservoir, heated to 160$^{\circ}$C under pressure, and re-injected.  The cycle is reversed during electricity generation.  A [3D reactive-transport study](https://mooseframework.org/modules/geochemistry/tests_and_examples/geotes_weber_tensleep.html) assessing potential precipitates in the engineering equipment, and mineralogy changes around the hot-water injection well, is one of the geochemistry module tests and examples.

The major ions in the Weber-Tensleep formation water are Cl$^{-}$, Na$^{+}$, SO$_{4}^{2-}$, HCO$_{3}^{-}$, K$^{+}$, Ca$^{2+}$ and HS$^{-}$, with a host of additional minor ions.  The pH is around 6.5, the ionic strength is around 1.7$\,$mol.kg$^{-1}$, and the temperature around 90$^{\circ}$C.  The observed mineralogy involves quartz (80\%), K-feldspar (8\%), Calcite (5\%), Siderite and Dolomite (each 2\%), Fe-chlorite and Illite (each 1\%) and a host of trace minerals.

The 3D MOOSE model of the hot-water injection involves coupling a model of the heat exchanger to a model of transport to the geochemical model of the reservoir.  These 3 models are loosely coupled using MOOSE's "multiapp" approach, using the geochemistry and PorousFlow modules.  The coupled modelling reveals that anhydrite is the main precipitate in the heat exchanger, that illite and kaolinite dissolve around the injection well, and that K-feldspar and quartz precipitate around the injection well.  \autoref{fig:weber_tensleep_eg} shows some 3D contours of the Weber-Tensleep aquifer.

![Temperature, porosity, pH and free volume of Quartz after 90 days of injection in the 3D coupled model.\label{fig:weber_tensleep_eg}](joss_paper_weber_tensleep.png)


# Acknowledgements - TODO

Funding for the development was provided by the Idaho National Laboratory - TODO check correct wording.  The authors thank the MOOSE framework team, past and present, for providing the MOOSE framework and auxiliary functionality (quality control, test harnesses, documentation scaffolds, build scripts, etc).

# References
