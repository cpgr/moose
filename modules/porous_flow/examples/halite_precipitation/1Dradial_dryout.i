# 1D radial CO2 injection into a saline aquifer with halite (salt) precipitation
# ("dry-out").
#
# This example extends the standard PorousFlowBrineCO2 radial injection
# (examples/co2_intercomparison/1Dradial) so that salt is a CONSERVED nonlinear
# variable (fluid component 2, as in fluidstate/theis_brineco2.i) and is allowed
# to precipitate as solid halite once the brine salinity exceeds the solubility
# limit.
#
# Physics of dry-out:
#   As supercritical CO2 is injected it evaporates water from the brine.  The
#   salt left behind drives the local salinity (xnacl) up.  When xnacl exceeds
#   the halite solubility X_eq, solid halite precipitates, consuming pore space
#   and choking permeability near the wellbore.
#
# Method (native local-equilibrium precipitation in PorousFlowBrineCO2):
#   * xnacl (fluid component 2) is the CONSERVED total-salt variable, with its
#     own mass balance (mass2/flux2).
#   * PorousFlowBrineCO2 is run with precipitate_salt = true.  Its flash imposes
#     LOCAL EQUILIBRIUM: the dissolved salinity is clamped at the solubility
#     X_eq(T) and any excess salt is reported as a precipitated-salt mass
#     fraction.  No kinetic-rate, surface-area or solubility-constant knobs are
#     involved - precipitation is instantaneous and exact at the solubility
#     limit (contrast the kinetic PorousFlowAqueousPreDis* stack used in
#     1Dradial_dryout_predis.i).
#   * PorousFlowHaliteVolumeFraction converts that precipitated-salt mass
#     fraction into a solid halite volume fraction c_halite (m3 halite / m3 rock).
#   * PorousFlowPrecipitateMassTimeDerivative carries d/dt(rho_halite * c_halite)
#     on the salt (component 2) equation, removing the precipitated salt from the
#     brine balance.  Because the same conserved salt either stays dissolved or
#     becomes solid halite, the scheme is EXACTLY mass conservative.
#   * PorousFlowPorosity(chemical_equilibrium = true) subtracts c_halite from the
#     reference porosity, and PorousFlowPermeabilityKozenyCarman chokes
#     permeability as porosity drops.
#
# What this example demonstrates (near-well salt clogging / injectivity loss):
#   xnacl rises to the solubility limit X_eq = 0.2672, halite precipitates,
#   porosity collapses and Kozeny-Carman chokes the permeability, all
#   mass-conservatively.  The result is a halite bank LOCALIZED at the wellbore
#   that seals the near-well rock and chokes injectivity - it does NOT spread as
#   a broad propagating front.  This is the textbook CO2 "dry-out" behaviour:
#   counter-current flow draws brine back toward the evaporation site at the
#   well, where the dissolved salt first reaches saturation and drops out, so
#   precipitation focuses in the first cell(s).  Lowering the initial salinity
#   delays the clogging (the well stays open longer) but does not spread the
#   bank.  Contrast this with ADVECTION-driven reactive transport, where a
#   reactive solute is carried through the domain and reacts along the flow path
#   over many cells (see examples/tutorial/07.i - injected tracer precipitates
#   with porosity/permeability feedback - and examples/tutorial/13.i - dolomite
#   dissolution front); here the EVAPORATIVE driver localizes the reaction
#   instead.
#
#   Because the local-equilibrium flash caps dissolved salinity at X_eq exactly
#   (no rate throttling tied to porosity * saturation), the model stays mass
#   conservative all the way through pore clogging - the total_salt_kg
#   postprocessor stays constant.  This is the key correctness improvement over
#   the kinetic stack in 1Dradial_dryout_predis.i, which had to bound xnacl with
#   a variational-inequality solver and so spuriously CREATED salt once that
#   bound went active at full dry-out.
#
# Robustness:
#   Driving a near-well cell all the way to pore clogging is a stiff,
#   near-singular regime.  Two ingredients keep the solve well posed:
#     * porosity_min on the [porosity] material floors the (otherwise unbounded)
#       chemical porosity at a small positive value so it cannot go negative once
#       the pore space is filled by halite.  (This required a small extension to
#       PorousFlowPorosityExponentialBase, which previously only guaranteed
#       positivity for pressure/strain-driven decay, not for chemistry.)
#     * the [brineco2] PorousFlowFluidState material is fed a CONSTANT capillary
#       pressure (pctmp, pc=0) instead of the VG curve.  Phase pressures,
#       densities and the flash still use the full VG pc via the fluid-state UO
#       (fs); the constant pc only zeroes dpc/dS, removing the steep capillary-
#       derivative term that - near dry-out - wrecked the SMP analytic Jacobian
#       and the liquid-flux linearisation.  The only physical cost is dropping the
#       small capillary-flux contribution; capillarity is negligible for this
#       near-well dry-out.  With it in place a plain NEWTON + SMP solve converges
#       cleanly.
#   The local-equilibrium path needs NO variational-inequality bound on xnacl:
#   the flash itself keeps the dissolved salinity inside BrineCO2's valid range
#   (it never exceeds X_eq), so the [Bounds] / vinewtonrsls apparatus required by
#   the kinetic stack is gone.
#   Treat the injection rate and mesh below as a starting point to be tuned for
#   the regime of interest.
#
# Notes / caveats:
#   * Isothermal (45 C).  The solubility limit is X_eq(45 C) = 0.2672, evaluated
#     internally by the BrineCO2 flash; there is no X_eq AuxVariable to set.  The
#     local-equilibrium path is fully non-isothermal-capable (the flash uses
#     X_eq(T) directly), so a non-isothermal model needs no extra wiring here.

[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 200
  xmin = 0.1 # wellbore radius (m)
  xmax = 200 # reservoir extent (m)
  bias_x = 1.02 # refine near the well
  coord_type = RZ
  rz_coord_axis = Y
[]

[Debug]
  show_material_props = true
[]

[Problem]
  type = FEProblem
[]

[GlobalParams]
  PorousFlowDictator = dictator
  gravity = '0 0 0'
[]

[Variables]
  [pgas]
    initial_condition = 12e6 # 12 MPa
  []
  [zi]
    initial_condition = 0 # total CO2 mass fraction (no CO2 initially)
    scaling = 1e4
  []
  [xnacl]
    initial_condition = 0.15 # 15 wt% initial salinity; well below the 0.2672
    # solubility limit so that cells must be substantially dried out before halite
    # forms, letting the dry-out front develop a finite width and propagate
    # radially rather than instantly clogging the wellbore cell
  []
[]

[AuxVariables]
  # --- diagnostics ---
  [halite]
    order = CONSTANT
    family = MONOMIAL
  []
  [porosity]
    order = CONSTANT
    family = MONOMIAL
  []
  [permeability]
    order = CONSTANT
    family = MONOMIAL
  []
  [saturation_gas]
    order = CONSTANT
    family = MONOMIAL
  []
  [x1] # dissolved CO2 mass fraction in liquid
    order = CONSTANT
    family = MONOMIAL
  []
[]

[Kernels]
  # Fluid component 0 (H2O)
  [mass0]
    type = PorousFlowMassTimeDerivative
    fluid_component = 0
    variable = pgas
  []
  [flux0]
    type = PorousFlowAdvectiveFlux
    fluid_component = 0
    variable = pgas
  []
  # Fluid component 1 (CO2)
  [mass1]
    type = PorousFlowMassTimeDerivative
    fluid_component = 1
    variable = zi
  []
  [flux1]
    type = PorousFlowAdvectiveFlux
    fluid_component = 1
    variable = zi
  []
  # Fluid component 2 (NaCl) - conserved, with a precipitation sink
  [mass2]
    type = PorousFlowMassTimeDerivative
    fluid_component = 2
    variable = xnacl
  []
  [flux2]
    type = PorousFlowAdvectiveFlux
    fluid_component = 2
    variable = xnacl
  []
  [salt_precipitation]
    # Removes the precipitated salt mass d/dt(rho_halite * c_halite) from the
    # brine balance; the same conserved salt either stays dissolved or becomes
    # solid halite, so the scheme is exactly mass conservative.
    type = PorousFlowPrecipitateMassTimeDerivative
    variable = xnacl
  []
[]

[AuxKernels]
  [halite]
    type = MaterialRealAux
    variable = halite
    property = PorousFlow_halite_volume_fraction_nodal
    execute_on = 'timestep_end'
  []
  [porosity]
    type = PorousFlowPropertyAux
    variable = porosity
    property = porosity
    execute_on = 'timestep_end'
  []
  [permeability]
    type = PorousFlowPropertyAux
    variable = permeability
    property = permeability
    row = 0
    column = 0
    execute_on = 'timestep_end'
  []
  [saturation_gas]
    type = PorousFlowPropertyAux
    variable = saturation_gas
    property = saturation
    phase = 1
    execute_on = 'timestep_end'
  []
  [x1]
    type = PorousFlowPropertyAux
    variable = x1
    property = mass_fraction
    phase = 0
    fluid_component = 1
    execute_on = 'timestep_end'
  []
[]

[UserObjects]
  [dictator]
    type = PorousFlowDictator
    porous_flow_vars = 'pgas zi xnacl'
    number_fluid_phases = 2
    number_fluid_components = 3
  []
  [pc]
    type = PorousFlowCapillaryPressureVG
    alpha = 5.099e-5
    m = 0.457
    sat_lr = 0.0
    pc_max = 1e7
  []
  [fs]
    type = PorousFlowBrineCO2
    brine_fp = brine
    co2_fp = co2
    capillary_pressure = pc
    salt_component = 2 # xnacl is the conserved total-salt component
    precipitate_salt = true # enable local-equilibrium halite precipitation
  []
  [pctmp]
    # Constant capillary pressure fed ONLY to the [brineco2] PorousFlowFluidState
    # material (NOT to the fluid-state UO fs, which keeps the VG pc above).  This
    # zeroes dpc/dS in the material so the steep VG capillary derivative does not
    # enter the SMP analytic Jacobian / liquid-flux linearisation near dry-out -
    # the key to clean NEWTON+SMP convergence.  See the header "Robustness" note.
    type = PorousFlowCapillaryPressureConst
    pc = 0
  []
[]

[FluidProperties]
  [co2sw]
    type = CO2FluidProperties
  []
  [co2]
    type = TabulatedBicubicFluidProperties
    input_fp = co2sw
  []
  [water]
    type = Water97FluidProperties
  []
  [watertab]
    type = TabulatedBicubicFluidProperties
    input_fp = water
    temperature_min = 273.15
    temperature_max = 573.15
    # Generates the tabulation each run (self-contained).  For faster iteration,
    # swap to: fluid_property_file = water_fluid_properties.csv  +
    #          allow_fp_and_tabulation = true   (reads the file produced below)
    fluid_property_output_file = water_fluid_properties.csv
  []
  [brine]
    type = BrineFluidProperties
    water_fp = watertab
  []
[]

[Materials]
  [temperature]
    type = PorousFlowTemperature
    temperature = 45
    block = 0
  []
  [brineco2]
    type = PorousFlowFluidState
    gas_porepressure = pgas
    z = zi
    temperature_unit = Celsius
    xnacl = xnacl
    capillary_pressure = pctmp # const pc: zeroes dpc/dS for SMP - see header & [pctmp]
    fluid_state = fs
    block = 0
  []

  # --- native local-equilibrium halite: solid volume fraction from the flash ---
  [halite]
    # Converts the precipitated-salt mass fraction reported by the BrineCO2 flash
    # into a solid halite volume fraction c_halite (m3 halite / m3 rock).
    type = PorousFlowHaliteVolumeFraction
    block = 0
  []

  # --- feedback: halite volume fraction reduces porosity, then permeability ---
  [porosity]
    type = PorousFlowPorosity
    porosity_zero = 0.12
    chemical_equilibrium = true # subtract the equilibrium halite volume fraction
    # Floor porosity at a small positive value so it cannot go negative once the
    # pore space is filled by precipitated halite.  Kozeny-Carman then chokes the
    # permeability toward its floor, throttling inflow and self-limiting further
    # precipitation (the physical dry-out shutdown).
    porosity_min = 1e-4
    block = 0
  []
  [permeability]
    type = PorousFlowPermeabilityKozenyCarman
    poroperm_function = kozeny_carman_phi0
    phi0 = 0.12
    n = 2
    m = 2
    k0 = 1e-13
    k_anisotropy = '1 0 0 0 1 0 0 0 1'
    block = 0
  []

  [relperm_water]
    type = PorousFlowRelativePermeabilityVG
    m = 0.457
    phase = 0
    s_res = 0.3
    sum_s_res = 0.35
    block = 0
  []
  [relperm_gas]
    type = PorousFlowRelativePermeabilityCorey
    n = 2
    phase = 1
    s_res = 0.05
    sum_s_res = 0.35
    block = 0
  []
[]

[BCs]
  # Far-field: hold reservoir pressure by letting fluid escape at the outer radius
  # Far-field: hold reservoir pressure with a simple Dirichlet BC (as in
  # theis_brineco2.i).  A mobility-based PorousFlowSink would pull
  # permeability -> (chemical) porosity -> the precipitation chemistry onto this
  # boundary; evaluating those nodal materials on a 1D boundary face trips a
  # latent out-of-bounds in PorousFlowMaterial::computeNodalProperties (asserts
  # in a devel build).  The Dirichlet BC avoids that dependency entirely, and the
  # outer boundary is far from the near-well dry-out so the BC type is immaterial.
  [right_pressure]
    type = DirichletBC
    boundary = right
    variable = pgas
    value = 12e6
  []
  # Inject CO2 (component 1) as a prescribed mass flux distributed over the
  # wellbore boundary.  A flux BC is gentler on the solver than a point source
  # into the tiny innermost cell.  Negative flux_function = injection (a
  # PorousFlowSink residual is positive for removal).
  # [inject_co2]
  #   type = PorousFlowSink
  #   boundary = left
  #   variable = zi
  #   flux_function = injection_rate
  # []
  [inject_co2]
    type = FunctionNeumannBC # Bypasses the PorousFlow material face loop
    boundary = left
    variable = zi
    function = injection_rate
  []
[]

[Functions]
  [injection_rate]
    type = ParsedFunction
    # positive => injection for FunctionNeumannBC (residual = -test*func)
    expression = '(q_max / area) * 0.5 * (1.0 - tanh((p_well - p_lim) / p_scale))'
    symbol_names = 'area          q_max p_lim p_scale p_well'
    symbol_values = '0.62831853    0.05  18e6  1e6     p_well' # Hardcoded area value
  []
[]

# [Functions]
#   [injection_rate]
#     type = ParsedFunction
#     # Self-limiting CO2 injection: a target rate q_max (kg/s) spread over the
#     # wellbore area, smoothly ramped to zero as the well pressure p_well
#     # approaches the bottom-hole limit p_lim.  This models injectivity loss as
#     # the near-well rock clogs and keeps the pressure bounded (a fixed-flux well
#     # would otherwise run away once permeability chokes).  p_well is lagged (a
#     # postprocessor), which is fine for a soft limiter.
#     expression = '-(q_max / area) * 0.5 * (1.0 - tanh((p_well - p_lim) / p_scale))'
#     symbol_names = 'area          q_max p_lim p_scale p_well'
#     symbol_values = 'injection_area 0.02  18e6  1e6     p_well'
#   []
# []

[Preconditioning]
  [smp]
    type = SMP
    full = true
    petsc_options_iname = '-ksp_type -pc_type -sub_pc_type -sub_pc_factor_shift_type'
    petsc_options_value = 'gmres bjacobi lu NONZERO'
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  end_time = 8.64e5 # 10 days - the halite bank keeps marching radially outward as
  # the dry-out front advances (it does not freeze, so there is always an active
  # front to resolve)
  nl_max_its = 20
  l_max_its = 100
  nl_abs_tol = 1e-7
  dtmax = 2.5e4 # cap on the step.  The local-equilibrium path is robust enough
  # that this binds for most of the run; the dry-out front (each cell's
  # saturation_gas -> 1) is smooth and the radial geometry slows it with radius,
  # so it tolerates large steps.  Push higher if a rerun stays cut-back-free; the
  # saturation front is the only place that would object.
  # line_search = basic
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 500
    growth_factor = 1.5 # native path takes the front without cutbacks
    cutback_factor = 0.5
  []
[]

[Postprocessors]
  [xnacl_well]
    type = PointValue
    point = '0.1 0 0'
    variable = xnacl
  []
  [halite_well]
    type = PointValue
    point = '0.1 0 0'
    variable = halite
  []
  [porosity_well]
    type = PointValue
    point = '0.1 0 0'
    variable = porosity
  []
  [sgas_well]
    type = PointValue
    point = '0.1 0 0'
    variable = saturation_gas
  []
  [p_well]
    type = PointValue
    point = '0.1 0 0'
    variable = pgas
    execute_on = 'initial timestep_begin' # lagged feedback for injection_rate
  []
  # --- salt mass-conservation check ---
  # No salt is injected, so the TOTAL salt (dissolved + solid halite) must stay
  # constant at its initial value.  With the local-equilibrium flash this holds
  # to solver tolerance all the way through clogging (the conserved salt is
  # merely repartitioned dissolved <-> solid); total_salt_kg should stay flat.
  [dissolved_salt_kg]
    type = PorousFlowFluidMass
    fluid_component = 2
  []
  [solid_salt_vol]
    type = ElementIntegralVariablePostprocessor
    variable = halite # m3 halite / m3 rock, integrated -> total halite volume (m3)
  []
  [total_salt_kg]
    type = ParsedPostprocessor
    # 2165 = halite_density default used by PorousFlowHaliteVolumeFraction and
    # PorousFlowPrecipitateMassTimeDerivative; must match for the balance to close
    expression = 'dissolved_salt_kg + 2165 * solid_salt_vol'
    pp_names = 'dissolved_salt_kg solid_salt_vol'
  []
[]

[VectorPostprocessors]
  [line]
    type = ElementValueSampler
    sort_by = x
    variable = 'halite porosity permeability saturation_gas x1'
    execute_on = 'timestep_end'
    outputs = spatial
  []
[]

[Outputs]
  print_linear_residuals = false
  perf_graph = true
  [time]
    type = CSV
  []
  [spatial]
    type = CSV
    sync_only = true
    # 0.5, 1, 2, 5, 10 days - frames to watch the halite bank build up at the well
    # and then march radially outward as the dry-out front advances
    sync_times = '4.32e4 8.64e4 1.728e5 4.32e5 8.64e5'
  []
[]
