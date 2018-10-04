# Two phase density-driven convection of dissolved CO2 in brine
#
# Initially, the model has a gas phase at the top with a saturation of 0.1
# (which corresponds to an initial value of zi = 0.096).
# Diffusion of the dissolved CO2
# component from the saturated liquid to the unsaturated liquid below reduces the
# amount of CO2 in the gas phase. As the density of the CO2-saturated brine is greater
# than the unsaturated brine, a gravitational instability arises and density-driven
# convection of CO2-rich fingers descend into the unsaturated brine.
#
# The instability is seeded by a random perturbation to the porosity field.
# Mesh adaptivity is used to refine the mesh as the fingers form.
inactive = 'MeshModifiers'
[GlobalParams]
  PorousFlowDictator = 'dictator'
  gravity = '0 -9.81 0'
[]

[Adaptivity]
  max_h_level = 2
  marker = marker
  [Indicators]
    [indicator]
      type = GradientJumpIndicator
      variable = zi
    []
  []
  [Markers]
    [marker]
      type = ErrorFractionMarker
      indicator = indicator
      refine = 0.8
    []
  []
[]

[Mesh]
  type = GeneratedMesh
  dim = 2
  ymax = 1
  xmax = 1
  ny = 20
  nx = 20
  bias_y = 1
  elem_type = QUAD4
[]

[AuxVariables]
  [xnacl]
    initial_condition = 0.01
  []
  [saturation_gas]
    order = FIRST
    family = MONOMIAL
  []
  [xco2l]
    order = FIRST
    family = MONOMIAL
  []
  [density_liquid]
    order = FIRST
    family = MONOMIAL
  []
  [porosity]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[AuxKernels]
  [saturation_gas]
    type = PorousFlowPropertyAux
    variable = saturation_gas
    property = saturation
    phase = 1
    execute_on = 'timestep_end'
  []
  [xco2l]
    type = PorousFlowPropertyAux
    variable = xco2l
    property = mass_fraction
    phase = 0
    fluid_component = 1
    execute_on = 'timestep_end'
  []
  [density_liquid]
    type = PorousFlowPropertyAux
    variable = density_liquid
    property = density
    phase = 0
    execute_on = 'timestep_end'
  []
[]

[Variables]
  [pgas]
    scaling = 1
  []
  [zi]
    scaling = 1e4
  []
[]

[ICs]
  [pressure]
    type = FunctionIC
    function = 10e6-9.81*1000*y
    variable = pgas
  []
  [zi]
    type = BoundingBoxIC
    variable = zi
    x1 = 0
    x2 = 1
    y1 = 0.9
    y2 = 1
    inside = 0.06
    outside = 0
  []
  [porosity]
    type = RandomIC
    legacy_generator = false
    variable = porosity
    min = 0.25
    max = 0.275
  []
[]

[Kernels]
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
  [diff0]
    type = PorousFlowDispersiveFlux
    fluid_component = 0
    variable = pgas
    disp_long = '0 0'
    disp_trans = '0 0'
  []
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
  [diff1]
    type = PorousFlowDispersiveFlux
    fluid_component = 1
    variable = zi
    disp_long = '0 0'
    disp_trans = '0 0'
  []
[]

[UserObjects]
  [dictator]
    type = PorousFlowDictator
    porous_flow_vars = 'pgas zi'
    number_fluid_phases = 2
    number_fluid_components = 2
  []
  [pc]
    type = PorousFlowCapillaryPressureConst
    pc = 0
    sat_lr = 0.1
    pc_max = 1e+06
  []
  [fs]
    type = PorousFlowBrineCO2
    brine_fp = brine
    co2_fp = co2
    capillary_pressure = pc
  []
[]

[Modules]
  [FluidProperties]
    [co2sw]
      type = CO2FluidProperties
    []
    [co2]
      type = TabulatedFluidProperties
      fp = co2sw
    []
    [brine]
      type = BrineFluidProperties
    []
  []
[]

[Materials]
  [temperature]
    type = PorousFlowTemperature
    temperature = '45'
  []
  [brineco2]
    type = PorousFlowFluidStateBrineCO2
    gas_porepressure = 'pgas'
    z = 'zi'
    temperature_unit = Celsius
    xnacl = 'xnacl'
    capillary_pressure = pc
    fluid_state = fs
  []
  [porosity]
    type = PorousFlowPorosityConst
    porosity = 'porosity'
  []
  [permeability]
    type = PorousFlowPermeabilityConst
    permeability = '4e-12 0 0 0 4e-12 0 0 0 4e-12'
  []
  [relperm_water]
    type = PorousFlowRelativePermeabilityCorey
    phase = 0
    n = 2
    s_res = 0
    sum_s_res = 0
  []
  [relperm_gas]
    type = PorousFlowRelativePermeabilityCorey
    phase = 1
    n = 2
    s_res = 0.1
    sum_s_res = 0.1
  []
  [diffusivity]
    type = PorousFlowDiffusivityConst
    diffusion_coeff = '2e-9 2e-9 0 0'
    tortuosity = '1 1'
  []
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
    petsc_options_iname = '-ksp_type -pc_type -sub_pc_type -sub_pc_factor_shift_type'
    petsc_options_value = 'gmres asm lu NONZERO'
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  end_time = 5e6
  nl_max_its = 25
  l_max_its = 100
  nl_abs_tol = 1e-6
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 10
    growth_factor = 1.25
    cutback_factor = 0.5
  []
[]

[Postprocessors]
  [total_co2_in_gas]
    type = PorousFlowFluidMass
    phase = '1'
    fluid_component = 1
  []
  [total_co2_in_liquid]
    type = PorousFlowFluidMass
    phase = '0'
    fluid_component = 1
  []
  [numdofs]
    type = NumDOFs
  []
[]

[Outputs]
  print_linear_residuals = false
  perf_graph = true
  exodus = true
  csv = true
[]

[BCs]
  [Periodic]
    [periodic]
      variable = 'pgas zi'
      auto_direction = 'x'
    []
  []
[]

[MeshModifiers]
  [cap]
    type = SubdomainBoundingBox
    block_id = 1
    top_right = '1 1 0'
    bottom_left = '0 0.98 0'
  []
[]
