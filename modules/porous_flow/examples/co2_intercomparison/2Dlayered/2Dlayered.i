# Intercomparison problem 7: CO2 Injection into a 2D layered brine formation
#
# From Pruess et al, Code intercomparison builds confidence in
# numerical simulation models for geologic disposal of CO2, Energy 29 (2004)

[GlobalParams]
  PorousFlowDictator = 'dictator'
  gravity = '0 -9.81 0'
[]

[Mesh]
  type = GeneratedMesh
  dim = 2
  ymax = 183
  xmax = 6000
  nx = 200
  ny = 61
  bias_x = 1.01
[]

# [Adaptivity]
#   max_h_level = 2
#   marker = marker
#   initial_marker = initial
#   initial_steps = 2
#   [Indicators]
#     [indicator]
#       type = GradientJumpIndicator
#       variable = zi
#     []
#   []
#   [Markers]
#     [marker]
#       type = ErrorFractionMarker
#       indicator = indicator
#       refine = 0.5
#     []
#     [initial]
#       type = BoxMarker
#       bottom_left = '0 15 0'
#       top_right = '30 30 0'
#       inside = REFINE
#       outside = DO_NOTHING
#     []
#   []
# []

[MeshModifiers]
  [shale1]
    type = SubdomainBoundingBox
    block_id = 1
    bottom_left = '0 51 0'
    top_right = '6000 54 0'
  []
  [shale2]
    type = SubdomainBoundingBox
    block_id = 1
    bottom_left = '0 84 0'
    top_right = '6000 87 0'
  []
  [shale3]
    type = SubdomainBoundingBox
    block_id = 1
    bottom_left = '0 117 0'
    top_right = '6000 120 0'
  []
  [shale4]
    type = SubdomainBoundingBox
    block_id = 1
    bottom_left = '0 150 0'
    top_right = '6000 153 0'
  []
  [renameblocks]
    type = RenameBlock
    old_block_id = '0 1'
    new_block_name = 'aquifer shale'
    depends_on = 'shale1 shale2 shale3 shale4'
  []
[]

[AuxVariables]
  [saturation_gas]
    order = FIRST
    family = MONOMIAL
  []
  [x1]
    order = FIRST
    family = MONOMIAL
  []
  [y0]
    order = FIRST
    family = MONOMIAL
  []
  [xnacl]
    initial_condition = 0.032
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
  [x1]
    type = PorousFlowPropertyAux
    variable = x1
    property = mass_fraction
    phase = 0
    fluid_component = 1
    execute_on = 'timestep_end'
  []
  [y0]
    type = PorousFlowPropertyAux
    variable = y0
    property = mass_fraction
    phase = 1
    fluid_component = 0
    execute_on = 'timestep_end'
  []
[]

[Variables]
  [pgas]
  []
  [zi]
    initial_condition = 0
  []
[]

[ICs]
  [pressure]
    type = FunctionIC
    function = 11e6-9.81*1000*y
    variable = pgas
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
[]

[UserObjects]
  [dictator]
    type = PorousFlowDictator
    porous_flow_vars = 'pgas zi'
    number_fluid_phases = 2
    number_fluid_components = 3
  []
  [pc_aquifer]
    type = PorousFlowCapillaryPressureVG
    alpha = 2.857e-4
    m = 0.4
    sat_lr = 0.2
    pc_max = 1e7
  []
  [pc_shale]
    type = PorousFlowCapillaryPressureVG
    alpha = 1.613e-5
    m = 0.4
    sat_lr = 0.2
    pc_max = 1e7
  []
  [fs_aquifer]
    type = PorousFlowBrineCO2
    brine_fp = brine
    co2_fp = co2
    capillary_pressure = pc_aquifer
  []
  [fs_shale]
    type = PorousFlowBrineCO2
    brine_fp = brine
    co2_fp = co2
    capillary_pressure = pc_shale
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
      save_file = false
    []
    [water]
      type = Water97FluidProperties
    []
    [watertab]
      type = TabulatedFluidProperties
      fp = water
      save_file = false
    []
    [brine]
      type = BrineFluidProperties
      water_fp = watertab
    []
  []
[]

[Materials]
  [temperature]
    type = PorousFlowTemperature
    temperature = '37'
  []
  [brineco2_aquifer]
    type = PorousFlowFluidState
    gas_porepressure = 'pgas'
    z = 'zi'
    temperature_unit = Celsius
    xnacl = 'xnacl'
    capillary_pressure = pc_aquifer
    fluid_state = fs_aquifer
    block = 'aquifer'
  []
  [brineco2_shale]
    type = PorousFlowFluidState
    gas_porepressure = 'pgas'
    z = 'zi'
    temperature_unit = Celsius
    xnacl = 'xnacl'
    capillary_pressure = pc_shale
    fluid_state = fs_shale
    block = 'shale'
  []
  [porosity_aquifer]
    type = PorousFlowPorosityConst
    porosity = '0.35'
    block = 'aquifer'
  []
  [permeability_aquifer]
    type = PorousFlowPermeabilityConst
    permeability = '3e-12 0 0 0 3e-12 0 0 0 3e-12'
    block = 'aquifer'
  []
  [porosity_shale]
    type = PorousFlowPorosityConst
    porosity = '0.1025'
    block = 'shale'
  []
  [permeability_shale]
    type = PorousFlowPermeabilityConst
    permeability = '1e-14 0 0 0 1e-14 0 0 0 1e-14'
    block = 'shale'
  []
  [relperm_water]
    type = PorousFlowRelativePermeabilityVG
    m = 0.4
    phase = 0
    s_res = 0.2
    sum_s_res = 0.25
  []
  [relperm_gas]
    type = PorousFlowRelativePermeabilityVG
    m = 0.4
    phase = 1
    s_res = 0.05
    sum_s_res = 0.25
  []
[]

[BCs]
  [rightwater]
    type = PorousFlowPiecewiseLinearSink
    boundary = 'right'
    variable = pgas
    use_mobility = true
    PorousFlowDictator = dictator
    fluid_phase = 0
    multipliers = '0 1e9'
    PT_shift = '11e6'
    pt_vals = '0 1e9'
    mass_fraction_component = 0
    use_relperm = true
  []
  [rightco2]
    type = PorousFlowPiecewiseLinearSink
    variable = zi
    boundary = 'right'
    use_mobility = true
    PorousFlowDictator = dictator
    fluid_phase = 1
    multipliers = '0 1e9'
    PT_shift = '11e6'
    pt_vals = '0 1e9'
    mass_fraction_component = 1
    use_relperm = true
  []
[]

# [DiracKernels]
#   [source]
#     type = PorousFlowSquarePulsePointSource
#     point = '0 22 0'
#     mass_flux = 0.1585
#     variable = zi
#   []
# []

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
  end_time = 6.31152e7
  nl_max_its = 25
  l_max_its = 100
  automatic_scaling = true
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 1000
    growth_factor = 1.5
    cutback_factor = 0.75
    optimal_iterations = 6
    iteration_window = 1
  []
[]

[Postprocessors]
  [massgas]
    type = PorousFlowFluidMass
    fluid_component = 1
  []
  [xco2l]
    type = PorousFlowFluidMass
    fluid_component = 1
    phase = '0'
  []
[]

[Outputs]
  print_linear_residuals = false
  perf_graph = true
  exodus = true
  csv = true
  sync_times = '2.592e6 3.15576e7 6.31152e7'
  checkpoint = true
[]
