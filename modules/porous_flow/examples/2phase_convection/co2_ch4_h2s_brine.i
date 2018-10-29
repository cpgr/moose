# Two phase density-driven convection of dissolved CO2, CH4 and H2S in brine
#
# Initially, the model has a gas phase at the top with a saturation of 0.1
# The instability is seeded by a random perturbation to the porosity field.
# Mesh adaptivity is used to refine the mesh as the fingers form.

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
      variable = zco2
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
    initial_condition = 0
  []
  [saturation_gas]
    order = FIRST
    family = MONOMIAL
  []
  [xco2l]
    order = FIRST
    family = MONOMIAL
  []
  [xch4l]
    order = FIRST
    family = MONOMIAL
  []
  [xh2sl]
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
  [xch4l]
    type = PorousFlowPropertyAux
    variable = xch4l
    property = mass_fraction
    phase = 0
    fluid_component = 2
    execute_on = 'timestep_end'
  []
  [xh2sl]
    type = PorousFlowPropertyAux
    variable = xh2sl
    property = mass_fraction
    phase = 0
    fluid_component = 3
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
  [zco2]
    scaling = 1e4
  []
  [zch4]
    scaling = 1e4
  []
  [zh2s]
    scaling = 1e4
  []
[]

[ICs]
  [pressure]
    type = FunctionIC
    function = 10e6-9.81*1000*y
    variable = pgas
  []
  [zco2]
    type = BoundingBoxIC
    variable = zco2
    x1 = 0
    x2 = 1
    y1 = 0.95
    y2 = 1
    inside = 0.1
    outside = 0
  []
  [zch4]
    type = BoundingBoxIC
    variable = zch4
    x1 = 0
    x2 = 1
    y1 = 0.95
    y2 = 1
    inside = 0.1
    outside = 0
  []
  [zh2s]
    type = BoundingBoxIC
    variable = zh2s
    x1 = 0
    x2 = 1
    y1 = 0.95
    y2 = 1
    inside = 0.01
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
    variable = zco2
  []
  [flux1]
    type = PorousFlowAdvectiveFlux
    fluid_component = 1
    variable = zco2
  []
  [diff1]
    type = PorousFlowDispersiveFlux
    fluid_component = 1
    variable = zco2
    disp_long = '0 0'
    disp_trans = '0 0'
  []
  [mass2]
    type = PorousFlowMassTimeDerivative
    fluid_component = 2
    variable = zch4
  []
  [flux2]
    type = PorousFlowAdvectiveFlux
    fluid_component = 2
    variable = zch4
  []
  [diff2]
    type = PorousFlowDispersiveFlux
    fluid_component = 2
    variable = zch4
    disp_long = '0 0'
    disp_trans = '0 0'
  []
  [mass3]
    type = PorousFlowMassTimeDerivative
    fluid_component = 3
    variable = zh2s
  []
  [flux3]
    type = PorousFlowAdvectiveFlux
    fluid_component = 3
    variable = zh2s
  []
  [diff3]
    type = PorousFlowDispersiveFlux
    fluid_component = 3
    variable = zh2s
    disp_long = '0 0'
    disp_trans = '0 0'
  []
[]

[UserObjects]
  [dictator]
    type = PorousFlowDictator
    porous_flow_vars = 'pgas zco2 zch4 zh2s'
    number_fluid_phases = 2
    number_fluid_components = 4
  []
  [pc]
    type = PorousFlowCapillaryPressureConst
    pc = 0
    sat_lr = 0.1
    pc_max = 1e+06
  []
  [fs]
    type = PorousFlowBrineCO2CH4H2S
    brine_fp = brine
    co2_fp = co2
    ch4_fp = ch4
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
    [ch4]
      type = MethaneFluidProperties
    []
  []
[]

[Materials]
  [temperature]
    type = PorousFlowTemperature
    temperature = '45'
  []
  [brineco2ch4h2s]
    type = PorousFlowFluidState
    gas_porepressure = 'pgas'
    z = 'zco2 zch4 zh2s'
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
    diffusion_coeff = '2e-8 2e-8 0 0 0 0 0 0'
    tortuosity = '1 1'
  []
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
    petsc_options = '-ksp_converged_reason -snes_converged_reason'
    petsc_options_iname = '-ksp_type -pc_type -sub_pc_type -sub_pc_factor_shift_type'
    petsc_options_value = 'gmres asm lu NONZERO'
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  end_time = 5e6
  nl_max_its = 20
  l_max_its = 100
  nl_abs_tol = 1e-6
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 10
    growth_factor = 1.5
    cutback_factor = 0.75
  []
[]

[Functions]
  [./flux]
    type = ParsedFunction
    vals = 'delta_xco2 dt'
    vars = 'dx dt'
    value = 'dx/dt'
  [../]
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
  [total_ch4_in_gas]
    type = PorousFlowFluidMass
    phase = '1'
    fluid_component = 2
  []
  [total_ch4_in_liquid]
    type = PorousFlowFluidMass
    phase = '0'
    fluid_component = 2
  []
  [total_h2s_in_gas]
    type = PorousFlowFluidMass
    phase = '1'
    fluid_component = 3
  []
  [total_h2s_in_liquid]
    type = PorousFlowFluidMass
    phase = '0'
    fluid_component = 3
  []
  [total_liquid]
    type = PorousFlowFluidMass
    phase = '0'
    fluid_component = 0
  []
  [delta_xco2]
    type = ChangeOverTimePostprocessor
    postprocessor = total_co2_in_liquid
  []
  [dt]
    type = TimestepSize
  []
  [lux]
    type = FunctionValuePostprocessor
    function = flux
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
  [top]
    type = PresetBC
    variable = pgas
    value = 10e6
    boundary = bottom
  []
  [Periodic]
    [periodic]
      variable = 'pgas zco2 zch4 zh2s'
      auto_direction = 'x'
    []
  []
[]
