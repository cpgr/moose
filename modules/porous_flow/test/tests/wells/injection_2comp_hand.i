# Two-component injection well, built by hand.
# The [Wells] action in injection_2comp_action.i must reproduce this exactly.
# Component 0's mass balance is on massfrac0, component 1's is on pp.

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 4
  ny = 4
[]

[GlobalParams]
  PorousFlowDictator = dictator
[]

[Functions]
  [rate_fn]
    type = ConstantFunction
    value = 1e-4
  []
[]

[Variables]
  [massfrac0]
  []
  [pp]
  []
[]

[ICs]
  [massfrac0]
    type = ConstantIC
    variable = massfrac0
    value = 0.2
  []
  [pp]
    type = ConstantIC
    variable = pp
    value = 20e6
  []
[]

[Kernels]
  [mass0]
    type = PorousFlowMassTimeDerivative
    fluid_component = 0
    variable = massfrac0
  []
  [mass1]
    type = PorousFlowMassTimeDerivative
    fluid_component = 1
    variable = pp
  []
  [flux0]
    type = PorousFlowAdvectiveFlux
    fluid_component = 0
    variable = massfrac0
    gravity = '0 0 0'
  []
  [flux1]
    type = PorousFlowAdvectiveFlux
    fluid_component = 1
    variable = pp
    gravity = '0 0 0'
  []
[]

[UserObjects]
  [dictator]
    type = PorousFlowDictator
    porous_flow_vars = 'massfrac0 pp'
    number_fluid_phases = 1
    number_fluid_components = 2
  []
[]

[FluidProperties]
  [simple_fluid]
    type = SimpleFluidProperties
    bulk_modulus = 2e9
    density0 = 1000
    thermal_expansion = 0
  []
[]

[Materials]
  [temperature]
    type = PorousFlowTemperature
  []
  [ppss]
    type = PorousFlow1PhaseFullySaturated
    porepressure = pp
  []
  [massfrac]
    type = PorousFlowMassFraction
    mass_fraction_vars = massfrac0
  []
  [simple_fluid]
    type = PorousFlowSingleComponentFluid
    fp = simple_fluid
    phase = 0
  []
  [porosity]
    type = PorousFlowPorosityConst
    porosity = 0.2
  []
  [permeability]
    type = PorousFlowPermeabilityConst
    permeability = '1e-13 0 0  0 1e-13 0  0 0 1e-13'
  []
  [relperm]
    type = PorousFlowRelativePermeabilityConst
    phase = 0
  []
[]

[DiracKernels]
  [inj1_source_0]
    type = PorousFlowPointSourceFromPostprocessor
    variable = massfrac0
    mass_flux = inj1_q_0
    point = '0.5 0.5 0'
  []
  [inj1_source_1]
    type = PorousFlowPointSourceFromPostprocessor
    variable = pp
    mass_flux = inj1_q_1
    point = '0.5 0.5 0'
  []
[]

[Postprocessors]
  [inj1_q_0]
    type = FunctionValuePostprocessor
    function = rate_fn
    scale_factor = 0.3
    execute_on = 'initial timestep_begin'
  []
  [inj1_q_1]
    type = FunctionValuePostprocessor
    function = rate_fn
    scale_factor = 0.7
    execute_on = 'initial timestep_begin'
  []
  [inj1_cum_0]
    type = TimeIntegratedPostprocessor
    value = inj1_q_0
  []
  [inj1_cum_1]
    type = TimeIntegratedPostprocessor
    value = inj1_q_1
  []
  [p_well]
    type = PointValue
    variable = pp
    point = '0.5 0.5 0'
    execute_on = 'initial timestep_end'
  []
  [x_well]
    type = PointValue
    variable = massfrac0
    point = '0.5 0.5 0'
    execute_on = 'initial timestep_end'
  []
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  solve_type = Newton
  dt = 100
  end_time = 1000
  nl_abs_tol = 1e-12
[]

[Outputs]
  csv = true
  execute_on = 'initial timestep_end'
  file_base = injection_2comp
[]
