# Test Brooks-Corey capillary pressure with heterogeneous entry pressure

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 2
  ny = 2
[]

[GlobalParams]
  PorousFlowDictator = dictator
[]

[Variables]
  [p0]
    initial_condition = 1e6
  []
  [s1]
    initial_condition = 0.1
  []
[]

[AuxVariables]
  [pe]
    family = MONOMIAL
    order = CONSTANT
  []
  [pc]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[AuxKernels]
  [pc]
    type = PorousFlowPropertyAux
    property = capillary_pressure
    variable = pc
  []
[]

[ICs]
  [pe]
    type = RandomIC
    variable = pe
    min = 1e5
    max = 5e5
  []
[]

[Kernels]
  [p0]
    type = Diffusion
    variable = p0
  []
  [s1]
    type = Diffusion
    variable = s1
  []
[]

[UserObjects]
  [dictator]
    type = PorousFlowDictator
    porous_flow_vars = 'p0 s1'
    number_fluid_phases = 2
    number_fluid_components = 2
  []
  [pc]
    type = PorousFlowCapillaryPressureBC
    lambda = 2
    log_extension = false
    pe = pe
    sat_lr = 0
  []
[]

[Materials]
  [temperature]
    type = PorousFlowTemperature
  []
  [ppss]
    type = PorousFlow2PhasePS
    phase0_porepressure = p0
    phase1_saturation = s1
    capillary_pressure = pc
  []
  [kr0]
    type = PorousFlowRelativePermeabilityVG
    phase = 0
    m = 0.5
  []
  [kr1]
    type = PorousFlowRelativePermeabilityCorey
    phase = 1
    n = 2
  []
[]

[Executioner]
  type = Steady
  solve_type = NEWTON
  nl_abs_tol = 1e-6
[]

[Outputs]
  exodus = true
  execute_on = timestep_end
[]
