# Injection well shut in by an event trigger rather than by a schedule.  The shut-in time is not
# known in advance, so it cannot be a sync time: the well shuts one timestep after the condition is
# first met, since the postprocessor is read from the end of the previous step.

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
  # Shut the well once the well-node pressure exceeds 20.2 MPa.  ParsedFunction resolves 'p_well'
  # as a Postprocessor, so the trigger depends on the state of the simulation, not on t.
  [too_hot_fn]
    type = ParsedFunction
    expression = 'pw > 20.2e6'
    symbol_names = 'pw'
    symbol_values = 'p_well'
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

[Wells]
  [inj1]
    type = PorousFlowInjectionWell
    point = '0.5 0.5 0'
    rate = rate_fn
    injected_mass_fractions = '0.3 0.7'
    porepressure = pp
    mass_fraction_vars = massfrac0
    shut_if = too_hot_fn
  []
[]

[Postprocessors]
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
  end_time = 2000
  nl_abs_tol = 1e-12
[]

[Outputs]
  csv = true
  execute_on = 'initial timestep_end'
  file_base = injection_shut_if
[]
