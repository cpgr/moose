# Production well with a scheduled shut-in, opening at t=250 and shutting at t=650.
# A producer's rate comes from its deliverability curve, so it is shut by gating the sink's flux
# through 'multiplying_var' rather than by zeroing a rate function.  While shut, the produced mass,
# rate and cumulative must all read zero rather than holding their last open-step values.

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 4
  ny = 4
[]

[GlobalParams]
  PorousFlowDictator = dictator
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
  [prod1]
    type = PorousFlowProductionWell
    point_file = prod1.bh
    line_length = 1.0
    p_or_t_vals = '19e6 21e6'
    fluxes = '0 1e-4'
    porepressure = pp
    mass_fraction_vars = massfrac0
    start_time = 250
    end_time = 650
  []
[]

[Postprocessors]
  [p_well]
    type = PointValue
    variable = pp
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
  file_base = production_scheduled
[]
