# Scheduled injection well: opens at t=250, shuts at t=650, on a dt=100 timestepper.
# start_time and end_time are inserted into the sync times, so the solve lands exactly on them.
# The injected mass must therefore be exactly rate * (650 - 250) = 1e-4 * 400 = 0.04 kg,
# split 0.3/0.7 between the two components: inj1_cum_0 = 0.012, inj1_cum_1 = 0.028.

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

[Wells]
  [inj1]
    type = PorousFlowInjectionWell
    point = '0.5 0.5 0'
    rate = rate_fn
    injected_mass_fractions = '0.3 0.7'
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
  file_base = injection_scheduled
[]
