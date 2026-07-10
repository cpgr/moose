# A nonisothermal well ('temperature' supplied) that omits 'fp'.
# The action must error rather than build an enthalpy source it cannot evaluate.

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
    value = 1e-3
  []
  [T_in_fn]
    type = PiecewiseLinear
    xy_data = '
      0    350
      1000 400'
  []
[]

[Variables]
  [pressure]
  []
  [temperature]
    scaling = 1e-6
  []
[]

[ICs]
  [pressure]
    type = ConstantIC
    variable = pressure
    value = 1e6
  []
  [temperature]
    type = ConstantIC
    variable = temperature
    value = 400
  []
[]

[Kernels]
  [mass0]
    type = PorousFlowMassTimeDerivative
    fluid_component = 0
    variable = pressure
  []
  [flux0]
    type = PorousFlowAdvectiveFlux
    fluid_component = 0
    variable = pressure
    gravity = '0 0 0'
  []
  [energy_dot]
    type = PorousFlowEnergyTimeDerivative
    variable = temperature
  []
  [heat_conduction]
    type = PorousFlowHeatConduction
    variable = temperature
  []
  [heat_advection]
    type = PorousFlowHeatAdvection
    variable = temperature
    gravity = '0 0 0'
  []
[]

[UserObjects]
  [dictator]
    type = PorousFlowDictator
    porous_flow_vars = 'pressure temperature'
    number_fluid_phases = 1
    number_fluid_components = 1
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
    temperature = temperature
  []
  [ppss]
    type = PorousFlow1PhaseFullySaturated
    porepressure = pressure
  []
  [massfrac]
    type = PorousFlowMassFraction
  []
  [simple_fluid]
    type = PorousFlowSingleComponentFluid
    fp = simple_fluid
    phase = 0
  []
  [relperm]
    type = PorousFlowRelativePermeabilityConst
    phase = 0
  []
  [rock_heat]
    type = PorousFlowMatrixInternalEnergy
    specific_heat_capacity = 830.0
    density = 2750
  []
  [thermal_conductivity]
    type = PorousFlowThermalConductivityIdeal
    dry_thermal_conductivity = '2.5 0 0  0 2.5 0  0 0 2.5'
  []
  [porosity]
    type = PorousFlowPorosityConst
    porosity = 0.1
  []
  [permeability]
    type = PorousFlowPermeabilityConst
    permeability = '1e-14 0 0  0 1e-14 0  0 0 1e-14'
  []
[]

[Wells]
  [inj1]
    type = PorousFlowInjectionWell
    point = '0.5 0.5 0'
    rate = rate_fn
    injected_mass_fractions = '1'
    porepressure = pressure
    temperature = temperature
    injection_temperature = T_in_fn
  []
[]

[Postprocessors]
  [p_well]
    type = PointValue
    variable = pressure
    point = '0.5 0.5 0'
    execute_on = 'initial timestep_end'
  []
  [T_well]
    type = PointValue
    variable = temperature
    point = '0.5 0.5 0'
    execute_on = 'initial timestep_end'
  []
  [total_heat]
    type = PorousFlowHeatEnergy
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
  file_base = injection_thermal
[]
