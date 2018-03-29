# Pressure pulse in 1D with 1 phase - transient
# using the PorousFlowFullySaturatedDarcyBase Kernel
[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 20
  xmin = 0
  xmax = 100
[]

[GlobalParams]
  PorousFlowDictator = dictator
[]

[Variables]
  [./pp]
    initial_condition = 2E6
  [../]
[]

[Kernels]
  [./mass0]
    type = PorousFlowMassTimeDerivative
    fluid_component = 0
    variable = pp
  [../]
  [./flux]
    type = PorousFlowFullySaturatedDarcyBase
    variable = pp
    gravity = '0 0 0'
  [../]
[]

[UserObjects]
  [./dictator]
    type = PorousFlowDictator
    porous_flow_vars = 'pp'
    number_fluid_phases = 1
    number_fluid_components = 1
  [../]
[]

[Modules]
  [./FluidProperties]
    [./simple_fluid]
      type = SimpleFluidProperties
      bulk_modulus = 2e9
      density0 = 1000
      thermal_expansion = 0
      viscosity = 1e-3
    [../]
  [../]
[]

[Materials]
  [./temperature]
    type = PorousFlowTemperature
    at_nodes = true
  [../]
  [./temperature_qp]
    type = PorousFlowTemperature
  [../]
  [./ppss]
    type = PorousFlow1PhaseFullySaturated
    at_nodes = true
    porepressure = pp
  [../]
  [./ppss_qp]
    type = PorousFlow1PhaseFullySaturated
    porepressure = pp
  [../]
  [./massfrac]
    type = PorousFlowMassFraction
    at_nodes = true
  [../]
  [./simple_fluid]
    type = PorousFlowSingleComponentFluid
    fp = simple_fluid
    phase = 0
    at_nodes = true
  [../]
  [./simple_fluid_qp]
    type = PorousFlowSingleComponentFluid
    fp = simple_fluid
    phase = 0
  [../]
  [./dens_all]
    type = PorousFlowJoiner
    at_nodes = true
    material_property = PorousFlow_fluid_phase_density_nodal
  [../]
  [./dens_all_at_quadpoints]
    type = PorousFlowJoiner
    material_property = PorousFlow_fluid_phase_density_qp
    at_nodes = false
  [../]
  [./visc_all]
    type = PorousFlowJoiner
    material_property = PorousFlow_viscosity_qp
  [../]
  [./porosity]
    type = PorousFlowPorosityConst
    at_nodes = true
    porosity = 0.1
  [../]
  [./permeability]
    type = PorousFlowPermeabilityConst
    permeability = '1E-15 0 0 0 1E-15 0 0 0 1E-15'
  [../]
  [./relperm]
    type = PorousFlowRelativePermeabilityCorey
    at_nodes = true
    n = 0
    phase = 0
  [../]
  [./relperm_all]
    type = PorousFlowJoiner
    at_nodes = true
    material_property = PorousFlow_relative_permeability_nodal
  [../]
[]

[BCs]
  [./left]
    type = DirichletBC
    boundary = left
    value = 3E6
    variable = pp
  [../]
[]

[Preconditioning]
  [./andy]
    type = SMP
    full = true
    petsc_options_iname = '-ksp_type -pc_type -snes_atol -snes_rtol -snes_max_it'
    petsc_options_value = 'bcgs bjacobi 1E-15 1E-20 10000'
  [../]
[]

[Executioner]
  type = Transient
  solve_type = Newton
  dt = 1E3
  end_time = 1E4
[]

[Postprocessors]
  [./p005]
    type = PointValue
    variable = pp
    point = '5 0 0'
    execute_on = 'initial timestep_end'
  [../]
  [./p015]
    type = PointValue
    variable = pp
    point = '15 0 0'
    execute_on = 'initial timestep_end'
  [../]
  [./p025]
    type = PointValue
    variable = pp
    point = '25 0 0'
    execute_on = 'initial timestep_end'
  [../]
  [./p035]
    type = PointValue
    variable = pp
    point = '35 0 0'
    execute_on = 'initial timestep_end'
  [../]
  [./p045]
    type = PointValue
    variable = pp
    point = '45 0 0'
    execute_on = 'initial timestep_end'
  [../]
  [./p055]
    type = PointValue
    variable = pp
    point = '55 0 0'
    execute_on = 'initial timestep_end'
  [../]
  [./p065]
    type = PointValue
    variable = pp
    point = '65 0 0'
    execute_on = 'initial timestep_end'
  [../]
  [./p075]
    type = PointValue
    variable = pp
    point = '75 0 0'
    execute_on = 'initial timestep_end'
  [../]
  [./p085]
    type = PointValue
    variable = pp
    point = '85 0 0'
    execute_on = 'initial timestep_end'
  [../]
  [./p095]
    type = PointValue
    variable = pp
    point = '95 0 0'
    execute_on = 'initial timestep_end'
  [../]
[]

[Outputs]
  file_base = pressure_pulse_1d_fully_saturated
  print_linear_residuals = false
  csv = true
[]
