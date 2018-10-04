# Tests correct calculation of properties in PorousFlowBrineCO2CH4H2S

[Mesh]
  type = GeneratedMesh
  dim = 2
[]

[GlobalParams]
  PorousFlowDictator = dictator
  temperature = 30
[]

[Variables]
  [./pgas]
    initial_condition = 20e6
  [../]
  [./zco2]
     initial_condition = 0.2
  [../]
  [./zch4]
     initial_condition = 0.1
  [../]
  [./zh2s]
     initial_condition = 0.05
  [../]
[]

[AuxVariables]
  [./xnacl]
    initial_condition = 0.1
  [../]
  [./pressure_gas]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./pressure_water]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./saturation_gas]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./saturation_water]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./density_water]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./density_gas]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./viscosity_water]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./viscosity_gas]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./x0_water]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./x0_gas]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./x1_water]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./x1_gas]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./x2_water]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./x2_gas]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./x3_water]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./x3_gas]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[AuxKernels]
  [./pressure_water]
    type = PorousFlowPropertyAux
    variable = pressure_water
    property = pressure
    phase = 0
    execute_on = timestep_end
  [../]
  [./pressure_gas]
    type = PorousFlowPropertyAux
    variable = pressure_gas
    property = pressure
    phase = 1
    execute_on = timestep_end
  [../]
  [./saturation_water]
    type = PorousFlowPropertyAux
    variable = saturation_water
    property = saturation
    phase = 0
    execute_on = timestep_end
  [../]
  [./saturation_gas]
    type = PorousFlowPropertyAux
    variable = saturation_gas
    property = saturation
    phase = 1
    execute_on = timestep_end
  [../]
  [./density_water]
    type = PorousFlowPropertyAux
    variable = density_water
    property = density
    phase = 0
    execute_on = timestep_end
  [../]
  [./density_gas]
    type = PorousFlowPropertyAux
    variable = density_gas
    property = density
    phase = 1
    execute_on = timestep_end
  [../]
  [./viscosity_water]
    type = PorousFlowPropertyAux
    variable = viscosity_water
    property = viscosity
    phase = 0
    execute_on = timestep_end
  [../]
  [./viscosity_gas]
    type = PorousFlowPropertyAux
    variable = viscosity_gas
    property = viscosity
    phase = 1
    execute_on = timestep_end
  [../]
  [./x1_water]
    type = PorousFlowPropertyAux
    variable = x1_water
    property = mass_fraction
    phase = 0
    fluid_component = 1
    execute_on = timestep_end
  [../]
  [./x1_gas]
    type = PorousFlowPropertyAux
    variable = x1_gas
    property = mass_fraction
    phase = 1
    fluid_component = 1
    execute_on = timestep_end
  [../]
  [./x0_water]
    type = PorousFlowPropertyAux
    variable = x0_water
    property = mass_fraction
    phase = 0
    fluid_component = 0
    execute_on = timestep_end
  [../]
  [./x0_gas]
    type = PorousFlowPropertyAux
    variable = x0_gas
    property = mass_fraction
    phase = 1
    fluid_component = 0
    execute_on = timestep_end
  [../]
  [./x2_water]
    type = PorousFlowPropertyAux
    variable = x2_water
    property = mass_fraction
    phase = 0
    fluid_component = 2
    execute_on = timestep_end
  [../]
  [./x2_gas]
    type = PorousFlowPropertyAux
    variable = x2_gas
    property = mass_fraction
    phase = 1
    fluid_component = 2
    execute_on = timestep_end
  [../]
  [./x3_water]
    type = PorousFlowPropertyAux
    variable = x3_water
    property = mass_fraction
    phase = 0
    fluid_component = 3
    execute_on = timestep_end
  [../]
  [./x3_gas]
    type = PorousFlowPropertyAux
    variable = x3_gas
    property = mass_fraction
    phase = 1
    fluid_component = 3
    execute_on = timestep_end
  [../]
[]

[Kernels]
  [./mass0]
    type = PorousFlowMassTimeDerivative
    variable = pgas
    fluid_component = 0
  [../]
  [./mass1]
    type = PorousFlowMassTimeDerivative
    variable = zco2
    fluid_component = 1
  [../]
  [./mass2]
    type = PorousFlowMassTimeDerivative
    variable = zch4
    fluid_component = 2
  [../]
  [./mass3]
    type = PorousFlowMassTimeDerivative
    variable = zh2s
    fluid_component = 3
  [../]
[]

[UserObjects]
  [./dictator]
    type = PorousFlowDictator
    porous_flow_vars = 'pgas zco2 zch4 zh2s'
    number_fluid_phases = 2
    number_fluid_components = 4
  [../]
  [./pc]
    type = PorousFlowCapillaryPressureConst
    pc = 0
  [../]
  [./fs]
    type = PorousFlowBrineCO2CH4H2S
    brine_fp = brine
    co2_fp = co2
    ch4_fp = ch4
    capillary_pressure = pc
  [../]
[]

[Modules]
  [./FluidProperties]
    [./co2]
      type = CO2FluidProperties
    [../]
    [./brine]
      type = BrineFluidProperties
    [../]
    [./ch4]
      type = MethaneFluidProperties
    [../]
  [../]
[]

[Materials]
  [./temperature]
    type = PorousFlowTemperature
  [../]
  [./brineco2ch4h2s]
    type = PorousFlowFluidState
    gas_porepressure = pgas
    z = 'zco2 zch4 zh2s'
    temperature_unit = Celsius
    xnacl = xnacl
    capillary_pressure = pc
    fluid_state = fs
  [../]
  [./permeability]
    type = PorousFlowPermeabilityConst
    permeability = '1e-12 0 0 0 1e-12 0 0 0 1e-12'
  [../]
  [./relperm0]
    type = PorousFlowRelativePermeabilityCorey
    n = 2
    phase = 0
  [../]
  [./relperm1]
    type = PorousFlowRelativePermeabilityCorey
    n = 3
    phase = 1
  [../]
  [./porosity]
    type = PorousFlowPorosityConst
    porosity = 0.1
  [../]
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  dt = 1
  end_time = 1
  nl_abs_tol = 1e-12
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
  [../]
[]

[Postprocessors]
  [./density_water]
    type = ElementIntegralVariablePostprocessor
    variable = density_water
  [../]
  [./density_gas]
    type = ElementIntegralVariablePostprocessor
    variable = density_gas
  [../]
  [./viscosity_water]
    type = ElementIntegralVariablePostprocessor
    variable = viscosity_water
  [../]
  [./viscosity_gas]
    type = ElementIntegralVariablePostprocessor
    variable = viscosity_gas
  [../]
  [./xco2]
    type = ElementIntegralVariablePostprocessor
    variable = x1_water
  [../]
  [./xh2o]
    type = ElementIntegralVariablePostprocessor
    variable = x0_water
  [../]
  [./xch4]
    type = ElementIntegralVariablePostprocessor
    variable = x2_water
  [../]
  [./xh2s]
    type = ElementIntegralVariablePostprocessor
    variable = x3_water
  [../]
  [./yco2]
    type = ElementIntegralVariablePostprocessor
    variable = x1_gas
  [../]
  [./yh2o]
    type = ElementIntegralVariablePostprocessor
    variable = x0_gas
  [../]
  [./ych4]
    type = ElementIntegralVariablePostprocessor
    variable = x2_gas
  [../]
  [./yh2s]
    type = ElementIntegralVariablePostprocessor
    variable = x3_gas
  [../]
  [./sg]
    type = ElementIntegralVariablePostprocessor
    variable = saturation_gas
  [../]
  [./sw]
    type = ElementIntegralVariablePostprocessor
    variable = saturation_water
  [../]
  [./pwater]
    type = ElementIntegralVariablePostprocessor
    variable = pressure_water
  [../]
  [./pgas]
    type = ElementIntegralVariablePostprocessor
    variable = pressure_gas
  [../]
  [./x0mass]
    type = PorousFlowFluidMass
    fluid_component = 0
    phase = '0 1'
  [../]
  [./x1mass]
    type = PorousFlowFluidMass
    fluid_component = 1
    phase = '0 1'
  [../]
[]

[Outputs]
  csv = true
  execute_on = 'TIMESTEP_END'
  perf_graph = false
[]
