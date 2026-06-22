# 1D CO2 dry-out test for native local-equilibrium halite precipitation in PorousFlowBrineCO2.
#
# A short brine-saturated column (salt mass fraction z_s = 0.2, just below the X_eq = 0.2672
# solubility at 45 C). Dry CO2 is injected at the left: it displaces and evaporates water from the
# inlet cells, concentrating the dissolved salt until z_s exceeds solubility, at which point solid
# halite precipitates there (the textbook CO2 "dry-out").  This is the spatial-front complement to
# the closed-cell halite_precipitation_cycle.i: here the saturation-driven crossing of z_s = X_eq
# happens cell by cell as the dry-out front advances from the inlet.
#
# Only CO2 crosses the boundary -- no salt -- so the TOTAL salt (dissolved + solid) must stay
# constant: total_salt holds while solid_salt rises and dissolved_salt falls.  Porosity is constant
# (the native precipitation path does not wire clogging), so the dry-out stays away from the
# near-singular pore-filling regime and the solve is robust.

[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 10
  xmin = 0
  xmax = 2
  bias_x = 1.3 # refine toward the injection inlet where dry-out occurs
[]

[GlobalParams]
  PorousFlowDictator = dictator
  gravity = '0 0 0'
[]

[Variables]
  [pgas]
    initial_condition = 12e6
  []
  [zi]
    initial_condition = 0
    scaling = 1e4
  []
  [xnacl]
    initial_condition = 0.2
  []
[]

[AuxVariables]
  [halite_vf]
    order = CONSTANT
    family = MONOMIAL
  []
  [saturation_gas]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[AuxKernels]
  [halite_vf]
    type = MaterialRealAux
    variable = halite_vf
    property = PorousFlow_halite_volume_fraction_nodal
    execute_on = 'TIMESTEP_END'
  []
  [saturation_gas]
    type = PorousFlowPropertyAux
    variable = saturation_gas
    property = saturation
    phase = 1
    execute_on = 'TIMESTEP_END'
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
  [mass2]
    type = PorousFlowMassTimeDerivative
    fluid_component = 2
    variable = xnacl
  []
  [flux2]
    type = PorousFlowAdvectiveFlux
    fluid_component = 2
    variable = xnacl
  []
  [precipitate2]
    type = PorousFlowPrecipitateMassTimeDerivative
    variable = xnacl
  []
[]

[UserObjects]
  [dictator]
    type = PorousFlowDictator
    porous_flow_vars = 'pgas zi xnacl'
    number_fluid_phases = 2
    number_fluid_components = 3
  []
  [pc]
    type = PorousFlowCapillaryPressureConst
    pc = 0
  []
  [fs]
    type = PorousFlowBrineCO2
    brine_fp = brine
    co2_fp = co2
    capillary_pressure = pc
    salt_component = 2
    precipitate_salt = true
  []
[]

[FluidProperties]
  [co2sw]
    type = CO2FluidProperties
  []
  [co2]
    type = TabulatedBicubicFluidProperties
    input_fp = co2sw
  []
  [water]
    type = Water97FluidProperties
  []
  [brine]
    type = BrineFluidProperties
    water_fp = water
  []
[]

[Materials]
  [temperature]
    type = PorousFlowTemperature
    temperature = 45
  []
  [brineco2]
    type = PorousFlowFluidState
    gas_porepressure = pgas
    z = zi
    temperature_unit = Celsius
    xnacl = xnacl
    capillary_pressure = pc
    fluid_state = fs
  []
  [halite]
    type = PorousFlowHaliteVolumeFraction
  []
  [porosity]
    type = PorousFlowPorosityConst
    porosity = 0.2
  []
  [permeability]
    type = PorousFlowPermeabilityConst
    permeability = '1e-13 0 0 0 1e-13 0 0 0 1e-13'
  []
  [relperm_water]
    type = PorousFlowRelativePermeabilityCorey
    n = 2
    phase = 0
    s_res = 0.2
    sum_s_res = 0.25
  []
  [relperm_gas]
    type = PorousFlowRelativePermeabilityCorey
    n = 2
    phase = 1
    s_res = 0.05
    sum_s_res = 0.25
  []
[]

[BCs]
  [right_pressure]
    type = DirichletBC
    boundary = right
    variable = pgas
    value = 12e6
  []
  [inject_co2]
    type = FunctionNeumannBC
    boundary = left
    variable = zi
    function = co2_rate
  []
[]

[Functions]
  [co2_rate]
    type = ParsedFunction
    expression = 2e-3
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
  solve_type = NEWTON
  end_time = 3e4
  nl_abs_tol = 1e-8
  dtmax = 2.5e3
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 100
    growth_factor = 1.2
    cutback_factor = 0.5
  []
[]

[Postprocessors]
  [dissolved_salt]
    type = PorousFlowFluidMass
    fluid_component = 2
    phase = '0 1'
  []
  [halite_volume]
    type = ElementIntegralVariablePostprocessor
    variable = halite_vf
  []
  [solid_salt]
    type = ParsedPostprocessor
    pp_names = 'halite_volume'
    expression = '2165*halite_volume'
  []
  [total_salt]
    type = ParsedPostprocessor
    pp_names = 'dissolved_salt solid_salt'
    expression = 'dissolved_salt + solid_salt'
  []
[]

[Outputs]
  csv = true
  execute_on = 'TIMESTEP_END'
[]
