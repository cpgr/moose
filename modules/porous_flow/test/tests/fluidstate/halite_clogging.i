# 1D CO2 dry-out with halite-driven porosity/permeability clogging in PorousFlowBrineCO2.
#
# This is halite_dryout.i with the precipitation-to-porosity feedback switched ON: the solid
# halite volume fraction reported by the native local-equilibrium path is fed into
# PorousFlowPorosity (chemical_equilibrium = true), which subtracts it from the reference porosity,
# and PorousFlowPermeabilityKozenyCarman then chokes the permeability as porosity drops.
#
# As dry CO2 evaporates water from the inlet cells, dissolved salt concentrates past the X_eq =
# 0.2672 solubility, halite precipitates and fills pore space: porosity_min (the inlet) falls from
# 0.2 and permeability_min collapses by orders of magnitude.  The porosity_min floor keeps porosity
# strictly positive so the solve stays out of the singular pore-filling regime.  No salt crosses the
# boundary, so total_salt (dissolved + solid) is conserved throughout, confirming the clogging
# feedback is mass-conservative.

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
  [porosity]
    order = CONSTANT
    family = MONOMIAL
  []
  [permeability]
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
  [porosity]
    type = PorousFlowPropertyAux
    variable = porosity
    property = porosity
    execute_on = 'TIMESTEP_END'
  []
  [permeability]
    type = PorousFlowPropertyAux
    variable = permeability
    property = permeability
    row = 0
    column = 0
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
    type = PorousFlowPorosity
    porosity_zero = 0.2
    chemical_equilibrium = true
    # halite (m^3/m^3) eats pore space; floor keeps porosity positive once pores fill
    porosity_min = 1e-4
  []
  [permeability]
    type = PorousFlowPermeabilityKozenyCarman
    poroperm_function = kozeny_carman_phi0
    k0 = 1e-13
    phi0 = 0.2
    m = 2
    n = 3
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
  [porosity_min]
    type = ElementExtremeValue
    variable = porosity
    value_type = min
  []
  [permeability_min]
    type = ElementExtremeValue
    variable = permeability
    value_type = min
  []
[]

[Outputs]
  csv = true
  execute_on = 'TIMESTEP_END'
[]
