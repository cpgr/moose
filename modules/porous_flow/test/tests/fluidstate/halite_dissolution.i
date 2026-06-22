# Halite dissolution test for native local-equilibrium precipitation in PorousFlowBrineCO2.
#
# A single-phase aqueous column (no CO2: zi = 0) starts OVERSATURATED everywhere: the total-salt
# variable z_s = 0.32 exceeds the X_eq = 0.2672 solubility at 45 C, so at t = 0 every cell holds
# dissolved salt at saturation PLUS the excess as solid halite. The oversaturated initial condition
# is made consistent by PorousFlowFluidState publishing the precipitated-salt mass fraction at t = 0
# and PorousFlowHaliteVolumeFraction seeding the initial solid from it, so the initial total salt
# equals the inventory z_s represents (nothing is lost at the first step).
#
# Pure water (component 0) is injected at the left and fluid leaves through advective outflow
# boundaries (PorousFlowOutflowBC) on both the water and salt components at the right. With no CO2
# present there is nothing to dry the cells out and counteract dissolution (the failure mode of CO2
# dry-out): the fresh water sweeps in from the inlet, lowers z_s below solubility cell by cell, and
# the solid halite redissolves; the released salt is carried out of the column with the throughflow.
# The diagnostic is the solid halite mass solid_salt, which falls steadily from its oversaturated
# initial value as the dissolution front advances - confirming that injected water dissolves halite
# cleanly once CO2 dry-out is removed.

[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 20
  xmin = 0
  xmax = 2
[]

[GlobalParams]
  PorousFlowDictator = dictator
  gravity = '0 0 0'
[]

[Variables]
  [pgas]
    initial_condition = 15e6
  []
  [zi]
    initial_condition = 0
    scaling = 1e4
  []
  [xnacl]
    initial_condition = 0.32 # z_s > X_eq(45 C) = 0.2672 -> starts oversaturated (solid present)
  []
[]

[AuxVariables]
  [halite_vf]
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
  [co2]
    type = CO2FluidProperties
  []
  [brine]
    type = BrineFluidProperties
  []
  [water]
    type = Water97FluidProperties
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
    porosity = 0.1
  []
  [permeability]
    type = PorousFlowPermeabilityConst
    permeability = '1e-12 0 0 0 1e-12 0 0 0 1e-12'
  []
  [relperm0]
    type = PorousFlowRelativePermeabilityCorey
    n = 2
    phase = 0
  []
  [relperm1]
    type = PorousFlowRelativePermeabilityCorey
    n = 3
    phase = 1
  []
[]

[BCs]
  # Inject pure water (component 0) at the inlet to dilute the salt and dissolve the halite
  [inject_water]
    type = FunctionNeumannBC
    boundary = left
    variable = pgas
    function = water_rate
  []
  # Advective outflow of water and salt at the right, so the dissolved salt is carried out of the
  # column rather than re-precipitating at a trapped boundary
  [outflow_water]
    type = PorousFlowOutflowBC
    boundary = right
    variable = pgas
    mass_fraction_component = 0
  []
  [outflow_salt]
    type = PorousFlowOutflowBC
    boundary = right
    variable = xnacl
    mass_fraction_component = 2
  []
[]

[Functions]
  [water_rate]
    type = ParsedFunction
    expression = 5e-3
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
  end_time = 1e4
  nl_abs_tol = 1e-9
  dtmax = 500
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 20
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
