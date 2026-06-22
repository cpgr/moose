# Phase appearance/disappearance test for local-equilibrium halite precipitation in
# PorousFlowBrineCO2 (the persistent total-salt variable z_s should give smooth precipitation
# onset and redissolution with no primary-variable switching).
#
# A single, closed cell (no flux) holds a single-phase aqueous brine (the CO2 is fully dissolved at
# 15 MPa, so the liquid mass fraction f_l = 1 and the saturation threshold is exactly the halite
# solubility X_eq(T)). The salt mass fraction variable z_s = 0.27 is held constant by the closed
# cell. The
# temperature is cycled 90 -> 40 -> 90 C: X_eq(90 C) = 0.2769 > z_s so the cell starts
# undersaturated with no solid halite, X_eq drops below z_s near 60 C so halite precipitates
# (onset), bottoms out at X_eq(40 C) = 0.2664, then the warm-up redissolves it.
#
# Verifies: total_salt is conserved throughout the crossing in BOTH directions, solid_salt rises
# from and returns to zero, and the Newton solve converges at the onset/redissolution transition
# (no stall). z_s never switches, so onset and redissolution are continuous.

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 1
  ny = 1
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
    initial_condition = 0.01
  []
  [xnacl]
    initial_condition = 0.27
  []
[]

[AuxVariables]
  [temp]
    initial_condition = 90
  []
  [halite_vf]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[Functions]
  [tfunc]
    # Cool 90 -> 40 C (t in [0,10]) then warm 40 -> 90 C (t in [10,20]); X_eq crosses z_s = 0.27
    # near 60 C on both legs, so halite precipitates then redissolves.
    type = ParsedFunction
    expression = 'if(t <= 10, 90 - 5*t, 5*t - 10)'
  []
[]

[AuxKernels]
  [temp_aux]
    type = FunctionAux
    variable = temp
    function = tfunc
    execute_on = 'INITIAL TIMESTEP_BEGIN'
  []
  [halite_vf]
    # Copy the (nodal) halite volume fraction into an elemental AuxVariable so it can be integrated
    type = MaterialRealAux
    variable = halite_vf
    property = PorousFlow_halite_volume_fraction_nodal
    execute_on = 'TIMESTEP_END'
  []
[]

[Kernels]
  [mass0]
    type = PorousFlowMassTimeDerivative
    variable = pgas
    fluid_component = 0
  []
  [mass1]
    type = PorousFlowMassTimeDerivative
    variable = zi
    fluid_component = 1
  []
  [mass2]
    type = PorousFlowMassTimeDerivative
    variable = xnacl
    fluid_component = 2
  []
  [precipitate2]
    type = PorousFlowPrecipitateMassTimeDerivative
    variable = xnacl
  []
  # Advective flux is identically zero in this single uniform closed cell; included so the qp
  # fluid-state material (needed by the conservation postprocessors) is created.
  [adv0]
    type = PorousFlowAdvectiveFlux
    variable = pgas
    fluid_component = 0
  []
  [adv1]
    type = PorousFlowAdvectiveFlux
    variable = zi
    fluid_component = 1
  []
  [adv2]
    type = PorousFlowAdvectiveFlux
    variable = xnacl
    fluid_component = 2
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
    temperature = temp
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
  [porosity]
    type = PorousFlowPorosityConst
    porosity = 0.1
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  dt = 0.5
  end_time = 20
  nl_abs_tol = 1e-12
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
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
