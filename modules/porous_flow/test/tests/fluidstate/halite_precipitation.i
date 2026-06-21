# Conservation test for local-equilibrium halite precipitation in PorousFlowBrineCO2.
#
# A single, closed cell (no flux) holds a salt-saturated brine-CO2 mixture. The salt mass fraction
# variable is the total salt z_s; PorousFlowMassTimeDerivative carries the dissolved salt and
# PorousFlowPrecipitateMassTimeDerivative carries the solid halite. The temperature is ramped, which
# shifts the halite solubility X_eq(T) and so moves salt between the dissolved and solid states.
# Because the cell is closed, the total salt mass (dissolved + solid) must remain constant: the
# postprocessor total_salt should not change, while solid_salt and dissolved_salt do.

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
    initial_condition = 5e6
  []
  [zi]
    initial_condition = 0.4
  []
  [xnacl]
    initial_condition = 0.3
  []
[]

[AuxVariables]
  [temp]
    initial_condition = 70
  []
  [halite_vf]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[Functions]
  [tfunc]
    # Cool from 70 C to 30 C; X_eq drops, so more halite precipitates
    type = ParsedFunction
    expression = '70 - 8*t'
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
  dt = 1
  end_time = 5
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
