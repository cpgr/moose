# Tests the Jacobian of the salt equation when PorousFlowBrineCO2 computes local-equilibrium
# halite precipitation. The salt mass fraction is a transported variable interpreted as the total
# salt z_s; the PorousFlowPrecipitateMassTimeDerivative kernel carries the solid-halite mass. All
# nodes are salt-saturated (z_s > X_eq) so the precipitation branch (and its derivatives) is active
# and smooth for the finite-difference comparison.

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 2
  ny = 2
[]

[GlobalParams]
  PorousFlowDictator = dictator
  gravity = '0 0 0'
[]

[Variables]
  [pgas]
  []
  [zi]
  []
  [xnacl]
  []
[]

[ICs]
  [pgas]
    type = RandomIC
    min = 1e6
    max = 4e6
    variable = pgas
    seed = 1
  []
  [z]
    type = RandomIC
    min = 0.2
    max = 0.8
    variable = zi
    seed = 2
  []
  [xnacl]
    type = RandomIC
    min = 0.28
    max = 0.33
    variable = xnacl
    seed = 3
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
    type = PorousFlowCapillaryPressureVG
    m = 0.5
    alpha = 1e1
    pc_max = 1e4
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
    temperature = 50
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
  end_time = 1
  nl_abs_tol = 1e-12
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]
