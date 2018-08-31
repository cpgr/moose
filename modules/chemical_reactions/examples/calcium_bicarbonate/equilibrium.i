# Example of reactive transport model with precipitation and dissolution.
# Calcium (ca2) and bicarbonate (hco3) reaction to form calcite (CaCO3).
# Models bicarbonate injection following calcium injection, so that a
# moving reaction front forms a calcite precipitation zone. As the front moves,
# the upstream side of the front continues to form calcite via precipitation,
# while at the downstream side, dissolution of the solid calcite occurs.
#
# The reaction network considered is as follows:
# Aqueous equilibrium reactions:
# a)  h+ + hco3- = CO2(aq),             Keq = 10^(6.341)
# b)  hco3- = h+ + CO23-,               Keq = 10^(-10.325)
# c)  ca2+ + hco3- = h+ + CaCO3(aq),    Keq = 10^(-7.009)
# d)  ca2+ + hco3- = cahco3+,           Keq = 10^(-0.653)
# e)  ca2+ = h+ + CaOh+,                Keq = 10^(-12.85)
# f)  - h+ = oh-,                       Keq = 10^(-13.991)
#
# Kinetic reactions
# g)  ca2+ + hco3- = h+ + CaCO3(s),     A = 0.461 m^2/L, k = 6.456542e-2 mol/m^2 s,
# Keq = 10^(1.8487)
#
# The primary chemical species are h+, hco3- and ca2+. The pressure gradient is fixed,
# and a conservative tracer is also included.
#
# This example is taken from:
# Guo et al, A parallel, fully coupled, fully implicit solution to reactive
# transport in porous media using the preconditioned Jacobian-Free Newton-Krylov
# Method, Advances in Water Resources, 53, 101-108 (2013).
inactive = 'BCs'
[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 1
  ny = 1
  xmax = 1
  ymax = 1
[]

[Variables]
  [tracer]
  []
  [ca2+]
  []
  [h+]
    scaling = 1e6
  []
  [co2_aq]
  []
[]

[AuxVariables]
  [pressure]
    order = FIRST
    family = LAGRANGE
  []
  [hco3-]
    initial_condition = 1e-7
  []
  [caco3_aq]
    initial_condition = 3e-12
  []
  [caco3_s]
    initial_condition = 1e-6
  []
  [cahco3+]
    initial_condition = 1e-9
  []
  [caoh+]
    initial_condition = 2e-11
  []
  [co32-]
    initial_condition = 5e-14
  []
  [oh-]
    initial_condition = 5e-11
  []
  [conductivity]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[ICs]
  inactive = 'pressure_ic'
  [pressure_ic]
    type = FunctionIC
    variable = pressure
    function = pic
  []
  [co2_aq_ic]
    type = ConstantIC
    variable = co2_aq
    value = 1.5e-4
  []
  [ca2_ic]
    type = ConstantIC
    variable = ca2+
    value = 0.03
  []
  [h_ic]
    type = ConstantIC
    variable = h+
    value = 2e-4
  []
  [tracer_ic]
    type = BoundingBoxIC
    variable = tracer
    x1 = 0.0
    y1 = 0.0
    x2 = 1.0e-10
    y2 = 1
    inside = 1.0
    outside = 0.0
  []
  [conductivity]
    type = RandomIC
    legacy_generator = false
    variable = conductivity
    min = 1e-4
    max = 1e-3
  []
[]

[Functions]
  [pic]
    type = ParsedFunction
    value = '3-1*x'
  []
[]

[ReactionNetwork]
  [AqueousEquilibriumReactions]
    primary_species = 'ca2+ co2_aq h+'
    secondary_species = 'hco3- co32- caco3_aq cahco3+ caoh+ oh-'
    pressure = 'pressure'
    reactions = '- h+ + co2_aq = hco3- -6.341,
                 co2_aq - 2h+ = co32- -16.666,
                 ca2+ + co2_aq - 2h+ = caco3_aq -13.35,
                 ca2+ + co2_aq - h+ = cahco3+ -6.994,
                 ca2+ - h+ = caoh+ -12.85,
                 - h+ = oh- -13.991'
  []
  [SolidKineticReactions]
    primary_species = 'ca2+ co2_aq h+'
    kin_reactions = 'ca2+ + co2_aq - 2h+ = caco3_s'
    secondary_species = 'caco3_s'
    log10_keq = '1.8487'
    reference_temperature = '298.15'
    system_temperature = '298.15'
    gas_constant = 8.314
    specific_reactive_surface_area = '4.61e-4'
    kinetic_rate_constant = '6.456542e-7'
    activation_energy = '1.5e4'
  []
[]

[Kernels]
  [tracer_ie]
    type = PrimaryTimeDerivative
    variable = tracer
  []
  [tracer_pd]
    type = PrimaryDiffusion
    variable = tracer
  []
  [tracer_conv]
    type = PrimaryConvection
    variable = tracer
    p = 'pressure'
  []
  [ca2+_ie]
    type = PrimaryTimeDerivative
    variable = ca2+
  []
  [ca2+_pd]
    type = PrimaryDiffusion
    variable = ca2+
  []
  [ca2+_conv]
    type = PrimaryConvection
    variable = ca2+
    p = 'pressure'
  []
  [h+_ie]
    type = PrimaryTimeDerivative
    variable = h+
  []
  [h+_pd]
    type = PrimaryDiffusion
    variable = h+
  []
  [h+_conv]
    type = PrimaryConvection
    variable = h+
    p = 'pressure'
  []
  [co2_aq_ie]
    type = PrimaryTimeDerivative
    variable = co2_aq
  []
  [co2_aq_pd]
    type = PrimaryDiffusion
    variable = co2_aq
  []
  [co2_aq_conv]
    type = PrimaryConvection
    variable = co2_aq
    p = 'pressure'
  []
[]

[BCs]
  # [./ca2+_left]
  # type = PresetBC
  # variable = ca2+
  # boundary = left
  # value = 0.03
  # [../]
  # [./h+_left]
  # type = DirichletBC
  # variable = h+
  # boundary = left
  # value = 1e-7
  # [../]
  [tracer_left]
    type = DirichletBC
    variable = tracer
    boundary = 'left'
    value = 1.0
  []
  [tracer_right]
    type = ChemicalOutFlowBC
    variable = tracer
    boundary = 'right'
  []
  [ca2+_right]
    type = ChemicalOutFlowBC
    variable = ca2+
    boundary = 'right'
  []
  [co2_aq_left]
    type = PresetBC
    variable = co2_aq
    boundary = 'left'
    value = 0.001
  []
  [co2_aq_right]
    type = ChemicalOutFlowBC
    variable = co2_aq
    boundary = 'right'
  []
  [h+_right]
    type = ChemicalOutFlowBC
    variable = h+
    boundary = 'right'
  []
[]

[Materials]
  [porous]
    type = GenericConstantMaterial
    prop_names = 'diffusivity  porosity'
    prop_values = '1e-7 0.2'
  []
  [conductivity]
    type = Conductivity
    conductivity = 'conductivity'
  []
[]

[Executioner]
  type = Transient
  solve_type = PJFNK
  l_max_its = 50
  l_tol = 1e-5
  nl_max_its = 10
  nl_rel_tol = 1e-5
  end_time = 1e6
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 10
  []
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

[Outputs]
  perf_graph = true
  exodus = true
[]
