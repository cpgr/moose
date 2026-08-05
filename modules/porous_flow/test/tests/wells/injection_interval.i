# An injection well completed over an interval rather than at a single point.
#
# 5 kg/s is injected over a 20 m completion interval.  The interval is described by line_base, so
# one completion is generated for each of the ten elements it passes through, and the [Wells] action
# builds a line source per fluid component rather than a point source.  The rate is prescribed, so
# the total injected must be 5 kg/s however many completions the interval is divided into: the
# cumulative reporters of the two components sum to 5000 kg after 1000 s.

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 10
  ny = 10
  xmax = 100
  ymax = 20
[]

[GlobalParams]
  PorousFlowDictator = dictator
[]

[Functions]
  [rate_fn]
    type = ConstantFunction
    value = 5
  []
[]

[Variables]
  [massfrac0]
  []
  [pp]
  []
[]

[ICs]
  [massfrac0]
    type = ConstantIC
    variable = massfrac0
    value = 0.2
  []
  [pp]
    type = ConstantIC
    variable = pp
    value = 20e6
  []
[]

[Wells]
  [inj1]
    type = PorousFlowInjectionWell
    # The base sits inside the first element and the interval ends inside the last, so that neither
    # end lies on a mesh face
    line_base = '0.1 5 0.5'
    line_length = 19
    line_direction = '0 1 0'
    rate = rate_fn
    injected_mass_fractions = '0.3 0.7'
    porepressure = pp
    mass_fraction_vars = massfrac0
  []
[]

[Kernels]
  [mass0]
    type = PorousFlowMassTimeDerivative
    fluid_component = 0
    variable = massfrac0
  []
  [mass1]
    type = PorousFlowMassTimeDerivative
    fluid_component = 1
    variable = pp
  []
  [flux0]
    type = PorousFlowAdvectiveFlux
    fluid_component = 0
    variable = massfrac0
    gravity = '0 0 0'
  []
  [flux1]
    type = PorousFlowAdvectiveFlux
    fluid_component = 1
    variable = pp
    gravity = '0 0 0'
  []
[]

[UserObjects]
  [dictator]
    type = PorousFlowDictator
    porous_flow_vars = 'massfrac0 pp'
    number_fluid_phases = 1
    number_fluid_components = 2
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
  []
  [ppss]
    type = PorousFlow1PhaseFullySaturated
    porepressure = pp
  []
  [massfrac]
    type = PorousFlowMassFraction
    mass_fraction_vars = massfrac0
  []
  [simple_fluid]
    type = PorousFlowSingleComponentFluid
    fp = simple_fluid
    phase = 0
  []
  [porosity]
    type = PorousFlowPorosityConst
    porosity = 0.2
  []
  [permeability]
    type = PorousFlowPermeabilityConst
    permeability = '1e-13 0 0  0 1e-13 0  0 0 1e-13'
  []
  [relperm]
    type = PorousFlowRelativePermeabilityConst
    phase = 0
  []
[]

[Postprocessors]
  # One mass per fluid component, so that the rise in each can be checked against the cumulative
  # reported by the well.  fluid_component defaults to 0, so it must be named to see component 1.
  [mass_0]
    type = PorousFlowFluidMass
    fluid_component = 0
    execute_on = 'initial timestep_end'
  []
  [mass_1]
    type = PorousFlowFluidMass
    fluid_component = 1
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
  nl_abs_tol = 1e-10
[]

[Outputs]
  csv = true
  execute_on = 'initial timestep_end'
  file_base = injection_interval
[]
