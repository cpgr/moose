# Test PorousFlowLineSourceFromPostprocessor.
#
# 0.1 kg/s is injected for 1000 s into a closed column, so the fluid mass must rise from its
# initial 2000 kg by exactly 100 kg however the line is resolved.  There is no flux kernel, as in
# frompps.i, so nothing but the source can change the mass and the test isolates the question of
# whether the prescribed total is delivered.
#
# The same gold is used whether the line carries four completions or one, which is the property
# that matters: the total is prescribed, so it must not depend on the number of Dirac points the
# line is divided into.

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 1
  ny = 10
  xmax = 1
  ymax = 10
[]

[GlobalParams]
  PorousFlowDictator = dictator
[]

[Functions]
  [mass_flux_fn]
    type = ConstantFunction
    value = 0.1
  []
[]

[Variables]
  [pp]
  []
[]

[ICs]
  [PressureIC]
    variable = pp
    type = ConstantIC
    value = 20e6
  []
[]

[Kernels]
  [mass0]
    type = PorousFlowMassTimeDerivative
    fluid_component = 0
    variable = pp
  []
[]

[UserObjects]
  [dictator]
    type = PorousFlowDictator
    porous_flow_vars = pp
    number_fluid_phases = 1
    number_fluid_components = 1
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
[]

[DiracKernels]
  [source]
    type = PorousFlowLineSourceFromPostprocessor
    variable = pp
    mass_flux = mass_flux_in
    point_file = line_source.bh
  []
[]

[Postprocessors]
  [total_mass]
    type = PorousFlowFluidMass
    execute_on = 'initial timestep_end'
  []
  [mass_flux_in]
    type = FunctionValuePostprocessor
    function = mass_flux_fn
    execute_on = 'initial timestep_begin'
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
  nl_abs_tol = 1e-14
  dt = 100
  end_time = 1000
[]

[Outputs]
  csv = true
  execute_on = 'initial timestep_end'
  file_base = line_source
[]
