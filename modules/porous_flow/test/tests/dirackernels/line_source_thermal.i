# Test PorousFlowLineEnthalpySourceFromPostprocessor alongside
# PorousFlowLineSourceFromPostprocessor.
#
# 0.1 kg/s of 350 K fluid is injected for 1000 s into a closed column, so the fluid mass rises by
# exactly 100 kg.  There are no flux or conduction kernels, so nothing but the two sources can
# change the mass or the energy.
#
# Unlike the mass, the energy added is NOT independent of how the line is resolved: each completion
# evaluates the enthalpy of the injected fluid at the pressure where it sits, which is the point of
# distributing the heat rather than depositing it all at one place.

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
  [T_in_fn]
    type = ConstantFunction
    value = 350
  []
[]

[Variables]
  [pp]
  []
  [temp]
  []
[]

[ICs]
  [PressureIC]
    variable = pp
    type = ConstantIC
    value = 20e6
  []
  [TempIC]
    variable = temp
    type = ConstantIC
    value = 300
  []
[]

[Kernels]
  [mass0]
    type = PorousFlowMassTimeDerivative
    fluid_component = 0
    variable = pp
  []
  [energy]
    type = PorousFlowEnergyTimeDerivative
    variable = temp
  []
[]

[UserObjects]
  [dictator]
    type = PorousFlowDictator
    porous_flow_vars = 'pp temp'
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
    cv = 4000
    cp = 4000
  []
[]

[Materials]
  [temperature]
    type = PorousFlowTemperature
    temperature = temp
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
  [rock_heat]
    type = PorousFlowMatrixInternalEnergy
    specific_heat_capacity = 830.0
    density = 2750
  []
[]

[DiracKernels]
  [source]
    type = PorousFlowLineSourceFromPostprocessor
    variable = pp
    mass_flux = mass_flux_in
    point_file = line_source.bh
  []
  [source_h]
    type = PorousFlowLineEnthalpySourceFromPostprocessor
    variable = temp
    mass_flux = mass_flux_in
    point_file = line_source.bh
    pressure = pp
    T_in = T_in
    fp = simple_fluid
  []
[]

[Postprocessors]
  [total_mass]
    type = PorousFlowFluidMass
    execute_on = 'initial timestep_end'
  []
  [total_heat]
    type = PorousFlowHeatEnergy
    # The fluid phase must be named, or only the energy of the rock skeleton is counted and the
    # heat carried in by the injected fluid is invisible
    phase = 0
    execute_on = 'initial timestep_end'
  []
  [mass_flux_in]
    type = FunctionValuePostprocessor
    function = mass_flux_fn
    execute_on = 'initial timestep_begin'
  []
  [T_in]
    type = FunctionValuePostprocessor
    function = T_in_fn
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
  nl_abs_tol = 1e-12
  dt = 100
  end_time = 1000
[]

[Outputs]
  csv = true
  execute_on = 'initial timestep_end'
  file_base = line_source_thermal
[]
