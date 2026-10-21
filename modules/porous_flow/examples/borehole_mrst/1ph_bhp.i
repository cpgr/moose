# Finite-element counterpart of the single-phase pressure-controlled case that the finite-volume
# Peaceman borehole was verified against MRST with (see compare.py for the reference).
#
# The same box and borehole as 1ph_rate.i, but the bottom-hole pressure is prescribed: 9.95 MPa
# until t = 50 s (production) and 10.05 MPa afterwards (injection).  PorousFlowPeacemanBorehole
# only lets fluid flow in the direction its 'character' allows, so character changes sign with the
# bottom-hole pressure.  Every completion flows the same way in both periods here, so this matches
# the finite-volume borehole, which allows flow in either direction at every completion.
[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 4
  ny = 4
  nz = 6
  xmin = -10
  xmax = 10
  ymin = -10
  ymax = 10
  zmin = -15
  zmax = 15
[]

[GlobalParams]
  PorousFlowDictator = dictator
  # Peaceman's finite-difference value.  compare.py also runs 0.1594, the value for rectangular
  # finite elements (Chen and Zhang, Int J Num Analysis and Modeling 3 (2008) 375-388, eqn 4.13)
  re_constant = 0.28
[]

[Variables]
  [pp]
  []
[]

[ICs]
  [pp]
    type = FunctionIC
    variable = pp
    function = '1E7 - 1E4 * z'
  []
[]

[Kernels]
  [mass0]
    type = PorousFlowMassTimeDerivative
    fluid_component = 0
    variable = pp
  []
  [flux0]
    type = PorousFlowAdvectiveFlux
    fluid_component = 0
    variable = pp
    gravity = '0 0 -10'
  []
[]

[DiracKernels]
  [bh]
    type = PorousFlowPeacemanBorehole
    variable = pp
    SumQuantityUO = bh_outflow
    point_file = borehole.bh
    fluid_phase = 0
    unit_weight = '0 0 -1E4'
    use_mobility = true
    character = character
    bottom_p_or_t = bottom_p
  []
[]

[UserObjects]
  [dictator]
    type = PorousFlowDictator
    porous_flow_vars = 'pp'
    number_fluid_phases = 1
    number_fluid_components = 1
  []
  [bh_outflow]
    type = PorousFlowSumQuantity
  []
[]

[FluidProperties]
  [water]
    type = SimpleFluidProperties
    bulk_modulus = 2E9
    viscosity = 1E-3
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
  [water]
    type = PorousFlowSingleComponentFluid
    fp = water
    phase = 0
  []
  [porosity]
    type = PorousFlowPorosityConst
    porosity = 0.1
  []
  [permeability]
    type = PorousFlowPermeabilityConst
    permeability = '1E-12 0 0 0 1E-12 0 0 0 1E-12'
  []
  [relperm]
    type = PorousFlowRelativePermeabilityConst
    phase = 0
  []
[]

[Functions]
  [bottom_p]
    type = ParsedFunction
    expression = 'if(t <= 50, 9.95E6, 1.005E7)'
  []
  [character]
    # production while the bottom-hole pressure is below the reservoir, injection once above it
    type = ParsedFunction
    expression = 'if(t <= 50, 1, -1)'
  []
[]

[Postprocessors]
  [p_bot]
    type = FunctionValuePostprocessor
    function = bottom_p
  []
  [mass]
    type = PorousFlowFluidMass
  []
  [mass_old]
    type = PorousFlowFluidMass
    execute_on = 'timestep_begin'
    outputs = none
  []
  [dt]
    type = TimestepSize
    outputs = none
  []
  [bh_rate]
    # mass rate leaving the porous medium, from the change in fluid mass
    type = ParsedPostprocessor
    expression = '(mass_old - mass) / dt'
    pp_names = 'mass mass_old dt'
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
  end_time = 100
  [TimeStepper]
    type = TimeSequenceStepper
    time_sequence = '5 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100'
  []
  nl_abs_tol = 1E-10
  nl_rel_tol = 1E-12
[]

[Outputs]
  csv = true
[]
