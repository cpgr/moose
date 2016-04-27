# Checking that gravity head is established in the transient situation when 0<=saturation<=1 (note the less-than-or-equal-to).
# 2phase (PP), 2components, vanGenuchten, constant fluid bulk-moduli for each phase, constant viscosity, constant permeability, Corey relative perm.
# A boundary condition enforces porepressures at the right boundary
# For better agreement with the analytical solution (ana_pp), just increase nx

[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 10
  xmin = -1
  xmax = 0
[]

[GlobalParams]
  PorousFlowDictator_UO = dictator
[]

[Variables]
  [./ppwater]
    initial_condition = 0
  [../]
  [./ppgas]
    initial_condition = 0.5
  [../]
[]

[AuxVariables]
  [./massfrac_ph0_sp0]
    initial_condition = 1
  [../]
  [./massfrac_ph1_sp0]
    initial_condition = 0
  [../]
[]

[BCs]
  [./ppwater]
    type = PresetBC
    boundary = right
    variable = ppwater
    value = 0
  [../]
  [./ppgas]
    type = PresetBC
    boundary = right
    variable = ppgas
    value = 0.5
  [../]
[]

[Kernels]
  [./mass0]
    type = PorousFlowMassTimeDerivative
    component_index = 0
    variable = ppwater
  [../]
  [./flux0]
    type = PorousFlowAdvectiveFlux
    component_index = 0
    variable = ppwater
    gravity = '-1 0 0'
  [../]
  [./mass1]
    type = PorousFlowMassTimeDerivative
    component_index = 1
    variable = ppgas
  [../]
  [./flux1]
    type = PorousFlowAdvectiveFlux
    component_index = 1
    variable = ppgas
    gravity = '-1 0 0'
  [../]
[]


[Functions]
  [./ana_ppwater]
    type = ParsedFunction
    vars = 'g B p0 rho0'
    vals = '1 2 pp_water_top 1'
    value = '-B*log(exp(-p0/B)+g*rho0*x/B)' # expected pp at base
  [../]
[]

[UserObjects]
  [./dictator]
    type = PorousFlowDictator
    porous_flow_vars = 'ppwater ppgas'
    number_fluid_phases = 2
    number_fluid_components = 2
  [../]
[]

[Materials]
  [./ppss]
    type = PorousFlowMaterial2PhasePP_VG
    phase0_porepressure = ppwater
    phase1_porepressure = ppgas
    m = 0.5
    al = 1
  [../]
  [./massfrac]
    type = PorousFlowMaterialMassFractionBuilder
    mass_fraction_vars = 'massfrac_ph0_sp0 massfrac_ph1_sp0'
  [../]
  [./dens0]
    type = PorousFlowMaterialDensityConstBulk
    density_P0 = 1
    bulk_modulus = 1.2
    phase = 0
  [../]
  [./dens1]
    type = PorousFlowMaterialDensityConstBulk
    density_P0 = 0.1
    bulk_modulus = 1
    phase = 1
  [../]
  [./dens_all]
    type = PorousFlowMaterialJoinerOld
    material_property = PorousFlow_fluid_phase_density
  [../]
  [./dens_qp_all]
    type = PorousFlowMaterialJoiner
    material_property = PorousFlow_fluid_phase_density_qp
    use_qps = true
  [../]
  [./porosity]
    type = PorousFlowMaterialPorosityConst
    porosity = 0.1
  [../]
  [./visc0]
    type = PorousFlowMaterialViscosityConst
    viscosity = 1
    phase = 0
  [../]
  [./visc1]
    type = PorousFlowMaterialViscosityConst
    viscosity = 0.5
    phase = 1
  [../]
  [./visc_all]
    type = PorousFlowMaterialJoiner
    material_property = PorousFlow_viscosity
  [../]
  [./permeability]
    type = PorousFlowMaterialPermeabilityConst
    permeability = '1 0 0  0 2 0  0 0 3'
  [../]
  [./relperm_water]
    type = PorousFlowMaterialRelativePermeabilityCorey
    n_j = 1
    phase = 0
  [../]
  [./relperm_gas]
    type = PorousFlowMaterialRelativePermeabilityCorey
    n_j = 1
    phase = 1
  [../]
  [./relperm_all]
    type = PorousFlowMaterialJoiner
    material_property = PorousFlow_relative_permeability
  [../]
[]

[Postprocessors]
  [./pp_water_top]
    type = PointValue
    variable = ppwater
    point = '0 0 0'
  [../]
  [./pp_water_base]
    type = PointValue
    variable = ppwater
    point = '-1 0 0'
  [../]
  [./pp_water_analytical]
    type = FunctionValuePostprocessor
    function = ana_ppwater
    point = '-1 0 0'
  [../]
  [./ppwater_00]
    type = PointValue
    variable = ppwater
    point = '0 0 0'
  [../]
  [./ppwater_01]
    type = PointValue
    variable = ppwater
    point = '-0.1 0 0'
  [../]
  [./ppwater_02]
    type = PointValue
    variable = ppwater
    point = '-0.2 0 0'
  [../]
  [./ppwater_03]
    type = PointValue
    variable = ppwater
    point = '-0.3 0 0'
  [../]
  [./ppwater_04]
    type = PointValue
    variable = ppwater
    point = '-0.4 0 0'
  [../]
  [./ppwater_05]
    type = PointValue
    variable = ppwater
    point = '-0.5 0 0'
  [../]
  [./ppwater_06]
    type = PointValue
    variable = ppwater
    point = '-0.6 0 0'
  [../]
  [./ppwater_07]
    type = PointValue
    variable = ppwater
    point = '-0.7 0 0'
  [../]
  [./ppwater_08]
    type = PointValue
    variable = ppwater
    point = '-0.8 0 0'
  [../]
  [./ppwater_09]
    type = PointValue
    variable = ppwater
    point = '-0.9 0 0'
  [../]
  [./ppwater_10]
    type = PointValue
    variable = ppwater
    point = '-1 0 0'
  [../]
  [./ppgas_00]
    type = PointValue
    variable = ppgas
    point = '0 0 0'
  [../]
  [./ppgas_01]
    type = PointValue
    variable = ppgas
    point = '-0.1 0 0'
  [../]
  [./ppgas_02]
    type = PointValue
    variable = ppgas
    point = '-0.2 0 0'
  [../]
  [./ppgas_03]
    type = PointValue
    variable = ppgas
    point = '-0.3 0 0'
  [../]
  [./ppgas_04]
    type = PointValue
    variable = ppgas
    point = '-0.4 0 0'
  [../]
  [./ppgas_05]
    type = PointValue
    variable = ppgas
    point = '-0.5 0 0'
  [../]
  [./ppgas_06]
    type = PointValue
    variable = ppgas
    point = '-0.6 0 0'
  [../]
  [./ppgas_07]
    type = PointValue
    variable = ppgas
    point = '-0.7 0 0'
  [../]
  [./ppgas_08]
    type = PointValue
    variable = ppgas
    point = '-0.8 0 0'
  [../]
  [./ppgas_09]
    type = PointValue
    variable = ppgas
    point = '-0.9 0 0'
  [../]
  [./ppgas_10]
    type = PointValue
    variable = ppgas
    point = '-1 0 0'
  [../]
[]

[Preconditioning]
  active = andy
  [./andy]
    type = SMP
    full = true
    petsc_options_iname = '-ksp_type -pc_type -snes_atol -snes_rtol -snes_max_it'
    petsc_options_value = 'bcgs bjacobi 1E-12 1E-10 10000'
  [../]
  [./check]
    type = SMP
    full = true
    #petsc_options = '-snes_test_display'
    petsc_options_iname = '-ksp_type -pc_type -snes_atol -snes_rtol -snes_max_it -snes_type'
    petsc_options_value = 'bcgs bjacobi 1E-12 1E-10 10000 test'
  [../]
[]


[Executioner]
  type = Transient
  solve_type = Newton
  [./TimeStepper]
    type = FunctionDT
    time_t = '1E-3 1E-2 1E-1 2E-1'
    time_dt = '1E-3 1E-2 0.2E-1 1E-1'
  [../]
  end_time = 1.0
[]

[Outputs]
  execute_on = 'initial final'
  file_base = grav02d
  [./csv]
    type = CSV
  [../]
  exodus = true
[]
