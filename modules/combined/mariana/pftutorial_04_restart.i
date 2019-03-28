# Darcy flow with heat advection and conduction, and elasticity
[Mesh]
  file = pftutorial_04_out.e
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
  PorousFlowDictator = dictator
  biot_coefficient = 1.0
[]

[Variables]
  [./porepressure]
    initial_from_file_var = porepressure
  [../]
  [./temperature]
    initial_from_file_var = temperature
    scaling = 1E-8
  [../]
  [./disp_x]
    initial_from_file_var = disp_x
    scaling = 1E-10
  [../]
  [./disp_y]
    initial_from_file_var = disp_y
    scaling = 1E-10
  [../]
  [./disp_z]
    initial_from_file_var = disp_z
    scaling = 1E-10
  [../]
[]

[PorousFlowBasicTHM]
  porepressure = porepressure
  temperature = temperature
  coupling_type = ThermoHydroMechanical
  gravity = '0 0 0'
  fp = the_simple_fluid
  thermal_eigenstrain_name = thermal_contribution
  use_displaced_mesh = false
[]

[BCs]
  [./constant_injection_porepressure]
    type = PresetBC
    variable = porepressure
    value = 1E6
    boundary = injection_area
  [../]
  [./constant_injection_temperature]
    type = PresetBC
    variable = temperature
    value = 313
    boundary = injection_area
  [../]

  [./roller_tmax]
    type = PresetBC
    variable = disp_x
    value = 0
    boundary = tmax
  [../]
  [./roller_tmin]
    type = PresetBC
    variable = disp_y
    value = 0
    boundary = tmin
  [../]
  [./roller_top_bottom]
    type = PresetBC
    variable = disp_z
    value = 0
    boundary = 'top bottom'
  [../]
  [./cavity_pressure_x]
    type = Pressure
    boundary = injection_area
    variable = disp_x
    component = 0
    factor = 1E6
    use_displaced_mesh = false
  [../]
  [./cavity_pressure_y]
    type = Pressure
    boundary = injection_area
    variable = disp_y
    component = 1
    factor = 1E6
    use_displaced_mesh = false
  [../]
[]

[AuxVariables]
  [./stress_rr]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./stress_pp]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./initial_stress_xx]
    family = MONOMIAL
    order = CONSTANT
    initial_from_file_var = stress_xx
  [../]
  [./initial_stress_xy]
    family = MONOMIAL
    order = CONSTANT
    initial_from_file_var = stress_xy
  [../]
  [./initial_stress_xz]
    family = MONOMIAL
    order = CONSTANT
    initial_from_file_var = stress_xz
  [../]
  [./initial_stress_yx]
    family = MONOMIAL
    order = CONSTANT
    initial_from_file_var = stress_yx
  [../]
  [./initial_stress_yy]
    family = MONOMIAL
    order = CONSTANT
    initial_from_file_var = stress_yy
  [../]
  [./initial_stress_yz]
    family = MONOMIAL
    order = CONSTANT
    initial_from_file_var = stress_yz
  [../]
  [./initial_stress_zx]
    family = MONOMIAL
    order = CONSTANT
    initial_from_file_var = stress_zx
  [../]
  [./initial_stress_zy]
    family = MONOMIAL
    order = CONSTANT
    initial_from_file_var = stress_zy
  [../]
  [./initial_stress_zz]
    family = MONOMIAL
    order = CONSTANT
    initial_from_file_var = stress_zz
  [../]
[]

[AuxKernels]
  [./stress_rr]
    type = RankTwoScalarAux
    rank_two_tensor = stress
    variable = stress_rr
    scalar_type = RadialStress
    point1 = '0 0 0'
    point2 = '0 0 1'
  [../]
  [./stress_pp]
    type = RankTwoScalarAux
    rank_two_tensor = stress
    variable = stress_pp
    scalar_type = HoopStress
    point1 = '0 0 0'
    point2 = '0 0 1'
  [../]
[]

[Modules]
  [./FluidProperties]
    [./the_simple_fluid]
      type = SimpleFluidProperties
      bulk_modulus = 2E9
      viscosity = 1.0E-3
      density0 = 1000.0
      thermal_expansion = 0.0002
      cp = 4194
      cv = 4186
      porepressure_coefficient = 0
    [../]
  [../]
[]

[Materials]
  [./porosity]
    type = PorousFlowPorosity
    porosity_zero = 0.1
    fluid = true
    thermal = true
    mechanical = true
    thermal_expansion_coeff = 24e-6
    solid_bulk = 7.14e9
    reference_temperature = 293
    reference_porepressure = 101325
  [../]
  [./biot_modulus]
    type = PorousFlowConstantBiotModulus
    solid_bulk_compliance = 2E-7
    fluid_bulk_modulus = 1E7
  [../]
  [./permeability_aquifer]
    type = PorousFlowPermeabilityConst
    block = aquifer
    permeability = '1E-14 0 0   0 1E-14 0   0 0 1E-14'
  [../]
  [./permeability_caps]
    type = PorousFlowPermeabilityConst
    block = caps
    permeability = '1E-15 0 0   0 1E-15 0   0 0 1E-16'
  [../]

  [./thermal_expansion]
    type = PorousFlowConstantThermalExpansionCoefficient
    drained_coefficient = 0.003
    fluid_coefficient = 0.0002
  [../]
  [./rock_internal_energy]
    type = PorousFlowMatrixInternalEnergy
    density = 2500.0
    specific_heat_capacity = 1200.0
  [../]
  [./thermal_conductivity]
    type = PorousFlowThermalConductivityIdeal
    dry_thermal_conductivity = '10 0 0  0 10 0  0 0 10'
    block = 'caps aquifer'
  [../]

  [./elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 5E9
    poissons_ratio = 0.15
  [../]
  [./strain]
    type = ComputeSmallStrain
    eigenstrain_names = 'initial_stress thermal_contribution'
  [../]
  [./thermal_contribution]
    type = ComputeThermalExpansionEigenstrain
    temperature = temperature
    thermal_expansion_coeff = 8e-6 # this is the linear thermal expansion coefficient
    eigenstrain_name = thermal_contribution
    stress_free_temperature = 293
  [../]
  [./stress]
    type = ComputeLinearElasticStress
  [../]
  [./initial_strain]
    type = ComputeEigenstrainFromInitialStress
    initial_stress = '1 1 1  1 1 1  1 1 1'
    initial_stress_aux = 'initial_stress_xx initial_stress_xy initial_stress_xz initial_stress_yx initial_stress_yy initial_stress_yz initial_stress_zx initial_stress_zy initial_stress_zz'
    eigenstrain_name = initial_stress
  [../]
[]

[Preconditioning]
  active = basic
  [./basic]
    type = SMP
    full = true
    petsc_options = '-ksp_diagonal_scale -ksp_diagonal_scale_fix'
    petsc_options_iname = '-pc_type -sub_pc_type -sub_pc_factor_shift_type -pc_asm_overlap'
    petsc_options_value = ' asm      lu           NONZERO                   2'
  [../]
  [./preferred_but_might_not_be_installed]
    type = SMP
    full = true
    petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
    petsc_options_value = ' lu       mumps'
  [../]
[]

[Executioner]
  type = Transient
  solve_type = Newton
  start_time = 5e5
  end_time = 1E6
  dt = 1E5
  nl_abs_tol = 1E-10
[]

[Outputs]
  exodus = true
[]
