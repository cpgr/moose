# Test BrineFluidProperties calculations of density, viscosity and thermal
# conductivity
#
# Experimental density values from Pitzer et al, "Thermodynamic properties
# of aqueous sodium chloride solution", Journal of Physical and Chemical
# Reference Data, 13, 1-102 (1984)
#
# Experimental viscosity values from Phillips et al, "Viscosity of NaCl and
# other solutions up to 350C and 50MPa pressures", LBL-11586 (1980)
#
# Thermal conductivity values from Ozbek and Phillips, "Thermal conductivity of
# aqueous NaCl solutions from 20C to 330C", LBL-9086 (1980)
#
#  --------------------------------------------------------------
#  Pressure (Mpa)                |   20      |    20     |   40
#  Temperature (C)               |   50      |   200     |  200
#  NaCl molality (mol/kg)        |    2      |     2     |    5
#  NaCl mass fraction (kg/kg)    |  0.1047   |  0.1047   |  0.2261
#  --------------------------------------------------------------
#  Expected values
#  --------------------------------------------------------------
#  Density (kg/m^3)              |  1068.52  |  959.27   |  1065.58
#  Viscosity (1e-6Pa.s)          |  679.8    |  180.0    |  263.1
#  Thermal conductivity (W/m/K)  |  0.630    |  0.649    |  0.633
#  --------------------------------------------------------------
#  Calculated values
#  --------------------------------------------------------------
#  Density (kg/m^3)              |  1067.18  |  958.68   |  1065.46
#  Viscosity (1e-6 Pa.s)         |  681.1    |  181.98    |  266.1
#  Thermal conductivity (W/m/K)  |  0.637    |   0.662    |  0.658
#  --------------------------------------------------------------
#
# All results are within expected accuracy

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 3
  ny = 1
  xmax = 3
[]

[Variables]
  [./dummy]
  [../]
[]

[AuxVariables]
  [./pressure]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./temperature]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./xnacl]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./density]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./viscosity]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./k]
    family = MONOMIAL
    order = CONSTANT
  [../]
[]

[Functions]
  [./pic]
    type = ParsedFunction
    value = 'if(x<2,20e6, 40e6)'
  [../]
  [./tic]
    type = ParsedFunction
    value = 'if(x<1, 323.15, 473.15)'
  [../]
  [./xic]
    type = ParsedFunction
    value = 'if(x<2,0.1047, 0.2261)'
  [../]
[]

[ICs]
  [./p_ic]
    type = FunctionIC
    function = pic
    variable = pressure
  [../]
  [./t_ic]
    type = FunctionIC
    function = tic
    variable = temperature
  [../]
  [./x_ic]
    type = FunctionIC
    function = xic
    variable = xnacl
  [../]
[]

[AuxKernels]
  [./density]
    type = MaterialRealAux
     variable = density
     property = density
  [../]
  [./viscosity]
    type = MaterialRealAux
     variable = viscosity
     property = viscosity
  [../]
  [./k]
    type = MaterialRealAux
     variable = k
     property = k
  [../]
[]

[Modules]
  [./FluidProperties]
    [./water]
      type = Water97FluidProperties
    [../]
    [./brine]
      type = BrineFluidProperties
      water_fp = water
    [../]
  []
[]

[Materials]
  [./fp_mat]
    type = BrineFluidPropertiesTestMaterial
    pressure = pressure
    temperature = temperature
    xnacl = xnacl
    fp = brine
  [../]
[]

[Kernels]
  [./diff]
    type = Diffusion
    variable = dummy
  [../]
[]

[Executioner]
  type = Steady
  solve_type = NEWTON
[]

[Postprocessors]
  [./density0]
    type = ElementalVariableValue
    variable = density
    elementid = 0
  [../]
  [./density1]
    type = ElementalVariableValue
    variable = density
    elementid = 1
  [../]
  [./density2]
    type = ElementalVariableValue
    variable = density
    elementid = 2
  [../]
  [./viscosity0]
    type = ElementalVariableValue
    variable = viscosity
    elementid = 0
  [../]
  [./viscosity1]
    type = ElementalVariableValue
    variable = viscosity
    elementid = 1
  [../]
  [./viscosity2]
    type = ElementalVariableValue
    variable = viscosity
    elementid = 2
  [../]
  [./k0]
    type = ElementalVariableValue
    variable = k
    elementid = 0
  [../]
  [./k1]
    type = ElementalVariableValue
    variable = k
    elementid = 1
  [../]
  [./k2]
    type = ElementalVariableValue
    variable = k
    elementid = 2
  [../]
[]

[Outputs]
  csv = true
[]
