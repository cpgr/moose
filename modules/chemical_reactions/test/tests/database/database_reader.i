[Mesh]
  type = GeneratedMesh
  dim = 2
[]

[Executioner]
  type = Transient
  end_time = 1
[]

[ReactionNetwork]
   primary_species = 'H+ Ca++ HCO3-'
  [./GeochemicalDatabase]
    filename = ../../../database/eq36.dat
    primary_species = 'H+ Ca++ HCO3-'
    secondary_species = 'OH- CaOH+ CO3-- CaHCO3+ CaCO3(aq) CO2(aq)'
    mineral_species = 'Calcite'
  [../]
[]

[ICs]
  [./H+]
    type = ConstantIC
    variable = H+
    value = 1e-7
  [../]
  [./Ca++]
    type = ConstantIC
    variable = Ca++
    value = 1e-7
  [../]
  [./HCO3-]
    type = ConstantIC
    variable = HCO3-
    value = 1e-7
  [../]
[]

[Materials]
  [./porous]
    type = GenericConstantMaterial
    prop_names = 'porosity'
    prop_values = '0.1'
  [../]
[]

[Outputs]
  exodus = true
[]
