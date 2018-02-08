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
    filename = sample.dat
    primary_species = 'H+ Ca++ HCO3-'
    secondary_species = 'OH- CaOH+'
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
