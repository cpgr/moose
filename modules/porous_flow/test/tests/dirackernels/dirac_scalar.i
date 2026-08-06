# A point sink whose strength is a scalar unknown, constrained to deliver a prescribed rate.
# The scalar has no ScalarKernel: its only equation is assembled by the DiracKernel.
#
# The constraint is flux_coefficient*(u(point) - bhp) = rate, so at convergence
# bhp = u(0.5) - rate/flux_coefficient = u(0.5) - 0.05
[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 4
  xmin = 0
  xmax = 1
[]

[Variables]
  [u]
  []
  [bhp]
    family = SCALAR
    order = FIRST
  []
[]

[Kernels]
  [diff]
    type = Diffusion
    variable = u
  []
[]

[DiracKernels]
  [source]
    type = ConstrainedPointSource
    variable = u
    scalar_variable = bhp
    point = '0.5 0 0'
    flux_coefficient = 2.0
    rate = 0.1
  []
[]

[BCs]
  [sides]
    type = DirichletBC
    variable = u
    boundary = 'left right'
    value = 1
  []
[]

[ScalarKernels]
  [bhp_null]
    # checkKernelCoverage builds its list of covered variables from Kernels, ScalarKernels, FV
    # objects and IntegratedBCs, but never from DiracKernels (idaholab/moose#9325), so 'bhp' looks
    # uncovered even though the DiracKernel above assembles its equation.  NullScalarKernel is the
    # framework's remedy for exactly this: it contributes a zero residual and, with the default
    # jacobian_fill of zero, a zero Jacobian, so it changes nothing that is solved.  Note this is
    # preferable to setting Problem/kernel_coverage_check=false, which would switch the check off
    # for every variable in the input rather than just this one.
    type = NullScalarKernel
    variable = bhp
  []
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Steady
  solve_type = NEWTON
[]

[Postprocessors]
  [bhp_value]
    type = ScalarVariable
    variable = bhp
  []
  [u_mid]
    type = PointValue
    variable = u
    point = '0.5 0 0'
  []
[]

[Outputs]
  csv = true
[]
