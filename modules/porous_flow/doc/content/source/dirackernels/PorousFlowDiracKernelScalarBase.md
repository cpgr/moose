# PorousFlowDiracKernelScalarBase

A base class for `DiracKernels` that, as well as contributing to their own field variable, assemble
the row of a coupled scalar variable. The scalar is named `scalar_variable` in the input file and
`kappa` in the source, following [KernelScalarBase](KernelScalarBase.md), of which this class is the
point-source analogue.

It exists so that a line sink can carry a bottom-hole pressure that is itself a Newton unknown: the
well then delivers a prescribed rate, and the pressure that achieves it is solved for implicitly
along with the reservoir equations, rather than chased by an outer controller.

Supplying no `scalar_variable` leaves the class behaving exactly as a plain `DiracKernel`, so it is
safe as the base of point sources that do not need a scalar.

## Why a scalar variable

A reservoir simulator normally eliminates the bottom-hole pressure rather than carrying it as an
unknown. The constraint is linear in it, so no iteration is needed to do so:

!equation
P_{bh} = \frac{\sum_c W_c \lambda_c P_c - q_{target}}{\sum_c W_c \lambda_c}

That route is attractive. It needs far less sparsity, since eliminating the pressure couples the
completions only to each other rather than adding a row and column spanning the whole system, and
the maximum-pressure limit becomes a clamp evaluated inside the residual rather than a constraint
that has to be swapped in and out between Newton iterations.

It is not used here because the sums run over *every* completion and are needed before *any*
completion's residual can be assembled, whereas point sources are assembled in a single pass over
the elements. Only the kernel can form `W*lambda`, and only during that pass. Computing the sums in
a separate object would mean a second implementation of the well index and mobility, and lagging
them by an iteration would give up the exact rate matching that motivated solving for the pressure
implicitly in the first place. Eliminating the pressure also couples every completion to every
other, and those entries lie outside the sparsity pattern, so it would additionally need a
relationship manager to ghost and preallocate them.

Carrying the pressure as a scalar unknown asks the framework for that cross-element coupling
instead. The dense row and column libMesh allocates for a scalar degree of freedom are the price of
it, which is around 8% of the matrix per well on a single-variable model, and roughly half that on
a two-component one. That fraction is independent of mesh size, so it is a limit on the number of
wells rather than on the size of the model, and elimination only becomes worth its extra machinery
somewhere beyond a few dozen wells.

## How the assembly works

The Dirac threading only ever calls `computeResidual()`, `computeJacobian()` and
`computeOffDiagJacobian()` on a `DiracKernel`. Its Jacobian loop walks
`Assembly::couplingEntries()`, which holds pairs of *field* variables, so
`computeOffDiagJacobianScalar()` is never reached for a `DiracKernel`. All four blocks are therefore
assembled from inside the three calls that are made, and no framework loop needs to change:

| Block | Assembled in | Derived-class hook |
| - | - | - |
| field-field | `computeJacobian` | `computeQpJacobian` |
| field-kappa | `computeJacobian` | `computeQpJacobianScalar` |
| kappa-field | `computeOffDiagJacobian` | `computeQpScalarOffDiagJacobian` |
| kappa-kappa | `computeJacobian` | `computeQpScalarJacobian` |

The kappa row of the residual is flushed to the global vector by `computeResidual()` itself. The
Dirac Jacobian loop flushes cached Jacobians once it has finished, but the residual loop has no such
step: it only calls `Assembly::addResidual`, which visits field variables alone. Flushing per element
does not double count, because `Assembly::addCachedResidualDirectly` clears the cache as it adds.

Where a well is expanded into one kernel per fluid component, every kernel assembles its own
field rows and its own field-kappa column, but only one may own the shared constraint. Set
`compute_scalar_residuals = false` on the others.

## Kernel coverage check

The scalar variable's only equation comes from a `DiracKernel`, and
`NonlinearSystemBase::checkKernelCoverage` builds its list of covered variables from `Kernels`,
`ScalarKernels`, finite volume objects and `IntegratedBCs`, never from `DiracKernels`
([idaholab/moose#9325](https://github.com/idaholab/moose/issues/9325)). The scalar therefore looks
uncovered. The check then calls `getSubdomainsForVar` on it, which only knows field variables, so
the resulting error is the misleading "Variable ... does not exist in this system"
([idaholab/moose#20674](https://github.com/idaholab/moose/issues/20674)); that second part happens
for any scalar variable lacking a `ScalarKernel`, whether or not a `DiracKernel` is involved.

Add a [NullScalarKernel](NullScalarKernel.md) on the scalar, which is the framework's own remedy for
a variable the coverage check cannot see:

```
[ScalarKernels]
  [bhp_null]
    type = NullScalarKernel
    variable = bhp
  []
[]
```

It contributes a zero residual, and with the default `jacobian_fill` of zero a zero Jacobian, so
nothing that is solved changes. Prefer it to `Problem/kernel_coverage_check = false`, which switches
the check off for every variable in the input rather than just the one it cannot see.

Should `checkKernelCoverage` ever learn about `DiracKernels`, the null kernel becomes unnecessary.
The mechanism for that already exists: `ResidualObject::additionalROVariables()` lets an object
declare variables it covers beyond its own, `IntegratedBCs` are already consulted through it, and
the finite volume Lagrange multiplier constraints do the same thing by hand a few lines earlier in
that function. Only a call to `_dirac_kernels.subdomainsCovered()` is missing.
