# Numerical evolution, neutrino fitting and validation

The numerical calculation connects model couplings and threshold-generated operators to low-energy neutrino masses and mixing. Its physical scope is determined by the source coefficients and interactions that actually enter the numerical solver, not by the apparent completeness of a file name or a successful optimisation status.

## Data passed between effective theories

A reproducible calculation specifies the representation, initial scale and couplings, physical masses, fermion and scalar matching scales, RG export provenance, operator basis, tensor shape/complexity, normalisation and unit system. Scalar and fermion blocks may change dimension at physical thresholds. Any numerical state conversion must preserve the intended flavour symmetry and scalar-basis convention.

The Weinberg coefficient is a complex symmetric \(3\times3\) matrix, not a general Hermitian matrix. In the documented normalisation, \(m_\nu=-v^2 C_5/2\). Check the resulting Takagi diagonalisation residual, mass-square differences and lepton-mixing extraction against the charged-lepton basis.

## Independent numerical checks

The relevant tests include:

1. Round-trip packing and unpacking of complex UV and scalar-EFT states, including optional permitted bilinear blocks.
2. Agreement between symbolic and numerical beta-function normalisations and derivative conventions.
3. Threshold handoff shape, units and source hash checks.
4. Perturbative-order and operator-inventory checks; numerical agreement cannot repair missing physics.
5. ODE convergence under tighter tolerance and sensitivity to scale choices.
6. Takagi reconstruction, neutrino mass dimensions and mixing angle conventions.
7. Comparison with saved benchmarks only when their source inputs and approximations match.

Synthetic unit tests establish interface behaviour and algebraic identities. A recorded local model run checks that particular numerical data and assumptions. Neither one certifies complete one-loop EFT matching for every representation.

## Fits and scientific interpretation

A Yukawa fit that reproduces chosen oscillation measurements is **not** a parameter-free prediction. The objective function, fitted variables, accepted/failed cases and fixed benchmark assumptions must be reported alongside the observables. Keep model-dependent fitting separate from transfers of a fixed Yukawa matrix between models.

Selected historical studies used conditional Wilson coefficients, restricted scalar interactions and frozen input scales. They remain useful controlled comparisons, but must not be presented as fully general, physically renormalised results without the missing matching and operator evolution being verified.

For the currently documented capability boundaries consult [implementation status](IMPLEMENTATION_STATUS.md). For the actual saved numerical values and input assumptions consult [results](RESULTS.md) and the individual manifests. For reproducing existing outputs consult [reproducibility](REPRODUCIBILITY.md).
