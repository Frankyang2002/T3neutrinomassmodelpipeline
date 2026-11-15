# Calculation guide

The purpose of the code is to connect the electroweak quantum numbers and interactions of a T3 radiative mass model to neutrino masses and mixing parameters. This guide describes that scientific calculation; it does not assert that every possible interaction has already been connected to a complete numerical result.

## Choose the model

The ordinary T3 catalogue contains 16 electroweak assignments from five families:

| Family | \((d_{S_1},d_{S_2},d_F)\) | \(\alpha\) values |
| --- | --- | --- |
| A | (1, 3, 2) | −4, −2, 0 |
| B | (2, 2, 1) | −3, −1, +1 |
| C | (2, 2, 3) | −3, −1, +1 |
| D | (3, 1, 2) | −2, 0, +2 |
| E | (3, 3, 2) | −4, −2, 0, +2 |

The project's hypercharge convention is \(Q=T_3+Y\), with
\[
Y(S_1)=\frac{\alpha}{2},\qquad
Y(S_2)=\frac{\alpha+2}{2},\qquad
Y(F)=\frac{\alpha+1}{2}.
\]

The new T3 fields are odd under the imposed \(\mathbb Z_2\) symmetry. The complete renormalisable theory must include every independent Lorentz-, gauge- and \(\mathbb Z_2\)-invariant term for the selected representations. A chosen benchmark may set allowed couplings to zero, but that choice must be reported.

## Connect high-energy physics to neutrino observables

**UV theory.** Construct the model's Lagrangian and evolve its couplings from their specified input scale. The heavy fermion and scalar mass matrices, Yukawa matrices and scalar couplings determine the boundary data.

**Fermion decoupling.** Integrate out *all heavy fermions*, matching onto an effective theory that still contains the heavy scalars. If there are distinct fermion masses, apply the required matching at the corresponding scales.

**Scalar effective theory.** Evolve its renormalisable interactions and induced operator coefficients. This running must use the operator basis and conventions appropriate to the selected model and perturbative accuracy.

**Scalar decoupling.** Integrate out *all heavy scalars*, including effects of scalar mixing where allowed. The remaining low-energy theory is SMEFT, with the Weinberg coefficient and any other operators retained at the specified order.

**Low-energy observables.** Evolve the Weinberg coefficient and form the symmetric Majorana mass matrix using the documented convention
\[
m_\nu=-\frac{v^2}{2}C_5.
\]
The neutrino masses follow from Takagi factorisation. PMNS parameters also require the consistent charged-lepton flavour basis.

## Inputs and outputs

A reproducible calculation identifies the model and hypercharge, all relevant coupling choices, fermion and scalar masses, input and matching scales, RGE and matching source files, solver tolerances and numerical conventions. Numerical studies can produce neutrino masses, \(\Delta m^2_{21}\), \(\Delta m^2_{31}\), PMNS mixing parameters and comparisons under controlled changes of model inputs.

An optimised Yukawa matrix is a fitted input. Agreement with oscillation data under that fit is not an independent prediction of the fitted observables.

## How to use the available programs

From the repository root:

```powershell
# Available combined calculation options.
python -X utf8 pipeline.py --help

# Model-catalogue symbolic calculation options.
python -X utf8 full16_symbolic.py --help

# Inspect an existing numerical scan configuration without performing a scan.
python -X utf8 -m Numerical.fitting.ScanCLI path\to\existing_scan.json --dry-run

# Core regression checks.
.\scripts\RunCoreTests.ps1
```

New symbolic calculations require Wolfram Language, Matchete and RGBeta; consult [Installation](INSTALLATION.md). Symbolic expressions, saved-source numerical calculations and complete physical matching must not be treated as interchangeable guarantees of accuracy.

## Evidence and limitations

The repository has useful working components for model construction, running, matching and fitting. The available source does **not** verify one fully general, complete-one-loop, fermion-first-to-SMEFT numerical fit for all 16 models. Additional interactions and flavour structures are particularly relevant to B(−1) and C(−1). Evaluate a result against the actual model configuration and source files, not solely its final fit residual or an output folder's name.

For detailed scientific conventions see [T3 models and EFT prescription](MODEL_AND_PIPELINE.md); for the precise boundary of supported calculations see [Capabilities and limitations](IMPLEMENTATION_STATUS.md). For saved studies see [Results overview](RESULTS_OVERVIEW.md) and [Reproducibility](REPRODUCIBILITY.md).
