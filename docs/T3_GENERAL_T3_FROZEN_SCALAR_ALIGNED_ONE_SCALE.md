# T3 B_-1/C_-1: frozen scalar diagonal-aligned direct C5 benchmark

## Why

The earlier hybrid reference-Yukawa evaluation uses scalar masses from a **different**
experimental quadratic benchmark. The checked frozen fitted-reference configurations
instead specify (at their saved base scale) `mS1Sq` and `mS2Sq`, and do not specify
the new gauge-invariant scalar mixing `mu2`. Mixing **cannot** be inferred from
`lambdaT3` or the pre-existing exploratory mass-basis input.

This batch deliberately evaluates a **scenario**, not a refit:

- frozen 3x3 complex Yukawa matrices and diagonal fermion masses from
  `Numerical/benchmarks/fitted_reference/..._alpha_m1.json`;
- frozen scalar **diagonal mass squared** inputs exactly, with independent
  complex mass mixing `mu2` in GeV^2; default `mu2 = 0` explicitly assumed;
- `lambdaT3` transferred directly; the two quartics used by the exact
  one-loop coefficient are retained at their frozen **zero** values only.
  Nonzero values fail closed, since their Matchete-to-RGBeta normalization
  has not been certified;
- the **original Matchete-produced one-loop ordered coefficient** is evaluated
  with the new *verified* mass eigenvalues and mixing matrix, not replaced
  by a postulated analytic T3 loop formula.

The scalar convention is

```
Phi=(S1, epsilon Conjugate[S2]); epsilon={{0,1},{-1,0}}
A = {{mS1Sq,mu2},{Conjugate[mu2],mS2Sq}}
U^dagger A U = DiagonalMatrix[{mX^2,mY^2}], mX <= mY
Phi = U.{X,Y}
```

For `mu2 = 0` and frozen `mS1Sq > mS2Sq`, the sorted eigenbasis is
`X = epsilon Conjugate[S2]`, `Y = S1` with `U={{0,1},{1,0}}`, not
`U=IdentityMatrix[2]`.

## Installation paths

This ZIP contains **full files**, not diffs:

| Relative project path | Action |
|---|---|
| `scripts/RunGeneralT3NumericalC5.py` | Full replacement; adds strictly validated scalar-basis override |
| `scripts/RunGeneralT3FrozenScalarC5.py` | New scenario runner |
| `tests/python/test_general_t3_frozen_scalar_c5.py` | New integration and fail-closed tests |
| `docs/T3_GENERAL_T3_FROZEN_SCALAR_ALIGNED_ONE_SCALE.md` | New documentation |

It does **not** overwrite the one-loop Matchete coefficient, mass-basis
validator, or frozen fitted reference configurations.

## Usage (PowerShell, repository root)

```powershell
Expand-Archive -LiteralPath "C:\path\to\T3_General_Frozen_Scalar_Aligned_C5.zip" -DestinationPath . -Force
python -m pytest -q tests\python\test_general_t3_frozen_scalar_c5.py tests\python\test_general_t3_reference_c5.py tests\python\test_general_t3_numerical_c5.py
python -m scripts.RunGeneralT3FrozenScalarC5 --all
```

Optional *independent* mixing scenario, replacing the `mu2 = 0` default:

```powershell
python -m scripts.RunGeneralT3FrozenScalarC5 --all --mu2-real-gev2 10000 --mu2-imag-gev2 0
```

Using any specified `mu2` does **not** imply it was fitted or fixed by
existing EFT1 RGEs. Tachyonic/degenerate masses and broken input/hash
certificates are rejected.

Output directory (does not overwrite hybrid comparisons):
`output/general_t3_numerical_c5_frozen_scalar_scenario/`.

The ordered flavour coefficient is evaluated for each diagonal heavy fermion
mass and the symmetric Weinberg coefficient uses `C5 = A5 + Transpose[A5]`.
`mnu=-v^2 C5/2` and its Takagi singular values use the same conventions as
`RunGeneralT3NumericalC5.py`.

## Interpretation and remaining work

This is a **one-scale simultaneous all-heavy** result. It has correct *frozen
scalar diagonals* and flavour Yukawa input but the quadratic mixing is an
**explicit scenario**, not a frozen fit parameter. The saved Yukawa/coupling
values are taken from `base_state`, **not** evolved to each threshold in the
general EFT1 pipeline. It is **not** sequential fermion-first matching,
not a replacement for the frozen fit, and does not certify the physical
Matchete/RGBeta quartic normalization for nonzero extra couplings. In
particular, the result is **not** one-loop-matched EFT1 at the low scale.

After the scenario comparison, the next physics task remains integrating the
matched eigenbasis coefficients into the fermion-first EFT1 matching and
running sequence, with consistent mixing-quadratic running and thresholds.
