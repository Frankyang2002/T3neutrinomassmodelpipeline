# Batch 7 — A_-4 native source scalar Hessian and SU(2) weight-basis reconciliation

**Status:** source-driven numerical audit, **not physical UV-to-SMEFT matching**.

## Motivation

Batch 5 directly evaluates the one-loop UV mass with the *source-audited SU(2)
weight-basis Yukawa vertices*. Batch 6 uses native Matchete source quartics to
construct the heavy-scalar Cartesian Hessian, then passes that matrix to the UV
loop evaluator. The two scalar bases need not be identical. The exact
Batch6/Batch5 `1/sqrt(2)` spectrum ratio is therefore **not evidence by itself**
for a scalar Lagrangian normalization change or a value of `N_GammaCC`.

The current repository's original generator analysis gives a nontrivial
triplet rotation with columns in the weight basis (native-component rows):

```
 U_3 = [ -i/sqrt(2), 0, +i/sqrt(2) ]
       [     0,      1,       0     ]
       [ +1/sqrt(2), 0, +1/sqrt(2) ]
```

This unitary matrix is **loaded afresh from the checked local audit source**;
it is not inserted as a hard-coded numerical transformation. The scalar-field
replacement is `S_native = U S_weight`, for `U = diag(1, U_3)`. In the complex
scalar quadratic convention

```
 V2 = S_native^dagger A_native S_native + Re(S_native^T B_native S_native)
 A_weight = U^dagger A_native U
 B_weight = U^T B_native U
```

Batch 7 performs this unitary conversion before using the *unchanged* direct
UV Yukawa vertices. Both `A` and possible holomorphic `B` are rotated; neither
matrix is merely transposed or conjugated by convention.

## Direct native potential differentiation (independent route)

The Batch 6 Hessian is computed from fully symmetrized `lambda_ABab` and
`Delta(M2)_ab = v^2/2 * lambda_[H0R,H0R,a,b]`. Batch 7 independently:

1. Reads the original `A_-4` Matchete scalar-quartic seed and verifies the
   exact four HHSS terms (`Index=2,3,4,5`) and their native couplings.
2. Parses each **original** `TermInputForm` and evaluates its *native CGs*
   using the existing trusted source registry/parsing utilities. It does
   **not** call the quartic tensor builder or the quartic factorial adapter.
3. Inserts `H+=(0)` and `H0=v/sqrt(2)` directly in `V=-L` and expands the two
   remaining scalar fields as `(S_R+i S_I)/sqrt(2)` (with conjugation where
   specified by `Bar[Field[...]]`).
4. Analytically differentiates each resulting quadratic monomial. The four
   term Hessians are added only once, retaining the explicit Hermitian
   conjugate of the mixed quartic.
5. Checks this direct-source Hessian against Batch 6. A disagreement is
   recorded in the release status; it is **not** hidden or overwritten.
6. Converts the native real Hessian to `A_native`, `B_native`, then to the
   generator-audited UV weight basis, and re-evaluates the UV neutrino mass
   using both existing independent and production loop evaluators.

The common unbroken scalar mass parameters are included exactly once, and
all fields share the same `v` from the pinned fitted benchmark. No new
numerical couplings, SU(2) CG coefficients, signs or Yukawa conventions are
introduced. `lambdaS2Adj` remains classified as a four-heavy-scalar quartic
with no tree-level Higgs-VEV mass shift at vanishing heavy-scalar background.

## Branches and outputs

Batch 7 retains Batch 6's **two unverified alternatives**, not a new fit:

- `adjoint_alias`: native `lambdaH2Inv2` assigned to RGBeta `lambdaH2Adj`.
- `ordinary_alias`: native `lambdaH2Inv2` assigned to RGBeta `lambdaH2`.

Both are evaluated at the three existing labelled scalar stages (`UV`, `F`,
`S`). In **all** stages, the same dense fitted heavy-fermion state at `muF`
is retained. `UV` and `S` are intentionally *hybrid scale diagnostics*. They
are not silently identified with a physical common-scale theory.

Each of six scenarios writes:

- `aligned_UV_input.json` — complete weight-basis complex scalar A/B and
  original direct-UV fermion/Yukawa vertices.
- `independent_direct_UV.json` and `production_direct_UV.json` — two UV loop
  evaluations with pole cancellation and agreement checks.
- `native_Hessian_interleaved_GeV2.json` and
  `native_source_HHSS_Hessian_shift_GeV2.json` — direct source contraction.
- `weight_basis_scalar_A_gev2.json` and
  `weight_basis_scalar_B_gev2.json` — complex matrices actually supplied to
  the direct-UV evaluator.

Root output: `output/general_observable/lo_certification_benchmark/batch7_native_basis_<run-id>/`,
containing `summary.json` and `report.md`.

Numerical gates are intentionally separated:

- **Source derivative vs Batch6**: relative Hessian disagreement is reported;
  the root `status` reflects whether all six agree. This test does not itself
  decide native-to-RGBeta alias certification.
- **Rotation**: unitarity, Hermiticity, symmetry and exact inverse rotation
  checks must pass.
- **Direct loop**: finite UV pole, independent-vs-production agreement and
  positive scalar mass spectrum must pass.
- **EFT comparison**: use the Batch4 conditional complex mass matrix, fitting
  at most one free complex *diagnostic* multiplier. It is not `GammaCC`.

## Installation and execution

From PowerShell in the full original checkout:

```powershell
$repo = "C:\thesistemp\T3neutrinomassmodelpipeline"
$zip = "$env:USERPROFILE\Downloads\T3_Batch7_Native_Scalar_Basis.zip"
if (-not (Test-Path $zip)) { throw "ZIP not found: $zip" }
Expand-Archive -LiteralPath $zip -DestinationPath $repo -Force
Set-Location $repo
$env:PYTHONPATH = (Get-Location).Path

python -X utf8 -m unittest discover -s tests -p 'test_lo_native_scalar_basis_batch7.py' -v

python -X utf8 scripts\RunLONativeScalarBasisBatch7.py --check-only `
  --batch4-run-id source01 --batch5-run-id source01 --batch6-run-id source01

python -X utf8 scripts\RunLONativeScalarBasisBatch7.py --run-id source01 `
  --batch4-run-id source01 --batch5-run-id source01 --batch6-run-id source01
```

A rerun must choose a fresh `--run-id source02`; no prior output is overwritten.

All entries in the ZIP are **complete new files** with exact project-relative
paths:

```
Numerical/lo_certification/benchmark/native_scalar_basis_batch7.py
Numerical/lo_certification/benchmark/README_BATCH7.md
scripts/RunLONativeScalarBasisBatch7.py
tests/test_lo_native_scalar_basis_batch7.py
```

The three preceding batch implementations and all original GitHub/repository
files remain unchanged and read-only.

## Verification and physical limits

The isolated synthetic tests cover direct original-source monomial
Hessian normalization, Hermitian-conjugate mixed terms, complex portal phases,
nontrivial unitary triplet basis changes, holomorphic B rotation and reverse
map, scalar spectra, required source inventory, false-certification guards,
output isolation and the grouped EFT order.

The full **source-backed** parser, actual Batch6 Hessian comparison and six
new UV loop evaluations require your complete local Windows checkout; these
have not been run in the packaging container.

Still **uncertified**: the absolute signed native `GammaCC`-to-two-component
Weyl normalization, Matchete-to-RGBeta quartic parameter aliases, hard one-loop
heavy-fermion matching, full EFT1 Wilson anomalous-dimension mixing and a full
physical one-loop UV/EFT matching comparison. The complete EFT order is:

`UV -> integrate out ALL F together at muF -> EFT1(S1,S2) -> integrate out BOTH S together at muS -> SMEFT`.
