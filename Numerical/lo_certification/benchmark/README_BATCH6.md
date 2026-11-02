# Batch 6 — `A_-4` native HHSS EWSB scalar Hessian and candidate alias branches

**Status: conditional diagnostics only. Not a physical neutrino-mass calculation.**

## What is computed

Batch 6 consumes the successful, **unaltered** `batch4_common_input_source01` and
`batch5_uv_eft_source01` outputs. It imports the existing Matchete source-derived
quartic parser (`Numerical/joint_scalar_matching/source_quartic_v4.py`), original
scalar seed, two independent existing UV one-loop numerical evaluators, and Batch 5
matrix comparison. The fitted Yukawas, mass-basis transformations, heavy masses,
and original source CGs all remain those of Batch 4 and 5.

The native source's real-scalar potential convention is

\[
 V_4 = \frac{1}{4!}\lambda_{ABab}h_Ah_Bs_as_b+\cdots,
 \quad V=-\mathcal L\text{ on source scalar quartics}.
\]

With `H=(H+,H0)`, `H0=(v+h+iG)/sqrt(2)` and the project real Higgs axes
`(H+_R,H+_I,H0_R,H0_I)`, the mass Hessian at zero heavy-scalar background is

\[
 \left(M_s^2\right)_{ab}
 =\operatorname{diag}(m_{S_1}^2,m_{S_1}^2,
 m_{S_2}^2,\ldots,m_{S_2}^2)_{ab}
 +\frac{v^2}{2}\lambda_{H0_R H0_R ab}.
\]

The fully sourced 8x8 Hessian in the interleaved real-field basis is inverted
algebraically into the UV evaluator's four-complex-field Hermitian `A` and
symmetric `B`, verifying that it reconstructs the original Cartesian Hessian.
The old provisional Batch 5 `scalar_A_gev2` and `scalar_B_gev2` are **replaced**,
not incrementally modified or double-counted. The Matchete-derived quartic
normalization is read directly. The *identification* of a Matchete coupling
with an RGBeta numerical parameter is **not** certified.

For `A_-4`, the source HHSS inventory contains:

| Matchete coupling | Candidate RGBeta parameter |
|---|---|
| `lambdaH1Inv1` | `lambdaH1` in both branches |
| `lambdaT3` | `lambdaT3` in both branches |
| `lambdaH2Inv2` | **`lambdaH2Adj` OR `lambdaH2`**, two explicit branches |

The source does not independently certify which `lambdaH2Inv2` identification
is physical. Both variants are retained. These two alternatives are **not**
two contributions to be added together. In particular, the `adjoint_alias`
branch does not include `lambdaH2` as an additional portal; `ordinary_alias`
does not include `lambdaH2Adj` as an additional portal.

The two Batch 5 projected nonzero operators are classified separately:

* `lambdaH2Adj`: **Higgs–scalar quartic**, a candidate HHSS portal term entering
  an EWSB scalar Hessian depending on the unresolved coupling identification.
* `lambdaS2Adj`: **pure-scalar quartic**. At zero heavy-scalar background it has
  no tree-level Higgs-VEV quadratic contribution. Its effects on the evolved
  running scalar parameters are already present in the input states; its
  separate full-theory higher-order contributions are not computed.

The input stages are `UV` (frozen scalar parameters), `F` (scalar parameters
at `mu_F`), and `S` (EFT1 scalar parameters at `mu_S`). All six combinations
(three stages times two alias branches) use the **same fermion sector at**
`mu_F=1e5 GeV`. `UV` and `S` are labelled cross-scale/hybrid diagnostics.
For each, the package computes an independent UV one-loop mass matrix, verifies
against the repository production UV evaluator, checks the UV divergence and
fermion/scalar hierarchy, and compares with the unchanged Batch 4 conditional
EFT mass matrix **after one arbitrary complex scaling**. The matrix-shape fit
is not a determination of the physical `GammaCC` factor.

## Installation

Extract `T3_Batch6_Native_Portal_Hessian.zip` at the repository root
`C:\thesistemp\T3neutrinomassmodelpipeline`. The ZIP contains these **complete,
new, project-relative files**:

- `Numerical/lo_certification/benchmark/native_portal_batch6.py`
- `Numerical/lo_certification/benchmark/README_BATCH6.md`
- `scripts/RunLONativePortalBatch6.py`
- `tests/test_lo_native_portal_batch6.py`

No original source, tests, model code, or existing Batch 4/5 output is replaced.

```powershell
$repo = 'C:\thesistemp\T3neutrinomassmodelpipeline'
$zip = "$env:USERPROFILE\Downloads\T3_Batch6_Native_Portal_Hessian.zip"
if (-not (Test-Path $zip)) { throw "ZIP not found: $zip" }
Expand-Archive -LiteralPath $zip -DestinationPath $repo -Force
Set-Location $repo
$env:PYTHONPATH = (Get-Location).Path

python -X utf8 -m unittest discover -s tests -p 'test_lo_native_portal_batch6.py' -v
python -X utf8 scripts\RunLONativePortalBatch6.py --check-only `
    --batch4-run-id source01 --batch5-run-id source01
python -X utf8 scripts\RunLONativePortalBatch6.py --run-id source01 `
    --batch4-run-id source01 --batch5-run-id source01
```

If an output directory already exists, choose `--run-id source02`. Outputs are
immutable, and all calculations are assembled before publishing the directory.

### Output locations

`output/general_observable/lo_certification_benchmark/batch6_portal_source01/`

- `summary.json` — six UV/EFT spectra and residuals, input hashes, native source
  channel list, alias assumptions, pure-scalar classification, diagnostics.
- `report.md` — concise stage-by-stage numerical comparison.
- `UV_adjoint_alias/`, `UV_ordinary_alias/`, `F_.../`, `S_.../` — each contains
  `input.json`, `independent.json`, `production.json`,
  `Hessian_real_interleaved_GeV2.json`, and
  `delta_Hessian_from_HHSS_GeV2.json`.

## Verification status

Synthetic regression suite: 19 passing tests, one skipped when the original
Windows checkout and Batch 4/5 reports are unavailable. These tests verify
Hessian factorial, Higgs index, real/complex field ordering, Hermitian `A`,
symmetric `B`, source branch selection, pure-scalar classification, fail-closed
omissions, stage separation and no physical-certification promotion.

**Not tested in this packaging environment:** the full Windows native source
parser, the six UV loop calculations, and compatibility with actual Batch 4/5
outputs. Run the commands above and inspect exact failures if any. The
script does not manufacture a numerical result before your local run.

No physical source-to-Weyl normalization, full one-loop hard-fermion threshold,
Matchete/RGBeta parameter normalization, full EFT1 Wilson RGE, or complete
one-loop UV–EFT matching is asserted. The logical EFT sequence is always:

`UV -> integrate out ALL heavy fermions at muF -> EFT1(S1,S2) -> integrate out BOTH scalar multiplets at one muS -> SMEFT`.
