# Batch v6 — all-16 exploratory joint-scalar benchmark

**Scope:** produce one conditional complex `3×3` LLSS-induced `C5` estimate per unknown Matchete `N_GammaCC`, with rough mass singular-value scales, for the original eight plus eight additional mixed/split scalar T3 models. Numerical results are **not physically certified** and are **not complete one-loop neutrino masses**.

## Existing source code reused

- `Numerical/experimental/T3SelfBilinearMassBasis.py`: six real self-bilinear splits, including spherical SU(2) metric, intertwiner `W` and phase conventions. The potential convention is `mu2 (Sdag)^T C (Sdag) + h.c.` with no extra `1/2`: `m_plus^2=m^2+2|mu2|`, `m_minus^2=m^2-2|mu2|`.
- `Numerical/experimental/GeneralT3MixedScalarMassBasis.py`: two mixed-doublet cases `B_-1`, `C_-1` with `Phi=(S1,epsilon conj(S2))`, `epsilon12=+1`, `U† [[m1²,mu²],[mu²*,m2²]] U=diag(mX²,mY²)`.
- `Numerical/joint_scalar_matching/numerical_benchmark_v5.py`: eight ordinary model numerical backgrounds and conditional matching.
- Source-resolved fermion release + v4 Matchete-native `HHSS` channels + existing 1 TeV renormalisable EFT1 snapshots. v2/v3 kernel/parity checks are still used.

Only **four new files** are installed. No existing source code is replaced. No GitHub changes.

## Required action

Run in PowerShell from project root after downloading ZIP:

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
$zip = Get-ChildItem "$HOME\Downloads" -Filter "T3_joint_scalar_all16_exploratory_v6*.zip" -File | Sort-Object LastWriteTime -Descending | Select-Object -First 1
if (-not $zip) { throw "Download the v6 ZIP into Downloads first." }
Expand-Archive -LiteralPath $zip.FullName -DestinationPath . -Force
python -m unittest discover -s tests -p "test_t3_joint_scalar_*py" -v
python -m scripts.RunT3JointScalarAll16V6 --allow-exploratory-map --allow-exploratory-special-bases --output-dir "output\general_observable\lo_joint_scalar_matching\exploratory_all16_v6"
```

Results:

```text
output/general_observable/lo_joint_scalar_matching/exploratory_all16_v6/
  summary.json
  all16_conditional_mass_comparison.csv
  numerical_backgrounds/<model>.json
  numerical_backgrounds/<model>.npz
  C5_per_unknown_NGamma_exploratory/<model>.json
```

The summary reports **actual** model counts. It does not claim all 16 before the local run confirms it.

## Source matching / scalar bases

All eight special source layouts were independently checked against the actual all-F-first Matchete release.

| Models | Boundary source scalar layout | Numerical quadratic parameter at `muS` |
|---|---|---|
| `A_+0`, `D_+0`, `E_+0` | `T3Scalar1Minus`, `T3Scalar1Plus`, `NewScalar2` | `muT3GeneralS1Sq` |
| `A_-2`, `D_-2`, `E_-2` | `NewScalar1`, `T3Scalar2Minus`, `T3Scalar2Plus` | `muT3GeneralS2Sq` |
| `B_-1`, `C_-1` | `NewScalarX`, `NewScalarY` | `muT3GeneralS1S2` |

The `mu2` parameter and both diagonal masses come **only from** `output/general_observable/all16_v1/<slug>/running.json` / `eft1_at_scalar_threshold` (`muS=1000 GeV`). Missing or non-positive squared masses block the model. The unbroken scalar mass eigenvalues are calculated, not replaced by unsplit masses.

The v4 native `HHSS` quartics are built in **original** `(NewScalar1,NewScalar2)` real component order. Their scalar legs are transformed into the eight special source-boundary field layouts by a real orthogonal matrix `R`, with

`source_original_real_fields = R @ boundary_special_real_fields`.

`lambda_boundary = R^T lambda_source R` in scalar indices. The fermion LLSS boundary is already in the special basis, so it is **not rotated a second time**. Both scalars are integrated jointly at the common `muS`.

## Explicit exploratory assumptions

1. The v5 candidate native quartic aliases (`lambdaH1Inv1` → `lambdaH1`, etc.) **may not have the physical coefficient normalization**, and are not independently mapped to RGBeta components.
2. The experimental spherical/component SU(2) mass basis is **assumed** to equal the original Matchete CG component basis. This is not verified.
3. The special scalar mass eigenfields used in the fermion-first source boundary at `muF` are **identified** with the eigenfields at `muS`; no separate running rotation of the source tensor is computed.
4. The signed `N_GammaCC_to_Weyl` is **unknown**. Results are `C5_per_N`, and estimated mass scales are `eV_per_abs_N` (not actual eV masses).
5. We include only the existing fixed-order LLSS-induced EFT1 log plus joint finite scalar bubble, **not** all hard heavy-fermion one-loop matching or the complete operator-mixing RGE.
6. The v5 scalar trajectory is illustrative. These numbers are not fitted physical values.

If full tensor contraction fails the existing SU(2) Weinberg projection, the affected model is **reported as an error** rather than silently being declared a Weinberg coefficient. This matters even for rough numerical estimates: a non-SU(2)-consistent amplitude cannot unambiguously be called `C5`.

## Mass estimate convention

`mnu_rough = (246.22 GeV)^2 * C5 / 2`. Singular values of the complex symmetric `3×3` `C5` are converted to eV with `1 GeV=1e9 eV`. These are **per absolute unknown factor** `|N_GammaCC|`. Only singular-value magnitudes, not PMNS fits, are calculated.

## Development checks

The local suite reported 55/55 tests passing, including seven new v6 tests. **Limitation:** The local validation container does not include the upstream `Numerical/experimental` package; the two small test-only local stand-ins replicate the inspected upstream formulas. These stand-ins are **not included in the ZIP**. The user's Windows repository will run the tests against the actual checked-in experimental modules. Synthetic tests verify the plumbing but not a source-complete physical prediction.
