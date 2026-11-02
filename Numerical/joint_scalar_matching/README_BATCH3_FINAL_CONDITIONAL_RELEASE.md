# Batch 3 — final T3 exploratory numerical release

**STOP after this batch.** This is the third and final batch of the agreed rough-numerical milestone. It does not begin a new certification programme.

## Exact scope

- Preserve Batch 2's complete, validated 16-model conditional `C5` table **verbatim**.
- Revalidate 16 matrices, the specified `mnu = (246.22 GeV)^2 C5/2` conversion, mass singular values, Ward residuals, and source conditionality.
- For the representative models `A_-4` (ordinary scalar basis) and `B_-1` (mixed basis), **rerun** the repository's existing `Numerical/direct_one_loop/core.py` on each model's existing `input/lo_direct_one_loop/*.json` frozen UV input. Verify the frozen benchmark SHA-256; independently check existing archived direct results if present.
- Present separate direct-UV and EFT matrices and mass singular values. Compute the exact direct-UV finite-B0 heavy-mass expansion tail *at the UV point*, without modifying the sequential EFT result.
- Explicitly record why a physical or point-matched UV-versus-EFT **relative difference is not presently defined**.

No generated file changes UV, EFT1, or scalar matching code. GitHub remains read-only.

## The sequence is unchanged

```
UV -> integrate out ALL BSM fermions at 100000 GeV
   -> EFT1 (BOTH BSM scalar multiplets present)
   -> integrate out BOTH BSM scalar multiplets together at 1000 GeV
   -> SMEFT
```

All SM fermions remain active. There are no individually separated BSM fermion or scalar thresholds.

## Why the two calculation outputs MUST NOT be subtracted as physics

The inspected source files show:

| Property | Direct UV (`core.py` input) | Sequential EFT (Batch 2) |
|---|---|---|
| Source | Existing frozen `Numerical/benchmarks/fitted_reference/T3_*.json`, with provisional UV component vertices | Illustrative running endpoint `output/general_observable/all16_v1` and source-resolved LLSS/HHSS tensors |
| Input reference | `base_state.mu_gev = 10000000` GeV | `muF = 100000` GeV; `muS = 1000` GeV |
| Mass treatment | Broken-phase scalar mixing and exact `B0(0;MF^2,mS^2)` | Unbroken illustrative EFT1 masses, scalar quartic log + simultaneous finite scalar bubble |
| Normalization | Direct source four-spinor-to-Weyl mapping **uncertified** | Unknown signed `N_GammaCC` |
| Couplings | Restricted frozen UV inputs (not all general interactions) | Provisional Matchete/RGBeta quartic aliases |
| Terms | Direct `p=0` one-loop mass from supplied vertices | LLSS-induced leading-log and finite scalar matching only; omits hard heavy-fermion one-loop and complete EFT1 Wilson mixing |

Thus the matrices use *different renormalized parameter points, scalar masses, and unverified convention maps*. A raw Frobenius matrix difference or eigenvalue ratio is mathematically calculable, but would be physically misleading. The release reports **`null`** for `relative_uv_vs_eft_matrix_difference`, `matrix_difference_gev`, and the cross-method singular-value ratio. It reports **both separate matrices and spectra** instead. No fit or unknown normalization is inserted.

`Numerical/direct_one_loop/comparison.py` has similar strict same-point safeguards; this batch has no valid same-point benchmark certificate to supply to it.

### Exact finite heavy-mass diagnostic (direct UV only)

For each direct-UV internal eigenstate, define `r=mS^2/MF^2` and use:

\[
 \bar B_0(0;M_F^2,m_S^2;\mu)
 =1-\log\frac{M_F^2}{\mu^2}+\frac{r\log r}{1-r}.
\]

At `r << 1`, the next finite term beyond the leading small-`r` B0 expansion is

\[
 \Delta\bar B_0=\frac{r^2\log r}{1-r}.
\]

The runner reweights the direct UV core's *individual* eigenstate loop terms by `Delta B0/B0` and sums them into one **complex symmetric** matrix. This is a same-point diagnostic of the finite heavy-mass tail in the **direct calculation only**. It is **not** a correction to Batch 2 EFT matching, nor a certified EFT truncation error. It fails closed if `r >= 0.1` or `B0` is near zero.

## Previously archived source values (context, not the new local run)

The inspected GitHub `lo_direct_one_loop/provisional` reference outputs list these direct singular values, **ascending**, in eV:

- `A_-4`: `[0.00047448293, 0.00758272810, 0.04360064404]`.
- `B_-1`: `[0.00452249250, 0.00980338343, 0.05044130936]`.

The corrected v6 reference summary lists these conditional sequential EFT singular values, **descending**, per `|N_GammaCC|`:

- `A_-4`: `[0.02590212217, 0.02158510158, 0.01942653400]`.
- `B_-1`: `[0.02369106692, 0.01974255556, 0.01776824742]`.

These values are **not paired same-point physics predictions**. The Batch 3 script reads and recomputes from your local files and writes its own values; it does not hard-code these reference numbers.

A later source audit at `lo_certification_direct_uv/charge_support_fix_01/all_16_certification_status.json` reports 16/16 provisional direct evaluations, **zero source-certified physical UV models**. Neither model is a certified physical source prediction.

## Full replacement files in this ZIP

- `scripts/BuildT3Batch3FinalRelease.py` — new finalizer, reads inputs, computes direct UV from existing engine, builds final release.
- `tests/test_t3_joint_scalar_batch3_final_release.py` — four focused tests of the new validation, B0 tail formula, honest noncomparability, and fail-closed release.
- `Numerical/joint_scalar_matching/README_BATCH3_FINAL_CONDITIONAL_RELEASE.md` — this file.

No existing source file is replaced. Everything is intended to be installed under the existing project root, preserving these relative paths.

## PowerShell from the project root

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive `
  -LiteralPath "$HOME\Downloads\T3_Batch3_Final_Conditional_All16_UV.zip" `
  -DestinationPath . -Force

python -m unittest discover -s tests `
  -p "test_t3_joint_scalar_batch3_final_release.py" -v

python -m scripts.BuildT3Batch3FinalRelease `
  --batch2-dir "output\general_observable\lo_joint_scalar_matching\conditional_all16_batch2_v1" `
  --uv-input-dir "input\lo_direct_one_loop" `
  --previous-direct-dir "output\general_observable\lo_direct_one_loop\provisional" `
  --output-dir "output\general_observable\lo_joint_scalar_matching\final_conditional_all16_batch3_v1"
```

No need to rerun Matchete, RGBeta, EFT1 or Batch 2. The script refuses to overwrite a non-empty destination.

## Final outputs

Under `output/general_observable/lo_joint_scalar_matching/final_conditional_all16_batch3_v1/`:

- `summary.json`: final all-16 records plus the two detailed direct-UV/EFT comparison records, statuses and assumptions.
- `all16_conditional_mass_summary.csv`: **byte-identical copy** of validated Batch 2 CSV.
- `direct_uv_vs_eft_A_m4.json`, `direct_uv_vs_eft_B_m1.json`: each model's two matrices, mass singular values, UV eigenmasses, hierarchy, exact finite-B0 heavy-mass tail, provenance, and blockers.
- `direct_uv_comparison.csv`: compact two-row cross-method status and separated spectra (no invented normalization).
- `validation.json`: prior user-reported tests, current file-level/matrix validations, both direct-UV numerical statuses, and non-certification flags.

Expected terminal line if the two existing direct inputs are intact:

```
Batch 3: all16=16/16; provisional direct UV=2/2; same-point physical comparisons=0
```

Any UV failure is recorded **per model** in JSON and CSV and causes exit code `2` (without fabricating a direct matrix). Even a `2/2` provisional run still has **zero** certified UV-vs-EFT same-point comparisons. This satisfies the agreed numerical release scope; full matching normalization and same-point physical matching are explicitly deferred.

**EXPLORATORY / CONDITIONAL — NOT A COMPLETE PHYSICAL ONE-LOOP NEUTRINO-MASS PREDICTION**
