# B(-1)/C(-1): fixed channel normalization for the independent leading-log audit

## Scope and result

This package implements the **fixed, no-fit** B/C-only candidate mapping

\[
 D_{\rm EFT1}^{\rm candidate}
 = \sqrt3\,D_{\rm native:\lambda T3}^{\rm signed}
 +2\,D_{\rm general:\lambda HHdagS2S2}^{\rm signed}
 +2\,D_{\rm general:\lambda HHdagS1barS1bar}^{\rm signed}.
\]

The two `general:` tensors have **already been negated** by the existing
`GeneralT3ConditionalBCQuarticMap`. The new positive factors `2` do **not**
reapply a sign change. The 14 nonmixed models, all RGE beta functions, the
fitted UV Yukawas, and the Matchete one-scale evaluator are unchanged.

Based on the previous three-channel B/C audit, the *independently extracted*
Matchete/EFT channel factors at the saved benchmark were:

| Model | topology | direct S2S2 | direct S1barS1bar |
|---|---:|---:|---:|
| B_-1 | 1.7324939051 | 2.0003861953 | 2.0005416290 |
| C_-1 | 1.7322974782 | 2.0003092034 | 2.0003553745 |

The fixed-factor recombination of the **earlier saved** matrices (independent
of this code change) gives:

| Model | fixed-factor relative matrix difference at MF | at 4 MF |
|---|---:|---:|
| B_-1 | 0.0259624% | 0.00162275% |
| C_-1 | 0.0170629% | 0.00106645% |

These are numerical checks, **not an independent analytical proof** of the
Matchete/RGBeta/CG/spinor normalization. The agreement improves as MF increases,
consistent with finite-mass corrections.

## Files (project-relative paths)

- `Numerical/running/GeneralT3ConditionalBCQuarticNormalization.py` — **new** fixed channel map
- `scripts/RunT3GeneralSequentialFitAblation.py` — **complete replacement** controller (other 14 models untouched)
- `scripts/CheckT3BCFixedNormalization.py` — **new** read-only check against previously saved signed outputs
- `tests/python/test_bc_fixed_quartic_normalization.py` — **new** regression tests
- `README_BC_FIXED_CHANNEL_NORMALIZATION.md` — **new** instructions

Prerequisites: the B/C source-exact backend and preceding signed direct-quartic
map must already be installed. This package does not replace either prerequisite.

## Install from PowerShell at repository root

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "$HOME\Downloads\T3_BC_Fixed_Channel_Normalization.zip" -DestinationPath . -Force

python -X utf8 -m pytest -q `
  tests/python/test_bc_fixed_quartic_normalization.py `
  tests/python/test_general_bc_conditional_direct_quartic_map.py `
  tests/python/test_general_mixed_exact_fermion_boundary.py
```

## Verify the fixed factors without any new running

Reads the existing local `general_seq_bc_signed_01` and `log_bc_signed_01`
reports. This is a direct fixed-factor comparison (no fitted scale in the
physics prediction):

```powershell
python -X utf8 -m scripts.CheckT3BCFixedNormalization `
  --case-run general_seq_bc_signed_01 `
  --reference-audit log_bc_signed_01
```

The default fail-closed limit is **0.1%** fixed-factor relative matrix difference
at each of MF, 2MF, and 4MF. If additional quartic channels contribute, the
check rejects the three-channel assumption instead of silently discarding them.

## Recompute the B/C-only ablation with the fixed map

```powershell
python -X utf8 -m scripts.RunT3GeneralSequentialFitAblation `
  --source-run thesis16_merged_v1 `
  --run-id general_seq_bc_normalized_01 `
  --models B_-1,C_-1 `
  --spinor-factor 1
```

This changes the independent LLSS-induced EFT1 contribution and its ablation.
The W=1 one-scale-anchored C5 remains the same **by construction**; it is not a
new independent full hard-matching prediction.

## Repeat the independent full Matchete common-log check

```powershell
python -X utf8 -m scripts.RunT3GeneralMatcheteLogCoefficientAudit `
  --source-run thesis16_merged_v1 `
  --case-run general_seq_bc_normalized_01 `
  --run-id log_bc_normalized_01 `
  --models B_-1,C_-1 `
  --hierarchy-factors 1,2,4

$p = 'output\studies\general_scalar_extended\general_log_coefficient_audit\log_bc_normalized_01'
Import-Csv "$p\summary.csv" | ForEach-Object {
  '{0,-5} ratio={1:F8} fit={2:F8} fixed_diff={3:P5} residual={4:P5}' -f `
    $_.model, `
    [double]$_.ratio_Matchete_over_EFT_at_MF, `
    [double]$_.best_fit_factor_real_at_MF, `
    [double]$_.relative_matrix_diff_factor1_at_MF, `
    [double]$_.relative_matrix_residual_best_factor_at_MF
}
```

Expected at `MF`: B_-1 raw fixed-factor matrix difference ~0.026%; C_-1
~0.017%; best-fit diagnostic factors close to 1.000258 and 1.000162.
These are **expected**, not asserted as actual new-run results until rerun.

Output locations:

- `output/studies/general_scalar_extended/general_sequential_rge_ablation/general_seq_bc_normalized_01/`
- `output/studies/general_scalar_extended/general_log_coefficient_audit/log_bc_normalized_01/`

## Physics limitations

- Fixed factors `sqrt(3), 2, 2` are supported independently by the three
  source-separated common-log matrices, but are **not yet analytically certified**
  as a physical Matchete-to-RGBeta coupling dictionary.
- The complex mixed-scalar basis is provisionally identified at `mu_F` and `mu_S`.
- `GammaCC` to physical Weyl normalization is not independently certified.
- Hard F-threshold matching and missing EFT1 operator mixing are not supplied.
- This is a successful numerical *LLSS scalar-pair leading-log consistency*
  test when the rerun passes, **not** a complete physical one-loop SMEFT match.
- Order is always UV -> integrate out ALL fermions -> EFT1 with BOTH scalars
  -> integrate out ALL scalars -> SMEFT.
