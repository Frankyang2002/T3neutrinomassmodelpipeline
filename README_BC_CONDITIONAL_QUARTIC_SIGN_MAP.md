# B(-1)/C(-1) conditional direct-quartic sign consistency check

## Scope

This is a **conditional numerical convention**, not a certified physical
Matchete/RGBeta coupling dictionary. It aligns the B/C EFT1 HHSS contraction
with the existing one-scale B/C evaluator's declared sign mapping:

    lambda_Matchete = -lambda_RGBeta

for exactly two direct HHSS couplings:

- `lambdaHHdagS2S2` / `general:lambdaHHdagS2S2`
- `lambdaHHdagS1barS1bar` / `general:lambdaHHdagS1barS1bar`

Only the two already-calculated quartic **tensor channels** are inverted. The
runner creates a separate immutable `TensorBundle` with those channels negated
and its total recomputed. Input scalar couplings, UV/EFT1 RGEs, heavy-fermion
source-exact Wilson tensors, mass eigenbasis, one-scale Matchete coefficient,
and the 14 other models are unchanged. Each numerical case records the
conditional provenance in its `quartic_provenance` field.

The earlier **read-only** recombination of the saved B/C channel outputs found
that these sign flips reduce the best-normalization matrix residuals from
77.3% to ~5.57% for B(-1) and from 95.8% to ~7.75% for C(-1). The remaining
best-fit factors (~1.775, ~1.812) are NOT inserted into physics. These are
expected diagnostics, not guarantees for the newly recomputed runs.

## Files (all paths relative to repository root)

- `Numerical/running/GeneralT3ConditionalBCQuarticMap.py` — **NEW**, standalone
  sign-map validator and pure tensor transformation.
- `scripts/RunT3GeneralSequentialFitAblation.py` — **FULL REPLACEMENT**. Only
  `_mixed` constructs a signed tensor bundle; `_nonmixed` is unchanged.
- `tests/python/test_general_bc_conditional_direct_quartic_map.py` — **NEW**
  tests that ensure sign exclusivity, no mutation, no double application,
  preservation of mass eigenbasis, and B/C-only controller integration.
- `README_BC_CONDITIONAL_QUARTIC_SIGN_MAP.md` — this file.

This package **requires** the previous B/C source-exact backend installation;
it does not overwrite or include that backend. No GitHub files are changed.

## Install and test (PowerShell; repository root)

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline

Expand-Archive `
  -LiteralPath "$HOME\Downloads\T3_BC_Conditional_Quartic_Sign_Map.zip" `
  -DestinationPath . -Force

python -X utf8 -m pytest -q `
  tests/python/test_general_bc_conditional_direct_quartic_map.py `
  tests/python/test_general_mixed_exact_fermion_boundary.py

python -X utf8 -m scripts.RunT3MixedExactSourceSmoke
```

## Recompute B(-1)/C(-1) only, preserving previous runs

```powershell
python -X utf8 -m scripts.RunT3GeneralSequentialFitAblation `
  --source-run thesis16_merged_v1 `
  --run-id general_seq_bc_signed_01 `
  --models B_-1,C_-1 `
  --spinor-factor 1

python -X utf8 -m scripts.RunT3GeneralMatcheteLogCoefficientAudit `
  --source-run thesis16_merged_v1 `
  --case-run general_seq_bc_signed_01 `
  --run-id log_bc_signed_01 `
  --models B_-1,C_-1 `
  --hierarchy-factors 1,2,4

Import-Csv `
  output\studies\general_scalar_extended\general_log_coefficient_audit\log_bc_signed_01\summary.csv |
  Select-Object model,status,best_fit_factor_real_at_MF,relative_matrix_residual_best_factor_at_MF |
  Format-Table -AutoSize
```

Results are written to:

- `output/studies/general_scalar_extended/general_sequential_rge_ablation/general_seq_bc_signed_01/`
- `output/studies/general_scalar_extended/general_log_coefficient_audit/log_bc_signed_01/`

The original `general_seq_05`, `general_seq_bc_exact_01`, `log_coefficient_01`,
and `log_bc_exact_01` outputs are left unchanged.

## Interpretation

Check `summary.csv` for B/C normalization and post-best-fit matrix residuals.
Any residual may reflect an additional scalar-basis/CG/mapping discrepancy.
The one-scale C5 anchor and GammaCC/Weyl normalization are not independently
certified, and hard matching is not independently calculated. Do not treat the
new Wilson-off ablations as physical predictions merely because the residual
has decreased.
