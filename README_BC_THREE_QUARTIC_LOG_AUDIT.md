# Independent B(-1)/C(-1) three-quartic leading-log diagnostic

## Scope

Read-only diagnostic using existing fitted B/C models, the **signed** B/C EFT1
ablation, and its independent common-log Matchete audit. No Yukawa refits, RGE
changes, spinor normalization redefinitions, or modifications to the 14
nonmixed models. The EFT sequence remains
UV → integrate out **all fermions** → EFT1 with **both scalars** → integrate out
**all scalars** → SMEFT.

Three independently sourced channels:

| RGBeta/EFT1 coupling | Matchete symbol | EFT1 channel |
|---|---|---|
| `lambdaT3` | `lambdaT3` | `native:lambdaT3` |
| `lambdaHHdagS2S2` | `lambdaT3ExtraQ020020I1` | `general:lambdaHHdagS2S2` |
| `lambdaHHdagS1barS1bar` | `lambdaT3ExtraQ020200I1` | `general:lambdaHHdagS1barS1bar` |

In the one-scale evaluator, the last two couplings obey the existing
**conditional** assumption `lambda_Matchete=-lambda_RGBeta`. The signed EFT1
run has the same two signs inverted at its scalar tensor construction step.
This is not certified as a physical RGBeta/Matchete coupling dictionary.

For each channel `k` the script evaluates the *existing* compiled one-scale
Matchete coefficient with only quartic `k` enabled (in a deep copy of the
already evolved EFT1 endpoint). It shifts every explicit
`Log[MF^2/Mscalar^2]` by `Log[2]`, leaving masses, Yukawas, U and all RGEs
unchanged. The derivative is

    D_M(k) = [C5_M(k; shifted logs)-C5_M(k; original logs)]/log(2).

The independent EFT1 input is read from `1111.per_general_interaction`,

    D_E(k) = C5_EFT1_leading_log(k)/[2 log(1e5/1e3)].

The best-fit **complex** multiplier `z=<D_E,D_M>/<D_E,D_E>` is calculated
only as a diagnostic, never inserted into the physics pipeline.

The script fails closed if all three Matchete pieces do not sum to the full
one-scale result, if the EFT1 channels do not sum to the stored EFT1 log, if
the original full coefficient disagrees with the signed run, or if the
independent full-log result disagrees with `log_bc_signed_01`.

## Install

Extract the ZIP at repository root `C:\thesistemp\T3neutrinomassmodelpipeline`.
The ZIP adds only:

- `Numerical/running/GeneralT3BCQuarticChannelLogAudit.py`
- `scripts/RunT3BCQuarticChannelLogAudit.py`
- `tests/python/test_bc_quartic_channel_log_audit.py`
- `README_BC_THREE_QUARTIC_LOG_AUDIT.md`

## PowerShell

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline

python -X utf8 -m pytest -q tests/python/test_bc_quartic_channel_log_audit.py

python -X utf8 -m scripts.RunT3BCQuarticChannelLogAudit --check-only

python -X utf8 -m scripts.RunT3BCQuarticChannelLogAudit `
    --source-run thesis16_merged_v1 `
    --case-run general_seq_bc_signed_01 `
    --reference-audit log_bc_signed_01 `
    --run-id bc_channel_log_01
```

The default tests the saved heavy masses. To include an optional MF hierarchy
scan, append `--hierarchy-factors 1,2,4`; this reevaluates the expression at
three mass scales and will take more time.

Results:

```
output/studies/general_scalar_extended/general_bc_quartic_channel_log_audit/bc_channel_log_01/
  summary.csv
  summary.json
  B_m1/channel_audit.json
  C_m1/channel_audit.json
```

The `summary.csv` contains 6 channel rows at the default MF scale. The
important columns are `best_fit_factor_real`,
`matrix_residual_after_factor` and `relative_difference_factor1`.
The per-model JSON files include the full 3x3 channel matrices, comparison
factors, and independent three-channel reconstruction checks.

## Limits

Only the logarithmic coefficient is compared. The result neither computes an
independent hard fermion matching residual nor certifies GammaCC/Weyl
normalization, the mixed scalar basis at MF, or the quartic map. The full
sequential one-loop neutrino mass is not claimed.
