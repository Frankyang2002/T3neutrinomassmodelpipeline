# Independent full-general Matchete vs EFT1 leading-log coefficient audit

**Applies to the current 16-representation T3 general-Lagrangian pipeline.** Read-only source inspection: GitHub `Frankyang2002/T3neutrinomassmodelpipeline`, 4 November 2026. This package adds new files and **does not modify any existing matching, RGE, numerical fitting, or ablation implementation**. GitHub is not modified.

## Question tested

Does the general fermion-first EFT1 scalar-pair mixing actually reproduce the coefficient of the mass-ratio logarithm in the full one-scale Matchete one-loop neutrino mass?

For matched fermion scale `muF=10^5 GeV` and scalar scale `muS=10^3 GeV`, define

`L = ln(muF/muS) = ln(100)`.

The current EFT1 calculation gives `C5_LLSS_log = log(muS/muF) * beta_C5/(16*pi^2)`.

The independent Matchete coefficient contains logs `Log[MF^2/MSi^2]`. Shift **only their arguments** by 2 while keeping all masses and couplings unchanged:

`Log[MF^2/MSi^2] -> Log[2*(MF^2/MSi^2)]`.

The exact one-loop coefficient of a common additive logarithmic shift is

`D_Matchete = (C5_shifted - C5_original)/ln(2)`.

The corresponding leading-order EFT prediction for the same normalization is

`D_EFT1 = C5_LLSS_log/(2*ln(muF/muS))`.

This compares **two separately calculated sources of the coefficient**: the saved full one-scale Matchete expression and the independently constructed fermion-tree LLSS × general HHSS tensor. It **does not** use `C5_hard = C5_one_scale - C5_log - C5_bubble` as validation; that identity is bookkeeping, not an independent test.

For an exactly matched normalization and sufficiently hierarchical masses, `D_Matchete ≈ D_EFT1`. The ratio may deviate at finite masses by terms of order `MS^2/MF^2` or because the conditional scalar/Yukawa basis mappings are not equivalent.

The report includes an **optimal complex scale** `N_fit = <D_EFT1,D_Matchete>/<D_EFT1,D_EFT1>`, and two separate residuals:

- `relative_matrix_difference_at_assumed_factor_1`: direct physical-convention check with the current assumed factor 1.
- `relative_matrix_residual_after_best_fit_factor`: mismatch that cannot be removed by a single complex factor (basis/channel mismatch).

**The fitted factor is never used to change the computed one-scale or EFT coefficients.** It is a diagnostic to assess sign, normalization and tensor direction.

## Implementation and input sources

| Batch | Source | Operation |
|---|---|---|
| Nonmixed 14 models | Saved `thesis16_merged_v1/<model>/result.json`, current `output/general_matchete_matching/<model>/c5_coefficient.txt` | Reconstruct 3 heavy mass eigenstates and all Yukawa channels at the actual fitted `muF`, use existing `EvaluateGeneralC5Channels.wl` on unchanged and log-shifted expression at fitted `muS` |
| Mixed B(-1), C(-1) | Same saved fit, current mixed Matchete ordered C5 source | Regenerate **only** U=E=1 via existing mixed UV/EFT1 running implementation and evaluate original and shifted sources with existing compiled ordered C5 evaluator |
| EFT1 coefficient | Saved `general_seq_05/<model>/cases.json` case `1111` | Read `c5_induced_wilson_log`, independently calculated from dynamic F-tree Wilson tensor and all general HHSS interactions |
| Hierarchy scan | MF × `{1,2,4}`, Yukawa and scalar sector held fixed | Test whether any finite-MF discrepancy shrinks toward the large hierarchy. For comparison scale EFT tree coefficient as `1/MF` |

All 16 current one-scale expressions were inspected: their logarithm arguments are exactly `MF^2/(scalar mass)^2`, including split real scalar masses and mixed X/Y eigenmasses. The extractor fails on any other logarithm, so an unrelated log is never inadvertently shifted.

The fitted original one-scale C5 is recalculated and checked against the stored all-running case `1111` before interpreting a discrepancy.

## Package paths

Files go directly to these relative project paths:

- `Numerical/running/GeneralT3MatcheteLogCoefficientAudit.py` — coefficient extractor, comparison matrices, fitted-endpoint adapters for all 16 models.
- `scripts/RunT3GeneralMatcheteLogCoefficientAudit.py` — CLI, preflight, all-16 batch, per-model JSON/CSV outputs.
- `tests/python/test_general_matchete_log_coefficient_audit.py` — unit tests and mock Wolfram driver contract.
- `README_GENERAL_MATCHE_LOG_AUDIT.md` — this document.

## Install and run (PowerShell)

From `C:\thesistemp\T3neutrinomassmodelpipeline`:

```powershell
Expand-Archive -LiteralPath "$HOME\Downloads\T3_General_Matchete_Log_Coefficient_Audit.zip" -DestinationPath . -Force
python -X utf8 -m pytest -q tests/python/test_general_matchete_log_coefficient_audit.py
```

Preflight, without creating output or running RGEs:

```powershell
python -X utf8 -m scripts.RunT3GeneralMatcheteLogCoefficientAudit --check-only
```

All models, single benchmark plus `MF` hierarchy scaling by 2 and 4:

```powershell
python -X utf8 -m scripts.RunT3GeneralMatcheteLogCoefficientAudit `
  --source-run thesis16_merged_v1 `
  --case-run general_seq_05 `
  --run-id log_coefficient_01 `
  --hierarchy-factors 1,2,4
```

To do the baseline comparison more economically:

```powershell
python -X utf8 -m scripts.RunT3GeneralMatcheteLogCoefficientAudit `
  --source-run thesis16_merged_v1 `
  --case-run general_seq_05 `
  --run-id log_baseline_01 `
  --hierarchy-factors 1
```

The 14 nonmixed models use the Wolfram executable `wolframscript` already used by the general matching pipeline; add `--wolframscript C:\Path\To\wolframscript.exe` if it is not on PATH. The two mixed models use mpmath and the existing compiled Matchete evaluator and **do not require Wolfram**:

```powershell
python -X utf8 -m scripts.RunT3GeneralMatcheteLogCoefficientAudit `
  --run-id log_mixed_01 --models B_-1,C_-1 --hierarchy-factors 1,2,4
```

**Outputs** (new files only):

```text
output/studies/general_scalar_extended/general_log_coefficient_audit/log_coefficient_01/
    summary.json
    summary.csv
    A_m4/log_audit.json
    ...
    C_m1/log_audit.json
```

In `summary.csv`, inspect `best_fit_factor_real_at_MF`, `relative_matrix_diff_factor1_at_MF`, `relative_matrix_residual_best_factor_at_MF`, and their hierarchy endpoints. If all matrix residuals shrink as `MF` rises, that supports finite-mass corrections to the leading EFT1 coefficient. If the overall factor persists away from 1, investigate the GammaCC/Weyl normalization or scalar potential sign. If the best-fit residual remains large, investigate extra-operator channels, different scalar bases or HHSS quartic maps rather than inserting a single normalization factor.

## Verification and qualification

The new extractor and control logic have passing synthetic local tests. There is **no local checkout of the full current repository or Mathematica executable in this environment**, so the numerical all-16 audit is not claimed executed here. Its result is generated on the user's machine from the committed Matchete expressions and saved fitted run data.

An independent analytic spot check of the **A(-4)** saved benchmark, using the explicit two-mass Matchete expression and fitted `muF/muS` matrices, gave a Matchete/EFT coefficient ratio `1.000197` and relative matrix difference `0.0198%` at finite mass ratio. Keeping the same fitted Yukawas and scalar masses while scaling MF by 2 and 4 gave ratios `1.0000493` and `1.0000123`, respectively, consistent with finite-mass corrections decreasing approximately as `1/MF^2` in this model. This check uses the one-scale source, not the ablation bookkeeping residual. The new tool is intended to reproduce and extend that check for every model.

The test is for the **scalar-pair leading-log coefficient**. It does not independently match the hard fermion loops, close every dimension-5 operator class, or prove the extra-quartic RGBeta ↔ Matchete basis map. Its output must not be called a full certified sequential EFT matching result.
