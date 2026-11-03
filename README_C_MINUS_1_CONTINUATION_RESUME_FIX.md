# General T3 sequential leading-log: C(-1) continuation and resume fix

Verified against GitHub `Frankyang2002/T3neutrinomassmodelpipeline` current `main`, `scripts/RunT3GeneralSequentialFitAblation.py` (Git blob `0e838204d93cc19b98c55379e08c41f0a9179b43`), and the saved `continue_general16_v1/C_m1` continuation artifacts.

## Replacement files (complete, not patches)

- `scripts/RunT3GeneralSequentialFitAblation.py`
- `tests/python/test_general_sequential_leading_log.py`

No existing physics machinery is altered. GitHub remains unchanged.

## Cause

Original B/C fitted `result.json` contains `extra_quartic_fraction` and `extra_mixed_scalar_bilinear_fraction` at top level. `continue_general16_v1/C_m1/result.json` instead contains `matching_and_rge_provenance.quartic_fraction` and `.mu2_fraction`, equal to 0.1 and 0.01, respectively, with `modified_scalar_boundary=false`. Previously, `_mixed` read only the former fields and raised `KeyError` for C(-1).

`mixed_fit_fractions` now supports both saved result layouts. It checks consistency if both appear and never inserts unsourced defaults. No refitting occurs.

## Resume the existing 15/16 study

From the repository root in PowerShell:

```powershell
Expand-Archive -LiteralPath "$HOME\Downloads\T3_General_Sequential_Cm1_Continuation_Resume_Fix.zip" -DestinationPath . -Force
python -X utf8 -m pytest -q tests/python/test_general_sequential_leading_log.py
python -X utf8 -m scripts.RunT3GeneralSequentialFitAblation `
    --source-run thesis16_merged_v1 `
    --run-id general_seq_05 `
    --spinor-factor 1 `
    --resume
```

`--resume` reuses existing 16-case `cases.json` files for completed models and retries only missing/previously blocked models, so C(-1) is the only model reevaluated in the reported `general_seq_05` run. It deletes a stale `error.json` only after C(-1) succeeds and then rewrites the all-model `summary.json` and `all16_cases.csv` based on the complete set. It checks the saved requested model list and spinor factor, and from new runs also checks the source campaign.

Expected if C(-1) succeeds: `Computed: 16` in the summary. Full numerical C(-1) trajectory requires local symbolic/numerical artifacts and has not been run in the standalone packaging environment.

## Verification

Local tests: 15 passed, including legacy/continuation metadata, missing and inconsistent values, CLI argument parsing, and a simulated 15-reuse + C(-1)-retry flow. Full 16-model integration has not been executed here.
