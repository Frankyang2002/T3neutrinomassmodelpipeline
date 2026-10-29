# General T3: WolframScript self-test argument fix

## Observed failure

The Batch 7 runner invoked `wolframscript -file scripts/CompareGeneralT3C5.wl --selftest`.
On Windows, `wolframscript` consumed the dash-leading argument rather than forwarding it to the Mathematica script. The script received zero positional arguments, printed its usage and exited with code 64. Consequently **no formal response or RGBeta/tensor check was run in that invocation**. The earlier 16-model matching and Batch 5 zero-extra proof are unaffected.

## Correction

Call `wolframscript -file scripts/CompareGeneralT3C5.wl SELFTEST` with a plain positional token. The Mathematica script now recognises the token and prints the received argument list if usage is wrong. A Python regression test guards the actual PowerShell invocation.

## Run

```powershell
Expand-Archive -LiteralPath .\T3_General_Scalar_Batch8.zip -DestinationPath . -Force
python -m pytest -q tests\python\test_general_scalar_trial.py tests\python\test_scalar_quadratic_prototype.py tests\python\test_general_scalar_c5_proof.py tests\python\test_general_rgbeta_coverage.py tests\python\test_general_t3_scalar_quartic_response_audit.py
wolframscript -file .\scripts\CompareGeneralT3C5.wl SELFTEST
powershell -ExecutionPolicy Bypass -File .\scripts\RunGeneralT3ScalarQuarticResponseAudit.ps1
```

The runner reuses your existing 16-model matching outputs and refreshes the formal response, then runs the RGBeta-source and Matchete-tensor checks. It does **not** enable a numerical general-T3 theory. If any step reports `ReviewRequired`, inspect the relevant JSON rather than assuming it passed.

## Changed files

- `scripts/RunGeneralT3ScalarQuarticResponseAudit.ps1` — full replacement; pass positional `SELFTEST`.
- `scripts/CompareGeneralT3C5.wl` — full replacement; recognise positional token and print usage diagnostics.
- `tests/python/test_general_t3_scalar_quartic_response_audit.py` — full replacement; check CLI contract.
- `docs/GENERAL_T3_WOLFRAMSCRIPT_SELFTEST_FIX.md` — new.

This archive is consolidated: it also contains all unchanged files from Batch 7.
