# Stream 2 endpoint audit — after 16-model run

This package ADDS a read-only independent endpoint checker for existing Stream 2 numerical outputs. It does not replace the ODE/ADM or modify any `wilson.json` or `.npz`. The check establishes arithmetic consistency, source-file hashes, retained exact Stream 1 identifiers, full `(45,45,ns,ns)` tensor, MF/MS projections, beta-channel sum, scalar/lepton coordinate slicing, and trajectory endpoints. It does not independently certify Matchete/RGBeta scalar/GammaCC mapping, the chosen one-loop ADM or the physical matching.

Files:

* `scripts/AuditT3RoughStream2.py` — audit executable
* `tests/test_t3_rough_stream2_audit.py` — four contract tamper/success tests
* `STREAM2_AUDIT_README.md` — this guide

PowerShell:

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
$zip = Join-Path $env:USERPROFILE 'Downloads\t3_stream2_endpoint_audit.zip'
Expand-Archive -LiteralPath $zip -DestinationPath . -Force
& .\.venv\Scripts\python.exe -m pytest -q tests\test_t3_rough_stream2_audit.py
& .\.venv\Scripts\python.exe -m scripts.AuditT3RoughStream2 --require-existing-engine
```

Output: `output/parallel_t3/rough_wilson/rough_fermion_first_v1/stream2_audit.json`.

The user's console showed 10 substantial positive `||C5(MS)||` values, four exact zeros and two ~1e-28. As a numerical decision, values ~1e-28 are below the previous solver's default reported significant-generation floor (10*atol=1e-20), not established physical generation. Successful all-16 runs show that the *new canonical Stream1 input directory* was present; the legacy saved Workstream B probe checks an entirely different directory and cannot be interpreted as saying the solver used B's native parameters directly.

If any model fails the audit, inspect `stream2_audit.json`. In particular, `input_files_sha256_rechecked` catches canonical inputs that changed after the runs. If the full beta backend is the portable formula mirror rather than `existing_repository_matrix_free_engine`, check the installed RGE imports. Only the latter passes `--require-existing-engine`.
