# T3 end-to-end software smoke test

This checks operational integration of the **existing production** `Numerical.fitting.ScanCLI` code. It does **not** modify any production source or certify physical accuracy.

## Installation

Extract the ZIP into the repository root `C:\thesistemp\T3neutrinomassmodelpipeline` with directory structure preserved.

## Commands (PowerShell)

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
$env:PYTHONUTF8 = '1'
python scripts/maintenance/SmokeTestFullPipeline.py
```

If no automatically discovered config validates, run with a known existing production scan JSON:

```powershell
python scripts/maintenance/SmokeTestFullPipeline.py --config path\to\existing_scan.json
```

Only one point is evaluated. An existing grid uses its first value in each dimension; a random scan uses one reproducible seeded draw. The repository's existing physics configuration, matching artifacts, and RGE payloads are reused unchanged.

## Output

`output/audit/end_to_end_smoke/smoke_result.json` — `PASS_END_TO_END_OPERATIONAL` only if one numerical point succeeds and a finite observable dictionary and output JSON are produced.

`output/audit/end_to_end_smoke/single_point_scan.json` — original ScanCLI output for the single point.

A failure can mean absent or invalid matching exports, a configuration error, or numerical integration failure. The JSON report records why; it must not be treated as a physical accuracy assessment.
