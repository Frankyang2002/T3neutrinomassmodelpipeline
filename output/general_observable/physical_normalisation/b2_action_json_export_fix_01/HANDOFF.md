# B2 internal-action native RawJSON transport correction

## Local user report

- `tests/test_physical_b2_rgbeta_internal_action*.py` plus prior B2 tests: 68 passed on Windows.
- Native RGBeta 1.2.0 loading succeeded.
- `scripts/RunPhysicalB2RGBetaInternalActionProbe.wl` printed `BLOCKED: native JSON export failed` and left a file that Python could not parse (`line 1 column 1`).
- Registry factors have **not** been independently measured by that failed run.

## Changed files (complete replacements)

1. `scripts/RunPhysicalB2RGBetaInternalActionProbe.wl`: preserve all five existing RGBeta representative builds and six registry probes. Add precise non-JSON leaf reporting (Wolfram `Missing` and symbolic expression candidates), explicitly invalidate native verification if a non-JSON value is present, serialize in memory using `ExportString[...,"RawJSON"]`, round-trip import, write a `.partial` file as UTF-8 text, verify, then rename. Retain `GroupInvBar` DownValues symbol fix. Refuse existing outputs and preexisting partials.
2. `Numerical/lo_certification/physical_normalisation/b2_rgbeta_internal_action.py`: complete replacement of previous Python validator, with explicit failures for zero-byte or malformed input JSON; all existing physical status checks remain.
3. `tests/test_physical_b2_rgbeta_internal_action_json_export.py`: 8 additional tests for new transport rules and failure handling.

## Local evidence and limits

- Local Python: 76 tests passed across all five B2 test files (68 existing + 8 new).
- External Wolfram Language kernel: `RawJSON` rejects `Missing[...]` / unresolved symbolic expressions; JSON preflight/normalizer and JSON text round-trip smoke-tested. This is **not** the user's installed RGBeta version or local source tree.
- The exact offending field in the user's first failed export is unknown because it never produced valid JSON.
- The replacement emits diagnostic JSON and reports offending fields if any non-JSON values remain. A diagnostic report has status `NativeInternalTensorDataNotJSONSafe`; the Python validator rejects this status and cannot convert it to a successful witness.
- No production files changed, no matching/running/EFT, no coupling map certified. Physical source certificates **0/16**.

## Windows run commands

```powershell
Set-Location "C:\thesistemp\T3neutrinomassmodelpipeline"
$zip = Join-Path $env:USERPROFILE "Downloads\T3_B2_Internal_Action_JSON_Export_Fix.zip"
Test-Path $zip
Expand-Archive -LiteralPath $zip -DestinationPath "." -Force
.\.venv\Scripts\python.exe -m pytest tests\test_physical_b2_rgbeta_internal_action_json_export.py tests\test_physical_b2_rgbeta_internal_action_probe_source.py tests\test_physical_b2_rgbeta_internal_action.py tests\test_physical_b2_rgbeta_ward_witness.py tests\test_physical_b2_gauge_coupling.py -q
$env:T3_EXPERIMENTAL_GENERAL_FERMIONS = "0"
$env:T3_EXPERIMENTAL_GENERAL_QUADRATICS = "0"
$env:T3_EXPERIMENTAL_GENERAL_QUARTICS = "0"
$env:T3_EXPERIMENTAL_FLAVOR_MASSES = "0"
wolframscript -file scripts\RunPhysicalB2RGBetaInternalActionProbe.wl "output\general_observable\physical_normalisation\b2_rgbeta_action_native_03\internal_action_registry.json"
.\.venv\Scripts\python.exe scripts\RunPhysicalB2RGBetaInternalAction.py --native "output\general_observable\physical_normalisation\b2_rgbeta_action_native_03\internal_action_registry.json" --output "output\general_observable\physical_normalisation\b2_action_local_03"
```

If the Wolfram script reports `B2 JSON issue: ...`, inspect the emitted `json_preflight_issues` and share the diagnostic JSON; do not claim the native factor checks succeeded. Reusing `_02` is discouraged because it may contain the failed/empty JSON.

## Physical caveat

RGBeta `Sqrt[2]`, `2` and `6` are internal tensor normalization multipliers; they do not establish physical Matchete–RGBeta action ratios or resolve SU(2) barred-index conventions or the B1 spinor phase. The B/C literal-index Ward mismatch remains conditional.
