# B2 batch 4: Native RGBeta Git provenance handoff fix

## Observed native blocker

`wolframscript` had successfully loaded RGBeta but its four Git subprocess
captures became `Missing["GitCommandFailed", ...]`: one Git HEAD and three
tracked clean-filter blob hashes. The RawJSON diagnostic correctly refused
native validation, and the Python validator correctly rejected a missing HEAD.
This did **not** establish an RGBeta coupling failure.

## Complete replacement and added files

- **Replace:** `scripts/RunPhysicalB2RGBetaInternalActionProbe.wl`
- **Add:** `scripts/RunPhysicalB2RGBetaInternalActionNative.py`
- **Add:** `tests/test_physical_b2_internal_git_handoff.py`
- **Add:** `output/general_observable/physical_normalisation/b2_action_git_handoff_fix_01/HANDOFF.md`
- **Add:** `output/general_observable/physical_normalisation/b2_action_git_handoff_fix_01/evidence.json`

No upstream source, RGBeta installation or previous audit output changes.

## Change

The new Python launcher uses `source_provenance(root)` from the existing
B2 Ward witness to require a matching Git root, six unmodified source Git
blobs across HEAD/index/clean-filter checkout, and unchanged RGBeta
operator source tokens. It passes four validated literal 40-hex strings
into the isolated Wolfram child process via `T3_B2_PROVENANCE_*`.
Wolfram no longer invokes Git, and blocks if these values are missing.

After the subprocess returns success the launcher immediately runs the
existing `native_snapshot_guard` and `validate_native_data`, requiring
30 native internal registry ratios. The separate existing Python evaluator
still expands the representation-level internal conventions to 16 models.

This is checkout consistency, **not** cryptographic attestation of installed
RGBeta execution or physical source equivalence.

## User commands (PowerShell, repo root)

```powershell
$zip = Join-Path $env:USERPROFILE "Downloads\T3_B2_Internal_Action_Git_Handoff_Fix.zip"
Expand-Archive -LiteralPath $zip -DestinationPath . -Force
.\.venv\Scripts\python.exe -m pytest `
 tests\test_physical_b2_internal_git_handoff.py `
 tests\test_physical_b2_rgbeta_internal_action_json_export.py `
 tests\test_physical_b2_rgbeta_internal_action_probe_source.py `
 tests\test_physical_b2_rgbeta_internal_action.py `
 tests\test_physical_b2_rgbeta_ward_witness.py `
 tests\test_physical_b2_gauge_coupling.py -q
$env:T3_EXPERIMENTAL_GENERAL_FERMIONS="0"
$env:T3_EXPERIMENTAL_GENERAL_QUADRATICS="0"
$env:T3_EXPERIMENTAL_GENERAL_QUARTICS="0"
$env:T3_EXPERIMENTAL_FLAVOR_MASSES="0"
.\.venv\Scripts\python.exe scripts\RunPhysicalB2RGBetaInternalActionNative.py `
 --output "output\general_observable\physical_normalisation\b2_rgbeta_action_native_04\internal_action_registry.json"
.\.venv\Scripts\python.exe scripts\RunPhysicalB2RGBetaInternalAction.py `
 --native "output\general_observable\physical_normalisation\b2_rgbeta_action_native_04\internal_action_registry.json" `
 --output "output\general_observable\physical_normalisation\b2_action_local_04"
```

Do **not** directly invoke `wolframscript` for this probe: native capture
requires the Python-verified Git metadata injection.

## Evidence/status

- Pure Python/isolated mocks + prior B2 tests: **84/84 passed** in build environment.
- Installed RGBeta native run: **not performed in build environment**.
- No physical action, Matchete-RGBeta signed factor, or numeric RGBeta
  component-axis certificate follows from these internal measurements.
- Physical certification: **0/16**, by design.
