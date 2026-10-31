# Workstream B — B/C(-1) mixed bilinear native RGBeta probe, v4

## Outcome from local run of v3

The user's WolframScript runtime successfully resolved `RGBeta`PackageScope`Lam` (`downvalues=1`) and built the B(-1) and C(-1) models, but both models exited with `Blocked: Live massive Lam does not contain added mixed mass coupling`. The v3 probe checked `FreeQ[live, muT3GeneralS1S2, ...]` **before** expanding any native sparse-array rules. This test can miss symbolic coefficients inside an opaque or incompletely traversable `TStructure`/`SparseArray` representation. It does not by itself establish that the native mass vertex is absent. That diagnosis is a plausible implementation failure, **not verified** without examining the live tensor.

## Implementation

The v4 full replacement of `scripts/ProbeT3BCMixedBilinearRGBeta.wl`:

1. Preserves the v3 fixes: literal RGBeta package-scoped `Lam`, ordinary Wolfram index symbols (not `_` patterns), and native-context Association key resolution.
2. Calls the installed native `RGBeta`PackageScope`Lam[t3BMa,t3BMb,t3BMc,t3BMd,True]`.
3. Extracts actual nested `SparseArray` objects and examines `ArrayRules` rather than relying on `FreeQ` applied to the outer `TStructure`.
4. Requires a nonzero count of native sparse-array rules explicitly mentioning `muT3GeneralS1S2` before reporting `LiveRGBetaMixedBilinearRegisteredAndLamDetected`.
5. On failure, retains a per-model `Blocked` JSON report with `native_lam_sparse_array_count`, `native_lam_total_sparse_array_rules_count`, `native_lam_excerpt_inputform`, and up to 12 native rule examples, plus source hashes. The verifier will continue to reject `Blocked`.
6. Does **not** assert RGBeta epsilon component-phase alignment or physical Matchete↔RGBeta coupling certification.

**Only these three project files are packaged:**

- `scripts/ProbeT3BCMixedBilinearRGBeta.wl` — full replacement
- `tests/test_t3_bc_mixed_bilinear_rgbeta_probe_contract.py` — full replacement
- `docs/WORKSTREAM_B_BC_MIXED_BILINEAR_RGBETA_PROBE_FIX_V4.md` — new report

All source code is limited to these audit files. No RGBeta, Matchete, matching, RGE or other workstream production code is modified.

## Installation and run (PowerShell)

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive `
  -LiteralPath "C:\path\to\Workstream_B_BC_MixedBilinear_RGBeta_Probe_Fix_v4.zip" `
  -DestinationPath . -Force

python -X utf8 -m pytest -q `
  tests\test_t3_bc_mixed_bilinear_map.py `
  tests\test_t3_bc_mixed_bilinear_rgbeta_probe_contract.py

wolframscript -file scripts\ProbeT3BCMixedBilinearRGBeta.wl `
  output\audits\workstream_b_bc_mixed_rgbeta_v4.json
Write-Host "Wolfram exit code: $LASTEXITCODE"
```

Inspect results **even when exit code is 7**:

```powershell
python -c "import json; d=json.load(open('output/audits/workstream_b_bc_mixed_rgbeta_v4.json')); [print(m,v.get('status'),v.get('error'), 'witnesses=',v.get('native_lam_mass_sparse_array_rules_count'), 'arrays=',v.get('native_lam_sparse_array_count'), 'rules=',v.get('native_lam_total_sparse_array_rules_count'), 'example=',v.get('native_lam_excerpt_inputform','')[:350]) for m,v in d['models'].items()]"
```

Only if **both** models report `LiveRGBetaMixedBilinearRegisteredAndLamDetected`, run:

```powershell
python -X utf8 -m scripts.VerifyT3BCMixedBilinearMap `
  --matchete output\audits\workstream_b_bc_mixed_matchete_v1.json `
  --rgbeta output\audits\workstream_b_bc_mixed_rgbeta_v4.json `
  --output output\audits\workstream_b_bc_mixed_bilinear_map_v4.json
```

No source hashes or expected physics conventions are changed in `VerifyT3BCMixedBilinearMap.py`.

## Verification and remaining limits

- Python source-interface/regression tests: **16 passed** when executed together with the prior mixed-bilinear Python tests in the local assembled test tree.
- Native WolframScript/RGBeta v4 evaluation: **not available in this environment**; requires user's installed RGBeta runtime.
- The first positive result proves native sparse-array symbol visibility, not the complete physical scalar action equivalence.
- Matchete↔RGBeta overall SU(2) field-basis and epsilon phase are **uncertified**.
- The Higgs mass RGE closure gap and Workstream A final basis certificate remain separate outstanding items.
