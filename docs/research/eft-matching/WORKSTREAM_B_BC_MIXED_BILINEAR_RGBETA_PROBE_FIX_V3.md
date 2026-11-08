# Workstream B — B/C(-1) live RGBeta mixed-bilinear probe, v3

## Observed Windows failure

The v2 probe resolved `RGBeta`PackageScope`Lam` correctly, then raised:

```
Rule::rhs: Pattern t3BM_a appears on the right-hand side of rule t3BM_a -> t3BM_a.
RGBeta mixed bilinear probe blocked
Wolfram exit code: 5
```

In Wolfram Language, `t3BM_a` is **pattern syntax**, not a scalar index symbol.
When the installed RGBeta `Lam` constructs its coupling/permutation rules it
encounters the pattern on the right-hand side of a rule. The top-level `Check`
then discards the whole report, explaining the missing v2 JSON file. This is
an error in the diagnostic probe, not physical evidence against the mass RGE.

## Replacement

Replace the **entire** file:

* `scripts/ProbeT3BCMixedBilinearRGBeta.wl`

with the included version. The native call uses the plain symbols:

```
RGBeta`PackageScope`Lam[t3BMa,t3BMb,t3BMc,t3BMd,True]
```

The standalone test file at
`tests/test_t3_bc_mixed_bilinear_rgbeta_probe_contract.py` is also a full
replacement. It includes regression and negative-control checks for the
pattern-name problem, as well as the earlier package-scoped symbol and
association-key issues. No Matchete, verifier, running, or matching source is
changed.

New runtime behavior: if one or both native calls emit Wolfram messages,
each result is serialized as `{"status":"Blocked",...}` and the probe exits
with code 7 **after writing the JSON** (provided export is possible). It never
marks a failed evaluation successful and does not certify physical mapping.

## Exact PowerShell commands

From `C:\thesistemp\T3neutrinomassmodelpipeline`:

```powershell
Expand-Archive -LiteralPath 'C:\path\to\Workstream_B_BC_MixedBilinear_RGBeta_Probe_Fix_v3.zip' -DestinationPath . -Force

python -X utf8 -m pytest -q tests\test_t3_bc_mixed_bilinear_map.py tests\test_t3_bc_mixed_bilinear_rgbeta_probe_contract.py

wolframscript -file scripts\ProbeT3BCMixedBilinearRGBeta.wl output\audits\workstream_b_bc_mixed_rgbeta_v3.json
Write-Host "Wolfram exit code: $LASTEXITCODE"

python -c "import json; p='output/audits/workstream_b_bc_mixed_rgbeta_v3.json'; x=json.load(open(p,encoding='utf-8')); print({m:(v.get('status'),v.get('error')) for m,v in x['models'].items()})"
```

Only if the native run returns exit code 0 and both statuses are
`LiveRGBetaMixedBilinearRegisteredAndLamDetected`, run:

```powershell
python -X utf8 -m scripts.VerifyT3BCMixedBilinearMap `
  --matchete output\audits\workstream_b_bc_mixed_matchete_v1.json `
  --rgbeta output\audits\workstream_b_bc_mixed_rgbeta_v3.json `
  --output output\audits\workstream_b_bc_mixed_bilinear_map_v3.json
```

The verifier intentionally fails closed on blocked native evidence.

## Certification boundary

Passing these probes checks native registered coupling metadata and that the
installed RGBeta mass-augmented `Lam` contains the mixed mass coupling.
The Python verifier checks an exact conditional same-component convention
map and algebraic mass-Hessian identities. Neither performs independent
kinetic/generator field-component basis alignment nor proves the complete
physical matching convention. Keep `physical_mapping_certified=false`.

## Local verification

`14 passed` in Python tests (existing 9 algebra tests plus 5 probe contract
tests). A Wolfram kernel was unavailable in the preparation environment;
v3 has not been live-run against the user's installed RGBeta.
