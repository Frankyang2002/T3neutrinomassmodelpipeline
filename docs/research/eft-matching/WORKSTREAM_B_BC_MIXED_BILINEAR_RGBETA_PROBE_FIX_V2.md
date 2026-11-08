# Workstream B: B/C(-1) native RGBeta mixed-mass probe fix, v2

## Observed local failure

The original probe printed `DownValues::sym` and wrote a diagnostic JSON, but the
Python verifier stopped on `B_-1: RGBeta native source contract incomplete`.
**Writing the JSON is not a passing native RGBeta probe.** The per-model
`status` and `error` must be inspected; in particular, a blocked row must not
be turned into a successful mapping.

Two source errors are fixed in the complete replacement of
`scripts/ProbeT3BCMixedBilinearRGBeta.wl`:

1. `DownValues[Symbol[lamName]]` did not evaluate its held argument to the
   package function symbol. The installed RGBeta 1.2.0 source declares
   `PackageScope["Lam"]`; use the literal qualified symbol
   `RGBeta\`PackageScope\`Lam` both for `DownValues` and native evaluation.
2. `info[Fields]`, `info[Invariant]`, etc. can address symbols from a different
   Wolfram context than the actual Association keys. This project already
   handles the issue for native quartics in `scripts/RGBetaNativeAccess.wl`.
   The new probe independently resolves a **unique existing native key** by
   its short name and fails if missing or ambiguous.

The RGBeta upstream sources read for this fix were `Package/RGBeta.m` and
`Package/FieldsAndCouplings.m` from `aethomsen/RGBeta`; `Lam` is the actual
mass-augmented `Lam[a,b,c,d,True]`, and `$scalarMasses` includes the native
`Fields`, `Invariant`, `Coupling`, `SelfConjugate`, and
`UniqueArrangements` records.

## Run on the existing Windows repository

Extract this ZIP at the project root. Only the RGBeta probe is replaced; the
original Matchete probe and Python verifier are preserved.

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline

python -X utf8 -m pytest -q `
  tests\test_t3_bc_mixed_bilinear_map.py `
  tests\test_t3_bc_mixed_bilinear_rgbeta_probe_contract.py

# NEW filename: original v1 report remains for provenance.
wolframscript -file scripts\ProbeT3BCMixedBilinearRGBeta.wl `
  output\audits\workstream_b_bc_mixed_rgbeta_v2.json

# Inspect the native status. Do not infer success from the mere JSON path.
python -c "import json; p='output/audits/workstream_b_bc_mixed_rgbeta_v2.json'; x=json.load(open(p)); print({m:(v.get('status'),v.get('error')) for m,v in x['models'].items()})"

python -X utf8 -m scripts.VerifyT3BCMixedBilinearMap `
  --matchete output\audits\workstream_b_bc_mixed_matchete_v1.json `
  --rgbeta output\audits\workstream_b_bc_mixed_rgbeta_v2.json `
  --output output\audits\workstream_b_bc_mixed_bilinear_map_v2.json
```

The verifier requires native evidence, checks source hashes and rejects
uncertified/malformed inputs. If the Matchete report itself has a blocked row,
inspect that row; do not modify the verifier to accept it.

## Scope

This fixes source/runtime interface bugs, **not** a physical normalization
ambiguity. A conditional mixed-bilinear action equality, even if the verifier
passes, is not a certificate for the common complex component basis or the
SU(2) field/generator intertwiner. Workstream A's final operator certificate
remains a prerequisite. Neither production RGEs nor matching are modified.
