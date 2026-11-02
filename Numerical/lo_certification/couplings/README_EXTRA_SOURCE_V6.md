# Workstream 2 V6 — live Matchete extra-quartic source completion

## Verified source hierarchy

**UV generator:** `Lagrangian/interactions/GeneralScalarOperators.wl` defines an opt-in set of independent *complex* extra quartic couplings `lambdaT3ExtraQ...I...`, registered with `SelfConjugate -> False`. Each term enters through `PlusHc[lambda[] product CG]`. The source calls the bosonic SU(2) tensor symmetrizer and checks `GeneralScalarVerifiedQuartics.wl`. The actual `BuildT3Lagrangian` path enters this generator through `T3Topology.wl`.

**EFT1 seeds:** `output/**/eft1_after_F_scalar_quartic_seed.json` are *not* the same as the UV generator. V5 found zero native terms matching its 24 extra-HHSS candidate operators in the inspected `output/full16_symbolic` seeds. V6 scans *all* local EFT1 quartic seed candidates (not just this directory), independently checking CG registries and field monomials.

**Experimental UV source exports:** `*.interaction.inputform.txt` produced by the repository's `Lagrangian/RunFullGeneralScalarBuild.wl` are read-only additional evidence of an accepted UV operator. They are not EFT1 matching boundary conditions.

**Live CG:** `scripts/RunLOExtraQuarticCGExport.wl` builds the opt-in renormalizable UV Lagrangian using the existing Matchete builder under `UsingFrontEnd`, requires `FullValidation`, and uses native `GetCGTensor` to export the actual registered extra CG coefficients in InputForm. Complex `I`, exact radicals and rationals are preserved; no float coefficients are accepted by the Python reader. This is a **source-only** live operation. There is no call to `Match` or physical threshold matching.

**RGBeta:** `RGE/general/generated/general_rgbeta_quartic_registry.json` and the prior V4 named-source projections provide 50 general quartic declarations, of which 24 are HHSS. Their names and tensor direction ordering do not by themselves certify equivalence to the Matchete `lambdaT3ExtraQ...` basis.

## Operator and threshold policy

Use the fixed two-threshold chain:

`UV -> integrate out all BSM heavy fermions together at mu_F -> EFT1 including both S1 and S2 -> integrate out both scalars together at mu_S -> SMEFT`.

V6 does **not** replace `WilsonTensorRGE`, `GeneralEFT1C5Mixing`, `IntegratedEFT1JointScalarBridge` or the joint scalar matcher. It exports no full numerical HHSS tensor and no C5, because mapping and running data are missing. Uncertified channels are **not** set to zero.

The same coordinate basis is required for comparing polynomial tensor components:

`V4 = lambda_abcd phi_a phi_b phi_c phi_d / 4!`.

The optional ratio probe only proves possible CG polynomial proportionality in a tested embedding. It does **not** include the verified RGBeta `AddQuartic` prefactor, `SelfConjugate` convention, the Matchete `PlusHc` factor, or the physical triplet/pseudoreal rotation. Its outputs must **never** be copied as physical couplings without separate independent validation.

## Files installed (new only; no shared source modified)

- `Numerical/lo_certification/couplings/extra_quartic_source_v6.py`
- `Numerical/lo_certification/couplings/extra_quartic_batch_v6.py`
- `Numerical/lo_certification/couplings/README_EXTRA_SOURCE_V6.md`
- `scripts/RunLOExtraQuarticSourceV6.py`
- `scripts/RunLOExtraQuarticCGExport.wl`
- `tests/test_lo_couplings_extra_quartic_v6.py`
- `tests/test_lo_couplings_extra_quartic_v6_integration.py`
- `output/general_observable/lo_certification_couplings/EXTRA_QUARTIC_SOURCE_SCHEMA_V6.json`

## PowerShell

Unzip into `C:\thesistemp\T3neutrinomassmodelpipeline`:

```powershell
Set-Location 'C:\thesistemp\T3neutrinomassmodelpipeline'
Expand-Archive -LiteralPath "$HOME\Downloads\T3_WS2_Extra_Quartic_Live_Source_V6.zip" -DestinationPath . -Force
python -m unittest discover -s tests -p 'test_lo_couplings_*.py' -v
```

**Stage 1 — always runnable: all-16 local source scan, no Wolfram required:**

```powershell
python -m scripts.RunLOExtraQuarticSourceV6 --repo-root . --run-name ws2_extra_source_v1
Get-Content 'output\general_observable\lo_certification_couplings\ws2_extra_source_v1\summary.json' -Raw
```

**Stage 2 — actual native extra CGs, requires your installed Wolfram + Matchete front end:**

```powershell
Get-Command wolframscript -ErrorAction SilentlyContinue
python -m scripts.RunLOExtraQuarticSourceV6 --repo-root . --run-name ws2_extra_live_v1 --run-live
Get-Content 'output\general_observable\lo_certification_couplings\ws2_extra_live_v1\summary.json' -Raw
```

Optionally add `--probe-orientations` to the live command to compare the *exported* Matchete CGs with V4 RGBeta invariant tensor projections. This symbolic step can be expensive, and only produces candidate exact tensor ratios. If the live stage cannot validate all its source exports, the CLI exits nonzero and records an error, without fabricating a physical certificate. For non-PATH Wolfram installations, use `--wolframscript 'C:\Program Files\Wolfram Research\...\wolframscript.exe'` with the actual local path.

The exporter disables optional general fermions, flavor masses and general quadratics **only for this UV scalar-source extraction**, not for the actual full EFT1 project. This does not change the required physical two-threshold integration sequence.

## Summary fields

- `all16_source_inventory_passed`: the 16-model static source inventory and scan executed successfully, **not** a physical EFT1 certificate.
- `all16_live_Matchete_extra_operator_exports_passed`: 50 native extra operators successfully returned by live Matchete, if requested.
- `live_native_extra_Matchete_CGs_verified`: CG tensors actually exported from `GetCGTensor` (implicit singlet channels have no explicit CG).
- `native_named_extra_HHSS_in_EFT1_seeds`: direct named matches in all scanned EFT1 seed variants. A matching sector alone does not count.
- `physically_certified_HHSS_coupling_maps`: zero unless an independent normalization/phase matching stage is implemented and validated.
- `physical_C5_calculated`: always false.

The next physics task after successful live extraction is to cross-certify **the exported Matchete complex operator** against **RGBeta `AddQuartic` including its exact HC/factor convention**, regenerate a grouped-fermion EFT1 seed containing all extra operators (in the correctly rotated scalar frame), and only then run their beta functions through both scales and feed the complete scalar input to the existing simultaneous scalar-matching solver. Do not treat the UV operator inventory as an EFT1 threshold integration.
