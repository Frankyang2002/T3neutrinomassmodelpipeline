# Workstream 2 / V4 — named RGBeta extra-HHSS projection audit

The preceding V3 `RunLOQuarticCompletion` succeeded for all 16 models. V4
continues that exact source audit, but **does not certify** the physical
Matchete↔RGBeta coefficient dictionary. It is deliberately not a replacement
for the existing beta-function, EFT1 mixing or joint-scalar solvers.

## Project paths (new files only)

- `Numerical/lo_certification/couplings/named_hhss_projection.py`
- `Numerical/lo_certification/couplings/named_hhss_batch.py`
- `Numerical/lo_certification/couplings/README_NAMED_HHSS_PROJECTION_V4.md`
- `scripts/RunLONamedHHSSProjection.py`
- `tests/test_lo_couplings_named_hhss.py`
- `tests/test_lo_couplings_named_hhss_integration.py`
- `output/general_observable/lo_certification_couplings/NAMED_HHSS_PROJECTION_SCHEMA_V4.json`

All generated per-run products are under
`output/general_observable/lo_certification_couplings/ws2_named_hhss_v1/`.
No existing `RGE/`, `Lagrangian/`, or `Numerical/joint_scalar_matching/`
files are edited. Existing V1/V2/V3 results are read only.

## Reproducible PowerShell commands

```powershell
Set-Location 'C:\thesistemp\T3neutrinomassmodelpipeline'
Expand-Archive -LiteralPath "$HOME\Downloads\T3_WS2_Named_HHSS_Projection_V4.zip" -DestinationPath . -Force
python -m unittest discover -s tests -p 'test_lo_couplings_*.py' -v
python -m scripts.RunLONamedHHSSProjection --repo-root . --run-name ws2_named_hhss_v1 --v3-run-name ws2_quartic_completion_v1
$j = Get-Content 'output\general_observable\lo_certification_couplings\ws2_named_hhss_v1\summary.json' -Raw | ConvertFrom-Json
$j | Select-Object all16_source_projections_passed,extra_HHSS_source_symbols_examined,valid_canonical_orientation_symbols,physical_named_maps_certified,physical_C5_calculated
$j.models | Format-Table model,status,source_projection_count,valid_canonical_orientation_symbols,orientation_ambiguous_symbols -AutoSize
```

If a model fails, inspect
`ws2_named_hhss_v1/<model_slug>/failure.json`. Re-run with a new run name
rather than overwriting old evidence. For a single-model diagnostic use
`--model B_-1` and a different `--run-name`.

## Actual source verification

1. Check the V3 source-reference symbol inventory and SHA-256 of the V3
   immutable output for each model.
2. Read `RGE/general/generated/general_rgbeta_quartic_registry.json` and
   the **actual generated** `RGE/general/generated/<model>_extra_quartics.wl`.
   For reused names in B/C, read the actual `RGE/running/rgbeta/T3RGBetaModel.wl`.
   Every `AddQuartic` field list must agree with the registry. Each new source
   `GroupInvariant` must match its registry expression literally apart from
   whitespace and surrounding parentheses.
3. Parse the declared `del`, `eps`, `delS2`, and
   `T3RGBetaS2TripleInvariant` vocabulary only; unsupported functions fail
   closed. The triple-triplet contraction is computed from the explicit
   six-epsilon definition present in `T3RGBetaModel.wl`, using the usual
   symmetric-pair normalization `D(1,01)=1/sqrt(2)` as a **declared reference
   convention** rather than a source-certified RGBeta normalization.
4. Generate the full exact complex index coefficient array, project each
   identical-boson pair, and test representation singlet membership by
   expressing it in the V3 canonical SU(2) invariant reference basis.
   The equality of the **entire tensor**, not a least-squares projection,
   determines whether it belongs to the invariant span.
5. Try pseudoreality duality maps (epsilon for doublets, charge-conjugation
   matrix for triplets), shared across identical fields. Record **all**
   compatible maps and exact coordinate coefficients, never arbitrarily
   select a single preferred map or phase.
6. Check availability of each exact extra coupling symbol independently
   in the local running snapshots at both `muF` and `muS`.
7. Read the native `output/full16_symbolic/.../eft1_after_F_scalar_quartic_seed.json`
   and check whether a named Matchete extra coupling is present. Finding a
   name does not certify its tensor CG normalization or hermitian-conjugate
   signs.

For a complex invariant P, the reference component convention from V3 is

    V4 = -c P - conjugate(c) conjugate(P)
       = lambda_abcd phi_a phi_b phi_c phi_d / 4!

The source RGBeta projection coefficients in V4 represent the **formal
`GroupInvariant` polynomial** before RGBeta's package-specific `AddQuartic`
multiplicities, conjugation handling, and Lagrangian-to-potential sign.
They are **not** a usable physical coefficient conversion matrix.

### Interpretation

- `candidate_orientation_projection_maps`: exact algebraic projections for
  possible index embeddings, with the names of the dual maps used. Some
  distinct embeddings can agree numerically while still being physically
  different conventions.
- `valid_canonical_orientation_symbols`: count of names with one or more
  exact source-to-reference projections; **not** number of physically
  normalized couplings.
- `physical_named_maps_certified`: remains zero by design in V4; it can only
  change following a source-verified Matchete CG and RGBeta runtime convention
  audit.
- `physical_C5_calculated`: false, because physical EFT1 couplings and the
  independent spinor/Majorana signed conventions remain unresolved.

Failures are written per model and make the overall batch return nonzero.
No unknown coupling value is inserted as zero, and no partial tensor is
submitted as a complete input to the existing seven-term EFT1 Wilson RGE
or the simultaneous scalar matching.

The physical threshold sequence remains exactly

    UV -> integrate out ALL heavy fermions at muF -> EFT1 with both S1 and S2
       -> integrate out BOTH scalar multiplets together at muS -> SMEFT.
