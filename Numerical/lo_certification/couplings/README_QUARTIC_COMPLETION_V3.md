# Workstream 2 — all-16 baseline HHSS portal completion and extra-invariant reference tensors

This is an additive, **conditional** quartic geometry batch. All files remain under Workstream 2 exclusive ownership; the existing physical solvers and sources are **read-only**. No new EFT stages.

## Implemented

1. Load the 16 actual source-hashed Matchete native HHSS quartic-channel caches and existing illustrative RGBeta running snapshots at both `mu_F=100000 GeV` and `mu_S=1000 GeV` (MSbar). Confirm local source parser drift through the existing structural recheck.
2. Construct singlet/adjoint RGBeta portal tensors **as exact real fourth derivatives** with `V4=-L4`, interleaved complex-to-real fields, and `V4=lambda_abcd phi_a phi_b phi_c phi_d/4!`.
3. For each scalar species, preserve each native Matchete portal column `N`, and append exactly the RGBeta directions needed to complete its span: `B=[N,R_complement]`. Solve the exact full-component identity `B X=R`. No least-squares fit and no alias presumed.
4. Form per-channel numerical portal tensors separately from `lambdaH1`, `lambdaH1Adj`, `lambdaH2`, `lambdaH2Adj` at both scales. Assert numerically that their sum equals the direct RGBeta tensor at that same scale. Missing couplings block the entire scale's numeric tensor rather than introducing zeros.
5. Construct the SU(2)-invariant boson-symmetrized **reference** tensor space for each extra HHSS field-content sector listed by the RGBeta candidate registry. Check rank, hypercharge, Ward generators, and identical-boson symmetries. Export exact canonical real/imag fourth derivatives of `V4=-cP-c*P*` for every reference direction. This does **not** identify a reference vector with any particular Matchete or RGBeta coupling.
6. Explicitly expose unresolved physical normalization, couplings, and scalar phases. No physical `C5`, no replacement of the existing seven-term EFT1 solver or joint scalar matcher.

The previous source-derived audit identified **14 missing native baseline portal directions** across 10 models: A (3), D (3), E (8), none for B/C. The prior candidate registry contains **24 independent extra HHSS entries**, structurally distinct from those 14 baseline portal directions. The new run should confirm these totals against the actual local files; a discrepancy is recorded or blocks completion. Missing two-scale portal values block *numerical* portal arrays for that scale while retaining symbolic geometry.

## Files to install (complete new files)

- `Numerical/lo_certification/couplings/portal_completion.py`
- `Numerical/lo_certification/couplings/extra_hhss_geometry.py`
- `Numerical/lo_certification/couplings/quartic_completion_batch.py`
- `Numerical/lo_certification/couplings/README_QUARTIC_COMPLETION_V3.md`
- `scripts/RunLOQuarticCompletion.py`
- `tests/test_lo_couplings_completion.py`
- `tests/test_lo_couplings_completion_integration.py`
- `output/general_observable/lo_certification_couplings/QUARTIC_COMPLETION_SCHEMA_V3.json`

This package **requires** the previously installed `T3_WS2_EFT1_Coupling_Batch_SourceDriftFix` and `T3_WS2_Portal_Normalization_V2` packages; it does not modify their files.

## PowerShell

From `C:\thesistemp\T3neutrinomassmodelpipeline`:

```powershell
Expand-Archive -LiteralPath "$HOME\Downloads\T3_WS2_Quartic_Completion_V3.zip" -DestinationPath . -Force
python -m unittest discover -s tests -p 'test_lo_couplings_*.py' -v
python -m scripts.RunLOQuarticCompletion --repo-root . --run-name ws2_quartic_completion_v1

$j = Get-Content 'output\general_observable\lo_certification_couplings\ws2_quartic_completion_v1\summary.json' -Raw | ConvertFrom-Json
$j | Select-Object models_processed, all16_model_audits_passed, baseline_portal_directions_completed, extra_HHSS_registry_entries, physical_C5_calculated
$j.models | Format-Table model, baseline_portal_directions_completed, extra_HHSS_registry_entries, physics_status -AutoSize
```

If the output directory already exists, choose a **new** `--run-name`. The runner never overwrites results. If model tests fail, do not use the corresponding conditional arrays as physically certified input.

## Generated output

`output/general_observable/lo_certification_couplings/ws2_quartic_completion_v1/summary.json`, and for each `<slug>` (e.g. `A_m4`, `B_m1`):

- `portal_and_extra_status.json`: exact change-of-basis coefficient matrix, relative source file paths and SHA-256, scale metadata, flags, blockers.
- `baseline_native_plus_rgbeta_portals_muF_muS.npz`: `muF_channel_tensors`, `muF_channel_labels`, `muF_portal_sum`, likewise `muS_*`. Each tensor has shape `(4,4,n,n)`; channel stack `(k,4,4,n,n)`. These are **portal-only** tensors in the original S1/S2 canonical real-component frame.
- `exact_extra_HHSS_reference_tensors.json`: per-monomial list of unresolved named RGBeta symbols, reference invariant rank and exact real/imag HHSS tensor components with zero-based `[H_a,H_b,S_c,S_d]` indices. Values describe **unit-strength algebraic reference invariants** only, not physical named couplings.

Both scales use the **same** quartic coordinate frame, not an asserted scalar mass-eigenbasis. For split-real and mixed B/C cases, the exact independent source-to-boundary rotation at `mu_F` and at `mu_S` is a separate unresolved certification. To contract with the existing RGE or joint scalar matching code, those rotations, the complete extra HHSS coefficients and the signed Matchete/Weyl normalization must be certified first. No unverified full solver input is created.

## Physics interpretation and blockers

The augmented portal algebra is exact in the declared component convention, but it does not establish absolute Matchete *versus RGBeta runtime* normalization or the phase convention of an extra CG. The native `lambdaT3` channel, 24 named additional HHSS couplings, extra quadratic RG-closed channels, complete SM Weyl Yukawa and gauge embedding, and the GammaCC/Weyl factor remain outstanding. The illustrative RG trajectory at `(mu_F,mu_S)` is *not* a certified two-scale physical EFT1 background.

The fixed matching sequence remains exactly:

`UV -> integrate out ALL heavy BSM fermions together at mu_F -> EFT1(H,S1,S2,SM fermions) -> integrate out BOTH scalar multiplets together at mu_S -> SMEFT`.

The next physical dependency is to determine the **named** Matchete ↔ RGBeta transformations for the extra reference monomial sectors by comparing exact source `GroupInvariant`/native CG components, including conjugate SU(2) intertwiners, and then obtain their actual running coefficients at both scales. Do not assign unit reference coefficients as numeric couplings.
