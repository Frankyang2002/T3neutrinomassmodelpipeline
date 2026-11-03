# B(-1)/C(-1) native mixed-source exact flavour adapter

## Scope

This is a source-exact **fermion-threshold LLSS tree** adapter, evaluated at the
fitted three-generation Yukawas, heavy Majorana mass eigenvalues, and scalar
mixing `U(muS)`. The archived Matchete PL source and native CG registry are the
starting point; **no** universal `Sym(y1 MF^-1 y2^T)` proxy is used.

- All 10 B(-1) and 12 C(-1) PL sources are handled separately.
- `y1*y1`, `y1*y2`, `y2*y2` and the distinct ordered external-lepton channels
  are retained; C(-1) bridge sources 23/24 are not deleted because they were
  zero at the old diagonal-flavour benchmark.
- Native source signs, `Bar[...]` conjugation, SU(2) CG and `X/Y` scalar
  components are evaluated through the existing `project_ordered_term` and
  commuting-spinor polynomial symmetrizer.
- A new reference check regenerates the archived **ordered** native-CG source
  at the original reference point before any fitted-source calculation.
- The 14 nonmixed model branches, their RGEs, and the two-scalar matching
  kernel remain unchanged.

**Uncertified:** The GammaCC/Weyl overall convention, exact physical
Matchete/RGBeta scalar-coupling mapping and independent hard matching are not
established. The `U(muS)` identification with the Matchete X/Y source basis
remains an exploratory fixed-basis assumption. The audit may still expose
scalar/CG convention differences. Never alter the normalization factor just to
force matrix agreement.

## Full replacement / new files (project-relative paths)

1. `Numerical/running/GeneralT3MixedExactFermionWilson.py` — **NEW** complete backend
2. `scripts/RunT3GeneralSequentialFitAblation.py` — **FULL REPLACEMENT** controller
3. `scripts/RunT3MixedExactSourceSmoke.py` — **NEW** source/CG smoke check
4. `tests/python/test_general_mixed_exact_fermion_boundary.py` — **NEW** regression tests

No code is written to GitHub. The ZIP does not include or change output files.

## PowerShell (repository root)

```powershell
Set-Location 'C:\thesistemp\T3neutrinomassmodelpipeline'
Expand-Archive -LiteralPath "$HOME\Downloads\T3_BC_Mixed_Source_Exact_Backend.zip" -DestinationPath . -Force

python -X utf8 -m pytest -q tests/python/test_general_mixed_exact_fermion_boundary.py
python -X utf8 -m scripts.RunT3MixedExactSourceSmoke

python -X utf8 -m scripts.RunT3GeneralSequentialFitAblation `
  --source-run thesis16_merged_v1 `
  --run-id general_seq_bc_exact_01 `
  --models B_-1,C_-1 `
  --spinor-factor 1 `
  --check-only

python -X utf8 -m scripts.RunT3GeneralSequentialFitAblation `
  --source-run thesis16_merged_v1 `
  --run-id general_seq_bc_exact_01 `
  --models B_-1,C_-1 `
  --spinor-factor 1

python -X utf8 -m scripts.RunT3GeneralMatcheteLogCoefficientAudit `
  --source-run thesis16_merged_v1 `
  --case-run general_seq_bc_exact_01 `
  --run-id log_bc_exact_01 `
  --models B_-1,C_-1 `
  --hierarchy-factors 1,2,4

Import-Csv 'output\studies\general_scalar_extended\general_log_coefficient_audit\log_bc_exact_01\summary.csv' |
  Select-Object model,status,ratio_Matchete_over_EFT_at_MF,best_fit_factor_real_at_MF,relative_matrix_residual_best_factor_at_MF,best_fit_factor_real_at_largest_MF |
  Format-List
```

If a run ID directory already exists, choose a new ID. The old
`general_seq_05` and `log_coefficient_01` remain unchanged.

## Checkpoints

- Smoke: B_-1 10 PL and C_-1 12 PL; archived native CG reconstruction
  residual below `5e-10`; dynamic flavour tensors have the expected shape.
- Campaign: B_-1 and C_-1 each compute 16 U/E/W/S switch combinations.
- Audit: compare each full Matchete log coefficient with the *new* EFT1
  source-exact coefficient. Inspect both best-fit factor and matrix residual.
  A residual below 1% is a goal, not a guaranteed result.
