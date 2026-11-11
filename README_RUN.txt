Optimised Yukawa thesis refresh (local-only)
===========================================
Project file: scripts/RefreshOptimisedYukawaThesisTables.py

Usage (PowerShell):
  cd C:\thesistemp\T3neutrinomassmodelpipeline
  Expand-Archive -Path "$HOME\Downloads\T3_Optimised_Yukawa_Thesis_Refresh.zip" -DestinationPath . -Force
  python -X utf8 -m scripts.RefreshOptimisedYukawaThesisTables --check-only
  python -X utf8 -m scripts.RefreshOptimisedYukawaThesisTables --run-id optimal_yukawa_thesis_refresh_v1 --allow-stale-symbolic-source

If the numerical source freshness gate is satisfied without compatibility mode,
omit --allow-stale-symbolic-source. Compatibility mode reuses archived symbolic
exports with explicit provenance; it does not certify them as current.

Generate four-axis stage-freeze heatmaps using a DIFFERENT, fresh run ID:
  python -X utf8 -m scripts.RefreshOptimisedYukawaThesisTables --run-id optimal_yukawa_thesis_heatmaps_v1 --full-cases --allow-stale-symbolic-source

If interrupted, repeat the same command with --resume appended.

Outputs:
  output/studies/general_scalar_extended/optimal_yukawa_thesis_refresh/<run-id>/
    manifest.json, results.json, summary.json, observables.csv
    tables/Dm2_masses.tex
    tables/Dm2_mixing.tex
    tables/Am4_masses.tex
    tables/Am4_mixing.tex
    heatmaps/*.png   [with --full-cases only]

Each target model retains its own other couplings, scalar masses, extra Yukawas.
The new donor Y1 and Y2 matrices are transferred identically to every model.
Source donor inputs:
  output/studies/general_scalar_extended/coupled_yukawa_least_squares/
    ls_A4_target001_v1/A_m4/{best.json,best_fitted_config.json,best_general_uv_overrides.json}
    ls_D2_target001_v1/D_m2/{best.json,best_fitted_config.json,best_general_uv_overrides.json}

Physics: UV -> integrate out all fermions -> EFT1 -> integrate out all scalars
-> SMEFT. Exploratory finite-bubble + coupled EFT1, with omitted independent
hard one-loop matching, cached/provisional quartic context. Not precision EFT.

Failure behavior:
  Source matrix integrity checked against all 3 donor files.
  Donor `1111` predictions and masses must reproduce saved LS fits.
  Invalid/missing cases are marked Blocked, never replaced by invented numbers.
  No historical inputs and no GitHub data are modified.
