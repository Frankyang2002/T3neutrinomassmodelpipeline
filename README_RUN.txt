T3 partial six-model diagnostic run: rendering fix
==================================================

Replace:
  scripts/RefreshOptimisedYukawaThesisTables.py

The UV/EFT1 physics is unchanged. This fix:
- supports partial --models runs with no A_-4 reference row;
- writes -- for unavailable A(-4)-normalised epsilon_C5;
- keeps D(-2)-normalised results when D_-2 was evaluated;
- removes the combined dual_amplitude.tex table (already excluded from thesis);
- writes summary.json from --render-only, including missing-reference diagnostics;
- returns success for a completed selected-model pilot without requiring an
  unselected A_-4 reference row.

From project root, use the exact same options that created the six-model run:

python -X utf8 -m scripts.RefreshOptimisedYukawaThesisTables `
  --run-id optimised_yukawa_zero_conj_six_v1 `
  --models A_-2,A_+0,D_-2,D_+0,E_-2,E_+0 `
  --full-cases `
  --zero-extra-yukawas `
  --allow-stale-symbolic-source `
  --render-only

No RGE rerun is needed: results.json already contains the completed model data.

NOTE: A(-4)-normalised residuals cannot be calculated in this six-model pilot
because the selected models exclude A_-4. They will be available in an all16
run, which includes the reference row.
