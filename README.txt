T3 Coupled EFT1 Yukawa least-squares refinement
===============================================

NEW FILE (no replacement of existing production files):
 scripts/RunT3CoupledYukawaLeastSquares.py

Requires scripts/RunT3CoupledYukawaOvernight.py, specifically the
FingerprintCompat edition installed earlier, which has the
symbolic_export_compatibility and setup_physics functions.

PowerShell, from project root with activated .venv:

  python -X utf8 -m scripts.RunT3CoupledYukawaLeastSquares --run-id ls_smoke_v1 --models D_+0 --from-run yukawa_overnight_compat_v2 --check-only

Small real numerical smoke test (~8 EFT evaluations):

  python -X utf8 -m scripts.RunT3CoupledYukawaLeastSquares --run-id ls_smoke_v1 --models D_+0 --from-run yukawa_overnight_compat_v2 --allow-stale-symbolic-source --active-counts 2 --max-evals 8 --stage-nfev 2

Preferred first campaign, starting from the already improved Yukawas:

  python -X utf8 -m scripts.RunT3CoupledYukawaLeastSquares --run-id ls_refine_four_v1 --models D_+0,C_+1,A_+0,D_-2 --from-run yukawa_overnight_compat_v2 --allow-stale-symbolic-source --max-evals 120 --max-hours 10

Then the remaining models (fresh run-id; default source paths remain unchanged):

  python -X utf8 -m scripts.RunT3CoupledYukawaLeastSquares --run-id ls_refine_remaining_v1 --models A_-4,A_-2,B_-3,B_-1,B_+1,C_-3,C_-1,D_+2,E_-4,E_-2,E_+0,E_+2 --from-run yukawa_overnight_compat_v2 --allow-stale-symbolic-source --max-evals 120 --max-hours 10

Output (new only):
 output/studies/general_scalar_extended/coupled_yukawa_least_squares/<run-id>/<model-tag>/
    settings.json, best.json, best_fitted_config.json,
    best_general_uv_overrides.json (nonmixed only), evaluations.jsonl,
    summary.json, symbolic_export_compatibility.json

The target is the whitened five-dimensional fit residual from
Numerical/fitting/OscillationFit.py evaluated by case 1111 of the existing
coupled pipeline. The script checks sum(residual**2) against the existing
coupled chi2 before accepting an improvement. It never uses the older
one-scale optimiser. Other Yukawa directions remain fixed during each
5/8/12-active-coordinate least-squares stage.

All source and overnight files are read-only. No remote GitHub changes.
The cached symbolic export override reuses historical RGEs, so numerical
results should not be attributed to newer edited symbolic source.

Static Python syntax checks passed. Full numerical execution requires the
user's local cached symbolic exports and numerical inputs and has not been
performed in this delivery environment.
