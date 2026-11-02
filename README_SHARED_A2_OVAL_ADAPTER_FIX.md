# Shared A(-2) general16 — oval renderer adapter correction

This is an **additive local fix**, not a change to the neutrino physics.

## Why the existing run produced index.html but no oval

The numerical campaign saved `A_m2/sensitivity_diagnostics.json`, while the original
`Numerical.plotting.BuildFittedModelEnvelopes` expects the same file under the
literal model name `A_-2/sensitivity_diagnostics.json`. The renderer therefore
skipped every successfully evaluated model. The fix stages **copies** of the
verified scan JSON using the renderer's expected folder names, retaining all
blocked models in the summary. No calculations are replayed for `--render-only`.

## Files
- `scripts/RunT3SharedA2General16Reports.py`: complete replacement.
- `tests/python/test_thesis_shared_a2_general16_reports.py`: complete replacement.

## Install and render seven already completed models immediately

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive "$HOME\Downloads\T3_Shared_A2_Oval_Adapter_Fix.zip" -DestinationPath . -Force
python -X utf8 -m pytest -q tests/python/test_thesis_shared_a2_general16_reports.py
python -X utf8 -m scripts.RunT3SharedA2General16Reports `
  --run-id shared_A2_general16_ovals_v1 `
  --render-only
Start-Process "output\studies\general_scalar_extended\shared_a2_general16_reports\shared_A2_general16_ovals_v1\sensitivity_oval.html"
```

## Generate nine missing general C5 matching expressions

```powershell
$missing = @('C_-3','C_+1','D_-2','D_+0','D_+2','E_-4','E_-2','E_+0','E_+2')
foreach ($model in $missing) {
  python -X utf8 -m scripts.RunT3GeneralMatcheteMatching --model $model --reuse-success
  if ($LASTEXITCODE -ne 0) { Write-Warning "Matchete failed for $model (inspect the model's general_matching.log)" }
}
```

The matching command requires the current installed mass-basis Wolfram runner and
Wolfram/Matchete. The six self-bilinear models are exploratory positive-real
bilinear one-scale matching, **not certified sequential EFT matching**.

## Resume unchanged 16-model sensitivity run

```powershell
python -X utf8 -m scripts.RunT3SharedA2General16Reports `
    --donor-run-id thesis16_general_c5_full_v1 `
    --run-id shared_A2_general16_ovals_v1 `
    --gamma-min 0.7 --gamma-max 1.4 --gamma-points 7 --resume
```

The already completed seven successes are reused, previously blocked models
are retried. Failed matching is recorded, not replaced with a restricted C5.

Expected files:
- `output/studies/general_scalar_extended/shared_a2_general16_reports/shared_A2_general16_ovals_v1/index.html`
- `.../sensitivity_oval.html`
- `.../comparison.csv`
- `.../gamma_points.csv`

Ovals show empirical T1/T2/T3 gamma sensitivity **not uncertainty/confidence**.
The current fit lacks saved scale-resolved trajectories so the old RGE HTML
renderer cannot be validly populated from these endpoint-only outputs.
