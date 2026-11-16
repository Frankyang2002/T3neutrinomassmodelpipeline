# T3 Zee-style sensitivity report — uniform +0.01 BSM scalar UV benchmark

## Exact plotting workflow found on GitHub

- Original/unmodified generator: `scripts/PlotT3GeneralZeeOvals.py`
- Original 16-model report: `Reports/output/thesis_general_5_3/zee_ovals_all16_dualview_v5/`
- Original gamma sensitivity: `scripts/RunT3GeneralZeeSensitivity.py`
- NEW inputs: `output/studies/thesis_qualitative/uniform_scalar_001/campaigns/uniform001_all16_reference_v1/`

The original report used **different** UV scalar boundaries and must not be
relabelled as the new benchmark. The preview image included in this package
was drawn from the **16 accepted NEW reference fits**; it reuses the original
plotting module's NuFIT/JUNO Gaussian ellipse helper, palette, and model order.
It intentionally contains **no T1/T2/T3 trajectories and no model-response
covariance ovals** because none can be inferred from single baseline points.

The figure `model_grid_uniform001_reference_only.png` is a reference-centre
preview, **not** a completed sensitivity plot. In the new accepted reference
fits, lambdaT3 = 0.1 for A–D and 0.05 for E, all other symmetry-allowed BSM
scalar quartics at UV are +0.01, and allowed conjugate Yukawas are zero at
UV. The dimensionful scalar bilinear policy is unchanged, as is SM lambdaH.

NuFIT and JUNO experimental ovals are **independent Gaussian approximations**
constructed from the source plotter. They are not official joint likelihoods.
The neutrino physics is still conditional one-scale C5 matching, not a
certified fully sequential UV -> all fermions -> all scalars -> SMEFT.

## Project-relative files added by this ZIP

- `scripts/RunT3Uniform001ZeeSensitivity.py`: computes genuine T1/T2/T3
  trajectories from the new accepted reference fits using the current physics
  evaluators; saves plotter-compatible `points.csv` and `summary.json`.
- `scripts/BuildT3Uniform001ZeeReport.py`: reference-only preview OR invokes
  the **unchanged** `scripts/PlotT3GeneralZeeOvals.py` on newly computed data.
- `tests/python/test_t3_uniform001_zee_report.py`: isolated input/provenance,
  missing-data, checkpoint, and report gating tests.
- `Reports/output/thesis_general_5_3/zee_uniform001_ovals_v1/`: initial
  report and actual saved preview image.

No existing code, matching, RGE, earlier reports, or GitHub files are edited.

## Exact PowerShell commands (from repository root)

Extract one ZIP preserving project structure:

```powershell
Expand-Archive -LiteralPath .\t3_uniform001_zee_report.zip -DestinationPath . -Force
```

Tests:

```powershell
python -X utf8 -m unittest discover -s tests/python -p "test_t3_uniform001_zee_report.py" -v
```

Recreate reference-only preview from **your local 16 accepted new reference
fits**, using the original plotting module:

```powershell
python -X utf8 -m scripts.BuildT3Uniform001ZeeReport --mode reference
```

Validate all 16 reference inputs without numerical running:

```powershell
python -X utf8 -m scripts.RunT3Uniform001ZeeSensitivity `
  --models all --preflight
```

Run one inexpensive 5-gamma point model to test the actual physics path:

```powershell
python -X utf8 -m scripts.RunT3Uniform001ZeeSensitivity `
  --models A_-4 --points 5 --run-id uniform001_zee_A_m4_pilot_v1
```

Plot that pilot with the unmodified original generator:

```powershell
python -X utf8 -m scripts.PlotT3GeneralZeeOvals `
  --source output/studies/thesis_qualitative/uniform_scalar_001/zee_sensitivity/uniform001_zee_A_m4_pilot_v1 `
  --models A_-4 `
  --output-dir Reports/output/thesis_general_5_3/zee_uniform001_ovals_v1/pilot_A_m4
```

**Full 16-model T1/T2/T3 scan:**

```powershell
python -X utf8 -m scripts.RunT3Uniform001ZeeSensitivity `
  --models all --points 11 --run-id uniform001_zee_all16_v1
```

For an interrupted run, resume using the same model/points/reference/run ID:

```powershell
python -X utf8 -m scripts.RunT3Uniform001ZeeSensitivity `
  --models all --points 11 --run-id uniform001_zee_all16_v1 --resume
```

Then render the **full actual Zee-style 16-panel overview and detail plots**
with the original plotter, rather than a new rendering implementation:

```powershell
python -X utf8 -m scripts.BuildT3Uniform001ZeeReport `
  --mode sensitivity --sensitivity-run-id uniform001_zee_all16_v1
```

Exact outputs:

- `output/studies/thesis_qualitative/uniform_scalar_001/zee_sensitivity/uniform001_zee_all16_v1/points.csv`
- `output/studies/thesis_qualitative/uniform_scalar_001/zee_sensitivity/uniform001_zee_all16_v1/summary.json`
- `Reports/output/thesis_general_5_3/zee_uniform001_ovals_v1/full_sensitivity/model_grid_zee_zoom_out.png`
- `Reports/output/thesis_general_5_3/zee_uniform001_ovals_v1/full_sensitivity/model_grid_zee_zoom_in.png`
- original plotter's per-model PNG, CSV diagnostic tables, and provenance JSON

The default grid has 11 gamma samples in each of T1, T2, and T3 for all 16
representations: **528 saved rows**, with identical gamma=1 baseline reused
for efficiency. Only **496 unique physical evaluations** are needed if all
points finish (31 per model). Each evaluation still runs the full conditional
UV/EFT1/general C5/low-energy physics; runtimes may be substantial.

No physical sensitivity run was performed in this packaging environment; only
reference-preview plotting and synthetic unit tests were executed. The full
16-model reference-fit data were inspected in the GitHub repository.
