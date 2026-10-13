Optimised T3 RGE — 16-model display grouping correction
=======================================================

Replace the file below in the existing LOCAL repository:
  Numerical/orchestration/RunFittedRGEDisplay.py

PowerShell (from the repository root, with .venv activated):
  python -m Numerical.orchestration.RunFittedRGEDisplay --run --no-figures
  Start-Process "Reports\output\display\optimised_rge_smallY_smallL.html"
  Start-Process "Reports\output\display\optimised_rge_parameters_smallY_smallL.html"

What was wrong:
  Original dashboard expects <study>/<model>/data/running_diagnostics.json.
  The fitted runner cached <study>/<model>/running_diagnostics.json.
  It consequently named every curve after the shared scenario smallY_smallL.

What changed:
  The runner builds the display from temporary model/data/running_diagnostics.json
  views of the completed trajectories, with short labels T3-A-m2 etc.
  It does NOT change the stored RGE files, rerun integrations or overwrite
  the original fixed-comparison outputs. The parameter page remains available.

Do not use --force: the 16 already completed fitted trajectories are reused.
