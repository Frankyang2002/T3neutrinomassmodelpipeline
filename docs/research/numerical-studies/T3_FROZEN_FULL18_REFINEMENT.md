# Frozen-UV 18-parameter extended refinement

Replace `scripts/RunT3FrozenEstablishedSobol.py` with the complete provided file. Preserve earlier output folders.

To refine the **A(+0) small-frozen** saved best fit (the script's first case in the model registry), starting from `frozen_escalating_v1` without repeating any Sobol points:

```powershell
python -u -X utf8 -m scripts.RunT3FrozenEstablishedSobol `
  --source-campaign overnight_six_v1_bulk48 `
  --resume-from-run frozen_escalating_v1 `
  --run-id frozen_full18_v1 `
  --models A_+0 `
  --refine-only `
  --seed-count 5 `
  --nfev-18 600 `
  --sensitivity-step 0.001 *> .\frozen_full18_v1.log
```

An accepted small-frozen fit proceeds to the large-frozen case; an unsuccessful case stops by default. Checkpoint: `output/studies/general_scalar_extended/bulk_48/frozen_full18_v1/small_uv_frozen_5e-9/A_p0/full18_refinement.json`.

To run Sobol followed by staged 5,8,12,18-coordinate refinement for all models and both frozen scenarios, omit `--refine-only`, `--models`, `--resume-from-run`, and choose a new `--run-id`. Set `--nfev-18 600`.

No physical formula, acceptance cut, or matching parity threshold is modified. Existing legacy-parity failures remain blockers. The 18-coordinate local fit can require substantial evaluations; success is not guaranteed. The current runner resumes per completed round/checkpoint, not an interrupted least-squares call.
