# Joint scalar-boundary and Yukawa Sobol exploration

This is a targeted candidate-level general-quartic search for `A_-2` in `small_uv_running_5e-9`. It reuses `RunT3GeneralScalar48Campaign`, Stage30/31 UV RGEs, native EFT1 RGBeta equations, and the existing model-specific C5 evaluator. It does not modify historical inputs or physics sources.

## Parameters

Sobol searches 18 real Yukawa coordinates, **two independently log-uniform UV scalar mass factors** multiplying the original mass *not squared*, and up to 8 real ordinary scalar quartic boundary values other than `lambdaT3`. The scenario topology coupling stays `5e-9`. Existing `extra-real` EFT1 generalized quartics remain at their chosen boundary value. Heavy fermion masses remain fixed.

Stage 1 evolves UV+EFT1 and screens the full candidate, without computing C5. Stage 2 runs the existing C5/low-energy neutrino evaluator for screen-passing candidates and performs bounded least-squares Yukawa fitting for up to `--local-seeds` per round. The original 5.86 chi-square and mass/perturbativity criteria remain unchanged. A result below 5.86 is **candidate accepted**, not full sequential-Wilson physical matching certification.

**Limitations:** The full-UV mass overshoot around 40 TeV may persist because it is dominated by UV-running terms insensitive to the initial scalar mass squared. This search may find zero passing candidates; no result is guaranteed. This does not include the independent A(-4) legacy factor-of-two C5 issue.

## Run (PowerShell)

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
python -X utf8 -m pytest -q tests\test_t3_joint_scalar_yukawa.py
python -u -X utf8 -m scripts.SearchT3JointScalarYukawa `
    --run-id joint_a2_fulluv_v1 `
    --beta-campaign overnight_six_v1_bulk48 `
    --sobol-points 128 --max-rounds 8 `
    --max-nfev 100 --local-seeds 3 `
    *> .\joint_a2_fulluv_v1.log
```

Resume with all identical arguments and `--resume`. State checkpoints are written after every completed round; incomplete rounds restart. `--max-rounds` must remain unchanged on resume. For a quick 32-point test, use a **different run ID** and `--sobol-points 32 --max-rounds 1 --max-nfev 0`.

Output: `output/studies/general_scalar_extended/bulk_48/<run-id>/state.json` and `best_fitted_config.json` if a mass-safe evaluated fit exists. A passing scalar screen is not the same as an accepted neutrino fit.
