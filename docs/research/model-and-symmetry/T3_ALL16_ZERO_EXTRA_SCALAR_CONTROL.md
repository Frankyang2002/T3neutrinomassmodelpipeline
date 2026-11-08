# All-16 zero-additional-scalar control

This is a **separate, non-destructive** diagnostic built against the current
`Frankyang2002/T3neutrinomassmodelpipeline` sources on 30 October 2026. Only
new files are installed. It uses the fitted `optimised_low`/`optimised_high`
**frozen-UV** configurations; it does not refit Yukawas or change GitHub.

## Exactly what is zero

For **each of the 16 ordinary T3 assignments**:

- Set to zero every additional scalar quartic registered in
  `RGE/general/generated/general_rgbeta_quartic_registry.json`, including both
  newly introduced RGBeta coordinates and the *reused* RGBeta coordinates
  which the registry identifies as **additional** operators.
- Set to zero the model-specific extra holomorphic quadratic scalar mass,
  wherever allowed (8/16 models).
- **Do not** set the original `lambdaT3` to zero. Do not set original scalar
  diagonal masses, baseline ordinary quartics, Yukawas, gauge couplings or
  fermion masses to zero.
- Keep the original fitted high/low Yukawas. The boundary at 100 TeV is the
  historical frozen fit, not the full 10 PeV -> 100 TeV UV running.

Two distinct computations share identical initial conditions:

| Run mode | Extra quartics / quadratic at 100 TeV | What happens during 100-to-1 TeV EFT1 running? |
|---|---|---|
| `held_zero` | Exactly 0 | Their beta functions are **projected to 0 in memory**, so they stay 0. Only an artificial restricted-potential control, *not the general scalar RG equations*. |
| `boundary_zero` | Exactly 0 | Evolve full RGBeta candidate beta system; the additional couplings **may be radiatively generated**. |

These answer different scientific questions. The first checks restricted-model
recovery at fixed matching prescription; the second tests whether setting a
zero **boundary condition** is RG stable.

## Install (PowerShell, repository root)

Place the ZIP at your checkout root (for example `C:\thesistemp\T3neutrinomassmodelpipeline`):

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath .\T3_All16_ZeroExtraScalar_Control.zip -DestinationPath . -Force
python -X utf8 -m pytest -q tests\test_t3_all_zero_extra_scalar_control.py
```

First run the read-only source/input preflight:

```powershell
python -X utf8 -m scripts.RunT3AllZeroExtraScalarControl `
    --run-id zero_extra_all16_low_v1 --benchmark optimised_low `
    --models all --modes held_zero,boundary_zero --plan
```

Then compute **all 16 models**, with both settings:

```powershell
python -u -X utf8 -m scripts.RunT3AllZeroExtraScalarControl `
    --run-id zero_extra_all16_low_v1 --benchmark optimised_low `
    --models all --modes held_zero,boundary_zero
```

View the summary (from the project root):

```powershell
$report = '.\output\studies\general_scalar_extended\zero_extra_all16\zero_extra_all16_low_v1\summary.csv'
Import-Csv $report | Format-Table model,mode,lambdaT3_input_abs,lambdaT3_endpoint_abs, `
    max_extra_quartic_endpoint,sin2_theta12,sin2_theta13,sin2_theta23 -AutoSize
```

Optional high-coupling historical fit: use `--benchmark optimised_high` and a **different** `--run-id`, e.g. `zero_extra_all16_high_v1`.

If interrupted, rerun with the **same complete arguments** and append
`--resume`. If inputs or the script changed, use a different run ID. A failed
cell is listed under `failures` in `summary.json`, not substituted with old
numerical data. This is a deliberately non-overwriting run.

## Output directory and diagnostic fields

`output/studies/general_scalar_extended/zero_extra_all16/<run-id>/`

- `control_settings.json`: exact configuration and script SHA256.
- `summary.csv`, `summary.json`: models/modes, original historical predictions,
  candidate predictions, differences, original and endpoint |lambdaT3|,
  maximum induced quartic, extra quadratic magnitude, and blocked cases.
- `<benchmark>/<mode>/<model_tag>/continuous_trajectory.json`: full candidate
  EFT1 scalar RGE trajectory and coupling values; separate endpoint sidecar.
- `<benchmark>/<mode>/neutrino_v1/<model_tag>_neutrino.json`: existing
  conditional matching and SM+Weinberg evaluation.
- `<benchmark>/held_zero/original_kernel_v1/<model_tag>_original_kernel.json`:
  original model-specific one-loop C5 kernel at exactly the same held-zero EFT1
  scalar endpoint, when evaluable; separate from candidate matching.

## Interpreting results correctly

1. At `held_zero`, all additional quartics and the extra quadratic **must**
   remain zero at 1 TeV (verified numerically, otherwise a run fails).
2. At `boundary_zero`, all are zero only at 100 TeV. Nonzero low-scale values
   demonstrate radiative generation by the candidate beta system.
3. Historical and new results can disagree even at held zero because the
   original project and candidate pipeline use **different C5 matching
   prescriptions**. For B(-1)/C(-1), the same-endpoint `original_kernel` test
   isolates this matching discrepancy.
4. All eight model assignments allowing a new holomorphic mass are included.
   When `boundary_zero` regenerates a self-bilinear, the existing A/D/E
   restricted neutrino C5 kernel **does not include its contribution**;
   this row is not a complete new-potential prediction.
5. The candidate Matchete/RGBeta physical coupling map and **full sequential
   Wilson/threshold matching remain uncertified**. These are controlled
   numerical diagnostics, not validated full-model predictions.
6. Zero extra boundary values do not establish technical naturalness. Large
   effects at illustrative `f=1` show conditional sensitivity but not a
   phenomenologically valid fit without a new fit and matching checks.

The output should be compared against the previously saved
`frozen_general64_v1` results **after** confirming the zero-extra controls
behave as described. Do not overwrite those original 64-case outputs.
