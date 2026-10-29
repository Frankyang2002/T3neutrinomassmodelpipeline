# Overnight general-scalar six-benchmark exploratory campaign

New files, no replacements:

- `scripts/RunT3GeneralScalarOvernight.py`
- `scripts/RunT3GeneralScalarOvernight.ps1`
- `scripts/CompareT3GeneralScalarSharedA2.py` (full replacement of earlier shared comparator)
- `tests/test_t3_general_scalar_overnight.py`
- `tests/test_t3_general_scalar_shared_a2.py` (full replacement of earlier tests)
- `docs/T3_GENERAL_SCALAR_OVERNIGHT.md`

Includes the full `scripts/CompareT3GeneralScalarSharedA2.py` and its tests, and requires the existing `scripts/RunT3GeneralScalar48Campaign.py`, fitted historical inputs, and original candidate matching sources.

## Work performed

1. Reuse per-model saved general-quartic EFT1 RGBeta exports from `bulk48_v1` without Wolfram regeneration. Source hashes are recorded. `bulk48_v1` is not modified.
2. Call the real existing `RunT3GeneralScalar48Campaign --phase fit --resume` with scrambled 18-dimensional Sobol and local least-squares fitting, three scenarios, all 16 models. The bulk script itself reuses specially validated B(-1)/C(-1) fits and some B(+1) accepted cases by design, and checkpoints case `fit.json` results. A new bulk run ID ensures the remaining cases are evaluated rather than merely using the existing saved campaign summary.
3. Evaluate shared same-scenario A(-2) Yukawas on all models using `CompareT3GeneralScalarSharedA2`, when the donor exists, producing 48 comparisons. `--allow-poor-a2` additionally allows numerical comparisons with a **clearly labelled** poor-fit donor; it does not turn the donor into a good fit.
4. Generate a sortable/filterable `index.html` and 96-row `comparison.csv`/`.json` with failures and conditional physics tags. These are summary comparisons, **not** scale-dependent neutrino RGE traces or T1/T2/T3 gamma sensitivity plots.

**Limitations:** General-quartic bulk campaign omits the additional holomorphic scalar quadratic mixing terms and full sequential Wilson matching. Mixed B/C cases are from conditional separately validated one-scale fits; some legacy matching sources are used. No guarantee Sobol finds a mass-safe A(-2) full-UV donor or acceptable fits. A non-mass-safe donor cannot be shared. The new `continuous_v1` complete-quadratic trajectory is not equivalent to the old bulk campaign fitting evaluator and is not silently substituted for it.

## Windows PowerShell, repository root

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
python -X utf8 -m pytest -q tests\test_t3_general_scalar_overnight.py tests\test_t3_general_scalar_shared_a2.py

# Confirm prerequisite beta exports and prepare new isolated campaign without doing numerics.
python -X utf8 -m scripts.RunT3GeneralScalarOvernight --run-id overnight_six_v1 --prepare-only

# Run all three scenarios, Sobol + local refinement + shared comparisons + sortable reports.
# The command keeps running until finished. In a separate PowerShell window:
powershell -NoProfile -ExecutionPolicy Bypass -File .\scripts\RunT3GeneralScalarOvernight.ps1 -RunId overnight_six_v1 -SobolPoints 128 -MaxNfev 100 -AllowPoorA2

# Reuse completed case checkpoints and shared results after interruption.
powershell -NoProfile -ExecutionPolicy Bypass -File .\scripts\RunT3GeneralScalarOvernight.ps1 -RunId overnight_six_v1 -SobolPoints 128 -MaxNfev 100 -AllowPoorA2 -Resume

# After successful completion:
Start-Process .\output\studies\general_scalar_extended\six_benchmarks\overnight_six_v1_six\index.html
Get-Content .\output\studies\general_scalar_extended\six_benchmarks\overnight_six_v1_overnight.log -Tail 50
```

Avoid calling the original campaign separately with `--phase all` after preparing beta exports: this would attempt expensive Wolfram regeneration and may collide with the copied files.

*For genuinely unattended Windows execution*, launch a separate PowerShell process rather than leaving the controlling terminal or PC in sleep mode. Ensure Windows is configured **not to sleep** and sufficient disk space and numerical dependencies are available. Logs are flushed after every printed line. Do not start two processes with the same run ID.

Results: `output/studies/general_scalar_extended/bulk_48/overnight_six_v1_bulk48/` and `output/studies/general_scalar_extended/six_benchmarks/overnight_six_v1_six/`. PowerShell and Python logs are under `six_benchmarks/`.
