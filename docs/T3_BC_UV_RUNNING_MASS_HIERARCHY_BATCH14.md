# Batch 14 — B(-1)/C(-1) UV-running mass hierarchy replay

## What this actually does

The saved fresh Batch 13 UV-running fits have \(\chi^2=18376.678022\) for B(-1) and \(20279.293144\) for C(-1). Their stored scalar mass-squared parameters at \(\mu_S=1\) TeV imply running scalar masses \(86.9,135.6\) TeV and \(135.1,167.3\) TeV, respectively. These values are much closer to the heavy fermions near 100 TeV than to the nominal 1 TeV scalar matching scale.

`DiagnoseT3BCUVScalarHierarchy.py` **replays, but does not modify**, the exact fitted 18 real Yukawas from `general_fresh_v1_B_m1_cont` and `general_fresh_v1_C_m1_cont`. It reuses the existing experimental native RGBeta UV and EFT1 beta expressions. It evaluates two different scalar UV boundaries on each *same* fitted Yukawa vector:

- `general_quartics`: precisely the quartic and quadratic-mixing boundary settings recorded in the saved exploratory fit (default 0.1 × |lambdaT3| and mu2=0).
- `restricted_zero`: all seven extra quartics and new mixing set to zero at the UV boundary, with the *identical fitted Yukawas and all other historical inputs*. It is **not** the historical restricted best fit.

This separates the already-present UV/RG mass evolution from changes caused by explicit extra-scalar boundary values. It does not evaluate or optimise C5.

The runner checks SHA-256 fingerprints from the saved source manifest, re-evolves UV 10^7 -> 10^5 GeV and EFT1 10^5 -> 10^3 GeV using the existing `DOP853` solver and candidate RGBeta source, and **requires** the replay to reproduce the stored final mass-squared values, mixing magnitude and heavy-fermion masses. Source changes or inconsistent outputs stop the calculation.

For each sampled running scale it records the two eigenvalues of the chosen 2x2 mixed-doublet mass matrix, the heavy fermion singular masses (above muF; below muF the **frozen muF reference values**), and ratios such as \(m_{S,\max}/\mu\), \(m_{S,\max}/M_{F,\min}\), and \(x_{\max}=m_{S,\max}^2/M_{F,\min}^2\). It reports approximate log-interpolated crossings of \(m_{S}(\mu)=\mu\) and \(m_{S,\max}(\mu)=M_{F,\min}(\mu)\). These are **diagnostics, not physical decoupling or pole-mass calculations**.

It also reports the implemented ODE derivatives \(dm_{S_i}^2/d\ln\mu\) at the endpoints of the two integration segments. Reproducing an endpoint verifies an implementation replay; it does not validate the physical beta functions, matching or threshold treatment.

## Full new project files

- `scripts/DiagnoseT3BCUVScalarHierarchy.py`
- `tests/python/test_t3_bc_uv_scalar_hierarchy.py`
- `docs/T3_BC_UV_RUNNING_MASS_HIERARCHY_BATCH14.md`

Nothing is replaced. The new output is only under `output/studies/general_scalar_extended/diagnostics/<new-run-id>/`.

## Windows PowerShell

Unpack the ZIP in `C:\thesistemp\T3neutrinomassmodelpipeline` (replace the placeholder ZIP path below), then run:

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "C:\path\to\T3_BC_UV_Mass_Hierarchy_Diagnostic_Batch14.zip" -DestinationPath . -Force
python -m pytest -q tests\python\test_t3_bc_uv_scalar_hierarchy.py
python -m scripts.DiagnoseT3BCUVScalarHierarchy --run-id uv_mass_hierarchy_v1
Get-Content -Encoding UTF8 output\studies\general_scalar_extended\diagnostics\uv_mass_hierarchy_v1\summary.md
```

Outputs:

```text
output/studies/general_scalar_extended/diagnostics/uv_mass_hierarchy_v1/
    diagnostic.json
    mass_trajectory_samples.csv
    summary.md
```

For visualisation elsewhere, the CSV can be loaded into any plotting tool. Keep the two branches distinct and avoid interpreting the 1 TeV endpoint mass parameter as the 1 TeV physical scalar mass.

If the runner reports an SHA mismatch, the saved fit refers to a different source revision. Do not bypass the check; inspect the changed source first. If it reports an endpoint mismatch, do not use the diagnosed threshold crossings until the discrepancy is resolved.

## Interpretation gate

- If the *restricted-zero* branch is also heavy, the mass growth is not primarily an effect of turning on the extra quartics at the boundary. Investigate the existing UV scalar-mass RG terms and hierarchy/threshold assumptions first.
- If the extra quartics materially change the masses, quantify this difference before interpreting fitted C5 results.
- If \(m_{S,\max}(\mu_F)\ge M_{F,\min}(\mu_F)\), there is no strict fermion-first ordering of the running mass parameters at the imposed 100 TeV matching scale. This is not a proof of physical threshold inversion; it flags the need for a consistent mass-dependent/inverted threshold analysis.
- In all cases, the current \(\chi^2\) values come from an exploratory **conditional one-scale** C5 calculation with 18 free real Yukawas. They are not physically certified fits.

## Batch 14 saved-fit wrapper correction (28 October 2026)

The first user run passed 21/21 original offline tests but stopped in
`fit_source()` with `KeyError: 'cfg'`, before any RGBeta integration.
The saved `*_exploratory_fitted_config.json` file is a raw configuration;
Batch 11's `seed_vector()` accepts a *record* containing that configuration
under the `cfg` key. The replacement diagnostic therefore calls
`b11.seed_vector({'cfg': cfg})` when validating the saved 18 real Yukawa
entries against `best_parameters`.

This replacement also supplies three regression tests exercising the
actual `fit_source()` input route for B(-1), C(-1) and a deliberately
inconsistent best-fit vector. All 24/24 diagnostic offline tests pass.
Nothing in the RG equations, masses, source hashes, fit or threshold
assumptions was changed.

From the project root, unpack the replacement ZIP and rerun:

```powershell
Expand-Archive -LiteralPath "C:\path\to\T3_BC_UV_Mass_Hierarchy_Diagnostic_Batch14_Fix.zip" -DestinationPath . -Force
python -m pytest -q tests\python\test_t3_bc_uv_scalar_hierarchy.py
python -m scripts.DiagnoseT3BCUVScalarHierarchy --run-id uv_mass_hierarchy_v1
Get-Content -Encoding UTF8 output\studies\general_scalar_extended\diagnostics\uv_mass_hierarchy_v1\summary.md
```

The failed pre-integration invocation should not have created the diagnostic
output directory. If that directory does exist, inspect it and choose a
fresh run ID; do not delete or overwrite a completed previous result.
