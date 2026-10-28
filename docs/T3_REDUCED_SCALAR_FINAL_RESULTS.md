# T3 reduced scalar-sector: final frozen-UV comparison and UV-running diagnostic

## Locked scope

**Primary numerical study:** 16 representations × two frozen-UV input scenarios, using the completed `reduced_frozen_overnight_v2` results. Compare individually fitted real Yukawa matrices (18 parameters) against the *newly optimised scenario-matched A(−2) reference*, without refitting the common reference in target models.

The two scenarios have topology inputs `lambdaT3=5e-9` and `0.5` at the fermion threshold `muF=1e5 GeV`. Class E uses half these benchmark coupling values. The scalar mass selection is retained at **0.5–5 TeV**, and the heavy-fermion hierarchy screen remains **rhoF <= 0.5**. Never silently enlarge either window.

**Supplementary study:** 16 historical UV-running inputs at `muUV=1e7 GeV`, running to `muF=1e5 GeV`; use saved formal nominal-scalar-EFT1 endpoint masses at `muS=1e3 GeV` only as mass/hierarchy diagnostics. These are *not* valid, physically rematched low-energy neutrino predictions when scalar masses grow far above 1 TeV or scalar/fermion decoupling order reverses. No scalar UV boundary retuning or further Yukawa optimisation is undertaken.

## Perturbative scalar approximation

This reporting stage uses the prior fixed-order relevance analysis in `docs/T3_REDUCED_SCALAR_POTENTIAL_CONCLUSION.md`. The original scalar masses, ordinary Hermitian self/portal quartics and T3 topology quartic are retained. All additional scalar quartic and quadratic couplings start at zero at the respective input scale. Additional quartics are evolved with available *candidate* one-loop RGEs, but radiatively generated B/C direct quartic insertions into one-loop C5 are omitted as formally higher fixed-loop order. The eight extra quadratic operators are not evolved, subject to the documented zero-boundary source audit limitations.

The existing model-specific C5 matching is hierarchical for twelve models, while four extreme-hypercharge models retain legacy one-scale matching. The physical Matchete/RGBeta scalar normalization and complete general sequential Wilson matching are not certified. NuFIT Gaussianised target chi2 <=5.86 is a **benchmark screening criterion**, not a model exclusion confidence level. Fitting 18 Yukawa parameters to five observables also does not establish parameter uniqueness or predictivity.

## Install (PowerShell)

Extract `T3_ReducedScalar_Final_Results.zip` into the local repository root, preserving the three project-relative paths:

- `scripts/BuildT3ReducedScalarFinalResults.py`
- `tests/python/test_t3_reduced_scalar_final_reporting.py`
- `docs/T3_REDUCED_SCALAR_FINAL_RESULTS.md`

No existing source or fitting results are replaced. The package only creates new reports in `output/studies/reduced_scalar_final_reports/<run-id>/`.

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath .\T3_ReducedScalar_Final_Results.zip -DestinationPath . -Force
python -m pytest -q tests\python\test_t3_reduced_scalar_final_reporting.py
python -u -m scripts.BuildT3ReducedScalarFinalResults --run-id final_reduced_v1
```

This reads existing local folders:

- `output/studies/reduced_scalar_sobol/reduced_frozen_overnight_v2/`
- `output/studies/reduced_scalar_sobol/reduced_48_v1/`

If you want only immediate table/endpoint figures without replaying UV trajectories:

```powershell
python -u -m scripts.BuildT3ReducedScalarFinalResults --run-id final_reduced_quick_v1 --skip-uv-replay
```

Do not use the same run ID twice; outputs are deliberately non-overwriting. The full command re-integrates **only** the existing one-loop UV RGEs for the 16 models, from `1e7` to `1e5` GeV, with their historical Yukawas and zero additional scalar boundaries. No new fits, scalar-threshold fitting or Matchete export is performed.

## Outputs

- `summary.md`: source-checked numerical result tables and scientific limitations.
- `manifest.json`: exact source hashes, acceptance counts, trajectory solver costs.
- `frozen_comparison.csv`: 32 individual fits and 32 shared A(−2) evaluations.
- `oscillation_observable_pulls.csv`: full fitted/shared five-observable comparison against the saved NuFIT reference (320 rows).
- `running_uv_formal_endpoints.csv`: 16 historical initial-input diagnostics, with fixed nominal-EFT1 scalar endpoints and rhoF.
- `frozen_chi2_comparison.png`: two-panel log-chi2 comparison, 16 models per scenario.
- `running_uv_formal_mass_endpoints.png`: original running-UV mass values relative to the retained 0.5–5 TeV window.
- `running_uv_trajectories.csv` and `running_uv_scalar_trajectories.png`: optional 16-model UV-only re-integrated scalar-mass trajectories from `1e7` to `1e5` GeV. These do not assert valid scalar EFT thresholds.

The frozen source campaign and all earlier studies are only read, not changed. All model-specific C5 results in the report are inherited from those saved runs; the reporter does not recompute neutrino masses.

## Thesis result to emphasize

The **32 independently fitted frozen-UV cases pass** the chosen five-observable acceptance screen, whereas **10/32** common A(−2) Yukawa evaluations pass (five per frozen scenario). Thus model-specific Yukawa adjustment substantially improves agreement with neutrino oscillation inputs. This does not imply that the other representations are experimentally excluded or that the fit is predictive.

The high-scale exploratory study shows that the chosen historical running-UV boundary drives scalar masses beyond the 0.5–5 TeV benchmark, invalidating the nominal scalar matching scale for those trajectories. The thesis does not claim that all possible high-scale T3 boundary conditions fail.


## UV trajectory replay: B/C Stage31 reused-quartic inventory fix

The original reporting script passed B(-1)/C(-1) through the historical bulk
UV-beta helper, which validates Stage31 using an empty reused-quartic list.
That helper was written for bulk fits excluding the special B/C models. The
saved B(-1) and C(-1) Stage31 exports actually contain six preexisting
additional quartics (plus one new quartic) per model. The report replay now
validates using the **complete** registry inventory; the reused quartic UV
inputs remain exactly zero, with their one-loop RG evolution retained.

This is a reporting-only fix. It does not change the RGEs, fit configurations,
fixed thresholds, coupling conventions, neutrino fits, or other production
scripts. A failed report run may have created a partial output directory, so
choose a **fresh run ID** rather than overwriting it:

```powershell
python -m pytest -q tests\python\test_t3_reduced_scalar_final_reporting.py
python -u -m scripts.BuildT3ReducedScalarFinalResults --run-id final_reduced_v2
```

If UV replay is unnecessary, `--skip-uv-replay` builds the frozen comparison,
observable pulls and saved formal running-endpoint plots without any new ODE
integration. Use a separate fresh run ID for that option.
