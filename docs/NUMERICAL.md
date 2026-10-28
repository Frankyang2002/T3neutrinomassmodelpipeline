# Numerical pipeline, fitting, diagnostics and reproducibility

**Scope:** ordinary T3 neutrino-mass studies, with emphasis on the completed reduced-scalar frozen-boundary campaign. This is the operational explanation; see [PHYSICS.md](PHYSICS.md) for the T3 EFT structure, [SCALAR_SECTOR.md](SCALAR_SECTOR.md) for the controlled scalar restriction and [RESULTS.md](RESULTS.md) for the numerical findings. Historical study instructions remain provenance, not competing current specifications.

**Source status (29 October 2026):** verified against the read-only repository commit `f21d15d2d38f4f313a0f2214bc52f90e1469dce2`, including `scripts/RunT3ReducedScalarSobolComparison.py`, `scripts/BuildT3ReducedScalarFinalResults.py`, `scripts/BuildReducedScalarExistingInteractiveReports.py`, the established plotting modules and saved final reporting artifacts. The installed local `README.md` and new `docs/SCALAR_SECTOR.md` are uncommitted documentation changes and are **not** replaced here. The project user's Windows checkout is authoritative if its file contents differ from the inspected commit.

## 1. Responsibility boundaries

| Layer | Source / responsibility |
|---|---|
| Model and input identity | `common/`, `model/`, `studies/`, `configs/` |
| Invariant Lagrangian and symbolic matching | `Lagrangian/` and external Matchete/Wolfram stages |
| Group factors, native beta exports, EFT and Weinberg matching | `RGE/group_factors/`, `RGE/running/`, `RGE/matching/` |
| Numerical state, complex-vector packing, beta evaluation | `Numerical/core/` |
| Scale evolution and threshold adapters | `Numerical/running/` |
| Neutrino objective, scanning and local optimization | `Numerical/fitting/` plus study-specific `scripts/` |
| Physical masses, mixing and oscillation observables | `physics/` |
| Diagnostic calculation, visualization and orchestration | `Numerical/diagnostics/`, `Numerical/plotting/`, `Numerical/orchestration/`, `Reports/` |
| Verification | `tests/`, `validation/`, saved native-source manifests |

Examples of core implementations include `Numerical/core/State.py`, `StateVector.py`, `RGBetaEvaluator.py` and `BetaVector.py`, and running components `UVRunner.py`, `FermionThresholdBoundary.py`, `IntermediateScalarRunner.py`, `ScalarThresholdBoundary.py` and `SMWeinbergEvolution.py`. **The numerical layer evaluates the source-mapped formulae; it should not silently substitute an independent C5 matching expression.** The original `pipeline.py` is the production entry point; individual `scripts/` commands may implement narrower experimental or reporting workflows.

## 2. Numerical conventions and EFT stages

Use the logarithmic scale variable and one-loop beta convention

\[
 t=\ln\mu,\qquad \frac{dX}{dt}=\frac{\beta_X^{(1)}}{16\pi^2}.
\]

A beta export supplying the one-loop numerator must be divided by `16*pi**2` **once**, not once in the exporter and again in the solver. Complex couplings and complex flavour matrices are packed into real solver vectors; preserve field order, shapes, conjugations and the charged-lepton basis prescription. Reference ODE settings recorded in the numerical overview include `rtol=1e-8`, `atol=1e-11` and approximately 120 recorded scale points for some workflows. These are choices of numerical controls, **not** a convergence proof for every study. Other purpose-built runners have their own settings, which are saved or defined in their source.

The documented standard hierarchy is

\[
 \mu_{\rm UV}=10^7\;\mathrm{GeV}>\mu_F=10^5\;\mathrm{GeV}>
 \mu_S=10^3\;\mathrm{GeV}>\mu_{\rm low}=10^2\;\mathrm{GeV}.
\]

The fermion-first architecture consists of UV evolution if enabled; integration out of the heavy fermion; scalar-only EFT1 evolution; integration out of the BSM scalars; and SM+Weinberg running. The physical coefficient is formed from the **existing hard and direct** contributions in the matching convention,

\[
 C_5(\mu_S)=C_{5,\mathrm{hard}}+\Delta C_{5,\mathrm{direct}},
 \qquad m_\nu=-\frac{v^2}{2}\,C_5.
\]

The matched Majorana matrix is symmetrised and Takagi-diagonalised with the existing conventions. LLSS running is separately useful diagnostically; a one-loop correction to LLSS subsequently inserted into a one-loop C5 graph has a different fixed-order counting and is not a new leading C5 term. A frozen RGE segment still **retains the threshold matching**, unless the individual experiment explicitly says otherwise. Scalar-first EFT ordering and arbitrary physically mixed scalar mass-eigenstate thresholds are not established production paths.

## 3. Separate study families: do not mix their input semantics

| Family | Cases and boundary | Interpretation |
|---|---|---|
| Original fixed comparison | 16 representations × (`smallY`/`largeY`) × (`smallL`/`largeL`) = 64; nominal Y = `0.005/0.5`, lambda reference = `0.1/1.0` | Numerical completed-grid study, **not** 64 neutrino fits |
| Historical scalar-oriented inputs | `small_uv_running_5e-9`, `small_uv_frozen_5e-9`, `large_uv_frozen_0p5` | Existing input scenarios used by later scalar studies; do not equate UV-frozen with UV-run |
| Primary reduced-scalar frozen results | 16 models × two frozen scenarios = 32, at `muF=1e5 GeV` | Independently optimised Yukawas and scenario-matched shared A(-2) comparison |
| Exploratory reduced-scalar running UV | Historical seeds at `muUV=1e7 GeV`, evolved to `muF` | Scalar mass/hierarchy diagnostic; no successful general low-energy matching for too-heavy scalars |
| Extended general-scalar 48-case campaign | 16 × three scalar-oriented scenarios, additional quartic running where available | **Conditional exploratory general-quartic** calculation, not physically certified full general-scalar EFT |

For the two **primary frozen** cases, `lambdaT3` is `5e-9` or `0.5` at `muF` in classes A–D and `2.5e-9` or `0.25` in class E. The class-E one-half factor is a benchmark initialisation convention, not an extra amplitude multiplicity. UV beta evolution is **disabled** for these frozen cases; scalar-only EFT1 and final SM+Weinberg evolution remain in the calculation according to the study runner.

The reduced-potential restriction keeps the scalar masses, ordinary Hermitian self/portal quartics and the original T3 topology coupling. Additional symmetry-allowed scalar quartics and bilinears have zero boundary coefficients. Candidate one-loop RGBeta can regenerate additional quartics, which are evolved where implemented; **radiatively generated direct B/C extra quartics are not inserted into the fixed-order leading-one-loop C5**, because their inclusion there is formally beyond that order. The eight extra quadratic couplings are not dynamically evolved in the reduced fit campaign. Their zero-boundary argument and limitations are derived in [SCALAR_SECTOR.md](SCALAR_SECTOR.md). Nothing here certifies the full physical Matchete/RGBeta scalar-normalisation map, complete extra-quartic Wilson matching, finite general thresholds or boundedness of the scalar potential.

## 4. Reduced-scalar Yukawa fitting: what was actually optimised

The completed frozen source is `output/studies/reduced_scalar_sobol/reduced_frozen_overnight_v2/`. Its owner is `scripts/RunT3ReducedScalarSobolComparison.py`; the historical `reduced_48_v1` run also supplies **exploratory running-UV diagnostic** inputs to the final reporter.

- The search varies two **real** 3×3 Yukawa matrices, i.e. 18 real parameters in the saved campaign, not 18 measured degrees of freedom. For each model/scenario it checks the seed, source/input identities, physical screens, saved beta inventory and model-specific matching artifact.
- A scrambled Sobol initial search and bounded SciPy `least_squares` refinement are used. Invalid or mass-unsafe trials are classified explicitly rather than assigned fabricated chi-square values. `--max-nfev` constrains optimizer evaluations but is **not** a cap on internal physics calls made for finite-difference Jacobians.
- `--warm-start-run` admits a separately validated previously accepted Yukawa vector; `--a2-as-start` can use the **newly optimised scenario-specific A(-2)** fit as an extra candidate. For the locked overnight run the command used `--a2-source new-fit --a2-as-start --warm-start-run reduced_frozen_a2_newfit_v1 --sobol-points 64 --starts 3 --max-nfev 80` with the two frozen cases. This describes provenance; it is **not a request to rerun** an already completed optimisation.
- In the shared-reference comparison, the same A(-2) Yukawa entries are evaluated unchanged at another model's representation, scalar inputs and C5 source. The target model is **not** independently refitted in that column. A shared-point failure is not model exclusion.
- All fit records retain objective values, screening outcomes, numerical evaluation counts, fitting sources, hashes and the distinction between acceptable input and successful low-energy evaluation.

The objective is an approximate diagonal Gaussian five-observable fit to `configs/nufit_6_1_ic24_no.json`:

\[
 \chi^2_{\rm approx}=\sum_{a=1}^{5}
 \left(\frac{O_a^{\rm predicted}-O_a^{\rm target}}{\sigma_a}\right)^2,
 \qquad\chi^2_{\rm approx}\leq 5.86\quad\text{(campaign screening)}.
\]

This is **not** the official NuFIT multidimensional likelihood, and 18 fitted parameters against five observables do not imply unique or independent predictions. For central values and actual acceptance counts, see [RESULTS.md](RESULTS.md).

The report uses existing model-specific leading fixed-order C5 sources: **12** models have hierarchical hard+direct matching; **four extreme-hypercharge models** retain explicitly nonuniform legacy one-scale sources. A cross-model comparison is therefore conditional on this recorded source heterogeneity.

## 5. Results and exact provenance paths

The source of truth for the final frozen comparison and running-UV diagnostic is:

```text
output/studies/reduced_scalar_sobol/reduced_frozen_overnight_v2/
    comparison.json
    comparison.csv
    <case>/<model>/fit.json

output/studies/reduced_scalar_sobol/reduced_48_v1/
    <running-case>/<model>/fit.json     # saved historical input/diagnostic records

output/studies/reduced_scalar_final_reports/final_reduced_v2/
    summary.md
    manifest.json
    frozen_comparison.csv
    oscillation_observable_pulls.csv
    running_uv_formal_endpoints.csv
    frozen_chi2_comparison.png
    running_uv_formal_mass_endpoints.png
    running_uv_trajectories.csv
    running_uv_scalar_trajectories.png
```

The latter two UV-trajectory files are present when UV replay was enabled (as reported in the saved `final_reduced_v2/summary.md`). The `manifest.json` fingerprints the source files, including the fitted inputs, registry and script, and records the acceptance counts and replay solver evaluations. **Preserve this manifest, its referenced inputs and the underlying matched C5 and beta-export artifacts** when cleaning the repository. Relative paths in the manifest refer to the project root. `frozen_comparison.csv` holds 32 × 2 chi-square evaluations; `oscillation_observable_pulls.csv` has 320 model/scenario/mode/observable rows. The reporter reads the saved fits; it does **not** refit or rematch C5.

`output/studies/reduced_scalar_final_reports/final_reduced_v2/` has already been generated and checked into the inspected repository. Do not repeat this exact run ID: the report writer refuses existing destinations. The read-only display adapter for the lambda-study report family uses these saved CSV/PNG inputs, not a new physics calculation.

## 6. Interactive RGE and sensitivity diagnostics

The two established rendering functions are `Numerical.plotting.InteractiveModelComparison.write_dashboard()` for RGE HTML, and `Numerical.plotting.BuildFittedModelEnvelopes.build()` for **descriptive sensitivity ovals**, not statistical confidence regions. The adapter `scripts.BuildReducedScalarExistingInteractiveReports` prepares the missing diagnostics from the completed frozen input states and then reuses these renderers unchanged.

The original frozen Sobol campaign saves **low-energy endpoints**, not intermediate gamma scans or 41-point SM+Weinberg trajectories. Therefore the adapter performs **forward diagnostic evolution, without new fits or new symbolic matching**:

- For RGE curves: scalar EFT1 running contributes to the matched boundary at `muS=1e3 GeV`; the plotted physical `C5` trajectory is the **SM+Weinberg** segment `1e3 -> 1e2 GeV`. There is no invented C5 curve above the scalar threshold and no UV-running curve in a frozen-UV study.
- For sensitivity: `T1/T2/T3` transformations vary the specified input Yukawa/quartic combinations, rerun the scalar EFT1 for changed lambdaT3, and check every gamma=1 central neutrino result against its saved fitted/shared reference.
- A large-gamma EFT1 solver error (for example `Required step size is less than spacing between numbers`) is recorded as a **NumericalFailure for that point**, not interpolated, extrapolated, relabelled as a poor oscillation fit or used to invalidate gamma=1. A nonpositive running scalar mass threshold is also a numerical failure. Unexpected source/metadata mismatches still fail closed.
- The original oval renderer requires at least **three successful gamma points in each T1/T2/T3 direction including gamma=1**. Insufficiently supported models are omitted from that dashboard, with their failed points and reasons kept in JSON. This is not exclusion of the physical representation. The early small-coupling optimised dashboard had 14 plotted models, two omitted; the exact coverage must be checked in saved report metadata for other configurations.

Source diagnostics and rendered reports have **different owners**:

```text
output/studies/reduced_scalar_interactive_diagnostics/reduced_frozen_overnight_v2/
    <case>/<mode>/models/<model>/running_diagnostics.json
    <case>/<mode>/sensitivity/<model>/sensitivity_diagnostics.json
    <case>/<mode>/sensitivity/comparison_summary.json

Reports/output/display_reduced_scalar/reduced_frozen_overnight_v2/
    <case>/reused/<mode>/interactive_rge.html
    <case>/reused/<mode>/sensitivity_ovals.html
```

Two frozen cases and two modes (`optimised`, `shared_a_m2`) yield up to eight principal HTML reports, subject to availability and valid-curve coverage. A one-model pilot writes in a separate `pilots/` tree and **must not** be mistaken for a complete 16-model dashboard. Saving the underlying diagnostic JSON is essential; the HTML is a generated presentation.

The general-quartic 48-study atlas (`Reports/output/general_scalar_48_final/`) is a **different scientific family**. Its frozen 48-row CSV contains endpoint records, not complete RGE trajectories or covariance scans, so never reuse the reduced/fitted curves as if they were general-quartic trajectory predictions.

## 7. Commands: inspection, verification, optional regeneration

From `C:\thesistemp\T3neutrinomassmodelpipeline` with `.venv` active:

```powershell
# Read already-computed accepted results; no integration or fitting.
Get-Content -Encoding UTF8 .\output\studies\reduced_scalar_final_reports\final_reduced_v2\summary.md
Import-Csv .\output\studies\reduced_scalar_final_reports\final_reduced_v2\frozen_comparison.csv |
    Format-Table model,case,chi2_individual,chi2_shared_A_m2 -AutoSize

# Current Python regression baseline, using Windows-independent UTF-8 mode.
python -X utf8 -m pytest -q --tb=line

# Report-only checks; no new physics fit.
python -X utf8 -m pytest -q tests/python/test_t3_reduced_scalar_final_reporting.py
python -X utf8 -m pytest -q tests/python/test_reduced_scalar_existing_interactive_reports.py

# Inspect commands before launching optional physics calculations.
python -m scripts.RunT3ReducedScalarSobolComparison --help
python -m scripts.BuildT3ReducedScalarFinalResults --help
python -m scripts.BuildReducedScalarExistingInteractiveReports --help
```

If a **new** copy of the final static report is needed, use a fresh ID and omit unnecessary UV replay:

```powershell
python -u -m scripts.BuildT3ReducedScalarFinalResults --run-id final_reduced_rebuild_v1 --skip-uv-replay
```

This command **creates a new report directory** from saved source files; it does not alter the existing final report. Remove `--skip-uv-replay` only when the 16-model historical candidate UV ODE replay is specifically desired. Creating *new* sensitivity source diagnostics is costlier than rendering stored HTML; if needed, begin with a single-model pilot:

```powershell
python -u -m scripts.BuildReducedScalarExistingInteractiveReports --case small_uv_frozen_5e-9 --mode optimised --model A_-2 --kind rge
```

Full gamma scans require additional forward numerical evaluations. Avoid `--force` unless the saved source-hash checks require recomputation. `python pipeline.py --full` may launch Matchete/RGBeta/Wolfram; it is **not** necessary for a documentation cleanup or to rebuild figures from saved data. `python pipeline.py --smoke` is a separate optional external symbolic regression when that toolchain is installed.

For the earlier fixed-grid displays, the existing saved-data commands remain separate: `python -m Numerical.diagnostics.ExtractFullComparisonResults`, `python -m Numerical.plotting.BuildDisplayResults --output Reports/output/display`, and `python -m Numerical.orchestration.RunFittedSensitivityComparison --help`. Do **not** treat the old fixed-grid display builder as a reduced-scalar or general-quartic report writer.

## 8. Validation tiers and cleanup boundaries

1. **Python regression:** `2144 passed, 10 skipped, 1 warning` in the user's Windows checkout at `f21d15d2` with the scalar documentation change. This verifies Python-level contracts; ten skips are not passes. The known warning is a `SyntaxWarning` for an invalid escape sequence in a test docstring. A documented `python -X utf8` run avoids Windows `cp932` test I/O failures.
2. **Saved-source provenance:** verify JSON/CSV completeness and SHA-256 fingerprints against the manifest; check C5 source per model, zero-extra boundary, positive scalar masses, ODE success, and fit/observable records. SHA-256 validates file identity, not physical completeness.
3. **Native symbolic regression:** any claim about freshly generated Matchete/RGBeta requires independent Wolfram/Matchete verification. Python tests and saved candidate outputs do not certify a full physical Matchete/RGBeta scalar mapping.
4. **Numerical stability:** compare tolerances and integrations where needed; a successful ODE or tiny fitting residual is not a proof of convergence, vacuum stability or model viability. Never convert a failed gamma point into a smooth figure by interpolation.

**Retain** source matched coefficients, candidate native beta exports, fitted Yukawa/seed files, per-case `fit.json`, the final report manifest/CSV, reduced interactive diagnostics JSON, and regression tests. **Generated HTML/PNG/PDF** may be archived or regenerated only after their corresponding saved source and renderer dependencies are accounted for. Historical batch scripts may still be imported by current reporting scripts; audit imports **before** retiring them. Do not delete historical sources simply because a short master document now describes their conclusions.

### Historical-source consolidation note

This document incorporates the operational material previously dispersed through `INFO_NUMERICAL.md`, `INFO_IMPLEMENTATION.md`, `INFO_TESTS.md`, old display READMEs, `T3_REDUCED_SCALAR_SOBOL_COMPARISON.md`, `T3_REDUCED_SCALAR_FINAL_RESULTS.md`, `T3_REDUCED_SCALAR_EXISTING_INTERACTIVE_REPORTS.md` and report rebuilding notes. These old notes remain in the repository until reference and artifact dependency audits permit safe external archiving. The source-code and data pathnames above, not archived prose, define reproducibility.
