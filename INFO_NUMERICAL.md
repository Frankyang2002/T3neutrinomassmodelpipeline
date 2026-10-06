# T3 Numerical Pipeline — File-by-File Guide

This document explains the numerical code after the structural refactor.

Current structure:

```text
Numerical/
├── core/
├── running/
├── fitting/
├── diagnostics/
├── plotting/
└── orchestration/
```

The separation is:

```text
core           What is the numerical state and how are beta values evaluated?
running        How is one benchmark evolved through the thresholds?
fitting        How are benchmark parameters scanned/fitted?
diagnostics    What detailed physical information should be retained?
plotting       How are stored results presented?
orchestration  How does pipeline.py invoke and persist the numerical workflow?
```

No matching formula is invented in the numerical layer. The authoritative matched $C_5$ remains owned by `RGE/matching/`.

---

# 1. `Numerical/core/`

## `State.py`

**Problem solved:** the ODE code needs one canonical, validated representation of all renormalisable UV parameters.

**What it does:**

- defines `T3Representation`;
- defines SM numerical data;
- defines ordinary/shared-scalar UV states;
- validates dimensions, scales, finite values and representation-dependent optional quartics.

It stores data only; it does not integrate RGEs.

## `StateVector.py`

**Problem solved:** SciPy integrates real vectors, while the T3 state contains real scalars, complex scalars and complex matrices.

**What it does:**

- defines a deterministic block layout;
- packs each complex quantity into real/imaginary components;
- unpacks a solver vector back to the canonical state;
- keeps $\mu$ as trajectory metadata rather than an ODE coordinate.

## `RGBetaEvaluator.py`

**Problem solved:** RGBeta exports Mathematica `InputForm` strings, but the Python ODE needs numerical beta values at every state.

**What it does:** tokenises/parses the supported exported expression grammar and evaluates traces, matrix products, bars/conjugation and scalar expressions directly on a numerical state.

It interprets the repository `report_betas` convention.

## `BetaVector.py`

**Problem solved:** evaluated beta dictionaries must be converted into the exact real-vector layout used by `StateVector.py`.

**What it does:**

- validates beta names;
- validates real/complex dimensions;
- applies the $1/(16\pi^2)$ loop factor exactly once;
- packs derivatives into the ODE vector.

## `__init__.py`

Package marker only.

---

# 2. `Numerical/running/` — multi-threshold benchmark evolution

## `UVRunner.py`

**Problem solved:** evolve the complete renormalisable UV T3 theory numerically from the UV scale to the first threshold.

**What it does:**

- validates requested scales;
- builds the ODE RHS from RGBeta beta expressions;
- integrates in $t=\ln\mu$ with `solve_ivp`;
- reconstructs validated states at requested save points;
- preserves the full trajectory, not only the endpoint.

For Majorana cases it maintains the expected mass-matrix structure numerically.

## `FermionThresholdBoundary.py`

**Problem solved:** after integrating out $F$, the numerical renormalisable state must be projected from the UV theory to the scalar-only theory.

**What it does:**

- extracts heavy-fermion threshold diagnostics;
- removes `MF`, `y1`, `y2`;
- preserves remaining SM/scalar renormalisable couplings;
- builds the ordinary or shared-scalar intermediate state.

It intentionally does not perform higher-dimensional matching; that is owned by the Matchete/RGE pipeline.

## `IntermediateScalarState.py`

**Problem solved:** the theory after $F$ is removed has different dynamical variables from the UV state.

**What it does:** defines canonical ordinary/shared scalar-only intermediate states.

Higher-dimensional Wilson coefficients are not stored here; this class represents only the renormalisable part evolved by RGBeta.

## `IntermediateScalarStateVector.py`

**Problem solved:** the scalar-only state also needs a deterministic real solver layout.

**What it does:** packs/unpacks the intermediate state using the same complex-to-real convention as the UV state.

## `IntermediateScalarRGBetaEvaluator.py`

**Problem solved:** the generic RGBeta expression parser needs a stage-specific environment and metadata checks for the scalar-only EFT.

**What it does:** constructs the intermediate-state expression environment, validates the RGBeta payload against the representation/stage and returns derivatives in the correct vector layout.

## `IntermediateScalarRunner.py`

**Problem solved:** evolve the renormalisable scalar-only EFT from the fermion threshold to the scalar threshold.

**What it does:** integrates the scalar-only ODE, retains requested save scales and returns the complete intermediate trajectory.

## `ScalarThresholdBoundary.py`

**Problem solved:** after the remaining scalars are integrated out, only the SM renormalisable state remains, but the final EFT must also receive the matched $C_5$.

**What it does:**

- extracts running scalar masses for diagnostics;
- projects the intermediate state to the final SM boundary;
- combines that state with an externally supplied symmetric $C_5$;
- builds the initial conditions for SM+Weinberg evolution.

## `T3Trajectory.py`

**Problem solved:** a complete benchmark crosses several different numerical theories and state types.

**What it does:**

```text
UV state
    -> UV running
    -> F boundary
    -> scalar-only running
    -> scalar boundary
    -> final SM state
```

and optionally continues into SM+Weinberg running when a physical $C_5$ is supplied.

It does not construct the matched $C_5$ itself.

## `FinalC5TrajectoryAdapter.py`

**Problem solved:** the authoritative final-C5 evaluator requires threshold-scale masses and Yukawas in the correct heavy-fermion mass basis.

**What it does:**

- reads the UV state at the fermion threshold;
- reads scalar quantities at the scalar threshold;
- leaves already-positive diagonal heavy masses unchanged;
- for Majorana $M_F$, performs a Takagi rotation and rotates both Yukawas;
- for vectorlike $M_F$, performs an SVD and rotates the two heavy indices appropriately;
- constructs the numerical config expected by the established final-C5 evaluator;
- evaluates the established matching result.

This file is an **adapter**, not a new matching formula.

## `SMWeinbergEvolution.py`

**Problem solved:** below the scalar threshold the project needs direct numerical evolution of the final SM+Weinberg EFT.

**What it does:**

- validates final-EFT initial conditions;
- packs/unpacks SM couplings and complex symmetric $C_5$;
- calls the analytic beta model from `RGE/running/weinberg/WeinbergRGE.py`;
- integrates with `solve_ivp`;
- returns the endpoint state.

It delegates $C_5\to m_\nu$ physics to `physics/NeutrinoMass.py`.

## `SMWeinbergStage.py`

**Problem solved:** final symbolic matched-$C_5$ artifacts must be numerically substituted and run through the final EFT while preserving the established JSON output contract.

**What it does:**

- parses/evaluates final-C5 JSON or symbolic matrices;
- builds numerical substitutions from config/threshold data;
- invokes `SMWeinbergEvolution`;
- writes machine-readable final-stage outputs.

## `WeinbergTrajectory.py`

**Problem solved:** diagnostics and plots need the entire final SM+Weinberg trajectory rather than only the low-energy endpoint.

**What it does:** integrates/saves $C_5(\mu)$ and SM states at requested scales.

Physical interpretation of saved points belongs to `physics/NeutrinoTrajectory.py`.

## `__init__.py`

Package marker only.

---

# 3. `Numerical/fitting/`

## `OscillationFit.py`

**Problem solved:** a model prediction must be compared to an external oscillation-data target without hard-coding one NuFIT release into model code.

**What it does:**

- defines the canonical prediction vector;
- loads target central values/covariance/errors from JSON;
- handles ordering consistency;
- computes correlated residuals and $\chi^2$.

Experimental values are data, not source constants.

## `ParameterScan.py`

**Problem solved:** parameter exploration needs deterministic grid/random machinery separate from the expensive T3 point evaluator.

**What it does:**

- defines scan parameter/result objects;
- generates grid or seeded-random points;
- evaluates each point;
- records failed points;
- supports progress/checkpoint callbacks;
- selects best/accepted points;
- serialises scan results.

The T3-aware evaluator requires an external matched-C5 builder.

## `OscillationOptimizer.py`

**Problem solved:** brute-force scans waste expensive full-pipeline evaluations far from promising regions.

**What it does:** provides a local optimisation workflow using sensitivity information and reduced active parameter sets to improve an existing seed/best point.

It retains full evaluation history for diagnostics.

## `SobolBenchmarkSearch.py`

**Problem solved:** each model needs a reproducible initial UV benchmark when no compatible cached benchmark exists.

**What it does:**

- defines the default 19-real-parameter ordinary-T3 fit profile;
- generates bounded Sobol samples;
- evaluates the full authoritative pipeline;
- selects the lowest-$\chi^2$ seed;
- optionally performs local refinement;
- caches/reuses model-compatible benchmark results.

## `BenchmarkSensitivity.py`

**Problem solved:** low-energy RGE sensitivity should be separated from the trivial leading prefactor scaling of the one-loop mass formula.

**What it does:** defines three controlled scaling directions that preserve the leading UV product $\lambda_{T3}(y_1y_2+y_2y_1)$ and evaluates residual low-energy change.

## `ScanCLI.py`

**Problem solved:** scans need a reproducible configuration-file-driven command line instead of ad-hoc scripts.

**What it does:** loads a scan JSON, resolves paths, validates bindings/representation/scales, mutates bound parameters in the base state, constructs the T3 evaluator and executes/dry-runs the scan.

## `__init__.py`

Package marker only.

---

# 4. `Numerical/diagnostics/`

## `BestFitDiagnostics.py`

**Problem solved:** optimiser JSON stores the best fit, but thesis analysis needs underlying physical matrices and a retained running trajectory.

**What it does:** reruns the authoritative best point and records low-energy masses, $C_5$, $m_\nu$, PMNS magnitudes, threshold diagnostics and scale-dependent running quantities.

## `RunningDiagnostics.py`

**Problem solved:** the thesis needs one structured JSON containing the most relevant scale dependence across all EFT regions.

**What it does:** samples UV, intermediate and final trajectories and stores coupling norms, neutrino masses/splittings/mixing angles and $C_5$ information.

## `IntermediateWeinbergDiagnostics.py`

**Problem solved:** the direct intermediate LLSS $\to O_5$ contribution should be visualised as it accumulates between thresholds without changing the fixed-order matching prescription.

**What it does:** samples the same symbolic direct-running expression used by the authoritative final-C5 result at intermediate scales.

## `FinalC5ContributionDiagnostics.py`

**Problem solved:** the final result should expose how much comes from the hard threshold and how much from direct intermediate running.

**What it does:** evaluates

$$
C_5^{\rm final}=C_5^{\rm hard}+C_5^{\rm direct}
$$

numerically on a trajectory and verifies the decomposition.

## `NumericalConfigCheck.py`

**Problem solved:** malformed or representation-incompatible generated configs should be caught without running Matchete/RGBeta/ODEs.

**What it does:** quickly constructs UV and post-$F$ scalar states for all ordinary neutral-compatible models and checks saved/comparison configs structurally.

## `__init__.py`

Package marker only.

---

# 5. `Numerical/plotting/`

These modules are presentation-only. They must not recompute or modify the physics.

## `RunningResultFigures.py`

**Problem solved:** `RunningDiagnostics` JSON needs a focused publication/thesis figure set.

**What it does:** generates UV coupling, direct intermediate-Weinberg, final-C5, mass-splitting, mixing-angle and sensitivity figures plus a README.

## `PlotOptimizerResults.py`

**Problem solved:** optimiser convergence/history is hard to assess from raw JSON.

**What it does:** plots $\chi^2$ trajectories, best-so-far convergence, pulls, sensitivity ranking and active-parameter trajectories and writes a summary.

## `PlotScanResults.py`

**Problem solved:** parameter-scan JSON needs standard quality/diagnostic plots.

**What it does:** generates best-$\chi^2$ progress, distributions, pulls, observable ratios and parameter-vs-$\chi^2$ plots.

## `ThesisResultFigures.py`

**Problem solved:** the final fitted benchmark needs compact thesis-ready figures and tables.

**What it does:** reads `BestFitDiagnostics` JSON and generates neutrino-mass running, $C_5$ running, mass/PMNS heatmaps, mass spectrum and CSV/Markdown/LaTeX tables.

It explicitly uses a headless Matplotlib backend.

## `InteractiveModelComparison.py`

**Problem solved:** comparing all model trajectories is cumbersome with static plots alone.

**What it does:** collects stored `running_diagnostics.json` files and builds a standalone interactive HTML dashboard with EFT-region shading and model toggles.

## `BuildDisplayResults.py`

**Problem solved:** supervisor-facing results should be compact without deleting complete reproducibility output.

**What it does:** creates `Reports/output/display/`, copying selected figures/dashboards and writing compact numerical/benchmark summaries.

## `__init__.py`

Package marker only.

---

# 6. `Numerical/orchestration/`

## `NumericalConfig.py`

**Problem solved:** model-specific numerical benchmark configs and resolved benchmark state need persistent, atomic bookkeeping independent of `pipeline.py`.

**What it does:** loads/writes config JSON, prepares a model-specific numerical config and freezes a resolved benchmark config after a successful run.

## `PipelineNumericalResults.py`

**Problem solved:** the central pipeline needs one call that executes the complete numerical benchmark/diagnostics/figure workflow using artifacts from the current `RunRecord`.

**What it does:**

- recognises the integrated numerical config schema;
- validates representation consistency;
- constructs the UV numerical state;
- resolves/uses automatic benchmarks;
- runs the multi-threshold trajectory;
- evaluates final $C_5$;
- continues to low-energy neutrino observables;
- runs configured sensitivity/diagnostics;
- writes numerical JSON;
- generates the standard figure set;
- records artifacts back onto the run.

It orchestrates existing physics modules; it does not define new matching or beta equations.

## `__init__.py`

Package marker only.

---

# 7. Default numerical scales and benchmark flow

The current default integrated numerical config uses approximately

```text
UV scale              1e7 GeV
fermion threshold     1e5 GeV
scalar threshold      1e3 GeV
low scale             1e2 GeV
```

The ordinary benchmark path is

```text
config
  -> canonical UV state
  -> UV RGBeta numerical running
  -> F threshold
  -> scalar-only numerical running
  -> scalar threshold values
  -> authoritative matched C5
  -> SM+Weinberg numerical running
  -> m_nu and oscillation observables
  -> chi2 / diagnostics / figures
```

---

# 8. Important heavy-mass-basis convention

At the fermion threshold:

### Already diagonal positive $M_F$

No rotation is applied.

### Majorana $M_F$

The complex-symmetric mass matrix is Takagi diagonalised,

$$
U^T M_F U=D,
$$

and both heavy Yukawa indices are rotated with the Takagi matrix.

### Vectorlike $M_F$

The mass matrix is bi-unitarily diagonalised,

$$
U_L^\dagger M_F U_R=D,
$$

with the two Yukawa structures transformed on the corresponding left/right heavy indices.

This basis handling is numerical preparation for the established matching result; it is not an extra threshold correction.

---

# 9. Overall numerical subsystem overview

```text
RGE analytic/exported beta functions
        |
        v
Numerical/core/
    parse/evaluate beta functions and define states
        |
        v
Numerical/running/
    evolve one physical benchmark through the thresholds
        |
        +----> physics/
        |       interpret C5 as masses/observables
        |
        v
Numerical/fitting/
    vary benchmark parameters and calculate chi2
        |
        v
Numerical/diagnostics/
    retain detailed physics information
        |
        v
Numerical/plotting/
    present stored results only
```

`Numerical/orchestration/` connects that chain to `pipeline.py`.

The central rule is:

> Numerical modules may evaluate and integrate established physics, but they should not become a second independent source of matching formulae or RGE equations.
