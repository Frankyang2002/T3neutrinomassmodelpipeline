# Research Numerical Fits And Benchmarks

> Historical research record. This document combines **original source
> Markdown** under original-path headings. Local Markdown links to moved
> files are updated; exact original bytes are in the external archive.
> This is not a new certification
> of physics, a rewritten conclusion, or a new computational workflow.
> Original bytes and SHA-256 values are stored in the external archive.

## Contents

- [One authoritative fitted benchmark per T3 model](#source-consolidate-fits-and-figures-md) — `docs/CONSOLIDATE_FITS_AND_FIGURES.md`
- [Fermion-start reports via established plot builders](#source-fermion-start-reuse-existing-reports-md) — `docs/FERMION_START_REUSE_EXISTING_REPORTS.md`
- [T3-B-m1 flavour alignment: read-only saved-matrix diagnostic](#source-t3b-flavour-alignment-md) — `docs/T3B_FLAVOUR_ALIGNMENT.md`
- [T3-B-m1: UV-frozen fitted flavour comparison](#source-t3b-uv-frozen-flavour-comparison-md) — `docs/T3B_UV_FROZEN_FLAVOUR_COMPARISON.md`
- [B(-1)/C(-1): conditional scalar-normalisation and mixed-mass C5 audit (Batch 10)](#source-t3-bc-conditional-mixed-c5-convention-audit-md) — `docs/T3_BC_CONDITIONAL_MIXED_C5_CONVENTION_AUDIT.md`
- [Batch 11 — exploratory B(-1)/C(-1) general-scalar Yukawa refit](#source-t3-bc-general-scalar-exploratory-refit-batch11-md) — `docs/T3_BC_GENERAL_SCALAR_EXPLORATORY_REFIT_BATCH11.md`
- [Batch 12 — continuation of the exploratory B(-1)/C(-1) general-scalar fits](#source-t3-bc-general-scalar-exploratory-refit-continuation-batch12-md) — `docs/T3_BC_GENERAL_SCALAR_EXPLORATORY_REFIT_CONTINUATION_BATCH12.md`
- [Batch 18 — Six-fit B/C general-scalar validation and consolidation](#source-t3-bc-general-scalar-six-fit-validation-batch18-md) — `docs/T3_BC_GENERAL_SCALAR_SIX_FIT_VALIDATION_BATCH18.md`
- [Batch 17 — constrained B(-1)/C(-1) UV-running Sobol and local optimisation](#source-t3-bc-uv-hierarchy-constrained-fit-batch17-md) — `docs/T3_BC_UV_HIERARCHY_CONSTRAINED_FIT_BATCH17.md`
- [Batch 16 — UV-running B(-1)/C(-1) hierarchy-safe neutrino normalization probes](#source-t3-bc-uv-neutrino-normalization-batch16-md) — `docs/T3_BC_UV_NEUTRINO_NORMALIZATION_BATCH16.md`
- [Batch 14 — B(-1)/C(-1) UV-running mass hierarchy replay](#source-t3-bc-uv-running-mass-hierarchy-batch14-md) — `docs/T3_BC_UV_RUNNING_MASS_HIERARCHY_BATCH14.md`
- [Batch 15 — B/C UV-running scalar-mass feasibility, before any more fits](#source-t3-bc-uv-scalar-mass-feasibility-batch15-md) — `docs/T3_BC_UV_SCALAR_MASS_FEASIBILITY_BATCH15.md`
- [Batch 20 — B(+1) native coupled general-scalar EFT1 beta extension](#source-t3-b-plus1-coupled-eft1-batch20-md) — `docs/T3_B_PLUS1_COUPLED_EFT1_BATCH20.md`
- [Batch 19 — First non-B(-1)/C(-1) general-scalar C5 source probe: B(+1)](#source-t3-b-plus1-general-c5-probe-batch19-md) — `docs/T3_B_PLUS1_GENERAL_C5_PROBE_BATCH19.md`
- [Batch 21 — B(+1) three-scenario conditional general-scalar numerical fits](#source-t3-b-plus1-three-general-scalar-studies-batch21-md) — `docs/T3_B_PLUS1_THREE_GENERAL_SCALAR_STUDIES_BATCH21.md`
- [A(-2) UV-running scalar mass-boundary feasibility pilot](#source-t3-reduced-uv-scalar-boundary-pilot-md) — `docs/T3_REDUCED_UV_SCALAR_BOUNDARY_PILOT.md`

---

<a id="source-consolidate-fits-and-figures-md"></a>

## Original document: `docs/CONSOLIDATE_FITS_AND_FIGURES.md`

**Original SHA-256:** `394e26ee59b4ab0465d0eff257d3a001e038c62d3b91550d4d611401367e2116`

**Status:** Historical record; check the current source and main
documentation before executing commands reproduced below.

<details>
<summary>Expand original: One authoritative fitted benchmark per T3 model</summary>

# One authoritative fitted benchmark per T3 model

Source snapshot verified: GitHub main `166ed4c57b8c32eb0bb55e5b0cc9e30b67e7bb64` (2026-10-20). This package does **not** change GitHub or any physics input values.

## Project paths (full replacements / additions)

- `tools/ConsolidateFittedBenchmarks.py`: new, preview and checksum-verify 16 winning fits; optionally archive three historical stages.
- `Numerical/orchestration/studies/fitted_rge/RunFittedSensitivityComparison.py`: full replacement, authoritative `best_candidate` reads `fitted_benchmarks` and verifies SHA256 against the manifest; original `fit_stage` is retained for provenance.
- `tools/RebuildCurrentReports.py`: full replacement, sends plotted figures into `Reports/output/display/full_studies/`, uses saved numerical data, leaves original HTML destinations.
- `tools/RedirectStudyFigures.py`: new display-only `Figure.savefig` redirection for current report rebuilding.
- `tools/MoveStudyFigures.py`: new preview-first relocation of existing PNGs from `full/` to `display/full_studies/`.

## Procedure (PowerShell at repository root)

```
cd C:\thesistemp\T3neutrinomassmodelpipeline

# Preview selections BEFORE updating code. Requires the ZIP files already extracted.
python -m tools.ConsolidateFittedBenchmarks
python -m tools.ConsolidateFittedBenchmarks --apply

# Confirm 16 copied configs + one manifest:
(Get-ChildItem .\Reports\output\full\fitted_benchmarks\smallY_smallL -Directory).Count
python -m tools.RebuildCurrentReports --only displays

# Inspect current and historical figures before moving:
python -m tools.MoveStudyFigures
python -m tools.MoveStudyFigures --apply

# Rebuild plots under display/; this does not rerun any RGEs:
python -m tools.RebuildCurrentReports --only plots

# Only when the above succeeds and there are no consumer scripts relying on old stage paths:
python -m tools.ConsolidateFittedBenchmarks --apply --archive-old

git status --short
```

**Important:** do not archive old stages before running the initial copy and inspecting the manifest. The archive is under `.cleanup_quarantine/fitted_stages_YYYYMMDDTHHMMSS/`, outside `Reports/output/full/` but inside the project. Keep the quarantine until future numerical runs have been checked.

The consolidation keeps bytes unchanged for `feasibility_summary.json` and `fitted_candidate_config.json`. Historical scan diagnostics still contain old absolute fit paths; the updated fitted sensitivity reader explicitly accepts those old path strings when reusing unchanged saved scans. They are provenance, not a reason to rewrite JSON.

The plot move only handles PNGs already inside `Reports/output/full/`. Afterward, the current `tools.RebuildCurrentReports` uses `display/full_studies/<study>/<scenario>/...` for its newly created figures. The direct numerical `--run` commands in individual study modules may still emit their original plot paths inside `full/`; rerun `tools.MoveStudyFigures --apply` after such runs, or use the report rebuild.

`Reports/output/analytical/` TeX/PDF reports and all unrelated numerical outputs are untouched.

## Validation

The migration tool and PNG move were tested with 16 synthetic model directories across staged fits, with checksum-matching copies and archived original directories. Python files were syntax-checked. No full user-specific Windows checkout was available, and the end-to-end report rebuild was not run against your saved data in this environment.

</details>

---

<a id="source-fermion-start-reuse-existing-reports-md"></a>

## Original document: `docs/FERMION_START_REUSE_EXISTING_REPORTS.md`

**Original SHA-256:** `8d6416e37830bff8cbf91ad9ca2a53517826610f2f82668ecbcb12da2d49063a`

**Status:** Historical record; check the current source and main
documentation before executing commands reproduced below.

<details>
<summary>Expand original: Fermion-start reports via established plot builders</summary>

# Fermion-start reports via established plot builders

This study does NOT replace `InteractiveModelComparison.py` or `BuildFittedModelEnvelopes.py`.
The adapter `Numerical/orchestration/BuildFermionStartExistingReports.py` prepares saved data in their input structures and then invokes these exact pre-existing HTML renderers.

## Numerical convention

`011`: input T3 UV-state couplings are specified at the fermion threshold mu_F;
UV evolution is frozen, fermion matching enabled, intermediate scalar EFT RGEs enabled,
scalar matching enabled, and final SM+Weinberg RGEs enabled.
**The 'UV' interval is not displayed as a physical running curve.**

Cases `5e-9`, `5e-10`; class E has lambdaT3 halved in the fitted configuration.
Modes `optimised`, `shared_a_m2`. Neither fitting nor symbolic RGEs are rerun.

## Run

```powershell
# Verify a single fitted model and established HTML builders first
python -m Numerical.orchestration.BuildFermionStartExistingReports `
  --case 5e-9 --mode optimised --model T3-A-m2 --kind all

# Generate all cases, models, and modes (about 2*2*16*(3*11+1) numerical evaluations)
python -m Numerical.orchestration.BuildFermionStartExistingReports `
  --case all --mode all --kind all
```

Files are saved under:

- `output/studies/fermion_start_lambda_comparison/<case>/reports_existing/<mode>/models/<model>/running_diagnostics.json`
- `output/studies/fermion_start_lambda_comparison/<case>/reports_existing/<mode>/sensitivity/<model>/sensitivity_diagnostics.json`
- `Reports/output/display_fermion_start/<case>/reused/<mode>/interactive_rge.html`
- `Reports/output/display_fermion_start/<case>/reused/<mode>/sensitivity_ovals.html`

The existing interactive RGE renderer consumes saved final-Weinberg curves. Its UV quantities are deliberately unavailable and will not plot. Scalar-only EFT1 state-vector samples and their exact block layout are saved for later plotting, without manufacturing an effective Weinberg coefficient above mu_S.

The existing fitted envelope renderer uses saved gamma directions T1/T2/T3. Their definitions are preserved from `Numerical/fitting/BenchmarkSensitivity.py`, but each evaluation uses the fermion-start 011 physics prescription instead of the previous full-UV evaluator. Ovals are empirical, descriptive envelopes, **not confidence contours**. Its NuFIT/JUNO overlays retain the renderer's existing definitions.

Default sensitivity gamma range is 0.1 to 10 with 11 log points; change to e.g. `--gamma-min .001 --gamma-max 1000 --points 25 --force` if the wider original scan is desired. The default is a practical pilot, **not identical** to the former fitted-sensitivity gamma range. `--force` recomputes all selected saved data.

### Checks / limits

- Each RGE curve's low-energy endpoint is checked against the corresponding `comparison` evaluation.
- A gamma=1 centre is checked against the saved benchmark for every sensitivity direction.
- Failed individual gamma points are retained as failed, not interpolated.
- Existing plotting templates are reused unchanged.
- The full `InteractiveModelComparison` legacy intermediate *direct C5 continuation* is omitted because the fermion-start evaluator does not save that separate matching diagnostic. Its other available plotting selectors still work, and no curve is invented.
- The script is syntax-checked, but not physics-executed against the user's local symbolic artifacts. Run the one-model command before the all-model command.

</details>

---

<a id="source-t3b-flavour-alignment-md"></a>

## Original document: `docs/T3B_FLAVOUR_ALIGNMENT.md`

**Original SHA-256:** `236522d17b96bc61dd2e9e95a3c9026b2e4423ccbc4840c10573d5ea09d88a09`

**Status:** Historical record; check the current source and main
documentation before executing commands reproduced below.

<details>
<summary>Expand original: T3-B-m1 flavour alignment: read-only saved-matrix diagnostic</summary>

# T3-B-m1 flavour alignment: read-only saved-matrix diagnostic

Checked against GitHub `main` tree `9abe36f583e3cb2d7a3d1ce709032af7bb357a7f` (24 October 2026). This is a **new independent analysis script**, not a replacement for existing physics code. It reads JSON already produced by the repository and **never executes UV running, intermediate running, Weinberg running, matching, fitting, or symbolic programs**.

## Install / execute (PowerShell; project root)

Review any local changes before extracting. Extract `T3B_FLAVOUR_ALIGNMENT.zip` into your local repository root.

```powershell
Expand-Archive -LiteralPath .\T3B_FLAVOUR_ALIGNMENT.zip -DestinationPath . -Force
python -m py_compile Numerical\analysis\AnalyseT3BFlavourAlignment.py

# Recalculate alignment from the complete, *local* previously saved JSON:
python -m Numerical.analysis.AnalyseT3BFlavourAlignment `
  --coupled output\t3b_fitted_common_scale_threshold.json `
  --control output\t3b_fitted_fixed_order_control.json `
  --low-energy output\t3b_fitted_low_energy_comparison.json `
  --output output\studies\t3b_flavour_alignment\fitted_reference_local

Start-Process output\studies\t3b_flavour_alignment\fitted_reference_local\index.html
```

Uses `numpy` and `matplotlib` already used by the project. It outputs `summary.json`, `index.html`, `c12_c5_alignment.png`, and `mass_basis_delta.png` in the selected output directory.

The ZIP also contains a precomputed example under `output/studies/t3b_flavour_alignment/fitted_reference_saved/` based on **selected numeric fields transcribed from saved GitHub files**. The three `source_excerpts/*.json` files are explicitly excerpts, not copies of the complete original diagnostics. The included HTML/PNGs and `summary.json` were generated with the same script from these excerpts; the normal command above recomputes from the complete local originals.

## Benchmark identity / important limitation

The example uses `Numerical/benchmarks/fitted_reference/T3_dS1_2_dS2_2_dF_1_alpha_m1.json`, where `lambdaT3=0.1`, and the saved `output/t3b_fitted_common_scale_threshold.json` coupled study. Its result is `2.4925034125%` at the scalar threshold. This is **not** the later local optimised benchmark reported as `2.29%`, nor the shared comparison reported as `2.56%`. The original saved example had UV running in the upstream calculation. The diagnostic itself performs no running, but this does not mean the input data came from a frozen-UV pipeline. Do not conflate these.

For your **exact** optimised 2.29% result, pass the corresponding `--coupled`, `--control`, and `--low-energy` files produced from the same benchmark and matching scheme. It will refuse missing or inconsistent files. For a UV-frozen study, all three source files must come from that UV-frozen evaluation; changing a CLI label is not sufficient.

The coupled treatment remains exploratory and does not establish complete three-generation operator-space closure or complete higher-order matching.

## Calculations

At the same physical flavour basis and scalar-threshold scale `mu_S=1000 GeV`:

- `fixed = hard + direct`, checked against the separate saved fixed-order control.
- `delta = coupled - fixed` and `||delta||F / ||fixed||F`.
- For reference `A` and candidate `B`, best complex scalar `alpha=Tr(A^dag B)/Tr(A^dag A)`, then shape residual `||B-alpha A||F/||B||F`.
- `||delta-alpha fixed||F/||fixed||F` distinguishes aligned versus flavour-changing correction.

At `mu_low=100 GeV`:

- Uses independently saved, evolved low-energy `C5` matrices; **does not** evolve them.
- For the real-symmetric saved benchmark, the eigenvectors of fixed `C5` are valid for the CP-conserving Takagi mixing directions up to diagonal phases. It checks `max||U|-saved |PMNS||` and refuses a mismatch larger than `1e-6`.
- Computes `U.T @ (C5_coupled-C5_fixed) @ U`, its off-diagonal relative norm, and signed-eigenvalue first-order rotation estimates for diagnostics. Always compare with the saved exact mixing angle differences.
- Refuses significantly complex cases instead of silently applying a real-symmetric diagonalisation to the wrong problem. For general complex Yukawa/charged-lepton mixing, use the project's full Takagi and charged-lepton-basis routines.

## Example result (from saved fitted reference, not the new UV-frozen 0.5 study)

| Quantity | Calculated result |
| --- | ---: |
| `||delta C5(MS)|| / ||C5_fixed(MS)||` | 0.0249250341 |
| shape residual `C12(MS)` vs `C5_fixed(MS)` | 0.0093153660 |
| shape residual `C12(MS)` vs `C5_direct(MS)` | 0.0082999120 |
| parallel factor of `delta` relative to fixed `C5` | -0.0249178355 |
| perpendicular part `||delta_perp|| / ||C5_fixed||` | 0.0005990010 |
| mass basis `||offdiag(U.T delta C5 U)|| / ||C5_fixed||` | 0.0003092792 |
| `|PMNS|` reconstruction max deviation | 3.34e-10 |

Saved exact low-energy angle differences: `theta12 = +0.0278269 deg`, `theta13 = +0.0089040 deg`, `theta23 = -0.0031850 deg`.

The data support a dominant rescaling and a much smaller rotation-generating matrix component. This is an observed alignment in the saved exploratory example, **not** a proof of exact proportionality or a validated high-order prediction.

</details>

---

<a id="source-t3b-uv-frozen-flavour-comparison-md"></a>

## Original document: `docs/T3B_UV_FROZEN_FLAVOUR_COMPARISON.md`

**Original SHA-256:** `79149b25c9dd93f7c841f5e0a16c8125ea08def03c1a998b3fc3865be2548cf5`

**Status:** Historical record; check the current source and main
documentation before executing commands reproduced below.

<details>
<summary>Expand original: T3-B-m1: UV-frozen fitted flavour comparison</summary>

# T3-B-m1: UV-frozen fitted flavour comparison

## What this runs

`Numerical/analysis/CompareFermionStartT3BFlavour.py` is an **independent, read-only** diagnostic of two existing, individually Sobol-optimised configurations:

- `output/studies/fermion_start_lambda_comparison/5e-9/fitted/T3_dS1_2_dS2_2_dF_1_alpha_m1.json`
- `output/studies/fermion_start_lambda_comparison/0p5/fitted/T3_dS1_2_dS2_2_dF_1_alpha_m1.json`

It also requires the corresponding `fits/<model>.json` records with `status=Success` and `chi2_fit<5.86`. **It will not substitute a prepared/unfitted 0.5 seed, silently refit parameters, or guess missing matrices.** The case label `0p5` denotes `lambdaT3=0.5` for class B.

The workflow is:

1. Load the existing fitted config and model-specific symbolic JSON artifacts. Require matching model ID, case, coupling, accepted fit and explicit fermion-start marker.
2. Impose saved initial couplings directly at `mu_F` via the repository's `frozen_result(..., uv=True)`; assert `uv.nfev=0`. Run actual scalar-only EFT1 renormalisable beta functions from `mu_F` to `mu_S` and collect sampled background.
3. Evaluate fixed-order scalar-threshold `C5 = C5_hard + C5_direct` using the existing symbolic hierarchical coefficient; cross-check with `FinalC5TrajectoryEvaluator`.
4. Check the existing *unfitted* heavy-log normalization of the direct contribution; refuse unexplained discrepancies.
5. Independently call the repository's **opt-in exploratory** full-flavour two-operator `C12`–`C5` leading-log system with its common-scale matching-log boundary. It does not overwrite the fixed-order `C5`.
6. Run the same final SM+Weinberg evolution from `mu_S` down to `mu_low` for fixed-order and exploratory matrices.
7. Quantify `C12`/`C5` alignment, the correction perpendicular to fixed `C5`, Yukawa SVD misalignment, individual heavy-column `C12` contributions and the transformed neutrino mass correction `U.T @ delta(mnu) @ U` in the charged-lepton mass basis.

All matrix projections use complex Frobenius products. Takagi diagonalisation is performed using the repository's full 3-generation complex-matrix routine, **not** a real-only eigenvalue shortcut. Ten-percent Yukawa-column variation tests are explicitly *tree `C12` texture-only*, not new low-energy neutrino predictions.

## Local execution — PowerShell

Unzip the package **at the repository root**:

```powershell
Expand-Archive -LiteralPath .\T3B_UV_FROZEN_FLAVOUR_COMPARISON.zip -DestinationPath . -Force
```

Run the already-fitted `5e-9` reference only (does not require a high-scalar fit):

```powershell
python -m Numerical.analysis.CompareFermionStartT3BFlavour --cases 5e-9
Start-Process "output\studies\t3b_uv_frozen_flavour_alignment\index.html"
```

Then, when `0p5/fitted/` and `0p5/fits/` exist from the previous high-scalar Sobol workflow, compare both:

```powershell
python -m Numerical.analysis.CompareFermionStartT3BFlavour --cases 5e-9,0p5
Start-Process "output\studies\t3b_uv_frozen_flavour_alignment\index.html"
```

If your model-specific artifact locations differ, pass the same `--artifact-map PATH.json` used by the existing fermion-start study. Adjust EFT1 sample count with `--samples 121` (default 121).

If the 0.5 point has **not** been fitted, first run the previously provided high-scalar workflow (do **not** use an unfitted rescaled seed as a completed benchmark):

```powershell
python -m Numerical.orchestration.RunHighScalarHalfFermionStart preflight --check-artifacts
python -m Numerical.orchestration.RunHighScalarHalfFermionStart prepare
python -m Numerical.orchestration.RunHighScalarHalfFermionStart fit --resume --model T3-B-m1
```

`fit --model T3-B-m1` uses the existing local workflow; it may require its Sobol search and EFT runs and is separate from this **no-refit diagnostic**.

## Output

Files are created locally under `output/studies/t3b_uv_frozen_flavour_alignment/`:

- `index.html`, `comparison.json`, `comparison.csv`: status and side-by-side metrics.
- `5e-9/summary.json`, `5e-9/index.html`, plus diagnostic PNGs.
- `0p5/summary.json`, `0p5/index.html`, plus diagnostic PNGs **only if the fitted input exists and the evaluation passes**.
- `saved_5e-9_yukawa_only/`: bundled independently verified *algebraic* input benchmark/figures, **no RGEs and no coupled `C5`**.

### Key metrics

For each case, the complex projection is

`alpha = vdot(C5_fixed, C5_coupled-C5_fixed)/vdot(C5_fixed,C5_fixed)`;

`eps_perp = ||deltaC5-alpha*C5_fixed||_F/||C5_fixed||_F`.

The neutrino mass-basis check is

`eps_mass = ||offdiag(U_fixed.T deltaMnu U_fixed)||_F/||Mnu_fixed||_F`.

The physical charged-lepton rotation and Takagi phase conventions are applied **first**. For complex symmetric Majorana matrices the transpose in the last expression, *not* Hermitian conjugation, is intentional.

### Interpretation boundaries

- **Verified from repository source:** `RunFermionStartLambdaStudy` freezes UV with `nfev=0` while running EFT1 and Weinberg; `T3BExploratoryCoupledBenchmark` explicitly has incomplete full operator closure; `FinalC5ContributionDiagnostics` splits fixed-order matching into hard and direct.
- **Generated by the new script:** optional exploratory coupled correction, its mass eigenbasis and plots, generated locally *only after the above inputs are found*. These are fresh diagnostics, not historical output files.
- **Not established:** a validated next-to-leading-order physical neutrino prediction or a perturbativity proof for `lambdaT3=0.5`.
- The GitHub-saved `5e-9` fit is an *available* UV-frozen baseline; it should not be mislabeled as the user's otherwise unavailable historical "2.29%" run. To analyze that different run, the corresponding local config, full scalar-threshold matrices and coupled matching scheme are needed.
- Missing `0p5` fitted files appear as `MissingFittedInput` in the comparison rather than being silently substituted.

Source inspected on GitHub `main` at commit `9abe36f583e3cb2d7a3d1ce709032af7bb357a7f` (24 October 2026). No source changes were sent to GitHub.

</details>

---

<a id="source-t3-bc-conditional-mixed-c5-convention-audit-md"></a>

## Original document: `docs/T3_BC_CONDITIONAL_MIXED_C5_CONVENTION_AUDIT.md`

**Original SHA-256:** `2cd43fa969971e9f1e84ba95cc620e1010df5672d081731e68f55635bed3c508`

**Status:** Historical record; check the current source and main
documentation before executing commands reproduced below.

<details>
<summary>Expand original: B(-1)/C(-1): conditional scalar-normalisation and mixed-mass C5 audit (Batch 10)</summary>

# B(-1)/C(-1): conditional scalar-normalisation and mixed-mass C5 audit (Batch 10)

## Scope

**New independent report only.** Read saved Batch 9 candidate RGBeta UV/EFT1 trajectories for the three historical benchmark cases. Re-evaluate the existing one-scale Matchete *ordered* C5 expression with the actual 1 TeV scalar mass matrix, alongside two explicit **conditional** Matchete/RGBeta coupling identifications. All Yukawa and heavy-fermion matrices remain frozen at the saved historical fitted input. No RGE is rerun; no fit, Matchete matching, historical result, or Batch 9 result is changed.

The script uses **both** `pilot_q0p1_c5_unmixed` and optionally the prior `pilot_q0p1_mu0p01`, with their SHA-256 provenance. The outputs are **not** a certified physical general-scalar neutrino prediction or a sequential two-threshold EFT result.

## Source-verified normalisation issue

From the actual Stage 23 report `output/general_t3_stage22/stage23_quartic_prefactor_candidates.json`:

- For both B(-1) and C(-1), `020020` identifies `lambdaHHdagS2S2` with `lambdaT3ExtraQ020020I1` at tensor ratio `R/M = 1`.
- For both, `020200` identifies `lambdaHHdagS1barS1bar` with `lambdaT3ExtraQ020200I1` at the **corrected-covariant candidate** tensor ratio `R/M = 1`.
- RGBeta's `AddQuartic` factor `6` is cancelled by its real-field/tensor symmetrisation and complex-field convention; it is **not** an extra physical scalar-potential factor.
- Stage 26 `output/general_t3_stage22/stage26_conditional_quartic_action.json` checks the source Lagrangian signs: Matchete uses `+m O + h.c.`, RGBeta uses `-r O + h.c.` under a common field convention.

Thus the source-supported **conditional** map for *these two one-dimensional blocks* is

\[
  m_1 = -r_1,\qquad m_2=-r_2
\]

rather than the earlier Batch 9 identity assumption `m_i = r_i`. Stage 23 and Stage 26 explicitly say that a **physical component-basis map is not certified**. This script fails if their stated source signs, relevant block matrices, total inventories, or uncertified status change. It does not invent a map for the other five B/C additional quartics, whose RG effects are retained but whose direct C5 insertions are not present in the saved coefficient.

## Scalar mass basis

At the 1 TeV RG endpoint, diagonalise

\[
 A=\begin{pmatrix}m_{S1}^2&\mu_2^2\\(\mu_2^2)^*&m_{S2}^2\end{pmatrix},\quad
 \Phi=(S_1,\epsilon S_2^*),\quad \Phi=U(X,Y),\quad U^\dagger A U=\operatorname{diag}(M_X^2,M_Y^2).
\]

Reuse `Numerical.experimental.GeneralT3MixedScalarMassBasis.MixedDoubletMassBasis.build` and feed **its** complex `U`, positive eigenmasses, the saved historical 3×3 Yukawas, `MF`, and the two quartic values into the **existing** ordered `B_m1_final` / `C_m1_final` Matchete one-loop expression. Symmetrise as `C5 = A5 + A5.T` (no complex conjugation), with `hbar=1/(16π²)`.

Unlike Batch 9's original preview, the candidate `C5_joint` here uses the **actual nonzero RG-generated mixing**. The script also provides unmixed controlled projections to isolate mixing effects. This re-evaluation is **not a new matching derivation**, and success does not establish the Matchete native mixed-field `CheckLagrangian`, `HeavyMassBasis`, one-loop EFT1 scalar threshold matching, or the physical signed map.

## Exact telescoping decomposition

For each of the two coupling maps, define, with `Y1,Y2,MF` held at their historical fitted input:

- `C0`: zero-extra endpoint diagonal masses and lambdaT3; zero direct quartics and zero mixing;
- `Cd`: general-quartic endpoint diagonal masses and lambdaT3; zero direct quartics and **mixing projected to zero**;
- `Cm`: same general-quartic endpoint with its actual mixing; zero direct quartics;
- `Cq`: same general-quartic endpoint, nonzero mapped direct quartics, **mixing projected to zero**;
- `Cj`: same general-quartic endpoint, nonzero mapped direct quartics and its actual mixing.

Then *exactly as complex 3×3 matrices*:

\[
 C_j-C_0=(C_d-C_0)+(C_m-C_d)+(C_q-C_d)+(C_j-C_m-C_q+C_d).
\]

The terms are, respectively, diagonal-background/`lambdaT3` effect, mixed-quadratic effect at zero direct quartics, direct-quartic effect at unmixed masses, and quartic–mixing interference. The displayed values are **individual Frobenius norms divided by** `||C0||_F`. Norms are *not additive*, although the complex matrix decomposition is. The report checks numerical closure.

A separate calculation at `C0`'s diagonal masses and `lambdaT3`, with positive identity-mapped quartic insertion and zero mixing, must reproduce the saved Batch 9 preview **exactly within numerical tolerance**. This guards against changing the historical definition of that diagnostic.

An optional `--mixing-run-id pilot_q0p1_mu0p01` compares the **separately evolved** explicit-mixing input scenario to the quartics-only scenario at the same mapping sign; this includes any accompanying RG mass/quartic changes, so it is *not* a pure fixed-background derivative with respect to `mu2`.

The script also records *singular values* of `-v² C5 / 2` in eV, for illustrative one-scale comparison only. They are neither low-energy RG-evolved fitted neutrino masses nor independent oscillation predictions.

## Installation paths

All new files:

- `scripts/AnalyzeT3BCMatchingConventions.py`
- `tests/python/test_t3_bc_matching_conventions.py`
- `docs/T3_BC_CONDITIONAL_MIXED_C5_CONVENTION_AUDIT.md`

All new results go to
`output/studies/general_scalar_extended/<historical_case>/<new_run_id>/summary.{json,md}`.

The three case directories are `small_uv_running_5e-9`, `small_uv_frozen_5e-9` and `large_uv_frozen_0p5`. Existing run-ID directories are never overwritten. The original `output/studies/lambda_5e-9_refit/` and `output/studies/fermion_start_lambda_comparison/` stay read-only.

## PowerShell — repository root

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline

Expand-Archive `
  -LiteralPath "C:\path\to\T3_BC_Conditional_Mixed_C5_Conventions_Batch10.zip" `
  -DestinationPath . -Force

python -m pytest -q tests\python\test_t3_bc_matching_conventions.py

# Re-evaluate existing quartic-only candidate RG endpoint, with its induced mu2:
python -m scripts.AnalyzeT3BCMatchingConventions `
  --all --source-run-id pilot_q0p1_c5_unmixed `
  --run-id matching_conventions_v1

# Optional: use the existing nonzero-input mu2 scenario as well (new run ID):
python -m scripts.AnalyzeT3BCMatchingConventions `
  --all --source-run-id pilot_q0p1_c5_unmixed `
  --mixing-run-id pilot_q0p1_mu0p01 `
  --run-id matching_with_mu_v1

Get-ChildItem output\studies\general_scalar_extended\*\matching_conventions_v1\summary.md |
  ForEach-Object { Get-Content -Encoding UTF8 $_.FullName }
```

If a stage-23/stage-26 JSON, fitted input, beta file, or saved Matchete ordered coefficient has been regenerated, the script **fails closed** instead of substituting a different reference. In that case review source compatibility before rerunning Batch 9 with a *new* run ID. If `matching_conventions_v1` already exists, change `--run-id` to another unused value; do **not** delete prior results.

## Result interpretation

The current calculations can quantify **conditional sensitivity** to an explicit normalisation convention and a nonzero scalar mass-mixing matrix. They **cannot** certify a unique physical nonzero-general-scalar `C5`. The native Matchete mass-eigenstate Lagrangian and the Matchete/RGBeta physical field identification remain unverified, and a complete sequential EFT1 matching/evolution of general-scalar operators has not been performed. Do not use the resulting percentage changes as experimental uncertainties or fitted neutrino predictions.

</details>

---

<a id="source-t3-bc-general-scalar-exploratory-refit-batch11-md"></a>

## Original document: `docs/T3_BC_GENERAL_SCALAR_EXPLORATORY_REFIT_BATCH11.md`

**Original SHA-256:** `7a015c6b473adb07c5a47344506f7bb6b616160e5ca824d446ff27e93b2e0fe5`

**Status:** Historical record; check the current source and main
documentation before executing commands reproduced below.

<details>
<summary>Expand original: Batch 11 — exploratory B(-1)/C(-1) general-scalar Yukawa refit</summary>

# Batch 11 — exploratory B(-1)/C(-1) general-scalar Yukawa refit

## Scope

**New independent numerical studies only.** `scripts/RunT3GeneralScalarBCExploratoryFit.py` reuses the verified historical fitted inputs and the **existing** Sobol / sensitivity-ranked (5/8/12 active-parameter) Yukawa optimiser. The scalar extra-coupling boundary is **fixed**, not fitted. There are 18 real Yukawa fit variables (`y1` and `y2`, three by three each).

Each **trial** (including its starting historical Yukawa point) executes:

1. Deep-copy the historically fitted model configuration and replace only its 18 real Yukawa entries in memory.
2. Set each of seven candidate general-scalar quartics to `--quartic-fraction * abs(lambdaT3)` at the appropriate boundary. Set the additional holomorphic quadratic mixing to `--mu2-fraction * min(mS1Sq,mS2Sq)` GeV². Defaults are `0.1` and `0`, respectively. These are **illustrative choices**.
3. For `small_uv_running_5e-9`, integrate candidate **UV** RGBeta from 10⁷ to 10⁵ GeV, retaining Yukawa and heavy-fermion mass running. For `small_uv_frozen_5e-9` and `large_uv_frozen_0p5`, **do not** integrate UV: apply the original parameter values directly at 10⁵ GeV.
4. Diagonalize the fermion mass matrix at 10⁵ GeV with the repository's existing Majorana Takagi / vectorlike SVD routines, rotating `y1,y2` consistently. Project to the scalar-only EFT1 and evolve the candidate RGBeta quartics and mass mixing 10⁵ → 10³ GeV.
5. Diagonalize the actual 2×2 mixed-doublet quadratic mass matrix, and evaluate the **saved** Matchete **one-scale ordered** coefficient using the RG endpoint. Insert the **conditional** source-supported mapping `lambda_M=-lambda_RGBeta` for the two additional direct-C5 quartics, without claiming certified physical field identification.
6. Reuse the existing SM+Weinberg integrator for 10³ → 10² GeV, charged-lepton basis change, neutrino mass conversion, NuFIT 6.1 (five-observable) approximate χ².

This is a *conditional one-scale* benchmark, **not** newly matched general-potential sequential EFT1 C5, physically certified cross-package coupling normalisation, or an independent neutrino prediction. No production numerical module, Matchete source, or RGE source is modified.

## Output segregation

Run IDs are mandatory and must be new. Each run writes only to:

`output/studies/general_scalar_extended/<case>/<run-id>/`

- `B_m1_exploratory_fit.json` or `C_m1_exploratory_fit.json`: numerical fit diagnostics and source hashes.
- `B_m1_exploratory_fitted_config.json` or `C_m1_exploratory_fitted_config.json`: complete **new** fitted input copy with updated Yukawas and explicit exploratory provenance. Inherited old study metadata refers only to historical provenance.
- `summary.json`, `summary.md`: run results and scope limitations.

For project completeness there is a single full new script, a new test module and this document. Previous `output/studies/lambda_5e-9_refit/`, `output/studies/fermion_start_lambda_comparison/`, Stage 36, Batch 7/8, Batch 9, and Batch 10 files are never opened for writing.

## Local PowerShell

Extract the ZIP at `C:\thesistemp\T3neutrinomassmodelpipeline`. Run from that directory in the project `.venv`:

```powershell
python -m pytest -q tests\python\test_t3_general_scalar_bc_exploratory_fit.py

# Verify exact historical seeds before optimisation (two independent runs).
python -m scripts.RunT3GeneralScalarBCExploratoryFit `
  --case small_uv_frozen_5e-9 --model B_-1 --run-id b_small_seed_v1
python -m scripts.RunT3GeneralScalarBCExploratoryFit `
  --case large_uv_frozen_0p5 --model C_-1 --run-id c_large_seed_v1

# After both probes succeed: separate local refits. Default: 8 Sobol points
# then the existing staged 5/8/12-parameter refinement if necessary.
python -m scripts.RunT3GeneralScalarBCExploratoryFit `
  --case small_uv_frozen_5e-9 --model B_-1 --fit `
  --sobol-points 8 --run-id b_small_fit_v1
python -m scripts.RunT3GeneralScalarBCExploratoryFit `
  --case large_uv_frozen_0p5 --model C_-1 --fit `
  --sobol-points 8 --run-id c_large_fit_v1

Get-Content -Encoding UTF8 output\studies\general_scalar_extended\small_uv_frozen_5e-9\b_small_fit_v1\summary.md
Get-Content -Encoding UTF8 output\studies\general_scalar_extended\large_uv_frozen_0p5\c_large_fit_v1\summary.md
```

`--fit` is opt-in; without it, the command evaluates the historical Yukawa seed in the *new* scalar scenario, still performing RG and low-energy evolution. Increase `--sobol-points` to 32 or 64 for a broader, more expensive search after the pilot. `--bound` (default `3.0`) is an absolute per-entry Yukawa bound; the B/C small-scalar historical Yukawa entries are O(1), so do **not** use the high-scalar fitted-study default `0.01` for the small-scalar cases.

Optional explicit scalar-mixing input: `--mu2-fraction 0.01` **and a new run ID**. The one-scale matrix handles its nonzero eigenbasis but does not certify a physical matching result.

## Validation status

Offline tests exercise the immutable output path, parameter shape, copied trial state, original seeded inputs, stage-interface calls, and a complex mixed-mass one-scale coefficient under the negative signed map. They do **not** substitute for full numerical integration against the user's local UV/EFT1 RGBeta beta files, source coefficient, or Matchete output. The first local two-seed probes are therefore a required numerical gate before interpreting any refit.

</details>

---

<a id="source-t3-bc-general-scalar-exploratory-refit-continuation-batch12-md"></a>

## Original document: `docs/T3_BC_GENERAL_SCALAR_EXPLORATORY_REFIT_CONTINUATION_BATCH12.md`

**Original SHA-256:** `324de211047cd86b4ee3ad4caa71a3ce8731701945da5937088707e630e4b511`

**Status:** Historical record; check the current source and main
documentation before executing commands reproduced below.

<details>
<summary>Expand original: Batch 12 — continuation of the exploratory B(-1)/C(-1) general-scalar fits</summary>

# Batch 12 — continuation of the exploratory B(-1)/C(-1) general-scalar fits

## Purpose

Continue the **already completed** Batch 11 Yukawa optimisations rather than repeat a fit from the historical restricted-potential seeds. The source fits provided by the user are:

| Case | Model | Previous run | Previous chi-squared | General scalar scenario |
|---|---|---|---:|---|
| `small_uv_frozen_5e-9` | `B_-1` | `b_small_fit_v1` | 118.57915264725226 | Seven additional quartics at 0.1\|lambdaT3\|, zero quadratic input |
| `large_uv_frozen_0p5` | `C_-1` | `c_large_fit_v1` | 116.30993286412242 | Same boundary rule, with lambdaT3=0.5 |

Both prior fits are marked `ExploratoryPoorFit`. The goal is chi-squared <= 5.86. The Batch 12 continuation is **not guaranteed to attain it**.

The historical fitted results, previous Batch 9–11 runs, Matchete symbolic coefficient, scalar parameters, and original fit target remain unchanged. All outputs are **new-only**.

## Actual numerical path

Every Yukawa trial uses the **existing Batch 11** `CandidateEvaluator` without altering its physics:

1. Set all 18 *real* Yukawa parameters on a fresh copy of the original historical configuration.
2. Integrate candidate UV RGBeta (when the original case enables UV running) and scalar-only EFT1 RGBeta with the same additional-quartic/quadratic boundary as the source Batch 11 run.
3. Determine heavy-fermion mass-basis Yukawas at `muF=1e5 GeV`, including Takagi/SVD when necessary.
4. Insert endpoint scalar mass mixing, diagonal scalar masses and two direct quartics into the saved **one-scale** Matchete coefficient with the Stage 23/26 **conditional** map `lambda_M = -lambda_RGBeta`.
5. Run SM+Weinberg from `muS=1e3 GeV` to `mu_low=100 GeV`; calculate the unchanged five-observable NuFIT chi-squared.

The optimiser also reuses the existing `Numerical.fitting.SobolBenchmarkSearch.refine_sobol_candidates` implementation, with **all 18 Yukawa coordinates** active simultaneously (instead of the previous staged 5/8/12). This is a numerical optimisation adjustment, not a new matching calculation.

### Important numerical change — local variable scaling

The previous best B(-1) Yukawas reach 1.54 while the previous best C(-1) Yukawas reach 3.09e-4. Using a common optimiser interval `[-3,3]` makes the dimensionless derivative steps excessively large for C(-1). Batch 12 uses individual finite local intervals, centred on **each saved best parameter**:

```
half_width_j = max(relative_window * abs(y_j),
                   floor_ratio * max_k(abs(y_k)), absolute_floor)
lower_j = max(-bound, y_j - half_width_j)
upper_j = min(+bound, y_j + half_width_j)
```

Defaults: `relative_window=0.5`, `floor_ratio=0.02`, `absolute_floor=1e-8`, `bound=3`. It computes the exact source seed once, then launches a single all-18 local least-squares stage with `max_nfev=60`. Optional additional local Sobol search is disabled by default (`--sobol-points 0`) to avoid unnecessary integrations. Set `--sobol-points 8` for a separate broader study only if needed.

## Install

Extract the ZIP at the project root:

`C:\thesistemp\T3neutrinomassmodelpipeline`

The ZIP has **only these three new files**:

- `scripts/ContinueT3GeneralScalarBCExploratoryFit.py`
- `tests/python/test_t3_general_scalar_bc_continuation.py`
- `docs/T3_BC_GENERAL_SCALAR_EXPLORATORY_REFIT_CONTINUATION_BATCH12.md`

The original Batch 11 runner remains untouched.

## PowerShell

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "C:\path\to\T3_BC_General_Scalar_Fit_Continuation_Batch12.zip" -DestinationPath . -Force
python -m pytest -q tests\python\test_t3_general_scalar_bc_continuation.py

# First continue B(-1), small scalar, UV frozen:
python -m scripts.ContinueT3GeneralScalarBCExploratoryFit `
  --case small_uv_frozen_5e-9 `
  --model B_-1 `
  --previous-run-id b_small_fit_v1 `
  --run-id b_small_fit_v2 `
  --max-nfev 60

# Then continue C(-1), large scalar, UV frozen:
python -m scripts.ContinueT3GeneralScalarBCExploratoryFit `
  --case large_uv_frozen_0p5 `
  --model C_-1 `
  --previous-run-id c_large_fit_v1 `
  --run-id c_large_fit_v2 `
  --max-nfev 60
```

Optional broader local search if a continuation fails the fit screen: use **another** new run ID, with `--sobol-points 8`; do not overwrite `v2`.

Read the two summaries:

```powershell
Get-Content output\studies\general_scalar_extended\small_uv_frozen_5e-9\b_small_fit_v2\summary.md
Get-Content output\studies\general_scalar_extended\large_uv_frozen_0p5\c_large_fit_v2\summary.md
```

### New outputs

Each run writes four files only to its new study directory:

- `<B_m1/C_m1>_exploratory_fit.json` — full results, parameter vectors, optimisation trace, source hashes and diagnostic status;
- `<B_m1/C_m1>_exploratory_fitted_config.json` — complete distinct fitted configuration, provenance and 18 updated real Yukawas;
- `summary.json`;
- `summary.md`.

The command rejects an existing destination directory, a source run not marked as a completed fit, a mismatch between the saved best 18-vector and previous fitted config, a changed source physics-input hash, a different UV regime or scalar boundary, or a nonreproducible previous fit chi-squared. It does not edit the source file.

## Verification and interpretation

- Offline unit tests: 16 passed in the build environment; includes mocked full 18-parameter continuation, immutable source checks, prior-run consistency, adaptive scaling, and parsing of both user-supplied saved fit JSONs (the latter test auto-skips on Windows where that uploaded analysis ZIP is not present).
- Both actual source runs were checked directly from the user-uploaded ZIP: 18 parameters, expected statuses, and fitted-config values matched the saved best point exactly.
- **Not run here:** full RGBeta/Matchete/SM-Weinberg numerical continuation on the user's Windows repository. The local runs determine whether the fit screen is reached.
- All results remain **exploratory, not physically certified**. The Matchete/RGBeta common component convention, new-operator normalisations and full sequential mixed-scalar matching remain unresolved; the fitted neutrino observables are not independent predictions.

</details>

---

<a id="source-t3-bc-general-scalar-six-fit-validation-batch18-md"></a>

## Original document: `docs/T3_BC_GENERAL_SCALAR_SIX_FIT_VALIDATION_BATCH18.md`

**Original SHA-256:** `6e945f67442e145a57ce6efd38eb1742f537bbfd642279b7caeae79454d7b6c3`

**Status:** Historical record; check the current source and main
documentation before executing commands reproduced below.

<details>
<summary>Expand original: Batch 18 — Six-fit B/C general-scalar validation and consolidation</summary>

# Batch 18 — Six-fit B/C general-scalar validation and consolidation

## Purpose

Consolidate the **six accepted conditional exploratory** B(-1)/C(-1) neutrino fits from Batches 11–13 and 17, check their provenance and saved observables, then independently replay their existing **candidate RGBeta** running and conditional **one-scale Matchete C5**. The intent is a numerical-physics benchmark audit, **not** a new EFT matching derivation.

The six input directories are **read only**:

| Study | Model | Input run ID | Format |
|---|---|---|---|
| `small_uv_running_5e-9` | `B_-1` | `b_uv_constrained_pilot_v1` | Batch 17 `fit.json` |
| `small_uv_running_5e-9` | `C_-1` | `c_uv_constrained_pilot_v1` | Batch 17 `fit.json` |
| `small_uv_frozen_5e-9` | `B_-1` | `b_small_fit_v2` | Batch 12 continuation |
| `small_uv_frozen_5e-9` | `C_-1` | `c_small_fit_v1` | Batch 11 fit |
| `large_uv_frozen_0p5` | `B_-1` | `b_large_fit_v1` | Batch 11 fit |
| `large_uv_frozen_0p5` | `C_-1` | `c_large_fit_v2` | Batch 12 continuation |

The UV-running solutions fit `lambdaT3` alongside the 18 real Yukawas; they are not fixed-`5e-9` repetitions of the original UV-running study. Frozen cases retain their original fixed scalar boundary.

## What is checked

**Saved-input checks** (both modes):

- Six exact identities, status/acceptance flags, `chi2 <= 5.86`.
- Fitted Yukawas and UV `lambdaT3` agree between the saved fit and its complete fitted configuration.
- Source input SHA-256 fingerprints match the files currently on disk. Any mismatch stops the run.
- Saved predictions reproduce the five-variable NuFIT chi-square with the currently saved `configs/nufit_6_1_ic24_no.json` target.
- Saved complex symmetric C5 matrices at 1 TeV and 100 GeV are finite and nonzero.
- Scalar 2x2 quadratic mass eigenvalues are positive; UV-running cases satisfy `rhoF <= 0.5` and 0.5–5 TeV at 1 TeV.
- Frozen-study `rhoF` estimates are formed from the **100 TeV scalar starting masses** and 100 TeV heavy masses, never mixing values evaluated at different scales.

**Default full replay** additionally:

- Recreate the fitted UV/EFT1 boundary using the *saved best-fit Yukawas*, not fresh seeds.
- UV-running B/C: candidate RGBeta 10^7 -> 10^5 GeV; frozen B/C: no UV RGE; all six: candidate EFT1 10^5 -> 10^3 GeV.
- Independently evaluate the same conditional saved one-scale C5 and run SM+Weinberg 10^3 -> 10^2 GeV, cross-checking the saved five observables and chi-square.
- Sample running scalar 2x2 eigenmasses, heavy-fermion ratios, and maximum raw scalar quartic-coordinate magnitude. For UV, also sample largest singular value of the running `y1` and `y2` matrices.
- Flag coordinate-level perturbativity when `max |lambda| > 4*pi` or `max singular(y1,y2) > sqrt(4*pi)`. **These are numerical screening conventions**, not a complete physical perturbative-unitarity or potential-boundedness analysis.

No production RGE, Matchete coefficient, fitted input, or historical result is edited. The runner fails **before writing output** when any saved-file check or replay comparison fails.

## Project paths (new files)

- `scripts/ValidateT3BCGeneralScalarSixFits.py`
- `tests/python/test_t3_bc_general_scalar_six_validation.py`
- `docs/T3_BC_GENERAL_SCALAR_SIX_FIT_VALIDATION_BATCH18.md`

## PowerShell commands

From `C:\thesistemp\T3neutrinomassmodelpipeline`, with `.venv` active:

```powershell
Expand-Archive -LiteralPath "C:\path\to\T3_BC_General_Scalar_Six_Fit_Validation_Batch18.zip" -DestinationPath . -Force
python -m pytest -q tests\python\test_t3_bc_general_scalar_six_validation.py

# Full RG and C5 replay, 10 saved samples per running interval, new-only output:
python -m scripts.ValidateT3BCGeneralScalarSixFits --run-id bc_six_validation_v1

Get-Content -Encoding UTF8 output\studies\general_scalar_extended\validation\bc_six_validation_v1\summary.md
```

For a fast *saved-file-only* source/fit check (not full trajectory validation), with a different output ID:

```powershell
python -m scripts.ValidateT3BCGeneralScalarSixFits --run-id bc_six_saved_v1 --saved-only
```

The outputs for the full run are:

- `output/studies/general_scalar_extended/validation/bc_six_validation_v1/report.json`
- `output/studies/general_scalar_extended/validation/bc_six_validation_v1/summary.md`
- `output/studies/general_scalar_extended/validation/bc_six_validation_v1/trajectory_samples.csv`

Use a **new** `--run-id` if an output directory already exists. The report never claims a full 16-model general-scalar fit. No new numerical fit is performed.

## Physical limits

The candidate signed Matchete/RGBeta quartic identification is **not** certified; the saved C5 is a one-scale insertion, not complete sequential 100 TeV fermion and 1 TeV scalar Wilson-operator matching. Scalar masses from 2x2 RG quadratic matrices are not pole masses. Sampling positive masses does not guarantee vacuum stability for the complete general potential. The 100 TeV/1 TeV common-scale approximation persists, particularly when a running scalar mass is several TeV. Extremely small `chi2` values are existence checks with many adjustable Yukawa parameters, not evidence of unique predictions.

</details>

---

<a id="source-t3-bc-uv-hierarchy-constrained-fit-batch17-md"></a>

## Original document: `docs/T3_BC_UV_HIERARCHY_CONSTRAINED_FIT_BATCH17.md`

**Original SHA-256:** `35049f1b085ffa9de9196bfcb74fb450f7b0cb7235824e54df5dd5da3cfb5f7c`

**Status:** Historical record; check the current source and main
documentation before executing commands reproduced below.

<details>
<summary>Expand original: Batch 17 — constrained B(-1)/C(-1) UV-running Sobol and local optimisation</summary>

# Batch 17 — constrained B(-1)/C(-1) UV-running Sobol and local optimisation

## What is actually implemented

A new **conditional exploratory** fitter reuses all unchanged Batch 11/14/16
interfaces and saved, input-fingerprint-verified historical and exploratory
source configurations. Each trial updates **all 18 real Yukawa coordinates**
and **one additional real coordinate** `log10_lambdaT3` at 10^7 GeV. The
seven candidate additional scalar quartics start at `0.1*abs(lambdaT3)`;
`mu2_fraction=0` at the UV boundary and its RG-generated value is retained.

On every trial:

1. Run the actual saved candidate RGBeta UV state (10^7 -> 10^5 GeV).
2. Require positive mixed-scalar running-mass eigenvalues and
   `max(mScalar(muF))/min(MF(muF)) <= 0.5`.
3. For a valid UV hierarchy only, run the existing EFT1 candidate scalar RG
   state (10^5 -> 10^3 GeV) and require both scalar masses in `[0.5,5]` TeV.
4. Only if both mass gates pass, use the saved one-scale Matchete C5 with
   **conditional** Stage23/26 quartic mapping, heavy mass-basis rotation,
   and SM+Weinberg running to 100 GeV to compute the five NuFIT observables.
5. Recover the whitened five-component residual from the exact experimental
   target covariance; independently verify its squared norm against the
   returned C5 chi-square.

No scalar-first or sequential EFT matching is generated. These are mass
**running-parameter** inequalities, not a pole-mass decoupling certificate.

## Search method

- Inputs are the **two existing Yukawa textures** (historical and previous
  poor fit) with Batch15-successful common scales: B `0.005,0.01,0.015,0.02`;
  C `0.005,0.01,0.015`. Each is only a seed, not the final fit.
- Use the existing SciPy scrambled-Sobol generator from
  `Numerical.fitting.SobolBenchmarkSearch`, sampling locally around **all 18**
  Yukawa entries and a two-decade window for `log10(lambdaT3)` (clipped to
  global bounds). The `lambdaT3=5e-9/f^2` relation sets only the seed centre.
- Keep screened candidates, rank by **actual recomputed** chi-square, and
  call `scipy.optimize.least_squares` in dimensionless unit coordinates for
  **all 19 variables** from up to two best seeds.
- Ineligible candidates have **no assigned NuFIT chi-square**. The local
  solver instead sees an explicit nonzero penalty residual for a failed mass
  gate. The penalty is only an optimisation aid and is never reported as
  an oscillation statistic. Numerical exceptions likewise cannot become
  accepted candidates.
- Screened acceptance requires both mass conditions and `chi2 <= 5.86`.
  This procedure does **not** guarantee that a screened fit will be found.

## Project files

- `scripts/FitT3BCUVHierarchyConstrained.py` — **new, complete** runner.
- `tests/python/test_t3_bc_uv_hierarchy_constrained_fit.py` — unit and
  mocked full-run tests.
- `docs/T3_BC_UV_HIERARCHY_CONSTRAINED_FIT_BATCH17.md` — this document.

Other project modules, including Batch 16 and all production RG equations,
are unchanged.

## PowerShell commands

Run in `C:\thesistemp\T3neutrinomassmodelpipeline` after extracting the ZIP
into the project root:

```powershell
python -m pytest -q tests\python\test_t3_bc_uv_hierarchy_constrained_fit.py

# First pilot: one texture, selected mass-safe seed factors, modest budget.
python -m scripts.FitT3BCUVHierarchyConstrained `
  --model B_-1 --run-id b_uv_constrained_pilot_v1 `
  --seed-source historical_fit --seed-factors "0.005,0.01,0.02" `
  --sobol-points 8 --seed-count 1 --max-nfev 20

python -m scripts.FitT3BCUVHierarchyConstrained `
  --model C_-1 --run-id c_uv_constrained_pilot_v1 `
  --seed-source historical_fit --seed-factors "0.005,0.01,0.015" `
  --sobol-points 8 --seed-count 1 --max-nfev 20

Get-Content -Encoding UTF8 output\studies\general_scalar_extended\small_uv_running_5e-9\b_uv_constrained_pilot_v1\summary.md
Get-Content -Encoding UTF8 output\studies\general_scalar_extended\small_uv_running_5e-9\c_uv_constrained_pilot_v1\summary.md
```

If the pilots remain poor, use **new run IDs** to widen the search. For
example, drop `--seed-source` and `--seed-factors`, use `--sobol-points 16`,
`--seed-count 2`, `--max-nfev 60` (significantly slower). Never reuse a run ID.

## Outputs and provenance

For each model, independent immutable directory:

`output/studies/general_scalar_extended/small_uv_running_5e-9/<run-id>/`

- `fit.json`: complete candidate search inventory counts, best screened fit,
  fitted Yukawa coordinates, UV topology coupling, predicted oscillation
  observables, running-mass diagnostics and source hashes.
- `fitted_config.json`: new complete input configuration **only if** a
  screened candidate was evaluated, including modified `lambdaT3` and
  Yukawas. Not a certified physical model configuration.
- `summary.md`: concise acceptance/limitations report.

Old campaigns, Batch 11–16 runs, historical fits and code are never overwritten.

## Scientific limitations

A successful Sobol+least-squares `chi2 <= 5.86` fits five observables using
19 free real coordinates, so it demonstrates candidate numerical fit
existence only, not predictivity. The physical Matchete–RGBeta map is not
certified. The existing one-scale coefficient does not include the complete
renormalised general-scalar sequential EFT matching. Positive running scalar
masses at endpoints do not prove vacuum stability or physical pole ordering.
This is a numerical feasibility benchmark, suitable for a thesis only with
explicit qualifications and caveats.

## CLI fix (Batch 17.1)

The first Batch 17 ZIP incorrectly referenced `a.scalar_min` and `a.scalar_max`
inside `main()` while argparse defines `--scalar-min-tev` and
`--scalar-max-tev`. This resulted in an `AttributeError` **before numerical
running or output creation**. The replacement CLI passes
`a.scalar_min_tev` and `a.scalar_max_tev` to `run_model`. Two new regression
tests execute the CLI dispatch with default and explicitly set scalar windows.
No physics functions, beta equations, matching, or fit parameters were changed.

</details>

---

<a id="source-t3-bc-uv-neutrino-normalization-batch16-md"></a>

## Original document: `docs/T3_BC_UV_NEUTRINO_NORMALIZATION_BATCH16.md`

**Original SHA-256:** `92e0ee42224fc14db25a7bb34e30e6dcceb5e6d98102099aed0803baf8f1d362`

**Status:** Historical record; check the current source and main
documentation before executing commands reproduced below.

<details>
<summary>Expand original: Batch 16 — UV-running B(-1)/C(-1) hierarchy-safe neutrino normalization probes</summary>

# Batch 16 — UV-running B(-1)/C(-1) hierarchy-safe neutrino normalization probes

## Purpose and verified source interfaces

Batch 15's smaller-Yukawa scan found seven *mass-screen-passing* candidate points, not neutrino oscillation fits:

- B(-1): common Yukawa factors `0.005,0.01,0.015,0.02`.
- C(-1): common Yukawa factors `0.005,0.01,0.015`.

With the original `lambdaT3=5e-9`, scaling both Yukawa matrices can strongly suppress the one-loop neutrino mass. Test an **illustrative normalization restoration** by varying

\[
(y_1,y_2)_{UV}=f(y_1,y_2)_{\rm seed},\qquad
\lambda_{T3,UV}=\frac{5\times10^{-9}}{f^2}\,k,
\]

where `k=0.1,1,10` by default. The candidate seven extra quartics at the UV boundary are still `0.1*abs(lambdaT3)`, so they scale with the topology coupling. The UV input scalar masses and heavy-fermion masses are unchanged.

**This is a probe, not a refit.** Holding each Yukawa matrix texture fixed means its neutrino mixing pattern and mass-ratio pattern may still be wrong. Scaling lambda principally changes neutrino mass normalization; optimizing 18 Yukawas and lambda is the subsequent step if needed. No first-principles relationship `lambda ∝ 1/f²` is asserted: it is a seed heuristic, and the RG equations are re-integrated for every point.

## Inputs and provenance

Reads original historical `small_uv_running_5e-9` files and source-hash-verified prior `general_fresh_v1_B_m1_cont`, `general_fresh_v1_C_m1_cont` outputs via Batch 14 `fit_source`. Two optional UV Yukawa texture sources are evaluated:

- `saved_poor_fit`: the saved general-scalar exploratory best Yukawa matrices;
- `historical_fit`: the original restricted-potential UV-running fit Yukawas.

Each trial changes copies of only the 18 real Yukawa entries and the `lambdaT3` UV boundary. `source['lambdaT3']` **and** `source['cfg']['base_state']['ordinary']['lambdaT3']` are changed together; this is essential because the existing Batch 9 boundary builder reads the former and Stage 31 state initialization reads the latter.

The established Batch 11 logic supplies candidate one-loop RGBeta UV/EFT1, the heavy mass-basis rotation, the saved one-scale ordered Matchete C5 with the **conditional** Stage23/26 minus-sign map, and SM+Weinberg running to 100 GeV. All five NuFIT observables and their approximate chi-square are recomputed for every eligible candidate.

Before evaluating C5, the runner requires positive 2×2 running scalar mass eigenvalues and both:

\[
\rho_F=\frac{m_{S,\max}(10^5\,\mathrm{GeV})}{M_{F,\min}(10^5\,\mathrm{GeV})}\le 0.5,
\qquad m_S(10^3\,\mathrm{GeV})\in[0.5,5]\,\mathrm{TeV}.
\]

Failures at the fermion threshold skip EFT1; failures at the scalar threshold skip the neutrino calculation. The one-scale C5 is never presented as a certified physical or sequential matched result. Running masses are not pole masses, and the scalar window is a numerical benchmark convention, not a validated decoupling prescription.

A candidate with chi-square <=5.86 is a conditional seed-level fit only. Because no optimization occurs, neither `chi2` nor any parameter is promised to improve.

## Files (new; no existing file replacements)

- `scripts/ProbeT3BCUVNeutrinoNormalization.py` — full standalone numerical probe runner.
- `tests/python/test_t3_bc_uv_neutrino_normalization.py` — 25 offline contract and mock-interface tests.
- `docs/T3_BC_UV_NEUTRINO_NORMALIZATION_BATCH16.md` — this file.

All inputs are read-only. The results are new-only under:

```
output/studies/general_scalar_extended/diagnostics/<new-run-id>/
  probe.json
  probe.csv
  summary.md
```

`probe.json` records source hashes, full trial predictions (including C5 matrices where computed), the threshold mass screen, and explicit noncertification flags. `probe.csv` is a compact flat table.

## PowerShell

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "C:\path\to\T3_BC_UV_Neutrino_Normalization_Batch16.zip" -DestinationPath . -Force
python -m pytest -q tests\python\test_t3_bc_uv_neutrino_normalization.py

# First: inexpensive 14-trial fixed-normalization probe (two Yukawa textures).
python -m scripts.ProbeT3BCUVNeutrinoNormalization `
  --run-id uv_neutrino_norm_k1_v1 `
  --lambda-multipliers "1"

Get-Content -Encoding UTF8 output\studies\general_scalar_extended\diagnostics\uv_neutrino_norm_k1_v1\summary.md

# If needed, a separate 42-trial lambda sensitivity run; never reuses an ID.
python -m scripts.ProbeT3BCUVNeutrinoNormalization `
  --run-id uv_neutrino_norm_grid_v1 `
  --lambda-multipliers "0.1,1,10"

Get-Content -Encoding UTF8 output\studies\general_scalar_extended\diagnostics\uv_neutrino_norm_grid_v1\summary.md
```

To run fewer trials, select `--model B_-1 --seed-source historical_fit`, or supply `--yukawa-factors "0.005,0.01"` with `--model`. Use a unique run ID for every new set.

## Validation boundary

Offline tests cover boundary consistency, seed ordering, mass-screen gating, skip-C5 on invalid hierarchy, successful UV/EFT1/C5/low-energy mock-interface sequence, input parsing and nonoverwrite guards. No live RGE/C5 calculation is claimed to have run in the offline staging tree: the full numerical services and saved parameter inputs require the user's Windows repository.

If the new candidate chi-square remains large, the next step is a **hierarchy-constrained 18-Yukawa plus lambda refinement** starting from the best mass-safe point, not another unconstrained historical-UV fit.

</details>

---

<a id="source-t3-bc-uv-running-mass-hierarchy-batch14-md"></a>

## Original document: `docs/T3_BC_UV_RUNNING_MASS_HIERARCHY_BATCH14.md`

**Original SHA-256:** `658b5959806ab8adee621189b650ead093ebff92a8d5ccfc03806759c3742588`

**Status:** Historical record; check the current source and main
documentation before executing commands reproduced below.

<details>
<summary>Expand original: Batch 14 — B(-1)/C(-1) UV-running mass hierarchy replay</summary>

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

</details>

---

<a id="source-t3-bc-uv-scalar-mass-feasibility-batch15-md"></a>

## Original document: `docs/T3_BC_UV_SCALAR_MASS_FEASIBILITY_BATCH15.md`

**Original SHA-256:** `93494a3ab1b31ec098063b4de4ec43fcf68013b0f481b361c8a44b1bdd6eb249`

**Status:** Historical record; check the current source and main
documentation before executing commands reproduced below.

<details>
<summary>Expand original: Batch 15 — B/C UV-running scalar-mass feasibility, before any more fits</summary>

# Batch 15 — B/C UV-running scalar-mass feasibility, before any more fits

## Scope

The latest Batch 14 live replay for B(-1) / C(-1) found

| Model | mS,max(100 TeV) | MF,min(100 TeV) | mS,max/MF,min | Former exploratory chi2 |
|---|---:|---:|---:|---:|
| B(-1) | 145.31 TeV | 84.395 TeV | 1.7218 | 18,376.7 |
| C(-1) | 180.64 TeV | 104.30 TeV | 1.7319 | 20,279.3 |

The same fitted Yukawa point and a zero-extra-scalar-quartic boundary give
relative scalar mass-squared changes of only ~10^-12, so adding those
small quartics at the boundary does not cause the hierarchy inversion.

**Verified source detail (GitHub main a05b3b7):** In the native projected
UV mass-squared beta expressions, the terms quadratic in the heavy-fermion
mass and Yukawas occur with coefficient **-8** on each of two `Tr[...]`
structures for B(-1), and **-12** for C(-1). The Yukawa-dependent `mS_i^2`
term has coefficient 2 (B) / 3 (C) multiplying `Tr[y_i.Trans[Bar[y_i]]]`.
The convention is `16*pi^2 dm_i^2/dln(mu)=beta^(1)`. Because running is
from 10^7 **down** to 10^5 GeV, negative MF²Y² beta contributions can
*increase* positive running scalar mass-squared. This is why scanning the
Yukawas' overall scale is more informative than changing small quartic inputs.

## What the new script genuinely does

**New source:** `scripts/StudyT3BCUVScalarMassFeasibility.py`.

1. Read saved Batch 13 / Batch 12 general-scalar UV-running best Yukawas from
   `general_fresh_v1_B_m1_cont` / `general_fresh_v1_C_m1_cont` using the
   corrected Batch 14 loader. Verify every saved source SHA-256.
2. Copy the historical initial state *in memory*. Apply a uniform factor
   `f_y` to all 18 real fitted Yukawas, and a positive factor `f_S` to both
   historical initial scalar masses: `mS_i^2(UV) -> f_S^2*mS_i^2(UV)`.
   No negative scalar mass input, new fermion masses, or modified beta is used.
   Each general quartic remains `0.1*abs(lambdaT3)` at the UV boundary,
   or the actual source-fit value if it differs.
3. Reintegrate native candidate UV RGBeta (`10^7 -> 10^5 GeV`) for each point,
   including running of the full Yukawas and MF. Diagonalise the candidate
   two-doublet 2x2 scalar quadratic matrix including its actual generated
   `mu2` and compare to the smallest heavy-fermion singular mass at `muF`.
4. **If and only if** the candidate has positive scalar eigenvalues and
   `mS,max(muF) < MF,min(muF)`, run the existing candidate scalar-only EFT1
   (`10^5 -> 10^3 GeV`), retaining the additional quartic and quadratic
   evolution. Otherwise mark the inversion and skip the invalid EFT1 stage.
5. Record two independent screening conditions:
   - conservative F-first mass ratio at muF: `mS,max/MF,min <= 0.5`;
   - both positive running scalar masses at muS inside `[0.5,5] TeV`.
   Both thresholds are adjustable *screening conventions*, not physical pole
   threshold prescriptions.
6. At the baseline (`f_y=f_S=1`) separately evaluate the **actual native UV**
   `mS1Sq` and `mS2Sq` beta expressions at the UV start and `muF`, substituting
   `(y1,y2) -> 0` and `MF -> 0` at the *same frozen state*. Inclusion–exclusion
   gives exact local contributions: no-(Y/MF), Yukawa-only at MF=0,
   MF-only at Y=0, and joint Y–MF dependence. They sum algebraically to the
   full beta. These are not unique diagrams and do not integrate by themselves
   to the net mass change. The baseline must reproduce the previous fitted
   MF(muF), and will independently verify the Batch 14 scalar mass at muF
   when `uv_mass_hierarchy_v1/diagnostic.json` is present.

**Not done:** no new C5 evaluation, no new chi-square minimisation, no
physical pole masses, no sequential thresholds, no new certification of the
Matchete/RGBeta coupling map. These are **potential starting points** for
future constrained Yukawa refits, not accepted neutrino fits. A uniform
Yukawa reduction generally changes the neutrino observables drastically.

If the scan has zero successes, it does not establish a no-go theorem: it
only excludes the finite positive input-factor grid under this candidate RGE.
For tests reporting numerical failures, inspect the saved reason rather than
interpreting these points as a physical exclusion.

## Files to extract (complete files, no patches)

- `scripts/StudyT3BCUVScalarMassFeasibility.py`
- `tests/python/test_t3_bc_uv_mass_feasibility.py`
- `docs/T3_BC_UV_SCALAR_MASS_FEASIBILITY_BATCH15.md`

The prior Batch 7, 9, 11, 12 and **fixed** Batch 14 scripts must already be
installed, and the existing source beta JSON exports must be available.

## PowerShell, project root

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "C:\path\to\T3_BC_UV_Scalar_Mass_Feasibility_Batch15.zip" -DestinationPath . -Force
python -m pytest -q tests\python\test_t3_bc_uv_mass_feasibility.py tests\python\test_t3_bc_uv_scalar_hierarchy.py

# Default: B and C, 4 Yukawa factors x 3 scalar-mass factors = 24 UV replays.
# Only positive, F-first ordered endpoints proceed into EFT1.
python -m scripts.StudyT3BCUVScalarMassFeasibility --run-id uv_mass_feasibility_v1

Get-Content -Encoding UTF8 output\studies\general_scalar_extended\diagnostics\uv_mass_feasibility_v1\summary.md
```

For a smaller B-only pilot with exact historical baseline included:

```powershell
python -m scripts.StudyT3BCUVScalarMassFeasibility `
  --model B_-1 --yukawa-factors "0.1,0.5,1" `
  --uv-scalar-factors "0.5,1" `
  --run-id uv_mass_feasibility_b_pilot_v1
```

## Outputs (new directories only)

`output/studies/general_scalar_extended/diagnostics/<new-run-id>/`

- `screen.json`: all candidates, source hashes, running-mass screens and
  local UV mass beta dependence, with no invented chi2 scores.
- `scan.csv`: one row per trial, ready for inspection/plotting.
- `summary.md`: threshold screen results and beta-sector table.

Existing historical files, previous fits and Batch 14 results are never written.
Run ID collisions are refused. The script does not edit GitHub.

## What to do with the output

If the screen finds well-separated, TeV-scale scalar candidates, pick the
candidate closest to the saved fit Yukawa direction and construct a **new**
constrained 18-real-Yukawa fit with those UV masses fixed or included in a
controlled optimisation. Require both `chi2 <= 5.86` and mass ordering on
recomputed trajectories. This requires a separate follow-up batch; do not
read the old chi2 as the new candidate's fit quality.

If no positive-UV-mass configuration reaches both screens, a more general
scan of Yukawa textures and/or heavy-fermion/scalar input masses might help,
but those are different benchmark definitions. If the running mass hierarchy
remains inverted near decoupling, a proper scalar-first or near-degenerate
matching construction is needed; more optimisation alone cannot cure the EFT.

</details>

---

<a id="source-t3-b-plus1-coupled-eft1-batch20-md"></a>

## Original document: `docs/T3_B_PLUS1_COUPLED_EFT1_BATCH20.md`

**Original SHA-256:** `c3f5bfafe1564d8d2ab76987acfc81b7c5246d38b1e194f0e53a725c8deb0a84`

**Status:** Historical record; check the current source and main
documentation before executing commands reproduced below.

<details>
<summary>Expand original: Batch 20 — B(+1) native coupled general-scalar EFT1 beta extension</summary>

# Batch 20 — B(+1) native coupled general-scalar EFT1 beta extension

## Why this batch

Batch 19's local Matchete matching probe succeeded: the single B(+1) new scalar
quartic `lambdaT3ExtraQ022000I1` is absent from the extracted **one-scale**
ordered C5 coefficient, and three numerical points did not detect any direct
one-scale response. This **does not remove** the additional quartic from the
renormalisable RG equations. Stage 31 includes the RGBeta vertex in the UV;
existing Stage 34 explicitly omitted its RG evolution **below** the fermion
threshold.

This batch fills precisely that gap, **for B(+1) only**:

- New independent native RGBeta EFT1 builder first calls
  `T3RGBetaBuildEFT1[2,2,1,1]` with F absent; registers the existing checked
  `RGE/general/generated/B_+1_extra_quartics.wl` interaction; projects the
  complete `BetaTerm[Quartic,1]` and checks the exact 18 baseline + 1 new
  quartic beta equations. The explicit quartic vertex is
  `lambdaT3GeneralQ022000I1` in RGBeta, corresponding *conditionally* to
  Matchete `lambdaT3ExtraQ022000I1`.
- Reads the checked **Stage 33** fermion-threshold state and its nonzero sidecar
  boundary, enforcing the 16-model summary and the saved Stage31 hashes via
  the current `IntermediateEFTRunning.verify_boundary` implementation.
- Uses the production InputForm beta evaluator and state/beta-vector conventions
  to integrate **all 18 baseline coordinates + one complex extra quartic**
  simultaneously from 100 TeV to 1 TeV, one-loop DOP853, rtol=1e-8,
  atol=1e-11. The additional coupling affects all baseline betas in which the
  native projected RGBeta result contains it; it is **not frozen**.
- Separately evolves the same system with the extra quartic initially zero.
  This control is required to remain zero. Where a matching original Stage34
  run exists with the same Stage33 input hash, the zero-control final scalar
  parameters must agree to <=3e-6 relative (the actual error is recorded).
- Writes immutable new output under `output/studies/general_scalar_extended/
  b_p1_eft1_running/<run-id>/` with full final scalar EFT1 state, new-quartic
  trajectory, differences from zero extra, UV source fingerprints, beta
  expressions, and local native Wolfram log.

## Installation

The ZIP has **four new full files**. No existing production or historical file
is replaced:

| Project path | Purpose |
|---|---|
| `scripts/ExportT3BPlus1EFT1GeneralOneLoop.wl` | Native B(+1) RGBeta 18+1 EFT1 beta export |
| `scripts/RunT3BPlus1GeneralEFT1.py` | Provenance checks + coupled numerical running and control |
| `tests/python/test_t3_b_plus1_general_eft1.py` | Offline interface, fail-closed, and synthetic ODE tests |
| `docs/T3_B_PLUS1_COUPLED_EFT1_BATCH20.md` | This file |

In PowerShell from your local repository root:

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive `
  -LiteralPath "C:\path\to\T3_B_Plus1_Coupled_EFT1_Batch20.zip" `
  -DestinationPath . -Force

python -m pytest -q tests\python\test_t3_b_plus1_general_eft1.py

python -m scripts.RunT3BPlus1GeneralEFT1 `
  --run-id b_p1_eft1_general_v1

Get-Content -Encoding UTF8 `
  output\studies\general_scalar_extended\b_p1_eft1_running\b_p1_eft1_general_v1\summary.md
```

Prerequisites (the already-generated files, do not rerun unless missing or stale):

- `output/general_t3_stage22/stage33_B_p1_fermion_threshold_handoff.json`
- `output/general_t3_stage22/stage33_16_model_fermion_threshold_summary.json`
- `output/general_t3_stage22/stage31_16_model_uv_summary.json`
- Stage 31 source and SHA-provenanced 16-model reports.
- Local native Wolfram `wolframscript` with RGBeta and the same generated
  `B_+1_extra_quartics.wl` vertex.

**If a run errors**, keep its `wolfram_export.log` and `native_eft1_betas.json`
(if present) and use a different run ID after fixing the cause. The runner
refuses to overwrite a previous run folder.

## Outputs

- `native_eft1_betas.json`: full native projected one-loop equations.
- `wolfram_export.log`: complete stdout/stderr from the native Wolfram process.
- `report.json`: full source hashes, complete final EFT1 state, trajectory,
  numerical differences, and explicit physics-scope flags.
- `summary.md`: small human-readable screening report.

Expected success marker:

`SuccessCoupledBPlus1GeneralEFT1Running`

If the native projector is singular or fails to include a previously existing
quartic beta, the export refuses to write a valid success report. The exporter
must not substitute unprojected `BetaTerm[symbol,1]` for missing class-level
quartics or supply zero for an absent beta.

## What this proves and what it does not

It establishes that B(+1)'s extra quartic can be evolved as a **coupled
one-loop renormalisable scalar-EFT1 coordinate** in the established
conditional Stage29/31 RGBeta basis, once the live export and integration
succeed. It retains the Stage33 no-finite-renormalisable threshold boundary
convention and the previously fitted illustrative UV state.

It does **not** yet calculate a neutrino mass, any of the five NuFIT residuals,
or any of the original three refits. No C5 or EFT Wilson coefficients are
matched or run in this batch. Neither the one-scale Matchete/RGBeta quartic
basis map nor the complete sequential EFT chain is certified. Positive
running quadratic scalar masses are not physical pole masses or a potential
stability certificate.

### Next numerical batch

Connect this coupled B(+1) EFT1 evolution to a **per-trial** fitted input
pipeline for three historical studies: small frozen UV, large frozen UV, and
small full UV. The fitter must not simply load the Stage29/31 saved benchmark
trajectory, because those values correspond to a single UV input rather than
the new Yukawa optimisation trials. Use B(+1)'s model-specific C5 matching
source, the established SM+Weinberg endpoint, and exact NuFIT target convention;
continue to mark one-scale matching as exploratory.

</details>

---

<a id="source-t3-b-plus1-general-c5-probe-batch19-md"></a>

## Original document: `docs/T3_B_PLUS1_GENERAL_C5_PROBE_BATCH19.md`

**Original SHA-256:** `87efaa439adad6e654b123031b192c1fb283e577408546279e45f011f7071f87`

**Status:** Historical record; check the current source and main
documentation before executing commands reproduced below.

<details>
<summary>Expand original: Batch 19 — First non-B(-1)/C(-1) general-scalar C5 source probe: B(+1)</summary>

# Batch 19 — First non-B(-1)/C(-1) general-scalar C5 source probe: B(+1)

## Scope and verified project interfaces

This experiment targets **B(+1)**, a distinct ordinary T3 model with one
additional complex quartic and **zero additional quadratic mixing operators**
according to `RGE/general/generated/general_rgbeta_quartic_registry.json`:

- RGBeta: `lambdaT3GeneralQ022000I1`;
- Matchete: `lambdaT3ExtraQ022000I1`;
- fields: `Hdag Hdag S1 S1 + h.c.`.

The existing `Lagrangian/interactions/GeneralScalarOperators.wl` enables that
quartic with environment variable `T3_EXPERIMENTAL_GENERAL_QUARTICS=1`, and
`Lagrangian/RunModel.wl` exports an ordered C5 in `c5_coefficient.txt` when
Matchete matching and operator extraction succeed.

This batch executes the **same existing** Matchete driver twice at EFT order 5,
one loop, with **all heavy particles removed at one common step**:

1. Restricted Lagrangian (`T3_EXPERIMENTAL_GENERAL_QUARTICS=0`).
2. Extended Lagrangian (`T3_EXPERIMENTAL_GENERAL_QUARTICS=1`).

It fails closed if the extended interaction is not registered or one of the
matching/extraction statuses is not `Success`. It records source SHA-256 hashes.
The second Wolfram script reads the two actual ordered C5 expressions and
compares them at three *numerical one-generation flavour spurion* points,
including the exact zero-extra-quartic reduction. It probes one nonzero extra
quartic value and records any numerical response. A very large dimensional
regulator variable is used for a finite-limit **sample**, not an MSbar-renormalised
C5 calculation. No neutrino mass fit is performed.

The resulting status distinguishes **source occurrence** from **numeric
sensitivity**. If the extra coupling is absent from the ordered C5 output,
that observation does not imply it has no effect through the UV/EFT1 RGEs.

## Files (new, full files, no replacements)

- `scripts/ProbeT3BPlus1GeneralC5.py`: new-only run coordinator and structural checks.
- `scripts/CompareT3BPlus1OneScaleC5.wl`: native Mathematica numerical ordered-C5 check.
- `tests/python/test_t3_b_plus1_general_c5_probe.py`: offline fail-closed regressions.
- `docs/T3_B_PLUS1_GENERAL_C5_PROBE_BATCH19.md`: this note.

## PowerShell, from project root

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "C:\path\to\T3_B_Plus1_General_C5_Probe_Batch19.zip" -DestinationPath . -Force
python -m pytest -q tests\python\test_t3_b_plus1_general_c5_probe.py
python -m scripts.ProbeT3BPlus1GeneralC5 --run-id b_p1_c5_probe_v1
Get-Content -Encoding UTF8 output\studies\general_scalar_extended\matching_probes\B_p1\b_p1_c5_probe_v1\summary.md
```

The native Matchete calculation can be expensive. If it fails, inspect
`baseline.log`, `general.log` or `numeric_audit.log` in the run directory.
A failed run retains its directory and must **not** be overwritten; rerun
with a new ID after correcting the specific problem.

Output:

`output/studies/general_scalar_extended/matching_probes/B_p1/<run-id>/`

with `baseline/`, `general/`, source certificates, log files, numeric audit,
`report.json` and `summary.md`. Historical data are read-only.

## Scientific limitations

- **Implemented:** symbolic restricted and extended Matchete one-scale matching,
  registered new quartic check, numeric ordered-C5 source comparison.
- **Not implemented:** physical Matchete–RGBeta quartic normalisation, full
  sequential EFT1 Wilson matching, general-scalar mixing for other models,
  new B(+1) UV/EFT1 parameter-to-C5 adapter, oscillation refitting.
- **Local verification pending:** the live Wolfram matching and numeric comparison
  must run successfully in the user's native Windows Matchete installation.

The B(+1) result is the first **source gate** toward extending the original
three numerical studies to the other 14 models, **not yet an accepted new
model-study fit**. The next numerical step depends on whether the new quartic
is present in, or drops out of, the actual one-loop C5 source.

## Batch 19 CLI correction and interrupted-run recovery

`RunModel.wl` identifies B(+1) in `comparison_summary.json` as
`ModelClass: T3-d2-d2-F1`, `Alpha: 1`, not the short CLI label `B`.
The probe now checks the actual representation and hypercharge values.
No physics equations or Mathematica matching source files were modified.

To reuse the **existing baseline match from the original failed run**, replace
the files in this archive and execute:

```powershell
python -m pytest -q tests\python\test_t3_b_plus1_general_c5_probe.py
python -m scripts.ProbeT3BPlus1GeneralC5 --run-id b_p1_c5_probe_v1 --resume
```

`--resume` validates the existing baseline output and only runs missing stages.
The original version did not save source SHA-256s at the start. For that
original interrupted run, the resume checks source modification times and
records `baseline_source_hash_certified=false`; timestamps are **not**
cryptographic proof. New runs save a `source_fingerprints.json` before matching
and verify it on resume. To avoid this limitation altogether, run with a fresh
`--run-id` without `--resume` (which repeats the Matchete baseline).

## Batch 19 correction 2: local Matchete model-class label

The user's live output reports `ModelClass: "T3-B"` while the saved
repository B(+1) reference reports `"T3-d2-d2-F1"`. The prior validator
accepted only the latter and incorrectly rejected two successful baseline
runs. The corrected validator accepts **only those two** class-label spellings,
while requiring **all** of `Alpha=1`, `Scalar1SU2=2`, `Scalar2SU2=2`,
`FermionSU2=1`, `Scalar1Hypercharge="1/2"`,
`Scalar2Hypercharge="3/2"`, and `FermionHypercharge="1"`. A class
label alone is never treated as sufficient evidence of B(+1).

The accepted label is also recorded in `report.json`. No physics code,
matching source, or coupling convention is changed.

Use the latest full-file replacement ZIP and **resume the already completed
baseline in `b_p1_c5_probe_v2`**, which has source fingerprints recorded
before its first native Matchete run:

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath 'C:\path\to\T3_B_Plus1_General_C5_Probe_Batch19_Fix2.zip' -DestinationPath . -Force
python -m pytest -q tests\python\test_t3_b_plus1_general_c5_probe.py
python -m scripts.ProbeT3BPlus1GeneralC5 --run-id b_p1_c5_probe_v2 --resume
Get-Content -Encoding UTF8 output\studies\general_scalar_extended\matching_probes\B_p1\b_p1_c5_probe_v2\summary.md
```

A legacy `b_p1_c5_probe_v1 --resume` is also supported if its timestamp
provenance check passes. Prefer `v2` because it has the full SHA-256 manifest.

The offline tests are mock-based; real Matchete and native numerical C5 checks
remain to be established by this local run.

## Batch 19 correction 3: malformed native Wolfram numeric JSON

The user's B(+1) `b_p1_c5_probe_v2 --resume` completed the restricted and
extended Matchete stages but failed while parsing `native_numeric_audit.json`
with Python `JSONDecodeError` (reported line 41, column 33). Without the actual
file bytes, the precise offending token is **not verified**. High-precision
Wolfram `*^`-style numeric notation is one plausible cause; the updated
comparison writer removes that possible incompatibility without changing any
of the underlying numerical C5 expressions or the high-precision comparisons.

- `scripts/CompareT3BPlus1OneScaleC5.wl` (full replacement): uses finite
  machine-precision real numbers for JSON diagnostics, normalises any Wolfram
  scientific `*^` exponent formatting to JSON `e`, and checks an independent
  Wolfram RawJSON import before writing the result. High-precision arithmetic
  for the three C5 comparisons is unchanged.
- `scripts/ProbeT3BPlus1GeneralC5.py` (full replacement): permits an updated
  **comparison script only** when resuming a source-fingerprinted run. Every
  production matching source must still have the original SHA-256. It reruns
  the *comparison only*, not the expensive Matchete matching stages, and
  archives previous `native_numeric_audit.json` and `numeric_audit.log` as
  `native_numeric_audit.before_json_fix.json` and
  `numeric_audit.before_json_fix.log`. Malformed JSON is never accepted as
  a successful matching result.
- `tests/python/test_t3_b_plus1_general_c5_probe.py`: regressions for corrupted
  numeric JSON, comparison-script-only changes, and rejecting changed matching
  source code. These are mocked tests, not a native Wolfram verification.

Run from the project root, after installing **this** full replacement archive:

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath 'C:\path\to\T3_B_Plus1_General_C5_Probe_Batch19_Fix3.zip' -DestinationPath . -Force
python -m pytest -q tests\python\test_t3_b_plus1_general_c5_probe.py
python -m scripts.ProbeT3BPlus1GeneralC5 --run-id b_p1_c5_probe_v2 --resume
Get-Content -Encoding UTF8 output\studies\general_scalar_extended\matching_probes\B_p1\b_p1_c5_probe_v2\summary.md
```

The new `report.json` distinguishes `baseline_source_hash_certified` (matching
sources), `full_original_source_manifest_unchanged` (false after the comparator
update), and `comparison_script_updated_since_match`. The previous malformed
JSON and log stay available for exact source inspection. If the corrected
report is still malformed, the runner fails and retains the new native audit.
Do not treat this as proof of a new B(+1) physical EFT prediction or a refit.

</details>

---

<a id="source-t3-b-plus1-three-general-scalar-studies-batch21-md"></a>

## Original document: `docs/T3_B_PLUS1_THREE_GENERAL_SCALAR_STUDIES_BATCH21.md`

**Original SHA-256:** `78daadce789f057c9b04b1992f71c304e6380d672c7cfe38234cce9f69f76b59`

**Status:** Historical record; check the current source and main
documentation before executing commands reproduced below.

<details>
<summary>Expand original: Batch 21 — B(+1) three-scenario conditional general-scalar numerical fits</summary>

# Batch 21 — B(+1) three-scenario conditional general-scalar numerical fits

## Scope

Runs exactly the original three study prescriptions for **B(+1)**, using the
original `lambdaT3` values and each scenario's historical fitted real Yukawa
matrices as seeds:

| `--case` | UV interval | Original lambdaT3 |
|---|---|---:|
| `small_uv_frozen_5e-9` | **frozen**, boundary imposed at 100 TeV | 5e-9 |
| `large_uv_frozen_0p5` | **frozen**, boundary imposed at 100 TeV | 0.5 |
| `small_uv_running_5e-9` | **coupled UV RGE**, 10^7 GeV to 100 TeV | 5e-9 |

After the fermion boundary the **native Batch20 RGBeta 18+1** system is
integrated from 100 TeV to 1 TeV. All three fit the 18 *real* Yukawa entries;
`lambdaT3` remains fixed. The new extra complex quartic is assigned a **real,
illustrative RGBeta boundary of -0.01 by default**, independently of the
small/large topology quartic. That is a numerical choice, not a derived physical
coupling map. Alternative: `--extra-real <value>` or `--extra-fraction 0.1`
(the latter uses `0.1 * abs(lambdaT3)` as the extra boundary instead).

The C5 evaluator uses **the model-specific existing**
`output/full/hypercharge/T3_B_alpha_p1/data/final_weinberg_coefficient.json`,
containing the hierarchical hard and direct logarithmic term already used by
Stage36. Batch19 showed no direct extra-quartic dependence in the one-loop
Matchete coefficient; **the new quartic affects C5 through RG running** in
this candidate approximation. A fast one-generation kernel is derived from
that existing symbolic source and checked against Stage36's full symbolic
`evaluate_final_weinberg_json` at the first mass-safe trial (relative error
must be < 3e-8). It is not a newly matched general-scalar Wilson result.

All three candidates evolve SM+Weinberg from 1 TeV to 100 GeV and use the
same five-observable NuFIT target and charged-lepton mass-basis fit as the
existing project. Refit cost is reduced in frozen UV by evaluating the scalar
EFT1 trajectory once; UV-running points independently rerun UV and EFT1.

## Install

Extract the archive into the project root. It adds:

- `scripts/FitT3BPlus1GeneralThreeStudies.py`
- `tests/python/test_t3_b_plus1_general_three_studies.py`
- `docs/T3_B_PLUS1_THREE_GENERAL_SCALAR_STUDIES_BATCH21.md`

It depends on **Batch20** being installed and on a successful, unchanged live
Batch20 output at
`output/studies/general_scalar_extended/b_p1_eft1_running/b_p1_eft1_general_v1/`.
The native beta SHA256 and original zero-extra control status are rechecked.
The Stage31 B(+1) coupled UV beta export and original three fitted source
configs must also exist. No Mathematica call occurs during the fit itself.

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
python -m pytest -q tests\python\test_t3_b_plus1_general_three_studies.py

# First run a cheap 3-scenario pilot (historical + 4 Sobol points, no local fit).
python -m scripts.FitT3BPlus1GeneralThreeStudies `
  --all --run-id b_p1_three_pilot_v1 --sobol-points 4 --max-nfev 0

# Full pilot optimisation. Reuse the same Batch20 betas, but a NEW run ID.
python -m scripts.FitT3BPlus1GeneralThreeStudies `
  --all --run-id b_p1_three_fit_v1 --sobol-points 8 `
  --max-nfev 25 --seed-count 1
```

Read outputs for every case:

```powershell
Get-Content -Encoding UTF8 output\studies\general_scalar_extended\small_uv_frozen_5e-9\b_p1_three_fit_v1\summary.md
Get-Content -Encoding UTF8 output\studies\general_scalar_extended\large_uv_frozen_0p5\b_p1_three_fit_v1\summary.md
Get-Content -Encoding UTF8 output\studies\general_scalar_extended\small_uv_running_5e-9\b_p1_three_fit_v1\summary.md
```

Each case writes a new immutable `fit.json`, `summary.md`, and (if a mass-safe
candidate exists) `fitted_config.json`. No historical input/output is changed.
The fitter prints the status and chi2 for each trial as it evaluates it.
When a case has no accepted fit, the report distinguishes `PoorFit` from
`NoMassSafeNeutrinoFit` and records up to 15 numerical or mass-screen failures.
Do not label these cases successful without checking their live summaries.

## Numerical screens

- `rhoF = max(sqrt(mS1Sq(muF)),sqrt(mS2Sq(muF)))/min(MF_singular(muF)) <= 0.5`.
- Positive scalar *running quadratic* masses at 1 TeV within 0.5–5 TeV.
- Yukawa maximum singular value at the fermion threshold <= sqrt(4 pi).
- Scalar coordinate magnitudes at the 1 TeV endpoint and the extra-quartic
  boundary <= 4 pi. This **does not prove** vacuum stability or a global
  perturbativity bound over every scale.
- Recomputed five-observable chi2 <= 5.86 constitutes a *candidate accepted*
  fit. The fit has 18 Yukawa entries for five observables, so it is not a
  prediction of the angles.

## Limitations

- B(+1) only; this does not add all remaining 42 results.
- Matchete–RGBeta general quartic coefficient normalization remains conditional.
- The reference Stage36 B(+1) hard-plus-direct C5 is reused; the extra quartic
  is not directly rematched at either threshold.
- No matched `C12`/Wilson running from the fermion threshold; grouped scalar
  matching is an approximation.
- The displayed masses are unbroken-phase running quadratic masses, not pole
  masses; vacuum stability and unitarity are not tested.

</details>

---

<a id="source-t3-reduced-uv-scalar-boundary-pilot-md"></a>

## Original document: `docs/T3_REDUCED_UV_SCALAR_BOUNDARY_PILOT.md`

**Original SHA-256:** `bedf945634b682184dc364847ba788e684f46c8d64d42a1647c2d182f6a2789a`

**Status:** Historical record; check the current source and main
documentation before executing commands reproduced below.

<details>
<summary>Expand original: A(-2) UV-running scalar mass-boundary feasibility pilot</summary>

# A(-2) UV-running scalar mass-boundary feasibility pilot

## Scope

The completed frozen-UV Sobol/A(-2) comparison achieved 32/32 accepted individual fits. The remaining original running-UV `small_uv_running_5e-9` scenario was screened out for all 16 models. The diagnostic records had scalar masses about 30–169 TeV at the **fixed scalar threshold** of 1 TeV, above the retained 0.5–5 TeV mass window; some models also failed the fermion-first hierarchy diagnostic.

This new script checks **one model, A(-2)** without modifying any earlier campaign files. It uses the source-verified candidate one-loop UV and EFT1 RGBeta equations of the current repository, including all five A(-2) additional quartic sidecar coefficients starting at **exactly zero**, and the existing fermion-threshold projection. No additional quadratic term is added. The topology coupling is fixed at `lambdaT3(1e7 GeV)=5e-9`.

The pilot varies **only** the two *UV scalar squared-mass input parameters* at `1e7 GeV`; all Yukawa matrices, fermion masses, other scalar quartics, thresholds and RGEs remain fixed. By default, the Yukawas come from the historical A(-2) running-low fitted configuration. An explicit second option injects the frozen-low overnight A(-2) Yukawas as a **hypothetical UV Yukawa input**, not a fitted UV solution.

It is a feasibility test, not a new neutrino fit, and it does **not** claim to rule out UV-running T3 models if no positive UV scalar boundary exists for these fixed Yukawas.

## Numerical calculation

At fixed Yukawas and other initial inputs, the one-loop beta functions for the two scalar mass-squared parameters form an inhomogeneous linear system, while the dimensionless RGEs and fermion mass evolution are independent of these scalar squared masses in the examined truncation. Their combined UV+EFT1 trajectory therefore defines a candidate affine map

\[
\begin{pmatrix} m_{S_1}^2 \\ m_{S_2}^2 \end{pmatrix}_{\mu_S=10^3\,\mathrm{GeV}}
= A\begin{pmatrix}m_{S_1}^2\\m_{S_2}^2\end{pmatrix}_{\Lambda=10^7\,\mathrm{GeV}}+b,
\]

with input/output mass squares expressed in **TeV\(^2\)**. The script does not assume a specific numerical `A,b`: it obtains them from the *real repository's numerical ODE solvers*. It runs the historical UV boundary and independent perturbations of `4 TeV^2` (default), then tests the affine formula against **two extra full RG trajectories** (zero UV masses and a simultaneous two-mass perturbation). The calculation fails closed if the error exceeds `0.05 TeV^2`.

It then:

1. Solves algebraically for UV masses that would give **2 TeV scalars at muS**, retaining negative `m^2` solutions *only as diagnostics*.
2. Uses a constrained linear program to determine whether any **nonnegative UV mass-squared inputs** can place *both* low-scale scalar masses in the fixed 0.5–5 TeV range.
3. Independently checks the more restrictive condition that UV scalar masses themselves lie in 0.5–5 TeV.
4. If feasible, reruns the **full UV+EFT1 integration** at each proposed positive UV input, and checks the original mass hierarchy, scalar threshold and perturbativity screens. An affine feasibility result does not automatically become a physically screened accepted input.

The original fixed matching scales remain:

- `muUV = 1e7 GeV`
- `muF = 1e5 GeV`
- `muS = 1e3 GeV`

The target 0.5–5 TeV band, `rhoF <= 0.5` and Yukawa/quartic perturbativity are **not relaxed**.

## Full files in this package

- `scripts/ProbeT3ReducedUVScalarBoundaries.py` — complete new pilot source
- `tests/python/test_t3_reduced_uv_scalar_boundary_pilot.py` — offline tests for the affine reconstruction, numerical validation, signed roots and constrained mass-window feasibility
- `docs/T3_REDUCED_UV_SCALAR_BOUNDARY_PILOT.md` — this procedure

All outputs are isolated under `output/studies/reduced_scalar_uv_boundary_pilot/<run-id>/`; each run ID is write-once and the script refuses an existing run directory. No previous frozen fit, UV campaign or GitHub source is modified.

## Windows PowerShell

Extract the ZIP into your local checkout:

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath ".\T3_Reduced_UV_Scalar_Boundary_Pilot.zip" -DestinationPath . -Force
python -m pytest -q tests\python\test_t3_reduced_uv_scalar_boundary_pilot.py
```

Run the baseline A(-2) pilot:

```powershell
python -u -m scripts.ProbeT3ReducedUVScalarBoundaries --run-id a_m2_historical_v1
Get-Content -Encoding UTF8 output\studies\reduced_scalar_uv_boundary_pilot\a_m2_historical_v1\summary.md
```

Optional independent Yukawa-boundary sensitivity check, using the A(-2) best fit from the **completed frozen-low overnight run**:

```powershell
python -u -m scripts.ProbeT3ReducedUVScalarBoundaries --run-id a_m2_frozen_yukawas_v1 --yukawa-source frozen-newfit
Get-Content -Encoding UTF8 output\studies\reduced_scalar_uv_boundary_pilot\a_m2_frozen_yukawas_v1\summary.md
```

Output per run:

- `result.json`: exact forward trajectories, scalar mass² affine map, validation errors, inferred signed target and LP feasibility, physically screened follow-up trials, and input source digests.
- `summary.md`: concise result and caveats.

The script reuses `scripts/RunT3GeneralScalar48Campaign.py`, `scripts/RunT3ReducedScalarSobolComparison.py`, the existing Stage31 UV RGEs and the prior native EFT1 beta export. If the native EFT1 export is not already present, the previous reduced-scalar runner may invoke local `wolframscript` to generate it in this new isolated output directory.

## Interpretation

If both constrained programs are infeasible, that means **no choice of the two specified positive UV squared-mass parameters** achieves two light scalars at the lower threshold **for the chosen fixed Yukawas and other inputs**, in the validated candidate one-loop affine approximation. This is not a proof against simultaneously changing Yukawas, heavy fermion masses, quartics, thresholds or allowing different physics.

If a negative UV mass-squared input would be necessary for a 2 TeV output, it violates the selected positive-UV benchmark assumption; the script will not use that point in a physical running calculation.

If positive UV mass-squared solutions exist, they must also pass the original scalar/fermion hierarchy and perturbativity checks before being considered a viable next-stage UV fitting boundary. No vacuum stability, positive masses *throughout* the whole RG trajectory, or independent physical Matchete/RGBeta coupling mapping is certified.

</details>

---
