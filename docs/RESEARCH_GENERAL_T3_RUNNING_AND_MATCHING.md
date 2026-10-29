# Research General T3 Running And Matching

> Historical research record. This document combines **original source
> Markdown** under original-path headings. Local Markdown links to moved
> files are updated; exact original bytes are in the external archive.
> This is not a new certification
> of physics, a rewritten conclusion, or a new computational workflow.
> Original bytes and SHA-256 values are stored in the external archive.

## Contents

- [General T3: first coupled one-loop UV trajectory (A_+0)](#source-general-t3-a0-coupled-uv-running-md) — `docs/GENERAL_T3_A0_COUPLED_UV_RUNNING.md`
- [General T3: Convention-labelled UV quartic benchmark handoff](#source-general-t3-conventional-uv-handoff-md) — `docs/GENERAL_T3_CONVENTIONAL_UV_HANDOFF.md`
- [General T3: Coupled UV → fermion-threshold renormalisable handoff](#source-general-t3-fermion-threshold-handoff-md) — `docs/GENERAL_T3_FERMION_THRESHOLD_HANDOFF.md`
- [General T3: supported scalar-only intermediate EFT running](#source-general-t3-intermediate-eft-running-md) — `docs/GENERAL_T3_INTERMEDIATE_EFT_RUNNING.md`
- [T3 general-scalar batch: eight quadratic operators and one-loop mass beta functions](#source-general-t3-quadratic-batch-notes-md) — `docs/GENERAL_T3_QUADRATIC_BATCH_NOTES.md`
- [General T3: scalar-threshold renormalisable handoff](#source-general-t3-scalar-threshold-handoff-md) — `docs/GENERAL_T3_SCALAR_THRESHOLD_HANDOFF.md`
- [General T3 scalar interaction trial — batched validation package](#source-general-t3-trial-md) — `docs/GENERAL_T3_TRIAL.md`

---

<a id="source-general-t3-a0-coupled-uv-running-md"></a>

## Original document: `docs/GENERAL_T3_A0_COUPLED_UV_RUNNING.md`

**Original SHA-256:** `abc29da224ff0a0e0dca016df3f399f26eaf6c5d0a0e0acbec5eac9069e5c5a9`

**Status:** Historical record; check the current source and main
documentation before executing commands reproduced below.

<details>
<summary>Expand original: General T3: first coupled one-loop UV trajectory (A_+0)</summary>

# General T3: first coupled one-loop UV trajectory (A_+0)

## Scope

This is a **local, non-production** first numerical run for model `A_+0`,
`(dS1,dS2,dF,alpha) = (1,3,2,0)`, over `1e7 -> 1e5 GeV`.
It adds exactly the **five** Stage 29 candidate extra quartics to the established
UV ODE state. The other 15 ordinary T3 models and EFT thresholds are not
integrated in this stage.

`RGE/general/generated/A_+0_extra_quartics.wl` is loaded only after the actual
`T3RGBetaBuild[1,3,2,0]` call, **before** evaluating any beta function.
The exporter calls production `T3RGBetaOneLoopBetas[]`, then evaluates
**RGBeta's class-level `BetaTerm[Quartic,1]`**. RGBeta projects the whole
quartic class onto independent couplings; the exporter replaces every
production base-quartic beta and assigns the five extra-quartic betas from
the same projected Association. This avoids relying on per-coupling
`BetaTerm[symbol,1]`, for which RGBeta documents potential mixed and
non-diagonalised results. `T3RGBetaConventionalReportBetas` handles
the gauge convention. A projection failure or key mismatch aborts the export.

The Python runner uses the existing modules:

- `Numerical.fitting.ScanCLI.build_uv_state_from_config` to read the frozen
  fitted-reference UV point for A_+0;
- `Numerical.core.StateVector.pack_uv_state` and `unpack_uv_state`;
- `Numerical.core.RGBetaEvaluator.evaluate_inputform_expression` and
  `state_environment`;
- `Numerical.core.BetaVector.beta_values_to_derivative`.

It creates a 10-real-coordinate **sidecar** for five complex extra quartics,
without modifying the standard production state. The coupled derivatives are

\[
16\pi^2\frac{dX}{d\ln\mu}=\beta_X^{(1)}
\]

for **all** baseline fields and all five new quartics simultaneously. Therefore
these are neither frozen-background equations nor an extra-only toy flow.
The `DOP853` integrator uses `rtol=1e-8`, `atol=1e-11`, and at most `0.35`
steps in `ln(mu)`. Saved points are logarithmically spaced.

## Inputs and provenance

- Stage 29 handoff: `output/general_t3_stage22/stage29_convention_uv_handoff.json`
  (verified 16-model, 50-slot, 46-block boundary, default 0.01 on the Matchete
  side; Stage 30 uses the **recorded actual RGBeta values**, not an
  assumed uniform minus sign).
- Frozen A_+0 benchmark:
  `Numerical/benchmarks/fitted_reference/T3_dS1_1_dS2_3_dF_2_alpha_p0.json`
  (source-fitted Yukawas, not re-fitted here).
- Live RGBeta beta export:
  `output/general_t3_stage22/stage30_A_p0_full_uv_betas.json`.

A production/source provenance hash is emitted for each input. The live export
must supply exactly the baseline scalar/matrix beta names
matching the production packed state, plus five new beta expressions. It
**fails closed** if names, model, convention, or beta evaluation fail. The
exact baseline count is determined from the local production state; the code
does not hard-code that count.

## PowerShell (from project root)

```powershell
Expand-Archive -LiteralPath .\T3_General_Scalar_Stage30_QuarticProjectionFix.zip `
    -DestinationPath . -Force

python -m pytest -q tests\python\test_general_t3_coupled_a0_uv_running.py
if ($LASTEXITCODE -ne 0) { throw "Stage 30 tests failed" }

wolframscript -file .\scripts\ExportGeneralT3A0UVOneLoop.wl `
    output/general_t3_stage22/stage30_A_p0_full_uv_betas.json
if ($LASTEXITCODE -ne 0) { throw "Stage 30 RGBeta export failed" }

python -m scripts.CoupledA0UVRunning
if ($LASTEXITCODE -ne 0) { throw "Stage 30 integration failed" }
```

Outputs:

- `output/general_t3_stage22/stage30_A_p0_full_uv_betas.json`
- `output/general_t3_stage22/stage30_A_p0_coupled_uv_running.json`

The second file contains UV and final values of each quartic, gauge couplings
at every saved point, ODE diagnostics, and provenance hashes.

## Interpretation and limits

**Verified from source:** RGBeta `T3RGBetaBuild`, `BetaTerm` and production
beta convention exist and the numerical evaluator is able to parse its usual
exported arithmetic. The generated A_+0 extension declares five new quartics.

**Tested locally:** 29 Stage 30 Python regression cases, including exact
inventory checks, rejected stale/mismatched evidence, and backward ODE
integration of an independently solvable equation.

**Requires your machine:** actual projected `BetaTerm[Quartic,1]` expressions may be
complex, computationally slow, or may contain syntax not supported by the
current numerical parser; the integration will stop with a concrete error,
not silently set missing terms to zero. This archive is **not** evidence of a
successful live RGBeta one-loop run until both commands finish locally.

**Still conditional:** Stage 27 is used as a provisional triplet component
convention. Stage 28 established RGBeta's abstract `S2` identities, not a
unique numerical component identification. The Matchete-to-RGBeta signed
coupling map therefore remains **physically uncertified**. No threshold
matching, neutrino-mass prediction, vacuum-stability claim or all-model scan
is performed at Stage 30.


## Stage 30 native class-level diagnostic refinement

`BetaTerm[Quartic,1]` can print Wolfram messages even if its returned
expression is a valid association. `Check[expr,$Failed]` rejects ANY message
and previously masked whether the native projector succeeded. The exporter
now calls `CheckAbort[BetaTerm[Quartic,1],$Aborted]`, prints the native
result head and message names if it fails, and still rejects results other
than a complete projected association. In particular an actual
`ProjectionToUnmixedBetas::singular` / `Abort[]` is NOT bypassed.

This change is diagnostic and fail-closed; a live successful exporter run
and numerical trajectory are still required.

</details>

---

<a id="source-general-t3-conventional-uv-handoff-md"></a>

## Original document: `docs/GENERAL_T3_CONVENTIONAL_UV_HANDOFF.md`

**Original SHA-256:** `392d00081a99ea06874a353ec66bf9fe86c0fe254318c03dff61774fed9e53be`

**Status:** Historical record; check the current source and main
documentation before executing commands reproduced below.

<details>
<summary>Expand original: General T3: Convention-labelled UV quartic benchmark handoff</summary>

# General T3: Convention-labelled UV quartic benchmark handoff

## Scope

Produce a **numerical initial-condition candidate**, **not** RGE evolution, for all
50 scalar-quartic slots (38 new + 12 registered structural couplings) across the
16 ordinary T3 models. The 46 exact conditional signed matrices from Stage 23
are used only after verifying the Stage 26–28 evidence. Nothing is written into
production `Numerical/`, model configuration files, or Wolfram sources.

Current source ownership was inspected on GitHub before creating this script:

- `RGE/general/generated/general_rgbeta_quartic_registry.json`: 38 new vertices,
  12 reused structural slots, 50 total, 16 models.
- `Numerical/core/State.py`: canonical `T3UVState` dataclass.
- `Numerical/core/StateVector.py`: real ODE packing for that state.
- `Numerical/core/RGBetaEvaluator.py`: numerical evaluator with explicit
  UV state symbol environment.
- `docs/PHYSICS.md`: canonical chosen UV scale of `1e7 GeV`,
  `16*pi**2*dX/dln(mu)=beta_X^(1)` and Class-E `lambdaT3` benchmark-only 1/2.

The present numerical ODE does **not** include the 38 newly generated quartics.
Stage 29 inspects local source with Python's AST to record coverage for each
layer without importing or changing the solver.

## Mathematical transformation

For each block the Stage 23 tensor geometry is

\[
R_j=\sum_i M_i A_{ij}.
\]

The conditional Lagrangian sign gives

\[
\lambda_R=-A^{-1}\lambda_M,
\qquad
A\lambda_R+\lambda_M=0.
\]

Stage 29 checks this identity with **exact SymPy rational/algebraic numbers**
for every block before creating floating-point candidate values. It chooses

\[
\lambda_{M,i}(\mu_{UV})=0.01
\]

for **each of the 50 additional scalar-quartic slots** by default. This is a
new explicit *illustrative initial-condition choice*, **not** the historical
production benchmark and **not** a fitted or measured value. For the familiar
B/C `110101` two-dimensional example,

\[
A=\begin{pmatrix}\tfrac12&1\\-\tfrac{\sqrt3}2&0\end{pmatrix},\quad
\lambda_R=\begin{pmatrix}\tfrac{2\sqrt3}3\\-1-\tfrac{\sqrt3}3\end{pmatrix}0.01.
\]

Both the exact expressions and floats are included. The `lambdaT3` topology
coupling is **not part of this 50-slot initialisation**. The Class-E factor 1/2
used for that topology coupling is **not** applied to the additional quartics.

The triplet convention is labelled explicitly as the *Stage 27 Matchete CG
reference basis, provisionally identified with abstract RGBeta `S2`*.
Stage 28 did not independently export a numerical RGBeta component basis or
physical barred-triplet intertwiner, so no physical map is certified.

## Run (PowerShell, repository root)

```powershell
Expand-Archive -LiteralPath .\T3_General_Scalar_Stage29_ConventionUVHandoff.zip `
    -DestinationPath . -Force

python -m pytest -q tests\python\test_general_t3_convention_uv_handoff.py

python -m scripts.UVConventionHandoff --extra-quartic 0.01
```

Expected terminal summary, **conditional on current local evidence passing**:

```
Stage29 convention-labelled numerical UV values: 16/16 models, 50/50 slots, 46/46 blocks
New slots needing production one-loop state/beta support: 38
Largest |lambda_RGBeta| at chosen UV point: <computed number>
Actual RGE integration performed: NO
Physical Matchete/RGBeta coupling mapping: UNCERTIFIED
```

Report written to
`output/general_t3_stage22/stage29_convention_uv_handoff.json`.

The default invocation re-runs Stage 23 and Stage 28 checks, with Stage 28
recursively checking Stages 27 -> 26 -> 25. The diagnostic compares the current
saved Stage 23/28 outputs and hashes all source/evidence inputs to prevent use
of stale mapping data. It must be run with the regenerated local Stage 22–28
JSON evidence, including the successful live Wolfram exports.

For a different **illustrative** coefficient, use for example
`--extra-quartic 0.001`. The value must be a finite decimal in `(0,1]`.

## Output structure

The single JSON file contains 16 model entries, each with:

- `matchete_quartics`: exact and numerical chosen input parameters;
- `rgbeta_quartics`: exact and numerical conditionally mapped values, classified
  as `new` or `reused` in the current general registry;
- `blocks`: the exact signed matrix used and an exact zero residual;
- `production_ode_symbol_coverage`: current local state/vector/evaluator
  coverage of the 38 new and 12 reused slots;
- SHA-256 fingerprints for source and upstream evidence.

`physical_coupling_mapping_certified=false`,
`one_loop_beta_functions_for_new_quartics_evaluated=false`, and
`rge_evolution_performed=false` are deliberate. These are **not** solver inputs
for the current production code until that solver has been explicitly extended.

## Next numerical work

Export one-loop RGBeta beta expressions for the 38 new quartic symbols for each
relevant model and extend the numerical UV state/vector/evaluator to include
those parameters, alongside the 12 reused couplings. Check complete beta-key
coverage and pack/unpack consistency before integrating from `1e7` to `1e5`
GeV. Keep these boundary values labelled convention-dependent until the global
barred-field component identification is independently resolved. Do not silently
set missing beta functions to zero.

</details>

---

<a id="source-general-t3-fermion-threshold-handoff-md"></a>

## Original document: `docs/GENERAL_T3_FERMION_THRESHOLD_HANDOFF.md`

**Original SHA-256:** `b47f65f968615b9225aa02080b39b429d88eb5ca2bff7a882028552b7fb9337c`

**Status:** Historical record; check the current source and main
documentation before executing commands reproduced below.

<details>
<summary>Expand original: General T3: Coupled UV → fermion-threshold renormalisable handoff</summary>

# General T3: Coupled UV → fermion-threshold renormalisable handoff

## Scope and verified source interfaces

- Stage 30: `scripts/CoupledA0UVRunning.py` runs `A_+0` with five new complex quartic coordinates.
- Stage 31: `scripts/CoupledUVRunning.py` runs the remaining 15 ordinary models, including B/C reused quartics and zero-extra cases.
- Both output **compact** RGE trajectories, not the complete evolved UV state; the compact JSON is **insufficient to reconstruct** `y1`, `y2`, `MF`, scalar quadratic masses or the complete SM Yukawa matrices at the threshold.
- Production `Numerical.running.FermionThresholdBoundary.fermion_threshold_masses` computes heavy-mass matrix singular values.
- Production `Numerical.running.FermionThresholdBoundary.build_intermediate_scalar_boundary` removes the heavy-fermion renormalisable coordinates and projects the conventional UV state onto the existing scalar-only intermediate state. This explicitly introduces **no finite renormalisable threshold correction** and performs **no Wilson coefficient matching**.

Stage 33 uses these existing interfaces. It does not change any Stage 30/31, production RGE or matching code.

## Implementation

For each model:

1. Validate the existing Stage 31 combined summary, its hash of the original numerical report, and all recorded SHA-256 source digests.
2. Reconstruct UV initial conditions from the same Stage 29 handoff, registry and frozen fitted-reference input as the successful original run.
3. Re-evaluate the same native projected one-loop RGBeta beta expressions through the same existing Python numerical evaluator, `DOP853` integration and saved-point grid.
4. Compare all compact saved points (`mu`, `gY`, `g2`, `g3`, `lambdaH`, all new quartics and all reused quartics) against the replay; fail closed if their values or the function-evaluation count differ.
5. Unpack and validate the **full evolved canonical UV state** at `mu_F = 1e5 GeV`.
6. Compute singular heavy fermion masses from `MF`; use the **existing** fermion boundary projection to construct the renormalisable `SM+S1+S2` state.
7. Preserve the five/six newly generated quartics in a separate explicit **unmatched** sidecar. Reused B/C quartics are already in the canonical EFT state and are numerically cross-checked through the projection.
8. Write complete JSON serialisations of the UV state and intermediate scalar state including complex Yukawa/heavy-mass matrices (`real`/`imag` component matrices).

**Important:** This is a **renormalisable state handoff**, not a full EFT matching calculation. The newly generated quartics are **not** mapped to or evolved in the existing intermediate EFT; this requires a separately verified intermediate-model extension. The Wilson coefficients (including `C5`, `C12`) are **not** matched here, nor is any low-energy neutrino result computed.

The Matchete–RGBeta cross-package quartic map remains **UNCERTIFIED** and the Stage 29 extra quartic choices are illustrative rather than fitted.

## Installation

Extract the supplied ZIP into the repository root, preserving paths:

- `scripts/FermionThresholdHandoff.py`
- `tests/python/test_general_t3_fermion_threshold_handoff.py`
- `docs/GENERAL_T3_FERMION_THRESHOLD_HANDOFF.md`

## Commands (PowerShell, from repository root)

```powershell
Expand-Archive -LiteralPath .\T3_General_Scalar_Stage33_FermionThresholdHandoff.zip -DestinationPath . -Force
python -m pytest -q tests\python\test_general_t3_fermion_threshold_handoff.py
python -m scripts.FermionThresholdHandoff --model A_+0
python -m scripts.FermionThresholdHandoff --all
```

The Stage 29 handoff, Stage 30 report and exports, Stage 31 reports and exports, and Stage 31 16-model summary **must already be present** in `output/general_t3_stage22`. The replay checks their original source hashes. If a source has changed, the script exits with an explicit failure rather than silently using a different integration.

## Outputs

- `output/general_t3_stage22/stage33_A_p0_fermion_threshold_handoff.json`
- `output/general_t3_stage22/stage33_A_m2_fermion_threshold_handoff.json` (analogous for all models)
- `output/general_t3_stage22/stage33_16_model_fermion_threshold_summary.json` after all 16 succeed

Each model JSON includes the full evolved `uv_threshold_state`, `intermediate_scalar_boundary`, heavy fermion singular masses, unmatched new-quartic sidecar, source SHA-256 digests and explicit physics limitations.

## Verification boundary

- Offline: Stage 30 + Stage 31 + Stage 33 Python regression tests pass in a minimal isolated environment; one test requiring the full registry is skipped there.
- Live Windows: Stage 33 replay and threshold projection are **not yet verified** until executed in the user's local repository with its RGBeta exports.
- No GitHub changes are performed.

</details>

---

<a id="source-general-t3-intermediate-eft-running-md"></a>

## Original document: `docs/GENERAL_T3_INTERMEDIATE_EFT_RUNNING.md`

**Original SHA-256:** `b92041ae75f0d64729dd2ff0933276548b22eb488debf80b1b36eaf070c38fe6`

**Status:** Historical record; check the current source and main
documentation before executing commands reproduced below.

<details>
<summary>Expand original: General T3: supported scalar-only intermediate EFT running</summary>

# General T3: supported scalar-only intermediate EFT running

## Scope

Follows the **successful** Stage 33 renormalisable threshold projections for all
16 ordinary T3 models. Computes the existing scalar-only SM+S1+S2 one-loop RGE
subsystem from `mu_F = 1e5 GeV` to `mu_S = 1e3 GeV`, using the full independent
quartic-class `RGBeta` projector. No modifications to the production numerical
models, ODE solvers, RGE definitions, or GitHub.

## Files (relative to repository root)

- `scripts/ExportGeneralT3EFT1OneLoop.wl`: independent EFT1 full
  quartic-class projected `report_betas` export, with standard RGBeta metadata.
- `scripts/IntermediateEFTRunning.py`: verifies Stage 33 handoff and
  SHA-256 lineage, deserialises the full scalar boundary, checks EFT1 inventory,
  evolves with `Numerical.running.IntermediateScalarRunner.run_intermediate_scalar_segment`,
  and writes one model report or the full summary.
- `scripts/RunGeneralT3IntermediateEFTRunning.ps1`: orchestrates Wolfram and Python per model.
- `tests/python/test_general_t3_intermediate_eft_running.py`: regression
  tests, including metadata and sidecar verification.

## Inputs

- `output/general_t3_stage22/stage33_<model>_fermion_threshold_handoff.json`
- `output/general_t3_stage22/stage33_16_model_fermion_threshold_summary.json`
- `output/general_t3_stage22/stage31_16_model_uv_summary.json`
- Original Stage 30/31 UV running report for each model
- Generated `stage34_<model>_eft1_full_betas.json`

The Stage 33 original UV report hash, Stage 31 summary hash, and Stage 33
combined-inventory assertions must agree. The EFT1 model representation,
`IntegratedField="F"`, active fields, beta convention, and **exact** ODE
coupling-key inventory must agree before the integration starts.

## Running (Windows PowerShell, repository root)

```powershell
Expand-Archive -LiteralPath .\T3_General_Scalar_Stage34_IntermediateEFTRunning.zip -DestinationPath . -Force
python -m pytest -q tests\python\test_general_t3_intermediate_eft_running.py
powershell.exe -NoProfile -ExecutionPolicy Bypass -File .\scripts\RunGeneralT3IntermediateEFTRunning.ps1 -Models 'A_+0'
powershell.exe -NoProfile -ExecutionPolicy Bypass -File .\scripts\RunGeneralT3IntermediateEFTRunning.ps1
```

## Outputs

- `output/general_t3_stage22/stage34_<model>_eft1_full_betas.json`
- `output/general_t3_stage22/stage34_<model>_intermediate_eft_running.json`
- `output/general_t3_stage22/stage34_16_model_intermediate_eft_summary.json`

The model report includes **full initial and final scalar-only numerical
states** (including complex Yukawa matrices), 24 logarithmically spaced
intermediate trajectory samples, scalar mass diagnostics at `mu_S`,
solver counters, and input hashes.

## Verified conventions and limitations

The scalar-only EFT state drops `F`, `y1`, `y2`, and `MF` at the common threshold.
The ordinary quartics, including the 12 reused B/C slots, participate in
running using the existing EFT1 theory. Stage 29's **38 new** quartics have
no representation in the existing intermediate state or native EFT1 beta
export. Their mu_F values are retained **only** in a separately labelled
sidecar; Stage 34 neither evolves nor matches them, and therefore describes
**supported/truncated renormalisable EFT running**, not complete general
T3 scalar quartic EFT evolution.

No one-loop finite renormalisable matching, no C5/C12 Wilson-coefficient matching
or running, no physical Matchete-RGBeta component-basis certificate, and no
neutrino-mass/oscillation calculation are asserted. The Stage 29 quartic UV
values remain illustrative, not refitted physical parameters. Fermion matching
is at the common reference `1e5 GeV`, though individual singular masses differ.

If the native projected EFT1 beta system is singular, incompatible with the
existing state, or produces divergent/non-finite values, the run fails closed.

**Verification boundary**: offline Python tests do not certify successful
live Wolfram exports and full 16-model ODE integration; those require running
this package in the user's local repository.

</details>

---

<a id="source-general-t3-quadratic-batch-notes-md"></a>

## Original document: `docs/GENERAL_T3_QUADRATIC_BATCH_NOTES.md`

**Original SHA-256:** `bd44aaa9fd2eb04a84b706c9c32f541230d69686c5e61b46c0666ad671a1c789`

**Status:** Historical record; check the current source and main
documentation before executing commands reproduced below.

<details>
<summary>Expand original: T3 general-scalar batch: eight quadratic operators and one-loop mass beta functions</summary>

# T3 general-scalar batch: eight quadratic operators and one-loop mass beta functions

## Scope / source checkpoint

Read-only source reference: `Frankyang2002/T3neutrinomassmodelpipeline`, public `main`, checked 2026-10-27, including:

- `docs/general_scalar_trial_operator_manifest.json`
- `Lagrangian/model/T3Fields.wl`, `Lagrangian/model/T3ModelCatalog.wl`
- `Lagrangian/interactions/GeneralScalarOperators.wl`, `GeneralScalarVerifiedQuartics.wl`
- `RGE/running/rgbeta/T3RGBetaModel.wl`, `RunGeneralT3RGBeta.wl`
- `Numerical/experimental/GeneralT3ScalarRGState.py`, `GeneralT3ScalarQuadraticMasses.py`
- `scripts/FermionThresholdHandoff.py`, `IntermediateEFTRunning.py`
- 16-model Stage29 registry and the Stage22–29 physical tensor-map documentation.

The user separately reports a successful 16-model **full quartic EFT1** computation. This batch does not rerun/rebuild the original matching or numerical infrastructure, and it does not overwrite Stage34–36.

## Verified independent operator inventory

For `Y(S1)=alpha/2`, `Y(S2)=(alpha+2)/2`, with complex Z2-odd scalars, the only additional dimension-two gauge singlets in the 16 ordinary models are:

| Model | Holomorphic bilinear, plus H.c. | Candidate mass symbol | Basis interpretation |
|---|---|---|---|
| `A_+0` | `(S1†)^2` | `muT3GeneralS1Sq` | singlet split |
| `A_-2` | `S2† · S2†` | `muT3GeneralS2Sq` | triplet split |
| `B_-1` | `eps(S1†,S2†)` | `muT3GeneralS1S2` | mixed doublets |
| `C_-1` | `eps(S1†,S2†)` | `muT3GeneralS1S2` | mixed doublets |
| `D_+0` | `S1† · S1†` | `muT3GeneralS1Sq` | triplet split |
| `D_-2` | `(S2†)^2` | `muT3GeneralS2Sq` | singlet split |
| `E_+0` | `S1† · S1†` | `muT3GeneralS1Sq` | triplet split |
| `E_-2` | `S2† · S2†` | `muT3GeneralS2Sq` | triplet split |

No extra bilinear appears in the other eight assignments.  **There are 8 model-specific complex slots, not 8 universal couplings in each model.**  In the real gauge-unbroken scalar component space, the bilinear contributes a symmetric holomorphic mass matrix `B`. The two B/C terms may be diagonalized after rewriting the conjugate scalar to the same electroweak representation. That is only a *basis choice*: Yukawa/quartic tensors must be transformed too, so omitting the interactions from the original declared field basis is not justified. The other six coefficients can have a removable phase, but their magnitude produces a physical splitting. No claim of a full theory basis-redundancy proof is made.

## What this batch actually implements

1. **Independent quadratic gauge-invariant inventory**: `scripts/AuditGeneralT3QuadraticBasis.py`, checks the original manifest without using it to generate the candidate list. Fail-closed 8/16 test.
2. **Matchete audit, experimental**: `Lagrangian/RunGeneralQuadraticMatcheteAudit.wl`, builds the unchanged ordinary Lagrangian, all experimental quartics and the one allowed quadratic (if any), and validates individually and collectively using `CheckLagrangian`. `scripts/RunGeneralT3QuadraticMatcheteAudit.py` orchestrates 16 models. Requires live Matchete front end.
3. **RGBeta one-loop mass betas, experimental**: `RGE/running/rgbeta/RunGeneralT3QuadraticRGBeta.wl`, adds each model's quadratic directly via the documented `AddScalarMass[..., SelfConjugate -> False]` API, alongside all Stage29 general quartics; obtains `BetaTerm` for the **new complex quadratic coefficient and the two diagonal masses**, plus the full projected quartic beta class, for `UV` or `EFT1`. The script does not modify baseline models.
4. **EFT1 coupled numerical mass/quartic RG flow, experimental**: `scripts/RunGeneralT3QuadraticEFT1.py` reconstructs the hashed Stage33 threshold boundary and Stage29 sidecars; evolves `mS1Sq`, `mS2Sq`, 38 new / 12 reused quartic slots across models, and **the 8 new complex mass coefficients** wherever allowed, subject to successful live beta export. It rejects missing beta keys and unsupported InputForm expressions. Optional `--also-export-uv` exports UV mass betas but does not integrate a generalized UV mass boundary.
5. **Full replacement Hessian library**: `Numerical/experimental/GeneralT3ScalarQuadraticMasses.py`, preserves the existing function signatures and uses independent analytic spectra and the direct complex potential as checks. The chosen spherical SU(2) basis, including the neutral triplet bilinear metric, is explicit.
6. **Offline tests**: `tests/python/test_general_t3_quadratic_full_batch.py` tests 16 zero-limit models, all 8 complex quadratic sectors with three complex phases, spectra, Hessian contractions, tachyon flags and manifest checks. No Wolfram numerical result is claimed without a local live run.

## Physical convention and important limitation

The RGBeta mass term uses the **potential** convention (RGBeta's Lagrangian scalar potential enters with a minus sign):

`V_extra = mu2 * (Sa†)^T C (Sb†) + h.c.`

No hidden `1/2` is inserted. In this convention, self-conjugate-direction branches are `m² ± 2|mu2|`, and the B/C two-doublet squared-mass branches are

`(m1²+m2²)/2 ± sqrt(((m1²-m2²)/2)² + |mu2|²)`.

**This does not certify the Matchete Lagrangian source sign or its triplet/Bar component convention.** The selected spherical SU(2) tensors are mathematical conventions. The Stage22–29 quartic mapping remains conditional. Matchete mass-parameter basis, mass eigenstate threshold matching, and the Wilson coefficients for the new interactions remain unresolved. RGBeta's API confirms `AddScalarMass` can declare a complex mass coupling, but its live execution and one-loop output must still be verified on the user's Wolfram installation.

## Installation

Extract ZIP preserving paths into:

`C:\thesistemp\T3neutrinomassmodelpipeline`

No existing production file is replaced. The *existing experimental* `Numerical/experimental/GeneralT3ScalarQuadraticMasses.py` **is intentionally replaced in full** (same exported function names, additional checks).

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath 'C:\path\to\T3_Quadratic_Generalisation_Batch.zip' -DestinationPath . -Force
python -m pytest -q tests\python\test_general_t3_quadratic_full_batch.py
python -m scripts.AuditGeneralT3QuadraticBasis
```

First *live RGBeta mass-beta + numerical EFT1* test:

```powershell
python -m scripts.RunGeneralT3QuadraticEFT1 --model B_-1 --mu2-fraction 0.01
```

Then all 16, optional UV mass beta exports:

```powershell
python -m scripts.RunGeneralT3QuadraticEFT1 --all --mu2-fraction 0.01 --also-export-uv
```

Independent *live Matchete* candidate-Lagrangian test:

```powershell
python -m scripts.RunGeneralT3QuadraticMatcheteAudit --model B_-1
python -m scripts.RunGeneralT3QuadraticMatcheteAudit --all
```

Output:

- `output/general_t3_quadratic_eft1/quadratic_basis_audit.json`
- `output/general_t3_quadratic_eft1/B_m1_eft1_quadratic_live_rgbeta.json`
- `output/general_t3_quadratic_eft1/B_m1_eft1_quadratic_running.json`
- `output/general_t3_quadratic_eft1/B_m1_uv_quadratic_live_rgbeta.json` (only with `--also-export-uv`)
- `output/general_t3_quadratic_eft1/all_models_quadratic_eft1_summary.json` (only after successful `--all`)
- `output/general_t3_quadratic_eft1/matchete/...` (only with Matchete audit)

The default `--mu2-fraction 0` is a **zero-extra-quadratic control**. With `--mu2-fraction 0.01`, the starting extra coupling in a model that permits one is `0.01 * min(mS1Sq,mS2Sq) * exp(i*phase)`; the other eight models have no such parameter. `--mu2-phase` is in radians.

## Verification gates and next physics work

- OFFLINE: syntax, 16-model operator count and Hessian spectra/direct complex potential checks.
- LIVE: Matchete `CheckLagrangian` must accept all eight alongside the 50 candidate quartic slots; RGBeta must return complete UV/EFT1 beta inventories for all 16 models; the extended solver must succeed across all 16 without missing keys.
- UNRESOLVED: independent Matchete/RGBeta physical real-field intertwiner, exact mixed-mass/barred-triplet normalization, generalized UV mass RG boundary, diagonalization at the scalar threshold, and new one-loop matching/Wilson mixing.

**A numerically successful EFT1 mass/quartic run is NOT proof of a physically complete general T3 neutrino mass model.** Production Stage34–36 results are not overwritten or reinterpreted.

</details>

---

<a id="source-general-t3-scalar-threshold-handoff-md"></a>

## Original document: `docs/GENERAL_T3_SCALAR_THRESHOLD_HANDOFF.md`

**Original SHA-256:** `31e512b0cfe354a49fe92db1fb40b73fb893df80ab2cfcd629d9671f6ff8a357`

**Status:** Historical record; check the current source and main
documentation before executing commands reproduced below.

<details>
<summary>Expand original: General T3: scalar-threshold renormalisable handoff</summary>

# General T3: scalar-threshold renormalisable handoff

## Scope

Read the **individual** successful Stage 34 scalar-only EFT1 trajectories and their 16-model summary; reconstruct the complete ordinary-T3 state at the scalar matching scale

\[
\mu_S=10^3\,\mathrm{GeV}.
\]

Use the **existing** `Numerical.running.ScalarThresholdBoundary` machinery, without changing any production code:

- `scalar_threshold_masses(state)` calculates the positive running scalar mass estimates at \(\mu_S\), from the two scalar mass-squared parameters.
- `build_final_sm_boundary(state, matching_scale_gev=1e3)` constructs the canonical renormalisable SM state by removing the extra scalar fields **without finite renormalisable threshold corrections**.

No one-loop ODE is run again in Stage 35. The underlying Stage 34 calculation has already reached the scalar threshold.

## Inputs and checks

Expected files, in `output/general_t3_stage22/`:

- `stage34_16_model_intermediate_eft_summary.json`
- `stage34_<model_tag>_intermediate_eft_running.json` for each of the 16 models.

Stage 35 verifies the individual report hash recorded in the summary, plus the individual report's Stage 34 input SHA-256 inventory. It verifies the Stage 34 status/metadata, 16-model inventory, the initial/final Higgs quartic summaries, the new/reused quartic counts, and that the saved trajectories start/end at the reported scales. The initial and final **all-coupling compact snapshots** are compared to reconstructed canonical states. The two scalar mass diagnostics are independently recalculated.

`IntermediateEFTRunning.decode_intermediate_state` requires the input scale \(\mu_F\). To avoid changing that established API, Stage 35 validates the final scale \(\mu_S\), uses the decoder on a temporary representation with \(\mu_F\) **only in scale metadata**, and sets the reconstructed state back to \(\mu_S\) with `dataclasses.replace`. No coupling value is altered by this conversion.

Finally, all seven SM coupling blocks (`gY`, `g2`, `g3`, `lambdaH`, `yu`, `yd`, `ye`) must be exactly unchanged by the threshold projection.

## Physics boundary

This is a **supported/truncated renormalisable** continuation:

- 12 reused quartic slots were evolved by the existing Stage 34 EFT1 equations and are then removed with the scalars at \(\mu_S\).
- The 38 newly generated Stage 29 quartic slots were **not** matched into / evolved by the existing EFT1 system. Their values exist only in a \(\mu_F\) sidecar; Stage 35 does not invent values at \(\mu_S\).
- Finite renormalisable threshold corrections, scalar Wilson-coefficient matching (including \(C_5\) and \(C_{12}\)), the final SM+Weinberg RGE and neutrino-mass predictions are **not** performed.
- The physical Matchete–RGBeta coupling map is not certified.

## Install

Unpack the delivered archive **in the project root**, replacing no existing production files:

```powershell
Expand-Archive -LiteralPath .\T3_General_Scalar_Stage35_ScalarThresholdHandoff.zip -DestinationPath . -Force
```

Test:

```powershell
python -m pytest -q tests\python\test_general_t3_scalar_threshold_handoff.py
```

Single model:

```powershell
python -m scripts.ScalarThresholdHandoff --model A_+0
```

All models and summary:

```powershell
python -m scripts.ScalarThresholdHandoff --all
```

Re-aggregate only:

```powershell
python -m scripts.ScalarThresholdHandoff --summary
```

## Outputs

- `output/general_t3_stage22/stage35_A_p0_scalar_threshold_handoff.json` (and the analogous 15 model files)
- `output/general_t3_stage22/stage35_16_model_scalar_threshold_summary.json`

Each per-model report contains the complete canonical intermediate scalar state at \(\mu_S\), the projected canonical SM boundary at \(\mu_S\), two running scalar mass estimates, unprojected Stage 29 sidecar at \(\mu_F\), new/reused quartic inventory, upstream SHA-256 records, and explicit scope flags.

Stage 35 completion **does not** mean that Weinberg-operator matching or low-energy neutrino observables are available. Those require the pre-existing matched Wilson-coefficient input and its own checks.

</details>

---

<a id="source-general-t3-trial-md"></a>

## Original document: `docs/GENERAL_T3_TRIAL.md`

**Original SHA-256:** `722a1d5846fdd2b0f795491b802bf0df9506a27543ae2e7b2a68aab2e78f5053`

**Status:** Historical record; check the current source and main
documentation before executing commands reproduced below.

<details>
<summary>Expand original: General T3 scalar interaction trial — batched validation package</summary>

# General T3 scalar interaction trial — batched validation package

**Source snapshot:** `Frankyang2002/T3neutrinomassmodelpipeline` main
`41194fcd78fc7fe3c9378c9f8ee41228ad6473e7` (24 Oct 2026).

**Scope: opt-in experimental extension. This is NOT yet a verified fully
interaction-general Matchete/RGBeta/numerical pipeline.** The working model's
ordinary interactions remain unchanged unless specifically enabled. The ZIP
only edits local files, never GitHub.

## What the user's first run proved, and what it did not

The first run of `RunGeneralScalarAudit.wl` generated all **seven** additional
complex B(-1) quartic candidates. `FrontEndObject::notavail` occurred at every
`CheckLagrangian` call, and every check returned `$Failed`. The previous
script called `UsingFrontEnd` only for package loading, then performed the
validation outside that scope. This was an **audit-runner execution bug**;
zero accepted candidates was not evidence that their SU(2) contractions were
invalid. A second issue was that even the combined baseline check happened
outside that scope, so it also failed. The corrected runner encloses its **full
field/model construction and all validation inside**
`UsingFrontEnd[T3TrialAuditRun[trialModel]]`, matching the pattern of the
normal `RunModel.wl` builder.

Wolfram's frontend must actually be available. If the service frontend cannot
start, the new audit reports `Status: FrontendUnavailable` and returns null
accepted/valid counts, *not* seven invalid operators. The `BaselineLagrangianValid`
check runs before new candidates are tested. The final combined check always
includes every extra quartic, even if individually rejected.

## Implemented (new/revised relative to first ZIP)

- **Full replacement** `Lagrangian/RunGeneralScalarAudit.wl` — corrected
  frontend scope; baseline-first and full combined checks; machine-readable
  status and per-candidate outcomes.
- **Full replacement** `Lagrangian/interactions/GeneralScalarOperators.wl`
  — SU(2) invariant enumeration and symmetrisation, hypercharge and Z2 selection,
  complex Matchete couplings, and **fail-fast** behavior when an invariant or
  coupling definition fails (previously a structure could be dropped silently).
- **Full replacement** `Lagrangian/interactions/T3Topology.wl` — same opt-in
  `T3_EXPERIMENTAL_GENERAL_QUARTICS=1` production hook from first package.
- `scripts/RunGeneralScalarTrial.ps1` — drives an audit of B(-1), B/C(-1),
  or all 16 ordinary model points. Can separately run baseline and extended
  Matchete matching into isolated directories, with logs and JSON status.
- `scripts/ReportGeneralScalarTrial.py` — checks observed audit JSON against
  the 16-model independently derived operator manifest.
- `scripts/InspectGeneralScalarMatching.py` — compares *symbolic textual*
  Weinberg artifacts from baseline and extended runs, reports occurrences
  of `lambdaT3Extra...` couplings. It **does not** prove numerical equality
  or physical completeness.
- `Numerical/experimental/ScalarQuadraticPrototype.py` — isolated
  diagonalisation of the B/C(-1) complex-doublet mass matrix and splitting
  of Y=0 complex singlet/triplet into two real multiplets. This illustrates
  the quadratic-mixing extension **without wiring invalid non-diagonal
  propagators into the present matching**.
- Static operator manifest, Python operator counter, standalone unit tests.

**Predicted scalar-sector additions**, counting independent couplings per
separate supported ordinary model: 50 complex quartic couplings over 16 model
points and 8 complex quadratic terms over the same set. These are not 50
couplings of a single model. B/C(-1) each need seven new complex quartic
contractions relative to the Matchete builder. Counts use the original
complex-scalar field definitions and the gauge plus Z2 symmetries.

## Explicitly not implemented

- The additional **quadratic terms in the live Matchete/RGBeta/numerical model**.
  This needs propagator/basis changes, mass-matrix RGEs, rotated Yukawas,
  thresholds and matching. They are catalogued and prototyped, not enabled.
- **RGBeta beta functions for the complete added quartic basis.** The
  B/C(-1) RGBeta model already has six special quartics in a different basis,
  but that does not mean the new Matchete couplings have been connected to
  them. No numeric state-vector or beta-file changes are included.
- Newly allowed fermion-sector terms at special quantum numbers, or a
  repaired physical shared-scalar Matchete matching treatment.
- Matchete/Mathematica tests in the build environment: Wolfram is not
  available here. All actual Matchete results require your machine.

## Installation (PowerShell, repository root)

If you have local uncommitted changes, save or commit them first. The ZIP
replaces `T3Topology.wl` and the previous trial files. Consider a local-only
branch before expanding the ZIP:

```powershell
# Optional: creates a local experiment branch, does not push to GitHub.
git switch -c experiment/general-t3
Expand-Archive -LiteralPath .\T3_General_Scalar_Batch2.zip -DestinationPath . -Force

python -m pytest -q tests\python\test_general_scalar_trial.py tests\python\test_scalar_quadratic_prototype.py
python tests\python\check_t3_scalar_quartics.py --output output\general_scalar_static_audit.json
```

## Run many models with one command

First run the corrected B/C(-1) Matchete validity checks, in separate
`output/general_t3_trial` subfolders, with full logs:

```powershell
powershell -ExecutionPolicy Bypass -File .\scripts\RunGeneralScalarTrial.ps1 -Scope BC -Phase Audit
python scripts\ReportGeneralScalarTrial.py --root output\general_t3_trial
```

If those pass, run **all 16** ordinary supported cases:

```powershell
powershell -ExecutionPolicy Bypass -File .\scripts\RunGeneralScalarTrial.ps1 -Scope All -Phase Audit
python scripts\ReportGeneralScalarTrial.py --root output\general_t3_trial --all --require-complete --output output\general_scalar_audit_report.json
```

Some models have zero extra couplings. They still must pass the baseline and
combined checks. A failed result exits nonzero and is retained as structured
JSON plus a log. `-Scope All -Phase Audit` runs each model in its own
`wolframscript` process so an individual failed session does not contaminate
other model registry states.

To test symbolic one-loop matching **separately** for B(-1), first ensure
its audit passed, then:

```powershell
powershell -ExecutionPolicy Bypass -File .\scripts\RunGeneralScalarTrial.ps1 -Scope B -Phase Matching
python scripts\InspectGeneralScalarMatching.py --root output\general_t3_trial --output output\general_matching_symbolic_report.json
```

This runs original (flag OFF) and extended (flag ON) matching in different
output folders and restores the previous environment flag afterward. Matching
can be expensive. To audit and then automatically match only passing model
points in one command, use `-Scope BC -Phase Both` or `-Scope All -Phase Both`.
**Caution: extended symbolic matching is not yet consistent with the present
RGBeta and numerical state. Do not feed these experimental matching artifacts
into production numerical results.**

## Quadratic matrix prototype only

```powershell
python -c "from Numerical.experimental.ScalarQuadraticPrototype import mixed_t3_bc_doublet_spectrum; s=mixed_t3_bc_doublet_spectrum(700**2,900**2,(200**2)+0j); print(s.mass_squared,s.mass)"
```

For B/C(-1) the gauge-basis matrix used is

\[
M_S^2 = \begin{pmatrix}m_{S_1}^2&-\mu_{12}^2\\
-(\mu_{12}^2)^*&m_{S_2}^2\end{pmatrix},
\quad \Phi=(\widetilde S_1,S_2)^T.
\]

This is an **unbroken-gauge-phase scalar spectrum prototype**. Quartic Higgs
VEV effects and genuine RG mass mixing are not included.

## How to interpret status

- `Success`: baseline and all individual extra operators and their combined
  Lagrangian passed Matchete `CheckLagrangian`.
- `FrontendUnavailable`: Mathematica could not provide a usable frontend;
  *no gauge-invariance conclusion follows*.
- `BaselineValidationFailed`: old/ordinary Lagrangian failed before any
  extra candidate was tested. Diagnose Matchete installation/baseline first.
- `ExtraQuarticGenerationFailed`: CG/tensor/coupling construction failed.
- `IndividualOperatorValidationFailed`: at least one new operator failed
  `CheckLagrangian`, names/results are in JSON.
- `CombinedLagrangianValidationFailed`: individual checks passed but the
  full combined expression did not.

## Reversion

With changes safely committed on an experiment branch, return to the previous
branch with `git switch <original-branch>`; no GitHub push is needed. If the
experiment files are uncommitted, use `git status --short` to inspect them
before reverting. In particular, `git restore --
Lagrangian/interactions/T3Topology.wl` discards **all** your uncommitted edits
to that file, not only edits from this ZIP.

## Batch 3: Windows PowerShell parser correction (26 October 2026)

The Batch 2 PowerShell orchestration file contained an interpolation error:
`"Skipping matching for $name: ..."` is not valid Windows PowerShell because
`$name:` is parsed as a scoped variable reference. Batch 3 replaces the
entire runner, using `${name}:`. This defect prevented *any* audits or
matching from starting, regardless of the previously passing Python tests.

Batch 3 also clarifies the no-audits-found result of
`ReportGeneralScalarTrial.py` and returns a nonzero exit code for that case.
The physics modules, symbolic operator catalogue, and mass prototype are
unchanged from Batch 2. This is a runner hotfix, not a newly verified
Matchete result.

Install the Batch 3 ZIP **over** Batch 2 from the repository root:

```powershell
Expand-Archive -LiteralPath .\T3_General_Scalar_Batch3_Fix.zip -DestinationPath . -Force
python -m pytest -q tests\python\test_general_scalar_trial.py tests\python\test_scalar_quadratic_prototype.py

# One large run: audit all 16 representations; if a point passes, compare
# baseline and extended symbolic one-loop matching for that point.
powershell -ExecutionPolicy Bypass -File .\scripts\RunGeneralScalarTrial.ps1 -Scope All -Phase Both

# Summarise all 16, including failed/unrun points.
python scripts\ReportGeneralScalarTrial.py --root output\general_t3_trial --all --output output\general_scalar_audit_report.json
python scripts\InspectGeneralScalarMatching.py --root output\general_t3_trial --output output\general_matching_symbolic_report.json
```

**Important:** `-Scope All -Phase Both` may initiate expensive one-loop symbolic
matching in up to 16 baseline and 16 extended runs, conditional on the audits
passing. If you prefer the faster physics validation first, choose
`-Scope BC -Phase Both` before running `-Scope All -Phase Both`.
A failed audit skips matching for that model and records its log.

## Batch 4: findings from the first successful 16-model Matchete run (26 October 2026)

User-provided Windows log showed **all 16 baseline and extended symbolic matching jobs exited with code 0** and all attempted extra quartics passed individual and combined Lagrangian checks. This is encouraging but not an independent end-to-end physics verification; tree matching in the sample emitted non-fatal `FrontEndObject::notavail` warnings.

The independent Python operator inventory and the observed Mathematica generator agreed for 12/16 models. The missing four were B/C at `alpha=-3,+1`, with one *valid* quartic invariant each:

- B/C(`-3`): `Hdag^2 S2dag^2`, counts `{0,2,0,0,0,2}`.
- B/C(`+1`): `Hdag^2 S1^2`, counts `{0,2,2,0,0,0}`.

This package now uses the independently enumerated **SU(2)-singlet verified monomial manifest** for the 16 production scalar representation points, and requires Matchete to reproduce the *number of independent contractions* of every listed structure. The original dynamic hypercharge enumeration is retained as an uncertified fallback for unsupported models. Crucially, a successful `CheckLagrangian` on zero generated operators does **not** count as completeness verification anymore. The direct cause of the previous four empty Mathematica lists was not reproduced locally (no Matchete installed here); the new approach avoids relying on those lists for supported cases.

The logs also showed *symbolic* `C5` dependence in B/C(`-1`) on precisely two new couplings:

- `lambdaT3ExtraQ020020I1` (`Hdag^2 S2^2`, up to Hermitian conjugation)
- `lambdaT3ExtraQ020200I1` (`Hdag^2 S1dag^2`, up to Hermitian conjugation)

Other tested points showed no newly named coupling in the **available extracted C5 artifacts**. This is **not** a general proof of zero effects in all EFT threshold configurations or at higher orders. No UV numerical input has yet assigned these new couplings nonzero values.

### One-command local retest of the four incomplete points

Install this full consolidated ZIP at your repository root and run:

```powershell
Expand-Archive -LiteralPath .\T3_General_Scalar_Batch4.zip -DestinationPath . -Force
python -m pytest -q tests\python\test_general_scalar_trial.py tests\python\test_scalar_quadratic_prototype.py
python scripts\GenerateGeneralScalarWLManifest.py --check
powershell -ExecutionPolicy Bypass -File .\scripts\RunGeneralScalarTrial.ps1 -Scope Problems -Phase Both
python scripts\ReportGeneralScalarTrial.py --root output\general_t3_trial --all --require-complete --output output\general_scalar_audit_report.json
python scripts\InspectGeneralScalarMatching.py --root output\general_t3_trial --require-complete --output output\general_matching_symbolic_report.json
```

`-Scope Problems` selects only B/C(-3,+1). The other 12 completed audit/matching artifacts in the same output folder are retained. A successful aggregate report should read **16 passed / 0 failed**, if the new contractions pass Mathematica validation. Run the full 16-model audit as a further check if desired:

```powershell
powershell -ExecutionPolicy Bypass -File .\scripts\RunGeneralScalarTrial.ps1 -Scope All -Phase Audit
```

The PowerShell runner now independently gates matching on the manifest count, even when the Mathematica output mistakenly says `Status=Success`. The symbolic report also checks the number of extra interactions named in the extended matching logs and records text equality of baseline/extended `c5_coefficient.txt` after whitespace removal. **Text equality is not an algebraic equivalence test.**

### Actual implementation boundary

This batch completes the **Matchete candidate generation plus count enforcement for the scalar quartic sector**, subject to the four new local checks. It does not claim to complete the physical UV theory. No active scalar quadratic-mixing propagators, missing Yukawa sectors, corresponding general RGBeta beta functions, RG-to-Matchete CG normalisation map, or numerical boundary/state/threshold propagation were implemented. These are substantial independent changes, and must not be treated as verified solely because symbolic matching succeeded.

## Batch 5: complete scalar-quartic symbolic matching has passed all 16 models

**User-supplied local run on 26 October 2026:** 16/16 Matchete operator
completeness audits passed; all original and extended one-loop matching runs
exited successfully. The text-level `C5` comparison reports new formal
couplings only in ordinary B/C at alpha=-1:

- `lambdaT3ExtraQ020020I1`: `(Hdag)^2 S2^2`, the same *field structure* as
  RGBeta's `lambdaHHdagS2S2`.
- `lambdaT3ExtraQ020200I1`: `(Hdag)^2 (S1dag)^2`, the same *field structure*
  as RGBeta's `lambdaHHdagS1barS1bar`.

The RGBeta structures are **not proven to have the same tensor normalisation**
(or phase conventions) as the dynamically normalised Matchete CGs. Do not
identify the numerical values of these pairs until this is checked.

The user's output also contains `FrontEndObject::notavail` messages during
some tree-matching steps, but every one-loop matching run completed and
extracted C5. This message alone does not invalidate a successful result;
formal regression tests are added below to verify the matching artifacts.

### Batch 5: ONE analysis command (reuses existing matching; no expensive reruns)

```powershell
powershell -ExecutionPolicy Bypass -File .\scripts\RunGeneralT3Analysis.ps1 -RequireProof
```

It runs four analyses: verifies all 16 operator counts, inspects textual C5
outputs, enumerates Matchete/RGBeta field-structure coverage, then uses Wolfram
to test formally that setting all **new** quartic couplings to zero in the
extended coefficient gives exactly the original coefficient.

The zero-extra check must either prove the identity or explicitly say the
symbolic reduction was not proven; a missing or timed-out result is not
considered a pass. If an identity cannot be shown, inspect its residual
rather than asserting it must vanish. This is an output-level perturbative
regression, **not** a proof that the extra couplings remain zero under RG flow.

Outputs relative to project root:

- `output/general_scalar_audit_report.json`
- `output/general_matching_symbolic_report.json`
- `output/general_rgbeta_coverage.json`
- `output/general_t3_trial/c5_algebraic_report.json`
- `output/general_t3_trial/<MODEL>/general_c5_algebraic/delta_c5_inputform.txt`
- `output/general_t3_trial/<MODEL>/general_c5_algebraic/delta_c5_atoms.json`
- `output/general_t3_physics_summary.json`
- `output/general_t3_physics_summary.csv`

The last two files are compact machine-readable and human-readable tables
for your thesis notes. All newly generated files go to `output/`.

### RGBeta coverage result and why numerical work remains separate

The 16-model manifest contains 50 **per-model** additional complex quartic
couplings. Counting gauge-invariant field structures, 12 of them have existing
B/C(-1) RGBeta counterparts (six in each of two models), and 38 do not.
Within each B/C(-1), one of the two independent mixed Higgs portal
contractions is absent in RGBeta.

This count is *structural* only; it does not verify CG tensor normalisation or
beta functions. `scripts/CheckGeneralT3RGBetaCoverage.py` explicitly labels
these distinctions and checks that the six known special RGBeta coupling
names are actually defined in the local source. This is the next needed
input to a fully general RGE state extension. The quadratic scalar mixing
terms are still prototypes, not a working generalized matching formalism.

## Batch 6: formal C5 response, Matchete tensor normalisation, and RG source terms

The complete 16/16 Matchete scalar-quartic audit and 16-model zero-extra C5
proof have passed **on the user's Windows environment**. Batch 6 does not
rerun those expensive matching calculations. It adds the next diagnostics;
**it is not a numerical neutrino-mass implementation**.

From the project root, after extracting the full Batch 6 ZIP:

```powershell
python -m pytest -q tests\python\test_general_scalar_trial.py tests\python\test_scalar_quadratic_prototype.py tests\python\test_general_scalar_c5_proof.py tests\python\test_general_rgbeta_coverage.py tests\python\test_general_t3_scalar_quartic_response_audit.py
powershell -ExecutionPolicy Bypass -File .\scripts\RunGeneralT3ScalarQuarticResponseAudit.ps1
```

This uses the **existing**
`output/general_t3_trial/c5_algebraic_report.json` and matching artifacts,
and the **existing** RGBeta JSON exports
`output/T3_B_alpha_m1/data/uv_rgbeta_rge.json` and
`output/T3_C_alpha_m1/data/uv_rgbeta_rge.json`.

Outputs:

- `output/general_c5_response.json`: validates the independent, linear
  coupling derivatives in the formal B/C(-1) \(\Delta C_5\); labels whether
  their coefficient expressions include poles and logarithms. A failure is
  never interpreted as a proven identity.
- `output/general_rgbeta_source_proof.json`: substitutes zero for all six
  special B/C(-1) RGBeta quartics and records the formal one-loop beta-function
  source for the two C5-active quartics. Demonstrates whether setting all
  these couplings to zero is RG-invariant.
- `output/general_t3_trial/matchete_active_quartic_tensors.json`: computes
  Matchete SU(2) tensors for the two C5-active quartics and compares them to
  a **defined** delta/epsilon symmetrised reference tensor. It reports the
  proportionality factor and an exact component-by-component residual.

**Normalisation limitation:** This last step determines only a Matchete-to-
canonical-index factor. In RGBeta the same operator is written using
`delS2` and RGBeta's own invariant conventions, and `AddQuartic` applies no
implicit combinatorial normalisation. A separate RGBeta-to-canonical check
and Lagrangian sign/factor check are required before identifying numerical
values of `lambdaT3ExtraQ020020I1` with `lambdaHHdagS2S2`, or
`lambdaT3ExtraQ020200I1` with `lambdaHHdagS1barS1bar`. Do **not** feed an
unverified transformation into the production neutrino-mass pipeline.

If the existing RGBeta output files are absent, run the script with `-SkipRGE`
and inspect the Matchete proof separately. If Matchete fails a canonical tensor
check, its JSON preserves the raw tensor to help diagnose SU(2) index
orientation rather than silently declaring a unit conversion.

The one-loop RGBeta normalisations and the missing 38 structurally independent
quartic couplings outside the existing RGBeta set remain unimplemented.
Quadratic scalar mass mixing remains catalogued but is not enabled in matching
or RG evolution. Shared-scalar matching is also a separate physical model.

## Stage 9 tensor exporter regression fix (26 October 2026)

See `docs/GENERAL_T3_MATCHETE_TENSOR_INITIALIZATION.md`. `LoadModel["SM"]` must run after
`ResetAll[]` in the standalone tensor exporter before BSM representations
are defined. For a targeted rerun only:

```powershell
powershell -ExecutionPolicy Bypass -File .\scripts\RunGeneralT3ScalarQuarticResponseAudit.ps1 -TensorOnly
```

This does not claim that RGBeta's scalar coupling convention has yet been mapped.

## Stage 10 (experimental numerical pilots)

See `docs/GENERAL_T3_MATCHING_RG_SOURCE_PILOTS.md` for the two independent numerical probes and their strict limitations. Run `scripts/RunGeneralT3MatchingRGSourcePilots.ps1` after the successful Stage 8/9 formal-C5, RGBeta and tensor diagnostics. Stage 10 adds no new production numerics or RGBeta-to-Matchete coupling identification.

## Stage 11 — finite scalar mass-ratio logarithms

See `docs/GENERAL_T3_C5_FINITE_MASS_LOGARITHMS.md`. The experimental Matchete hard-C5
pilot now retains finite physical mass-ratio logarithms after applying
`μbar2 = MF^2` (instead of refusing any surviving `Log`). The toy input now
supplies positive `MS1`, `MS2` in GeV. Unrecognized symbolic scales and
unsubtracted UV poles still block evaluation and generate a diagnostic JSON.
The PowerShell runner always attempts all four independent B/C probes before
reporting failure. No production RGE or matching code was changed.

## Stage 12 extra-only mass-matrix/phase diagnostics

After the successful Stage 11 numerical pilots, run:

```powershell
powershell -ExecutionPolicy Bypass -File .\scripts\RunGeneralT3ExtraC5PhaseResponse.ps1
```

This reuses the stored two-model Matchete `delta_c5` matrices and reports the additional contribution to `m_nu = -(v^2/2) C5`, the three extra-only Takagi singular values and a 16-point quartic phase scan. No further symbolic matching is performed. Full details and limitations: `docs/GENERAL_T3_EXTRA_C5_PHASE_RESPONSE.md`.

## Stage 13: coupled RGBeta UV numerical evolution and SU(2) tensor mapping

See `docs/GENERAL_T3_COUPLED_UV_RUNNING.md`. After the Stage 9 tensor audit and
Stage 8 RGBeta-source proof have succeeded, run:

```powershell
python -m pytest -q tests\python\test_general_t3_coupled_u_v_r_g_e.py
powershell -ExecutionPolicy Bypass -File .\scripts\RunGeneralT3CoupledUVRGE.ps1
```

This uses real one-loop coupled RGBeta beta functions for B/C(-1) **above the
heavy-fermion thresholds**, including the six complex scalar quartics already
registered in RGBeta. It is an isolated experiment, not a replacement for the
production hierarchical numerical neutrino-mass pipeline. Group tensor ratios
are checked, but no unverified equality of complex Matchete and RGBeta
couplings is asserted.


## Stage 14 (conditional UV-induced C5 sensitivity)

See `docs/GENERAL_T3_CONDITIONAL_UV_TO_C5.md`. After Stage 11 and Stage 13
outputs exist, use `python -m pytest -q tests\python\test_general_t3_conditional_u_v_to_c5.py`
then `powershell -ExecutionPolicy Bypass -File .\scripts\RunGeneralT3ConditionalUVToC5.ps1`.
This uses hypothetical sign/conjugation maps, **not** a proven Matchete/RGBeta
Lagrangian coupling conversion, and never modifies production files. One-loop
RGE-generated quartics inserted into one-loop C5 are formally of order hbar^2.

## Stage 15 — source-convention map for B/C(-1)

Stage 14's 16 hypothetical sign/conjugation scenarios are narrowed to the
source-derived convention \(\lambda_{\mathrm{Matchete}}=-\lambda_{\mathrm{RGBeta}}\)
for the two C5-active B/C(-1) quartics *in the same named scalar basis*.
This follows the RGBeta MBL quartic sign (Thomsen, EPJC 81 (2021) 408,
Eq. 2.3), the experimental Matchete positive `PlusHc` interactions, and the
previously verified +1 SU(2) tensor factors. The source checker fails closed
on field-word, tensor, or interaction-definition changes; the numerical
calculation must reproduce Stage 14's `negative/negative` case.

```powershell
python -m pytest -q tests\python\test_general_t3_quartic_convention_c5.py
powershell -ExecutionPolicy Bypass -File .\scripts\RunGeneralT3QuarticConventionC5.ps1
```

Results are written under `output/general_t3_stage15/`; no production modules
are modified. **Still not threshold-consistent matching or a complete two-loop
calculation**. See `docs/GENERAL_T3_QUARTIC_CONVENTION_C5.md`.

</details>

---
