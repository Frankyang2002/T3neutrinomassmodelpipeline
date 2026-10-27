# Stage 30: first coupled one-loop UV trajectory (A_+0)

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

python -m pytest -q tests\python\test_general_t3_stage30_coupled_a0_uv_running.py
if ($LASTEXITCODE -ne 0) { throw "Stage 30 tests failed" }

wolframscript -file .\scripts\ExportGeneralT3Stage30A0OneLoop.wl `
    output/general_t3_stage22/stage30_A_p0_full_uv_betas.json
if ($LASTEXITCODE -ne 0) { throw "Stage 30 RGBeta export failed" }

python -m scripts.Stage30CoupledA0UVRunning
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
