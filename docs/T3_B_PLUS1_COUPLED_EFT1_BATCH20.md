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
