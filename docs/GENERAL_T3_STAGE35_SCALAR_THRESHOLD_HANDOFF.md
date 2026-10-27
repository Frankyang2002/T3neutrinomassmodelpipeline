# Stage 35 — scalar-threshold renormalisable handoff

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

`Stage34IntermediateEFTRunning.decode_intermediate_state` requires the input scale \(\mu_F\). To avoid changing that established API, Stage 35 validates the final scale \(\mu_S\), uses the decoder on a temporary representation with \(\mu_F\) **only in scale metadata**, and sets the reconstructed state back to \(\mu_S\) with `dataclasses.replace`. No coupling value is altered by this conversion.

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
python -m pytest -q tests\python\test_general_t3_stage35_scalar_threshold_handoff.py
```

Single model:

```powershell
python -m scripts.Stage35ScalarThresholdHandoff --model A_+0
```

All models and summary:

```powershell
python -m scripts.Stage35ScalarThresholdHandoff --all
```

Re-aggregate only:

```powershell
python -m scripts.Stage35ScalarThresholdHandoff --summary
```

## Outputs

- `output/general_t3_stage22/stage35_A_p0_scalar_threshold_handoff.json` (and the analogous 15 model files)
- `output/general_t3_stage22/stage35_16_model_scalar_threshold_summary.json`

Each per-model report contains the complete canonical intermediate scalar state at \(\mu_S\), the projected canonical SM boundary at \(\mu_S\), two running scalar mass estimates, unprojected Stage 29 sidecar at \(\mu_F\), new/reused quartic inventory, upstream SHA-256 records, and explicit scope flags.

Stage 35 completion **does not** mean that Weinberg-operator matching or low-energy neutrino observables are available. Those require the pre-existing matched Wilson-coefficient input and its own checks.
