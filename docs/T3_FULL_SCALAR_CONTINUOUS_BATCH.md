# Continuous full-scalar candidate RG batch

**Repository:** `Frankyang2002/T3neutrinomassmodelpipeline` (local checkout only).  
**Source inspected:** GitHub `main`, commit `357e82cf2bacda023cc933e26c882895fd60e918` (29 October 2026). Check local differences before replacing any files.

## Files to add (no replacements)

- `scripts/RunT3FullScalarContinuous.py`
- `tests/test_t3_full_scalar_continuous.py`
- `docs/T3_FULL_SCALAR_CONTINUOUS_BATCH.md` (this file)

**Do not overwrite the root `README.md`.** The previous package accidentally did so.

## What was implemented

One continuous **candidate-basis** RG trajectory for each ordinary T3 model:

1. Read the Stage29 additional-quartic UV boundary at `muUV = 10^7 GeV`, the established fitted UV state, and the model-specific **full** RGBeta quadratic+quartic UV beta export.
2. Add one nonzero complex holomorphic quadratic coefficient where the symmetry allows it (8 of the 16 models). Set its UV magnitude to `--mu2-fraction * min(mS1Sq,mS2Sq)` with phase `--mu2-phase` (radians). The other eight models have no such coefficient.
3. Integrate *all ordinary and additional quartic beta functions, both diagonal scalar mass beta functions, Yukawa/gauge terms, and the allowed holomorphic mass beta together* from `10^7` to `10^5 GeV`. The holomorphic mass is also put in the evaluation environment for all other beta functions, so feedback is not discarded.
4. At `10^5 GeV`, use the existing renormalisable heavy-fermion threshold boundary for the **evolved UV state**, transport the additional quartics and holomorphic mass continuously, and deliberately omit uncomputed finite and Wilson threshold terms.
5. Integrate the full coupled scalar-only EFT1 system to `10^3 GeV`. Reconstruct all component scalar mass-squared eigenvalues and record any tachyonic endpoint without silently treating it as a physical spectrum.
6. For `B_-1` and `C_-1` with positive scalar masses, call the existing *simultaneous all-heavy one-scale* Matchete expression evaluator using the new EFT1 endpoint. Its frozen Yukawa and conditional two-quartic sign-map assumptions remain unchanged. This is not a sequential threshold result.

Every model gets its own resumable file. History under `output/general_t3_stage22/`, `output/general_t3_quadratic_eft1/`, and the reduced reference benchmarks is **read-only**. New results live under `output/studies/general_scalar_extended/continuous_full_scalar/<run-id>/`.

## Prerequisites

From the local repository root, ensure the existing saved artifacts are available:

- `output/general_t3_stage22/stage29_convention_uv_handoff.json`
- `output/general_t3_quadratic_eft1/<MODEL>_uv_quadratic_live_rgbeta.json`
- `output/general_t3_quadratic_eft1/<MODEL>_eft1_quadratic_live_rgbeta.json`

`<MODEL>` uses e.g. `B_m1`, `A_p0`. There should be **16 UV + 16 EFT1** quadratic-aware beta files. These are distinct from your previously validated `output/general_t3_stage20/*_rgbeta.json` files. The GitHub repository snapshot contains the quadratic-aware files, but their presence and integrity on the **local** machine must be checked.

For optional B/C C5, the previously used `scripts/RunT3QuadraticEndpointC5.py` and its Matchete one-scale inputs must be present. If these are absent, the scalar UV/EFT1 integration still produces outputs; C5 is recorded as blocked.

## Run in PowerShell

```powershell
# From the T3neutrinomassmodelpipeline root.
python -X utf8 -m pytest -q tests/test_t3_full_scalar_continuous.py

# All 16 models, including nonzero allowed quadratic coefficients:
python -m scripts.RunT3FullScalarContinuous --all --run-id continuous_v1 --mu2-fraction 0.01

# Resume after Ctrl+C (same command plus --resume; completed models are skipped):
python -m scripts.RunT3FullScalarContinuous --all --run-id continuous_v1 --mu2-fraction 0.01 --resume

# Or run one model to get started:
python -m scripts.RunT3FullScalarContinuous --model B_-1 --run-id b_m1_v1 --mu2-fraction 0.01
```

To run the 16-model RG **without** optional C5:

```powershell
python -m scripts.RunT3FullScalarContinuous --all --run-id continuous_rg_v1 --mu2-fraction 0.01 --skip-c5
```

A second run with different inputs **must use a new run ID**. Saved checkpoints use SHA-256 of source inputs and refuse resume after mismatched inputs. The script saves the summary after every completed model. `Ctrl+C` between model completion messages is the cleanest stopping point.

## Results

```
output/studies/general_scalar_extended/continuous_full_scalar/continuous_v1/
    summary.json
    A_p0/
        continuous_trajectory.json
        A_p0_eft1_quadratic_running.json
    B_m1/
        continuous_trajectory.json
        B_m1_eft1_quadratic_running.json
        one_scale_c5/B_m1_final_endpoint_one_scale.json
    ...
```

The summary records the number of successful models, scalar endpoint eigenvalues, and any B/C one-scale masses. The dedicated endpoint report records the **complete** scalar state including new quartics and a quadratic coefficient, with all status flags explicitly uncertified.

## Scientific limits — not implemented by this batch

- It **does not** replace the production Matchete full scalar Lagrangian or rewrite the original reduced-potential benchmarks.
- It **does not** compute one-loop fermion-threshold `C12`/`C5` matching, full scalar EFT Wilson running, or a physical sequential scalar-threshold `C5` for the general potential.
- It **does not** derive the mass-eigenfield Matchete one-loop coefficient for the six A/D/E self-bilinear model slots. A tree-level field rotation alone is insufficient to recover the correct neutrino C5.
- The cross-package RGBeta-to-Matchete component and coupling basis is not fully certified. For the optional B/C one-scale result only, the existing source-level two-quartic sign mapping is reused.
- The optional B/C one-scale C5 uses **frozen-reference Yukawa matrices**, not the UV-evolved Yukawas; accordingly its numerical masses are conditional scenarios, **not physical neutrino mass predictions**.
- The fermion-threshold treatment transports the *renormalisable* scalar couplings. Finite threshold corrections and Wilson matching remain unimplemented.

These are substantive remaining integration requirements. The script delivers a continuous 16-model *candidate one-loop scalar RG* calculation, **not** the complete physical full-scalar T3 neutrino-mass pipeline.

## Verification performed in preparation

- Python source/test syntax compiled successfully.
- 10 focused tests passed in an isolated dependency-stub harness. This verifies argument contracts, beta inventories, both-stage derivative assembly/backreaction, ODE integration, checkpointing, and path safety. It **does not** execute Wolfram/RGBeta, the full checkout, or any real T3 RG trajectory. Run the tests on the actual checkout and inspect `summary.json` before interpreting numerical outputs.
