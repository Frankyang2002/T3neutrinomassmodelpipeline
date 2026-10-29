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
