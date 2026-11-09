# Workstream A v28–v30 — ONE common scalar matching scale after all fermions

## Actual source implementation (not an inferred beta)

The code imports the existing `scripts.RunGeneralT3QuadraticEFT1` and `scripts.RunT3BCFullScalarMassEFT1` one-loop RGBeta numerical evaluators from the installed repository. It never fabricates a scalar beta coefficient. It reconstructs the 1e5 GeV stage33 candidate EFT1 renormalisable boundary, evolves **both** S1 and S2 until a single common scale, and retains their complete scalar-eigenstate basis.

- **`five_mass`**: requires `output/audits/workstream_b_bc_full_mass/{B,C}_m1_eft1_full_mass_v1.json` plus `output/general_t3_quadratic_eft1/{B,C}_m1_eft1_quadratic_live_rgbeta.json`, stage33 and the original registry. Uses the existing 5-real-mass coupled RGBeta beta evaluator (Higgs mH^2, mS1^2, mS2^2, complex mu2), plus all dimensionless beta functions.
- **`quadratic_only`**: fallback if full Higgs mass-betas are not exported. Requires the earlier `output/general_t3_quadratic_eft1/{B,C}_m1_eft1_quadratic_running.json` and corresponding native EFT1 beta exports. Replays its original initial boundary; **does not claim Higgs-mass running** and explicitly checks its stored 1 TeV quadratic endpoint.
- `auto` selects five_mass when a model's verified full mass beta file exists, otherwise quadratic_only. Backend choice is separately recorded for each model.

The self-consistent common scale solves

\[
 \log\mu_S-\tfrac14\left[\log m_X^2(\mu_S)+\log m_Y^2(\mu_S)\right]=0,
 \quad\mu_S=\sqrt{m_X(\mu_S)m_Y(\mu_S)}.
\]

Independent alternative: `upper` sets \(\mu_S=m_Y(\mu_S)\), **still integrating out both scalars simultaneously**, with no intermediate single-scalar EFT. Both roots are calculated and stored. The geometric prescription minimises the larger absolute mass logarithm at a fixed spectrum. The `muS` running state is taken from the *dense ODE solution*, not interpolated from a saved 1000 GeV endpoint. The 1000 GeV endpoint is replayed only for parity checks in the reduced backend.

The all-fermion-first requirement is respected as a *renormalisable boundary projection* from Stage33; that boundary did **not** match LLSS/C5 Wilson coefficients. These reports are not physical predictions, do not run Wilson operators, and do not claim any Matchete↔RGBeta normalization. The saved v21–v27 numerical B/C mass benchmark is **not** assumed to equal this new independent initial-data choice.

## Install paths

1. `Numerical/parallel_t3/wilson/bc_joint_scalar_scale_v28.py`
2. `Numerical/parallel_t3/wilson/bc_joint_scalar_scale_audit_v29.py`
3. `Numerical/parallel_t3/wilson/bc_fermion_first_common_handoff_v30.py`
4. `Numerical/parallel_t3/wilson/BC_JOINT_SCALAR_EFT1_V28_V30.md`
5. `scripts/RunT3ParallelABCJointScaleV28V30.py`
6. `scripts/RunT3ParallelABCJointScaleV28V30.ps1`
7. `tests/test_t3_parallel_a_bc_joint_scale_v28_v30.py`

## Commands

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
$zip = "$env:USERPROFILE\Downloads\t3_workstream_A_bc_joint_scalar_v28_v30.zip"
Expand-Archive -LiteralPath $zip -DestinationPath . -Force
powershell -ExecutionPolicy Bypass -File scripts\RunT3ParallelABCJointScaleV28V30.ps1
```

Optional explicit reduced mode: `powershell -ExecutionPolicy Bypass -File scripts\RunT3ParallelABCJointScaleV28V30.ps1 -Backend quadratic_only`.

A missing original native beta source is a **hard block**, not a reason to use invented coefficients. Generate source beta evidence by the existing Workstream B procedure and rerun. Existing v28/v29/v30 outputs are independently reverified and never overwritten; use a new output directory to change the benchmark parameters or scheme.

## Outputs under `output/parallel_t3/wilson/`

- `bc_joint_scalar_candidate_state_v28_B_m1.json`
- `bc_joint_scalar_candidate_state_v28_C_m1.json`
- `bc_joint_scalar_candidate_readiness_v29.json`
- `bc_fermion_first_joint_scalar_handoff_v30.json`

All have `physical_Wilson_evolution_authorized=false`. **NO Match[]**, no LLSS→C5 evolution, no scalar hard matching, no SMEFT coefficient. The EFT order is **UV → integrate out ALL fermions → integrate out ALL scalars → SMEFT**.
