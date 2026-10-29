# General T3: grouped fermion-first threshold continuation (experimental)

## Implemented

1. Run the existing *complete currently coded* B/C(alpha=-1) one-loop UV RGBeta system from `mu_UV` to a common `mu_F`, including both C5-active and four other existing RGBeta complex quartics.
2. Apply `build_intermediate_scalar_boundary` at `mu_F`: remove F and both Yukawas; transport the renormalisable quartics continuously (tree-level matching; no finite threshold correction).
3. Run the existing scalar-only **EFT1 RGBeta** one-loop system to `mu_S`, including all six complex quartics.
4. Use the existing scalar-threshold projection to SM.
5. With `final_weinberg_coefficient.json` supplied, use the **previously established baseline** hierarchical C5 evaluator and run the SM+Weinberg matrix to `mu_low`. This baseline DOES NOT contain the new quartic corrections to the Wilson matching.
6. Save all running extra quartic values at UV, F boundary and S boundary, plus baseline C5/mass matrices, and a machine-readable blocked flag for extra Wilson matching.

## NOT implemented

**The generated quartics do NOT yet contribute to the physical C5 at the scalar threshold.** The missing pieces are matching of the relevant higher-dimensional LLSS Wilson tensors at F, their EFT1 anomalous dimensions with extra quartic insertions, and the associated scalar matching and pole cancellation at S. The one-step Matchete hard response is not a substitute. No unsupported extra C5 or extra low-energy neutrino mass is output. The already existing `final_weinberg_coefficient.json` is a baseline coefficient, not an extended-model result. A complete O(hbar^2) calculation also needs other two-loop effects.

**The threshold prescription is grouped**, not individual thresholds for the three different MF and two different MS. The default benchmark enforces the strict ordering `mu_UV > max(MF) > min(MF) > mu_F > max(MS) > min(MS) > mu_S > mu_low` approximately, rechecking running masses. Grouped matching is an approximation until matched finite threshold corrections and sequential heavy-generation integration are supplied. The pre-existing Stage 11–15 benchmark does **not** have this ordering (MF1=1500 GeV < MS2=1780 GeV). Its numerical results are unchanged; this batch creates a separate benchmark with MF=[5100,6000,7000] GeV.

The benchmark values are illustrative, not measured couplings or a physics fit. The source-level `lambda_Matchete=-lambda_RGBeta` relation still assumes an identical scalar-field basis. The production numerical source files are untouched.

## Run

In the repository root:

```powershell
Expand-Archive -LiteralPath .\T3_General_Scalar_Batch16.zip -DestinationPath . -Force
python -m pytest -q tests\python\test_general_t3_fermion_first_thresholds.py
powershell -ExecutionPolicy Bypass -File .\scripts\RunGeneralT3FermionFirstThresholds.ps1
```

For the two renormalisable segments only (no baseline matching):

```powershell
powershell -ExecutionPolicy Bypass -File .\scripts\RunGeneralT3FermionFirstThresholds.ps1 -RenormalisableOnly
```

Requires existing `output/T3_{B,C}_alpha_m1/data/{uv_rgbeta_rge.json,eft1_rgbeta_rge.json,final_weinberg_coefficient.json}`. Outputs `output/general_t3_stage16/{B_m1,C_m1}_hierarchical_trajectory.json` and `stage16_summary.json`.
