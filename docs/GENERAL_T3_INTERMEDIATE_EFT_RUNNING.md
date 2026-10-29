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
