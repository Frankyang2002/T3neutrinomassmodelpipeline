# Stage 9: standalone Matchete tensor exporter initialisation

## Verified input

User's Batch 8 local run: 42 Python tests passed. The formal C5 response passed,
including the Yukawa-family classification for B/C(-1). The one-loop RGBeta
source check passed. The tensor exporter then emitted
`DefineRepresentation::unkwngroup: SU2L has not been defined as a group` and
reported `GenerationFailed` for its four active structures.

**Cause:** `ExportGeneralT3TensorConvention.wl` called `ResetAll[]` and then
built BSM SU(2) tensors without first calling `LoadModel["SM"]`. The generator's
`EnsureBSMRepresentation` registers a new representation under Matchete's
`SU2L` group, which does not exist until the SM model has been loaded.
`BuildT3Lagrangian` already follows the correct order: `ResetAll[]`,
`LoadModel["SM"]`, then BSM field/representation generation.

## Implemented

- Full replacement `scripts/ExportGeneralT3TensorConvention.wl`: reloads the
  SM group registry immediately after every Matchete reset, *inside*
  `UsingFrontEnd`, and distinguishes a failed SM load from a tensor mismatch.
- Full replacement `scripts/RunGeneralT3Stage6.ps1`: optional `-TensorOnly`
  for retrying the failed stage without repeating proven C5 and RGBeta work.
- Full replacement `tests/python/test_general_t3_stage6.py`: adds regression
  checks for the prerequisite initialization and targeted runner mode.

This corrects the `SU2L` registry error but does **not** claim the resulting
Matchete tensor is proportional to the canonical index tensor; that comparison
must still be run using the user's Matchete installation. Even a successful
canonical comparison does **not** by itself establish the Matchete-to-RGBeta
normalisation, signs or scalar potential convention.

## Run (PowerShell, repository root)

```powershell
Expand-Archive -LiteralPath .\T3_General_Scalar_Batch9.zip -DestinationPath . -Force
python -m pytest -q tests\python\test_general_scalar_trial.py tests\python\test_scalar_quadratic_prototype.py tests\python\test_general_scalar_c5_proof.py tests\python\test_general_rgbeta_coverage.py tests\python\test_general_t3_stage6.py
powershell -ExecutionPolicy Bypass -File .\scripts\RunGeneralT3Stage6.ps1 -TensorOnly
```

If the tensor stage passes, the existing Stage 8 formal-response and RGBeta
reports can be retained, or run the complete Stage 6 script to refresh them:

```powershell
powershell -ExecutionPolicy Bypass -File .\scripts\RunGeneralT3Stage6.ps1
```

Expected tensor JSON:
`output/general_t3_trial/matchete_active_quartic_tensors.json`.
The field `rgbeta_mapping_verified` remains `false` regardless of the tensor
check result. General RGBeta RGE and numerical C5 are still *not* activated.
