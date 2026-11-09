# Workstream A: all-16 native Matchete portal exporter — B/C(−1) mixed-mass exception

## Observed failure

`BuildT3Lagrangian` returns `Status = "RequiresMixedScalarMassDiagonalization"`
for B(−1) and C(−1) when the additional quadratic is enabled. Matchete's
`CheckLagrangian` rejects the unrotated heavy-scalar mass matrix, but the
builder retains `LInt`, the registered portal CGs, the accepted interaction
inventory, and explicit `MatchingMassBasisReady = False`.

The previous portal exporter required `Status = "Success"` for every model.
It additionally wrapped the whole builder (and entire run) in `Check[...]`,
which turns expected `HeavyMassBasis` warning messages into `$Failed` and
silently discards the useful symbolic builder association.

## Narrow change

* Retain full build and mass-basis validation for all models other than
  B(−1)/C(−1).
* Permit only the exact B(−1)/C(−1) status
  `RequiresMixedScalarMassDiagonalization`, and only with the verified
  additional-quadratic inventory, positive blocked mass term, absent
  `FullValidation`, and false `MatchingMassBasisReady`.
* Inspect only the native portal CGs in the accepted symbolic `LInt` and
  `AllowedInteractions`. Each tensor is still retrieved from native
  `GetCGTensor`, checked numerically and against the corrected in-`PlusHc`
  Hermitian phase. The A(−4) `I` requirement is unchanged.
* Preserve each model's true native build status in the JSON. The downstream
  Python audit checks this status and carries the mixed-mass blocker into
  its result. Any other failed source status is rejected.
* Remove `Check[...]` around the entire builder and model sweep, because
  these scope-wide `Check` expressions reject the **expected** warning.
  All critical data are independently checked with status and tensor guards.

## Scope

This allows **tensor evidence extraction only**. It does **not** diagonalize
the B/C(−1) scalar mass matrices, authorize `Match[]`, or certify physical
threshold inputs. The RGBeta adjoint generator orientation and native real
scalar coordinate map remain unresolved.

No change to `ScalarPotential.wl`, `LagrangianBuilder.wl`, or the RGBeta source.
No remote GitHub operations are performed.

## Commands from the repository root (PowerShell)

```powershell
$zip = "C:\path\to\t3_workstream_A_portal_v5_mixed_mass_fix.zip"
Expand-Archive -Path $zip -DestinationPath . -Force

.\.venv\Scripts\python.exe -m unittest discover -s tests -p "test_t3_parallel_a*.py" -v

wolframscript -file scripts\RunT3ParallelAExportPortalMatchete.wl `
  output\parallel_t3\wilson\matchete_portals_native_v5.json

Test-Path output\parallel_t3\wilson\matchete_portals_native_v5.json

.\.venv\Scripts\python.exe -m scripts.RunT3ParallelAPortalAudit `
  --matchete output\parallel_t3\wilson\matchete_portals_native_v5.json `
  --rgbeta output\parallel_t3\wilson\rgbeta_portals_native_v5.json `
  --project-root . `
  --output output\parallel_t3\wilson\portal_dictionary_conditional_v5.json
```

The Matchete exporter refuses to overwrite an existing output. The RGBeta
exporter has already run successfully and need not run again.

## Validation

Local Python unit suite: 53 tests passed using previously uploaded real
50-slot tensor evidence; Wolfram source passed lexical bracket/string checks.
Live native Matchete execution of this fix is pending on the user's system.
