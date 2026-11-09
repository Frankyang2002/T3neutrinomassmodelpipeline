# Workstream A v13 — native RGBeta `Lam`/projector failure forensics

## Observed Windows status

- v10 smoke EFT1 registration succeeded for A(-4), B(-1), E(-2): 38 quartic registrations, 29 real-quartic diagonal attempts.
- v10 saved the failed native witness `rgbeta_quartic_projectors_native_v10_smoke.json` and exited 7 (`Blocked_NativeRGBetaQuarticProjectorResponses`). The individual numeric/symbolic failure reasons were **not included in the console output**.
- v11 (`nonportal_scalar_invariants_conditional_v11.json`) and v12 (`topology_scalar_action_conditional_v12.json`) independently revalidated and remain reusable.
- A numerical scalar projector normalization certificate **has not** been established. The mere existence of an attempted native `Lam` contraction is not enough.

## Source inspection

Public project HEAD `2bc2e4b516776212645370105cc63495a3c47348` (9 Nov 2026), read-only. Original v10 exporter and installed RGBeta source use actual `$quartics`, `Projector` and `Lam` definitions. Upstream `RGBeta/Package/FieldsAndCouplings.m` defines both `UpdateProjectors` and `Lam`, including scalar conjugate entries, unique permutation averaging and the factor `6*Matrix[coupling]`. These expressions are symbolic; no numerical native S2 axes are available.

**Do not infer the cause of the failed response until the saved per-coupling rows are examined.** Potential categories include unresolved group structures, an actual non-unit diagonal, a singular response block or a timeout; the current console does not distinguish these.

## This batch: full replacement files, no changes to production physics

1. `scripts/RunT3ParallelAQuarticForensicsV13.py`: checks original native witness SHA256 against local source and installed RGBeta entrypoint. Validates 3 exact model inventories, 29 real diagonal rows, portal blocks and their response values; enumerates failures. Produces an immutable diagnostic named by original smoke SHA256.
2. `scripts/RunT3ParallelAProbeQuarticResponseV13.wl`: source-checked native differential re-probe on the first failed target/source pair per model, comparing **the original isolated `$quartics` method** with **full native `$quartics` then zeroing unrelated coupling symbols**. Both methods evaluate installed `Projector` and `Lam`; neither injects numeric `tGen` components. Time limits and concise native `InputForm` excerpts are included. This is diagnostic only.
3. `scripts/RunT3ParallelACompareQuarticReprobeV13.py`: independent classification of the two native routes; refuses mismatched forensic SHA or physics-certification flags and prints a small human-readable comparison.
4. `scripts/RunT3ParallelAScalarBatchV10V12.ps1`: complete resumable runner, retains v11/v12 source revalidation and isolated environment flags; now distinguishes a **blocked existing** v10 smoke file from a passed file. Preserves old witness; triggers v13 diagnostics; stops before full v10/Wilson. It intentionally returns code **7** while native v10 remains blocked.
5. `tests/test_t3_parallel_a_v10_failure_forensics_v13.py`: regression tests for source hashes, installed package identity, original output immutability, missing/malformed model evidence, multiple failure classes, both native routes, physics veto and runner guards.
6. `Numerical/parallel_t3/wilson/WORKSTREAM_A_V13_NATIVE_QUARTIC_FAILURE_TRIAGE.md`: this guide.

## PowerShell (project root)

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
$zip = "$env:USERPROFILE\Downloads\t3_workstream_A_v13_native_quartic_forensics.zip"
Expand-Archive -LiteralPath $zip -DestinationPath . -Force

.\.venv\Scripts\python.exe -m unittest discover `
  -s tests -p 'test_t3_parallel_a_v10_failure_forensics_v13.py' -v
if ($LASTEXITCODE -ne 0) { throw 'v13 forensic tests failed' }

powershell -ExecutionPolicy Bypass -File `
  scripts\RunT3ParallelAScalarBatchV10V12.ps1 -NativeOnly
```

An exit code of 7 means **source-guarded native evidence remains blocked**, not that the forensic classification itself failed. Paste the printed per-coupling diagnosis for the next targeted fix; the original JSON and v11/v12 outputs remain unchanged.

### Expected diagnostic files (names depend on original smoke SHA256)

- `output/parallel_t3/wilson/rgbeta_quartic_forensics_v13_<SHA12>.json`
- `output/parallel_t3/wilson/rgbeta_quartic_reprobe_v13_<SHA12>.json`

If native re-probe itself fails, forensic Python output is still available and should be used first. Do not delete or overwrite the old blocked smoke file.

## Certification boundary

v13 does **not** solve scalar-normalization defects, the full Matchete non-portal dictionary, B(-1)/C(-1) scalar mixing, signed GammaCC -> Weyl factors or the one-loop ADM. It does **not** authorize Wilson evolution or scalar matching. EFT sequence remains `UV -> integrate out ALL fermions -> integrate out ALL scalars -> SMEFT`.
