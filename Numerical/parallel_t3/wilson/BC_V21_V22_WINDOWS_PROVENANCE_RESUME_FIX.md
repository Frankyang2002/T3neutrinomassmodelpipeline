# Workstream A v21/v22 — Windows v19 provenance/resume repair

## Observed failure

The user's v21 full-repository integration test failed at
`v21.audit -> recorded.get(rel)==digest(f)` for
`output/parallel_t3/wilson/matchete_topology_native_v18.json`.
The failure occurred **before** numerical mass-basis rotation.

The source exporter v19 (`Numerical/parallel_t3/wilson/complex_scalar_bridge_v19.py`)
records input paths with `str(p.relative_to(root))`, a **backslash-separated
string on Windows**. The original v21 reader instead requested
`f.relative_to(root).as_posix()` (forward slashes). Dictionary lookup returned
`None` even for unchanged file bytes. This is a reproducible cross-platform
provenance-reader defect; the user's actual original SHA256 values were not
available in the preparation environment, so any additional changes remain
subject to the live revalidation.

## Replacement behavior

- Read v19's seven original `input_sha256` entries, normalize only their
  **path keys** (Windows backslash / POSIX slash), and compute **fresh
  SHA256 hashes of the actual original source bytes**. Saved JSON is never
  changed.
- Reject absent/incomplete seven-file inventories, noncanonical paths,
  absolute paths, traversal, duplicate path aliases, changed content, or
  missing native v18 topology / B/C source entries. The mismatch error
  distinguishes a missing key from a real SHA256 discrepancy.
- Recheck v20's saved v19 hash and B/C saved quadratic source hash chains.
- The optional read-only preflight prints the validated source counts and
  saved scalar mass eigenvalues before starting v21.
- Preserve v18–v20 JSON outputs and v21/v22 resume refusal to overwrite
  altered results. Do not regenerate v19 merely to change slash conventions.

## Source and scope

All paths are relative to repository root. Only the following files are
changed or added by this replacement batch:

1. `Numerical/parallel_t3/wilson/bc_mass_eigen_native_bridge_v21.py` — full replacement.
2. `scripts/RunT3ParallelABCMassEigenNativeV21.py` — full replacement CLI.
3. `scripts/RunT3ParallelABCMassEigenBatchV21V22.ps1` — full replacement runner.
4. `scripts/RunT3ParallelABCMassProvenanceV22.py` — new preflight CLI.
5. `tests/test_t3_parallel_a_bc_mass_bridge_v21_v22.py` — full replacement with Windows, POSIX, hash-tampering and simulated end-to-end tests.
6. `Numerical/parallel_t3/wilson/BC_V21_V22_WINDOWS_PROVENANCE_RESUME_FIX.md` — this documentation.

No Matchete source, RGBeta model, source Lagrangian, saved witness, physical
Wilson data, or GitHub file is modified by these files.

## Run on Windows (PowerShell)

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
$zip = "$env:USERPROFILE\Downloads\t3_workstream_A_bc_v21_v22_windows_hash_fix.zip"
Expand-Archive -LiteralPath $zip -DestinationPath . -Force
powershell -ExecutionPolicy Bypass -File scripts\RunT3ParallelABCMassEigenBatchV21V22.ps1
if ($LASTEXITCODE -ne 0) { throw 'v21/v22 batch failed' }
```

To diagnose the provenance chain alone without redoing tensor actions:

```powershell
.\.venv\Scripts\python.exe -m scripts.RunT3ParallelABCMassProvenanceV22 --project-root .
```

Expected new outputs after the full audited run:

- `output/parallel_t3/wilson/bc_native_mass_eigen_complex_conditional_v21.json`
- `output/parallel_t3/wilson/bc_quartic_mass_readiness_conditional_v22.json`

These outputs are **conditional** source-action / mass-transport witnesses.
No native Matchete redefined eigenfield matching, independent numerical RGBeta
S2 frame, signed GammaCC/Weyl normalization, physical one-loop ADM, or physical
Wilson evolution is certified. EFT order remains UV → all fermions out → all
heavy scalars out → SMEFT.
