# All-16 source-mapped rough numerical pipeline — batch v2

## Scope

This package is a **full-file replacement batch** for the existing v1 files, plus
one test and one Windows runner. It fixes the six blockers reported by the
first numerical batch. The existing one-scale scalar matching and SMEFT RGEs
are not rewritten. Workstream B source files are read, never changed.

The calculation retains the required ordering:

`UV -> integrate out ALL F at MF=100000 GeV -> EFT1 with BOTH scalars + full coupled dim-5 C -> integrate out BOTH scalars at MS=1000 GeV -> SMEFT -> neutrino masses`.

- **A(-2), D(+0), E(-2), E(+0):** The Matchete tensor exporter
  `Numerical/running/GeneralT3OrderedTensor.py` packs neutral self-bilinear
  scalar components as `Minus[1..d]`, then `Plus[1..d]`, *not* the interleaved
  `Re[1], Im[1], ...` order. The source-to-canonical Wilson transformation
  now composes that exact permutation with the existing triplet CG unitary.
- **B(-1), C(-1):** Their saved fermion source is an ordered mixed-X/Y tensor
  JSON, not a regular `eft1_after_F_wilson_seed.json`. The new code checks
  its original Workstream B provenance and uses the source's same-run MF
  `mixing_U_muS` matrix and the explicitly implemented
  `(S1, epsilon * conjugate(S2)) = U @ (X,Y)` basis relation to construct an
  orthogonal real 8x8 source-to-canonical tensor map.
- **Self-bilinear triplet scalar masses:** The old rough adapter represented
  the holomorphic bilinear as a Cartesian `I_3` in a spherical spin-1
  component frame. v2 converts it with the original triplet CG's
  `U.T @ U` invariant metric, preserving the eigenvalues and the complete
  real-component mixing. The corrected mass is written to this batch's own
  Stream 3 input and recorded with its new SHA256.
- **All other 10 cases:** Original completed outputs are kept unchanged on
  `--resume`. Their old tensor definitions are not replaced or projected.
- **All 16:** The v1 corrected conditional `lambdaT3` and `lambda12Cross`
  component tensors and the existing DOP853, one common MS matching,
  Weinberg running, and neutrino observable calculators are reused.

This is **conditional rough numerical implementation**, not a certified
full 1-loop physics prediction. In particular, the candidate ADM and the
signed GammaCC-to-Weyl and RGBeta coefficient normalisations have not been
independently derived in this batch. Unmatched hard C5 is deliberately
chosen to be zero for the rough benchmark; it is not asserted physically zero.
The input provenance, coordinate maps, checks, and status are written to JSON.

## Contents and intended paths (relative to project root)

| Project path | Action |
|---|---|
| `Numerical/parallel_t3/wilson/source_audited_all16.py` | Full replacement — six source maps and triplet mass basis |
| `scripts/RunT3SourceAuditedAll16.py` | Full replacement — scalar mass export, resumable all-16 runner |
| `tests/test_t3_source_audited_all16.py` | Existing v1 test file (unchanged, included for reproducibility) |
| `tests/test_t3_source_audited_all16_v2.py` | New six-case regressions against source B tensors |
| `scripts/RunT3SourceAuditedAll16V2.ps1` | New single-command Windows test + resume/fresh runner |
| `README_SOURCE_AUDITED_ALL16.md` | Replacement instructions (this file) |

No other project files need modification.

## Installation

From PowerShell:

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "$env:USERPROFILE\Downloads\t3_source_audited_all16_rough_v2.zip" -DestinationPath . -Force
```

### Recommended: resume only the six blocked models

```powershell
powershell -ExecutionPolicy Bypass -File scripts\RunT3SourceAuditedAll16V2.ps1
```

This reruns blocked cases and prints all 16 results. Run with `-WithTests`
to execute the source-map regression tests first. The original 10 completed models and their output directories
are untouched. `--resume` automatically detects the six blocked cases in the
old `output/parallel_t3/source_audited_all16_v1/summary.json`.

Equivalent Python commands:

```powershell
.\.venv\Scripts\python.exe -m pytest -q tests\test_t3_source_audited_all16.py tests\test_t3_source_audited_all16_v2.py
.\.venv\Scripts\python.exe scripts\RunT3SourceAuditedAll16.py --resume
.\.venv\Scripts\python.exe scripts\RunT3SourceAuditedAll16.py --show
```

### Fresh, fully reproducible all-16 recomputation

```powershell
powershell -ExecutionPolicy Bypass -File scripts\RunT3SourceAuditedAll16V2.ps1 -Fresh
```

This recreates the owned output folder for all 16. The underlying saved
Workstream B source files are read-only. To run just preflight/boundary checks:

```powershell
powershell -ExecutionPolicy Bypass -File scripts\RunT3SourceAuditedAll16V2.ps1 -BoundaryOnly
```

To rerun one model directly, use a separate fresh output directory or
`--overwrite --model E_-2` (the latter replaces the entire owned batch folder;
it is **not** equivalent to `--resume`).

## Outputs

- Batch summary:
  `output/parallel_t3/source_audited_all16_v1/summary.json`
- Source-mapped input copies:
  `output/parallel_t3/rough_inputs/source_audited_all16_v1/a/<slug>/`
- Complete canonical common-scalar matching backgrounds:
  `output/parallel_t3/rough_inputs/source_audited_all16_v1/c/<slug>/`
- Full EFT1 evolved Wilson tensor:
  `output/parallel_t3/source_audited_all16_v1/wilson/<slug>/wilson_at_MS.npz`
- Scalar-matching neutrino observables:
  `output/parallel_t3/source_audited_all16_v1/results/<slug>/observables.json`

Every completed model retains its provenance manifest and numerical masses;
failed models have a written `blocker` and are not filled with synthetic
numbers. The runner exits nonzero if any requested model remains blocked.

## Source basis references

- `Numerical/running/GeneralT3OrderedTensor.py`, `scalar_basis` and real scalar
  expansion (tracked source; grouped real fields for self-bilinear models).
- `Lagrangian/RunGeneralT3MassEigenMatchete.wl`, `T3MERotateInteraction`
  (native X/Y source field transformation and epsilon bridge).
- `Numerical/parallel_t3/threshold/models.py`, `_mixed_rotation` and
  `from_parameters` (mass mixing and complete real scalar Hessian).
- `Numerical/parallel_t3/background/fermion.py`,
  `evaluate_mixed_native_candidate` (source of the same-run X/Y tensor).
- `Numerical/parallel_t3/wilson/matchete_a4_map_v8.py` (A(-4) source CG anchor).

Source inspected: public repository `main`, commit
`9e50ca77e64e9d5fb37079f325bc2b8cab5112e9`. Local user execution is
necessary for the actual all-16 end-to-end numerical outcome.
