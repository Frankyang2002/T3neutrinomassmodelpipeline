# T3 neutrino mass project — architecture and research history

This guide consolidates the historical stage notes into a navigable record.
The original documents are preserved in the external cleanup archive.
This is an inventory and workflow guide, **not** a new physics certification.

## Verified production workflow (one-point smoke test, 2026-11-01)

The ordinary T3-B, alpha=-1 benchmark completed the current production numerical path:

1. Supply a representation, SM state, and heavy-sector UV benchmark.
2. Load the exported UV and intermediate RGBeta beta functions.
3. Integrate the UV RGEs to the heavy-fermion threshold.
4. Run the scalar-only intermediate EFT1 renormalisable state.
5. Evaluate the supplied Matchete-derived final Weinberg coefficient using `FinalC5TrajectoryAdapter`.
6. Run SM + Weinberg evolution to the low scale.
7. Construct the neutrino mass matrix, Takagi diagonalise, extract oscillation observables.

The prior smoke result reported `PASS_END_TO_END_OPERATIONAL` with finite oscillation outputs.
This demonstrates software operation for that benchmark, not physics validation for every model.

## Where to look first

- `pipeline.py`: overall project entry point (consult its CLI before invoking).
- `Numerical/orchestration/PipelineNumericalResults.py`: integrated numerical diagnostics.
- `Numerical/fitting/ScanCLI.py`: supported fixed/scan configuration entry point.
- `Numerical/fitting/ParameterScan.py`: per-point numerical workflow.
- `Numerical/running/`: numerical RGE and EFT evolution.
- `Numerical/core/`: state representations and beta evaluation.
- `physics/NeutrinoObservables.py`: Takagi and mixing observables.
- `Lagrangian/` and `RGE/`: matching, model construction and beta generation.
- `configs/`: benchmark and study configuration files.
- `tests/`: regression and physics checks.

## Research progression (conceptual chronology)

1. **Model setup:** T3 representations, group invariants, field content and Lagrangians.
2. **RG and matching exports:** Matchete and RGBeta source and coefficient artifacts.
3. **Numerical model:** UV parameters, state vectors, coupled running.
4. **EFT evolution:** heavy fermion threshold, scalar EFT, Weinberg coefficient.
5. **Phenomenology:** neutrino masses, mass splittings, mixing and scans.
6. **General extensions:** additional quartics/Yukawas and conditional one-step C5 studies.
7. **Operational integration:** single-point production smoke test passed.

## Important scope distinction

- **Production path:** ordinary split-scalar numerical scan, with a working one-point end-to-end run.
- **Experimental extensions:** additional couplings and hypothetical RGBeta-to-Matchete maps; not automatically production-certified.
- **Known restrictions:** a fully certified general coupling-basis map and sequential thresholds are not established for every extension.

## Historical source notes preserved by this cleanup

The archive manifest gives the original and archived location for every item.
The headings below are extracted from local notes when this cleanup runs.

### BATCH1_README.md

Original: `BATCH1_README.md`

Sections:
- T3 Batch 1: executable completion assessment
- Install
- Run (PowerShell, from the repository root)
- Source-level conclusions already supported
- Batch 2 priority
- Limitations

### BATCH2_README.md

Original: `BATCH2_README.md`

Sections:
- Batch 2 — conditional numerical result runner

### README.txt

Original: `README.txt`

Historical stage note (no Markdown headings detected).

### README_OnePromptT3Integration.md

Original: `README_OnePromptT3Integration.md`

Sections:
- T3 — one-prompt local general-Yukawa RGE integration
- Files (install ZIP at repository root)
- Exact PowerShell commands
- Pure Python structural tests before changing production source
- Guarded local full-file replacement and backup
- Reverify the corrected charged doublet Dirac mass
- Five model smoke tests, then all sixteen general-scalar RGBeta UV models
- All six extra-Yukawa RGBeta-only numerical one-step benchmarks
- Output paths
- Exact implementation and certification boundaries
- Undoing the local source change

### README_T3CompletionBatch.md

Original: `README_T3CompletionBatch.md`

Sections:
- T3 one-batch completion: full coupled general-Yukawa **UV-only** numerical integration
- Source status
- Files in this ZIP (complete new files, paths relative to repository root)
- Installation and full local execution (PowerShell)
- Uses your existing 16 *actual* Wolfram RGBeta exports.
- Includes an RHS preflight for every model, then integrates all couplings.
- Outputs
- What exactly is implemented
- What is NOT physically certified or performed

### README_T3EFT1WilsonTensor.md

Original: `README_T3EFT1WilsonTensor.md`

Sections:
- T3 EFT1 heavy-generation Wilson-tensor bridge (local experiment)
- Purpose and scope
- Project paths in this ZIP
- PowerShell (repo root)
- Isolate the exact 4 original A(0) tree monomials first.
- Next: attempt scalar tree matching of contracted EFT1.
- Next: scalar one-loop [B], still WITHOUT inherited hard [C] or EFT1 RGEs.
- Physical limitations — do not promote status

### README_T3ExtendedC5ExpressionReduction.md

Original: `README_T3ExtendedC5ExpressionReduction.md`

Sections:
- Native extended C5 coefficient reduction — source-derived, conditional numerical bridge
- Exact source findings
- Project paths
- PowerShell commands (local repository only)
- Numeric, parser and sign/index/loop negative controls, using your own real native files:
- Algebraically verify all six native matched expressions (no fresh Wolfram needed):
- Optional: one deliberately synthetic, conditional Matchete-basis 3x3 example:
- About the example
- Physical blockers retained

### README_T3ExtendedFlavorIndexFix.md

Original: `README_T3ExtendedFlavorIndexFix.md`

Sections:
- T3 native extended C5 — Matchete flavor-index fix
- Native failure established on Windows
- Implemented
- Files in ZIP (relative to repository root)
- PowerShell
- Small native validation: two flavors satisfy the installed Matchete API.
- All six with physically relevant three-heavy-generation input, only
- after the small native pilot succeeds.
- Actual test scope

### README_T3ExtendedLowEnergyBatch.md

Original: `README_T3ExtendedLowEnergyBatch.md`

Sections:
- T3 extended C5 — coordinated large-batch low-energy study
- Scope and physical status
- Install (Windows PowerShell)
- First test A(0), one endpoint and four phases
- Then all six models, three endpoint scales and four phases
- New files / intended installation paths
- Outputs
- Status distinctions

### README_T3ExtendedMatcheteLoadFix.md

Original: `README_T3ExtendedMatcheteLoadFix.md`

Sections:
- Extended Matchete C5 native load-order correction
- What failed
- Root cause and fix
- Files and project paths
- PowerShell
- Source-level test only (does not execute Matchete):
- Real live Matchete load + A(0) pilot:
- Only if pilot returns 0:

### README_T3ExtendedPhysicalMatching.md

Original: `README_T3ExtendedPhysicalMatching.md`

Sections:
- T3 general-Yukawa physical matching batch
- Scope and prior inputs
- New files / exact project paths
- PowerShell
- Check source and existing 16/16 RGE report without running Wolfram
- Recommended initial live native pilot for A(0), plus the numerical proxy
- All six real native matches and the six numerical proxies. Results survive
- a failed model; the command exits 4 if any model is blocked.
- Output locations
- Physics and interpretation
- Optional standalone numerical C5 check
- Tests performed by package author here

### README_T3HierarchicalExtraYukawa.md

Original: `README_T3HierarchicalExtraYukawa.md`

Sections:
- T3 extra-Yukawa hierarchical Matchete matching pilot
- Scope
- What this does not do
- Project paths (all complete new files)
- PowerShell (repo root)
- First substantial live stage: integrate out F, with seeded/zero action controls.
- Only if the F-only native pilot succeeded, test the second threshold.
- When the pilot succeeds, run the entire 6-model family.
- These are two separate large Wolfram batches; choose one or both.
- Main reports
- Physical meaning and follow-up
- Current validation

### README_T3RGAssistedExtendedC5.md

Original: `README_T3RGAssistedExtendedC5.md`

Sections:
- T3 Workstream A — UV-to-extended-C5 conditional numerical bridge
- Purpose and exact physics scope
- Implemented equations
- Files in ZIP (new files only)
- Local PowerShell
- Tests performed before release
- What genuinely remains

### README_T3Stage2IndexDiagnostic.md

Original: `README_T3Stage2IndexDiagnostic.md`

Sections:
- T3 hierarchical stage-2 index diagnosis and controlled common-mass experiment
- Cause identified from the actual `T3_Ap0_stage2_failure.zip`
- What this package adds
- Commands (PowerShell, repository root)
- CONTROLLED EXPERIMENT: three degenerate heavy flavours (one common MF).
- Do not confuse it with the previous three *independent* MF[r] benchmark.
- Status and scientific scope

### README_THESIS_COMPARISONS.txt

Original: `README_THESIS_COMPARISONS.txt`

Historical stage note (no Markdown headings detected).

### SMOKE_TEST_README.md

Original: `SMOKE_TEST_README.md`

Sections:
- T3 end-to-end software smoke test
- Installation
- Commands (PowerShell)
- Output

## Cleanup conventions

- `output/` and `Reports/output/` are **copied** to external archive, not removed: they may contain required runtime artifacts.
- Old stage README files are moved after being indexed above.
- Unreferenced probe source files and backup copies are moved to an external archive, not irreversibly deleted.
- Referenced probes stay in the project until their dependencies can be evaluated.
- No primary physics source, tests or configs are changed by this cleanup.

## Restore

Run `python scripts/maintenance/CleanRepository.py --restore "PATH_TO_ARCHIVE"` from the repository root.
Restoration refuses to overwrite modified or pre-existing paths.

