# Documentation consolidation plan (read-only)

- HEAD: `46b4f59915956d7568f91ede06e9162b82934429`
- Root docs Markdown: **120**
- Can consolidate and archive: **76**
- Retained or blocked: **44**
- New topic volumes: **7**
- Expected docs Markdown count after apply: **52**
- Text files scanned: **1092**; read errors: **0**

## Volume allocations

- `RESEARCH_GENERAL_T3_TENSORS_AND_QUARTICS.md`: 10 source documents
- `RESEARCH_GENERAL_T3_RUNNING_AND_MATCHING.md`: 7 source documents
- `RESEARCH_MATCHETE_AND_MASS_EIGEN.md`: 17 source documents
- `RESEARCH_EFT1_THRESHOLDS_AND_WEINBERG.md`: 17 source documents
- `RESEARCH_NUMERICAL_FITS_AND_BENCHMARKS.md`: 16 source documents
- `RESEARCH_GENERAL_SCALAR_AND_PILOTS.md`: 3 source documents
- `REPOSITORY_MAINTENANCE_HISTORY.md`: 6 source documents

## Retained and blocked

- `docs/FERMION_START_STUDY.md` — RETAIN_CANONICAL
- `docs/GENERAL_SCALAR_48_FINAL_REPORTING.md` — RETAIN_CANONICAL
- `docs/GENERAL_T3_16_MODEL_UV_RUNNING.md` — RETAIN_REFERENCED; source: tests/python/test_descriptive_research_file_names.py, tools/RenameNumberedResearchFiles.py
- `docs/GENERAL_T3_C5_FINITE_MASS_LOGARITHMS.md` — RETAIN_REFERENCED; source: tools/RenameNumberedResearchFiles.py
- `docs/GENERAL_T3_C5_LINEAR_RESPONSE_CORRECTION.md` — RETAIN_REFERENCED; source: tools/RenameNumberedResearchFiles.py
- `docs/GENERAL_T3_CONDITIONAL_UV_TO_C5.md` — RETAIN_REFERENCED; source: tools/RenameNumberedResearchFiles.py
- `docs/GENERAL_T3_COUPLED_UV_RUNNING.md` — RETAIN_REFERENCED; source: tools/RenameNumberedResearchFiles.py
- `docs/GENERAL_T3_EFT1_SCALAR_PAIR_MIXING.md` — RETAIN_REFERENCED; source: tests/python/test_descriptive_research_file_names.py, tools/RenameNumberedResearchFiles.py
- `docs/GENERAL_T3_EFT1_WILSON_SOURCE_TRANSPORT.md` — RETAIN_REFERENCED; source: tools/RenameNumberedResearchFiles.py
- `docs/GENERAL_T3_EXTRA_C5_PHASE_RESPONSE.md` — RETAIN_REFERENCED; source: tools/RenameNumberedResearchFiles.py
- `docs/GENERAL_T3_FERMION_FIRST_THRESHOLDS.md` — RETAIN_REFERENCED; source: tools/RenameNumberedResearchFiles.py
- `docs/GENERAL_T3_FINAL_C5_AND_WEINBERG_RUNNING.md` — RETAIN_REFERENCED; source: tools/RenameNumberedResearchFiles.py
- `docs/GENERAL_T3_GENERAL_QUARTIC_RGBETA_PROGRAMME.md` — RETAIN_REFERENCED; source: tools/RenameNumberedResearchFiles.py
- `docs/GENERAL_T3_MATCHETE_TENSOR_INITIALIZATION.md` — RETAIN_REFERENCED; source: tools/RenameNumberedResearchFiles.py
- `docs/GENERAL_T3_MATCHING_RG_SOURCE_PILOTS.md` — RETAIN_REFERENCED; source: tools/RenameNumberedResearchFiles.py
- `docs/GENERAL_T3_QUARTIC_CONVENTION_C5.md` — RETAIN_REFERENCED; source: tools/RenameNumberedResearchFiles.py
- `docs/GENERAL_T3_RGBETA_EXPORT_COVERAGE_AUDIT.md` — RETAIN_REFERENCED; source: tools/RenameNumberedResearchFiles.py
- `docs/GENERAL_T3_SCALAR_QUARTIC_VALIDATION.md` — RETAIN_REFERENCED; source: tools/RenameNumberedResearchFiles.py
- `docs/GENERAL_T3_SCALAR_TENSOR_CERTIFICATION.md` — RETAIN_REFERENCED; source: tools/RenameNumberedResearchFiles.py
- `docs/GENERAL_T3_TENSOR_110101_LOCAL_CORRECTION.md` — RETAIN_REFERENCED; source: tools/RenameNumberedResearchFiles.py
- `docs/GENERAL_T3_TENSOR_110101_SU2_COVARIANCE.md` — RETAIN_REFERENCED; source: tests/python/test_descriptive_research_file_names.py, tools/RenameNumberedResearchFiles.py
- `docs/GENERAL_T3_TENSOR_ALIGNED_GEOMETRY.md` — RETAIN_REFERENCED; source: tools/RenameNumberedResearchFiles.py
- `docs/GENERAL_T3_TENSOR_EXPLICIT_GEOMETRY.md` — RETAIN_REFERENCED; source: tools/RenameNumberedResearchFiles.py
- `docs/GENERAL_T3_WEINBERG_SOURCE_TRANSPORT.md` — RETAIN_REFERENCED; source: tools/RenameNumberedResearchFiles.py
- `docs/GENERAL_T3_WOLFRAMSCRIPT_SELFTEST_FIX.md` — RETAIN_REFERENCED; source: tools/RenameNumberedResearchFiles.py
- `docs/HIGH_SCALAR_HALF_FERMION_START.md` — RETAIN_REFERENCED; source: Numerical/orchestration/RunHighScalarHalfFermionStart.py
- `docs/HISTORICAL_DOCS_CONSOLIDATION_GUIDE.md` — RETAIN_CANONICAL
- `docs/NUMERICAL.md` — RETAIN_CANONICAL
- `docs/ONE_TIME_ARTIFACT_CLEANUP.md` — RETAIN_REFERENCED; source: tools/AuditOneTimeArtifacts.py
- `docs/PHYSICS.md` — RETAIN_CANONICAL
- `docs/REGENERATE_CURRENT_REPORTS.md` — RETAIN_CANONICAL
- `docs/REPOSITORY_CLEANUP_AUDIT.md` — RETAIN_REFERENCED; source: tools/AuditDocumentationRetirement.py
- `docs/RESULTS.md` — RETAIN_CANONICAL
- `docs/SCALAR_SECTOR.md` — RETAIN_CANONICAL
- `docs/SCRIPT_DEPENDENCY_REFINEMENT.md` — RETAIN_REFERENCED; source: tools/RefineScriptDependencyAudit.py
- `docs/T3_GENERAL_SCALAR_STUDY_SEPARATION.md` — RETAIN_CANONICAL
- `docs/T3_REDUCED_SCALAR_POTENTIAL_CONCLUSION.md` — RETAIN_REFERENCED; source: scripts/RunT3ReducedScalarSobolComparison.py, tools/AuditDocumentationRetirement.py
- `docs/T3_THESIS_PLAN.md` — RETAIN_CANONICAL
- `docs/T3_ZERO_BOUNDARY_LOW_ENERGY_BATCH8.md` — RETAIN_REFERENCED; source: tools/AuditDocumentationRetirement.py
- `docs/T3B_YUKAWA_BASIS.md` — RETAIN_CANONICAL
- `docs/TEST_BASELINE_FIX.md` — RETAIN_REFERENCED; source: SHA256_MANIFEST.json
- `docs/THESIS_EXISTING_COMPARISONS_PLACEMENT.md` — RETAIN_CANONICAL
- `docs/THESIS_FIGURE_CAPTIONS.md` — RETAIN_CANONICAL
- `docs/THESIS_FIGURE_PLACEMENT.md` — RETAIN_CANONICAL

## Safety

Archives are created OUTSIDE the repository with original bytes and SHA-256 hashes.
The planner checks current tracked status and references. No physics source or result is edited.
Do not use --apply before reviewing inventory.csv and blocking reasons.
