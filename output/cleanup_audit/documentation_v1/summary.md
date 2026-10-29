# Historical documentation reference audit

- Checkout HEAD: `fe0f511104cce8d7d914d114e3fd884eb47c6a95`
- Text files scanned: **10013**
- Candidates present: **14 / 14**
- Safe for external archival: **11**
- Blocked by retained references: **3**
- Source read errors: **0**
- Authoritative document fingerprint match: **3/3**

APPLIED: the safe set was copied to the external archive with SHA-256 verification; original active copies were then removed.
Archival additionally requires all three installed master documents to match the delivered consolidated versions.
Historical source mentions in authoritative documentation are retained as provenance, not live links.

| Candidate historical document | Classification | Live blockers | Replacement |
|---|---|---:|---|
| `docs/T3_SCALAR_POWER_COUNTING.md` | SAFE | 0 | `docs/SCALAR_SECTOR.md` |
| `docs/T3_SCALAR_VERTEX_BATCH2.md` | SAFE | 0 | `docs/SCALAR_SECTOR.md` |
| `docs/T3_SCALAR_QUADRATIC_BATCH3.md` | SAFE | 0 | `docs/SCALAR_SECTOR.md` |
| `docs/T3_SCALAR_RG_CLOSURE_BATCH4.md` | SAFE | 0 | `docs/SCALAR_SECTOR.md` |
| `docs/T3_SCALAR_REDUCED_RG_BATCH5.md` | SAFE | 0 | `docs/SCALAR_SECTOR.md` |
| `docs/T3_LEADING_SCALAR_POTENTIAL.md` | SAFE | 0 | `docs/SCALAR_SECTOR.md` |
| `docs/T3_REDUCED_SCALAR_POTENTIAL_CONCLUSION.md` | BLOCKED | 115 | `docs/SCALAR_SECTOR.md` |
| `docs/T3_ZERO_EXTRA_BOUNDARY_BENCHMARK_BATCH7.md` | SAFE | 0 | `docs/SCALAR_SECTOR.md` |
| `docs/T3_ZERO_BOUNDARY_LOW_ENERGY_BATCH8.md` | BLOCKED | 4 | `docs/SCALAR_SECTOR.md` |
| `docs/T3_GENERAL_SCALAR_STUDY_SEPARATION.md` | BLOCKED | 3 | `docs/SCALAR_SECTOR.md` |
| `docs/T3_REDUCED_SCALAR_SOBOL_COMPARISON.md` | SAFE | 0 | `docs/NUMERICAL.md` |
| `docs/T3_REDUCED_SCALAR_FINAL_RESULTS.md` | SAFE | 0 | `docs/RESULTS.md` |
| `docs/T3_REDUCED_SCALAR_EXISTING_INTERACTIVE_REPORTS.md` | SAFE | 0 | `docs/NUMERICAL.md` |
| `docs/T3_REDUCED_SCALAR_LAMBDA_STYLE_REPORTS.md` | SAFE | 0 | `docs/NUMERICAL.md` |

## What prevents archiving

### `docs/T3_REDUCED_SCALAR_POTENTIAL_CONCLUSION.md`

- `docs/T3_GENERAL_SCALAR_STUDY_SEPARATION.md:27` — blocking_reference: `Refer to docs/T3_REDUCED_SCALAR_POTENTIAL_CONCLUSION.md (prior Batch 8 package) and output/scalar_relevance/final_leading_scalar_potential.md for`
- `docs/T3_ZERO_BOUNDARY_LOW_ENERGY_BATCH8.md:20` — blocking_reference: `- docs/T3_REDUCED_SCALAR_POTENTIAL_CONCLUSION.md — permanent restricted-potential conclusion`
- `output/studies/reduced_scalar_sobol/reduced_48_v1/large_uv_frozen_0p5/A_m2/fit.json:1749` — blocking_reference: `"zero_boundary_relevance_argument": "docs/T3_REDUCED_SCALAR_POTENTIAL_CONCLUSION.md",`
- `output/studies/reduced_scalar_sobol/reduced_48_v1/large_uv_frozen_0p5/A_m4/fit.json:1743` — blocking_reference: `"zero_boundary_relevance_argument": "docs/T3_REDUCED_SCALAR_POTENTIAL_CONCLUSION.md",`
- `output/studies/reduced_scalar_sobol/reduced_48_v1/large_uv_frozen_0p5/A_p0/fit.json:1749` — blocking_reference: `"zero_boundary_relevance_argument": "docs/T3_REDUCED_SCALAR_POTENTIAL_CONCLUSION.md",`
- `output/studies/reduced_scalar_sobol/reduced_48_v1/large_uv_frozen_0p5/B_m1/fit.json:1752` — blocking_reference: `"zero_boundary_relevance_argument": "docs/T3_REDUCED_SCALAR_POTENTIAL_CONCLUSION.md",`
- `output/studies/reduced_scalar_sobol/reduced_48_v1/large_uv_frozen_0p5/B_m3/fit.json:1745` — blocking_reference: `"zero_boundary_relevance_argument": "docs/T3_REDUCED_SCALAR_POTENTIAL_CONCLUSION.md",`
- `output/studies/reduced_scalar_sobol/reduced_48_v1/large_uv_frozen_0p5/B_p1/fit.json:1745` — blocking_reference: `"zero_boundary_relevance_argument": "docs/T3_REDUCED_SCALAR_POTENTIAL_CONCLUSION.md",`
- `output/studies/reduced_scalar_sobol/reduced_48_v1/large_uv_frozen_0p5/C_m1/fit.json:1752` — blocking_reference: `"zero_boundary_relevance_argument": "docs/T3_REDUCED_SCALAR_POTENTIAL_CONCLUSION.md",`
- `output/studies/reduced_scalar_sobol/reduced_48_v1/large_uv_frozen_0p5/C_m3/fit.json:1745` — blocking_reference: `"zero_boundary_relevance_argument": "docs/T3_REDUCED_SCALAR_POTENTIAL_CONCLUSION.md",`
- `output/studies/reduced_scalar_sobol/reduced_48_v1/large_uv_frozen_0p5/C_p1/fit.json:1745` — blocking_reference: `"zero_boundary_relevance_argument": "docs/T3_REDUCED_SCALAR_POTENTIAL_CONCLUSION.md",`
- `output/studies/reduced_scalar_sobol/reduced_48_v1/large_uv_frozen_0p5/D_m2/fit.json:1749` — blocking_reference: `"zero_boundary_relevance_argument": "docs/T3_REDUCED_SCALAR_POTENTIAL_CONCLUSION.md",`
- `output/studies/reduced_scalar_sobol/reduced_48_v1/large_uv_frozen_0p5/D_p0/fit.json:1749` — blocking_reference: `"zero_boundary_relevance_argument": "docs/T3_REDUCED_SCALAR_POTENTIAL_CONCLUSION.md",`
- `output/studies/reduced_scalar_sobol/reduced_48_v1/large_uv_frozen_0p5/D_p2/fit.json:1743` — blocking_reference: `"zero_boundary_relevance_argument": "docs/T3_REDUCED_SCALAR_POTENTIAL_CONCLUSION.md",`
- `output/studies/reduced_scalar_sobol/reduced_48_v1/large_uv_frozen_0p5/E_m2/fit.json:1750` — blocking_reference: `"zero_boundary_relevance_argument": "docs/T3_REDUCED_SCALAR_POTENTIAL_CONCLUSION.md",`
- `output/studies/reduced_scalar_sobol/reduced_48_v1/large_uv_frozen_0p5/E_m4/fit.json:1743` — blocking_reference: `"zero_boundary_relevance_argument": "docs/T3_REDUCED_SCALAR_POTENTIAL_CONCLUSION.md",`
- `output/studies/reduced_scalar_sobol/reduced_48_v1/large_uv_frozen_0p5/E_p0/fit.json:1750` — blocking_reference: `"zero_boundary_relevance_argument": "docs/T3_REDUCED_SCALAR_POTENTIAL_CONCLUSION.md",`
- `output/studies/reduced_scalar_sobol/reduced_48_v1/large_uv_frozen_0p5/E_p2/fit.json:1743` — blocking_reference: `"zero_boundary_relevance_argument": "docs/T3_REDUCED_SCALAR_POTENTIAL_CONCLUSION.md",`
- `output/studies/reduced_scalar_sobol/reduced_48_v1/small_uv_frozen_5e-9/A_m2/fit.json:1749` — blocking_reference: `"zero_boundary_relevance_argument": "docs/T3_REDUCED_SCALAR_POTENTIAL_CONCLUSION.md",`
- `output/studies/reduced_scalar_sobol/reduced_48_v1/small_uv_frozen_5e-9/A_m4/fit.json:1682` — blocking_reference: `"zero_boundary_relevance_argument": "docs/T3_REDUCED_SCALAR_POTENTIAL_CONCLUSION.md",`
- `output/studies/reduced_scalar_sobol/reduced_48_v1/small_uv_frozen_5e-9/A_p0/fit.json:1749` — blocking_reference: `"zero_boundary_relevance_argument": "docs/T3_REDUCED_SCALAR_POTENTIAL_CONCLUSION.md",`
- `output/studies/reduced_scalar_sobol/reduced_48_v1/small_uv_frozen_5e-9/B_m1/fit.json:1752` — blocking_reference: `"zero_boundary_relevance_argument": "docs/T3_REDUCED_SCALAR_POTENTIAL_CONCLUSION.md",`
- `output/studies/reduced_scalar_sobol/reduced_48_v1/small_uv_frozen_5e-9/B_m3/fit.json:1745` — blocking_reference: `"zero_boundary_relevance_argument": "docs/T3_REDUCED_SCALAR_POTENTIAL_CONCLUSION.md",`
- `output/studies/reduced_scalar_sobol/reduced_48_v1/small_uv_frozen_5e-9/B_p1/fit.json:1745` — blocking_reference: `"zero_boundary_relevance_argument": "docs/T3_REDUCED_SCALAR_POTENTIAL_CONCLUSION.md",`
- `output/studies/reduced_scalar_sobol/reduced_48_v1/small_uv_frozen_5e-9/C_m1/fit.json:1752` — blocking_reference: `"zero_boundary_relevance_argument": "docs/T3_REDUCED_SCALAR_POTENTIAL_CONCLUSION.md",`
- `output/studies/reduced_scalar_sobol/reduced_48_v1/small_uv_frozen_5e-9/C_m3/fit.json:1745` — blocking_reference: `"zero_boundary_relevance_argument": "docs/T3_REDUCED_SCALAR_POTENTIAL_CONCLUSION.md",`
- `output/studies/reduced_scalar_sobol/reduced_48_v1/small_uv_frozen_5e-9/C_p1/fit.json:1622` — blocking_reference: `"zero_boundary_relevance_argument": "docs/T3_REDUCED_SCALAR_POTENTIAL_CONCLUSION.md",`
- `output/studies/reduced_scalar_sobol/reduced_48_v1/small_uv_frozen_5e-9/D_m2/fit.json:1749` — blocking_reference: `"zero_boundary_relevance_argument": "docs/T3_REDUCED_SCALAR_POTENTIAL_CONCLUSION.md",`
- `output/studies/reduced_scalar_sobol/reduced_48_v1/small_uv_frozen_5e-9/D_p0/fit.json:1749` — blocking_reference: `"zero_boundary_relevance_argument": "docs/T3_REDUCED_SCALAR_POTENTIAL_CONCLUSION.md",`
- `output/studies/reduced_scalar_sobol/reduced_48_v1/small_uv_frozen_5e-9/D_p2/fit.json:1743` — blocking_reference: `"zero_boundary_relevance_argument": "docs/T3_REDUCED_SCALAR_POTENTIAL_CONCLUSION.md",`
- ... plus 85 additional references in `references.csv`.

### `docs/T3_ZERO_BOUNDARY_LOW_ENERGY_BATCH8.md`

- `docs/T3_REDUCED_SCALAR_POTENTIAL_CONCLUSION.md:62` — blocking_reference: `The next **numerical continuation** reuses the existing model-specific hierarchical fixed-order \(C_5\) and SM+Weinberg running (Stage 36), but e`
- `docs/T3_ZERO_BOUNDARY_LOW_ENERGY_BATCH8.md:19` — blocking_reference: `- tests/python/test_t3_zero_boundary_low_energy_batch8.py — offline tests`
- `docs/T3_ZERO_BOUNDARY_LOW_ENERGY_BATCH8.md:21` — blocking_reference: `- docs/T3_ZERO_BOUNDARY_LOW_ENERGY_BATCH8.md — this implementation and interpretation guide`
- `docs/T3_ZERO_BOUNDARY_LOW_ENERGY_BATCH8.md:30` — blocking_reference: `python -m pytest -q tests\python\test_t3_zero_boundary_low_energy_batch8.py`

### `docs/T3_GENERAL_SCALAR_STUDY_SEPARATION.md`

- `output/studies/general_scalar_extended/README.md:19` — blocking_reference: `See docs/T3_GENERAL_SCALAR_STUDY_SEPARATION.md.`
- `output/studies/scalar_potential_boundary_audit/restricted_input_audit.md:24` — blocking_reference: `Store new results only under output/studies/general_scalar_extended/ with distinct run identifiers, input hashes, RG/matching scope and original-`
- `scripts/AuditT3SavedStudyScalarBoundaries.py:306` — blocking_reference: `'Store new results only under output/studies/general_scalar_extended/ with distinct run identifiers, input hashes, RG/matching scope and original`

## Safe procedure

1. Review `audit.json`, `references.csv` and this report, including blocked files.
2. Do not archive any blocked file until references are updated in full replacement files.
3. To archive only the safe group, explicitly run `--apply --confirm ARCHIVE_SAFE_ONLY`.
4. This copies every selected original to the external archive, verifies SHA-256, then removes its active copy.
5. Rerun `python -X utf8 -m pytest -q --tb=line`; do not remove fitting data or source code.

