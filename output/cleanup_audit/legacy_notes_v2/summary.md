# Legacy root README notes: local read-only audit

- HEAD: `b521201eeff037661553949981191138df49691a`
- Text source files scanned: **1083**
- Read errors: **0**

| File | Classification | References |
|---|---|---:|
| `README_C12_C5.txt` | `ARCHIVE_CANDIDATE` | 0 |
| `README_C12_GroupFactors_PDF.txt` | `ARCHIVE_CANDIDATE` | 0 |
| `README_C12_PDF_GroupFactor_Repair.txt` | `ARCHIVE_CANDIDATE` | 0 |
| `README_C12_Representation_GroupFactors_Existing_HTML.txt` | `ARCHIVE_CANDIDATE` | 0 |
| `README_NEXT_STEP.txt` | `ARCHIVE_CANDIDATE` | 0 |
| `README_T3BOperatorClosure.txt` | `ARCHIVE_CANDIDATE` | 0 |
| `README_T3B_CommonScale_Thresholds.txt` | `ARCHIVE_CANDIDATE` | 0 |
| `README_T3B_GroupFactors_HTML.txt` | `ARCHIVE_CANDIDATE` | 0 |
| `README_T3B_Integrated_OptIn.txt` | `ARCHIVE_CANDIDATE` | 0 |
| `README_T3B_Low_Energy_Comparison.txt` | `ARCHIVE_CANDIDATE` | 0 |
| `README_T3B_Perturbative_Audit_Fix.txt` | `ARCHIVE_CANDIDATE` | 0 |
| `README_T3B_Perturbative_Sidecar.txt` | `ARCHIVE_CANDIDATE` | 0 |
| `README_T3_C12_C5_Mixing_GroupFactors_PDF.txt` | `ARCHIVE_CANDIDATE` | 0 |
| `README_T3_C12_PDF_Table_Fix.txt` | `ARCHIVE_CANDIDATE` | 0 |
| `README_T3_C12_Reverse_Mixing_PDF.txt` | `ARCHIVE_CANDIDATE` | 0 |
| `README_THESIS_COMPARISONS.txt` | `RETAIN_THESIS_PROVENANCE` | 8 |

## Reference details

### `README_THESIS_COMPARISONS.txt`
- `docs/LEGACY_NOTES_CLEANUP.md:9` — `'README_THESIS_COMPARISONS.txt' is protected: it describes thesis comparison`
- `tools/AuditLegacyNotes.py:26` — `PROTECTED_NOTE_NAMES = {'README_THESIS_COMPARISONS.txt'}`
- `tools/AuditLegacyNotes.py:166` — `'- 'README_THESIS_COMPARISONS.txt' (thesis asset provenance).',`
- `tools/AuditLegacyNotes.py:204` — `'retained_unchanged': ['README_THESIS_COMPARISONS.txt',`
- `tests/python/test_legacy_notes_cleanup.py:31` — `(root/'README_THESIS_COMPARISONS.txt').write_text('Thesis provenance\n')`
- `tests/python/test_legacy_notes_cleanup.py:43` — `'README_THESIS_COMPARISONS.txt']`
- `tests/python/test_legacy_notes_cleanup.py:57` — `assert rows['README_THESIS_COMPARISONS.txt']['classification'] == 'RETAIN_THESIS_PROVENANCE'`
- `tests/python/test_legacy_notes_cleanup.py:72` — `assert (repo/'README_THESIS_COMPARISONS.txt').exists()`

## Explicitly retained / manual review

- Root `README.txt`, main `README.md`, `SHA256_MANIFEST.json`.
- `README_THESIS_COMPARISONS.txt` (thesis asset provenance).
- All three scientific-output ZIP archives.
- Four legacy PS1 utilities: they are scripts, not one-time console logs.
- Both existing `report_cleanup_archive` directories and their manifests.
- Historical `output/` and `Reports/` numerical/research artifacts.

No files are changed in dry-run mode. Apply only archives unreferenced, clean, tracked notes.
The external archive records each original path and SHA-256.
