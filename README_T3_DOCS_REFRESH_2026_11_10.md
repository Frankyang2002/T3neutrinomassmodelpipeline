# T3 documentation status refresh — 10 November 2026

This is a **documentation-only** package. It does **not** install a full-general EFT implementation, alter numerical output, change physics code or modify GitHub.

## Files included (complete files, relative to the local project root)

- `docs/IMPLEMENTATION_STATUS.md` — current code-based scope and what is actually executable.
- `docs/T3_A4_D2_GENERAL_PIPELINE_FEASIBILITY_2026_11_10.md` — two-model, full-general feasibility assessment.
- `docs/T3_A4_D2_FULL_GENERAL_IMPLEMENTATION_BATCHES.md` — proposed four major implementation batches with numerical acceptance criteria.
- `README_T3_THESIS_TWO_MODEL_PILOT.md` — restricted Sobol-v3 pilot correctly distinguished from future full-general fitting.
- `scripts/UpdateT3DocumentationStatus.py` — local, idempotent updater for the older overview pages.
- `tests/test_t3_documentation_status.py` — isolated updater tests.
- `README_T3_DOCS_REFRESH_2026_11_10.md` — this file.

The updater preserves the **full historical content and numerical tables** in old `README.md`, `docs/NUMERICAL.md`, `docs/RESULTS.md`, `docs/SCALAR_SECTOR.md`, `docs/MODEL_AND_PIPELINE.md`, `docs/PROJECT_GUIDE.md`, `docs/README.md`, `docs/T3_THESIS_PLAN.md`, `docs/T3_COMPLETE_GENERAL_SCALAR_SYMBOLIC.md`, `docs/T3_FULL_SCALAR_CONTINUOUS_BATCH.md`, `Numerical/parallel_t3/wilson/README.md`, and `README_T3_FINAL_REPRO_BATCH8.md`. It inserts a single visible current-scope note and changes a few dated, misleading `primary`/`authoritative` scope labels. It does **not** delete the old studies or claim they have been rerun.

## Local installation (PowerShell)

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "$HOME\Downloads\t3_documentation_refresh_2026_11_10.zip" -DestinationPath . -Force
.\.venv\Scripts\python.exe scripts\UpdateT3DocumentationStatus.py --apply
.\.venv\Scripts\python.exe scripts\UpdateT3DocumentationStatus.py --check
.\.venv\Scripts\python.exe -m pytest -q tests\test_t3_documentation_status.py
```

**No GitHub changes.** The script writes only to local Markdown files listed above; it does not modify any existing code/physics results. It rejects an incomplete status-block marker rather than overwrite uncertain content. Missing listed Markdown files are reported instead of silently invented.

## Source status

The audit references the read-only public `main` at `478c658c70a397cff7f04d397e6f7dd448179d01` (10 November 2026). The current Windows checkout may differ from GitHub; old documents not overwritten by this ZIP are modified **in place** by the updater so any additional local paragraphs remain intact. The existing v3 Sobol script continues to use restricted conditional matching.
