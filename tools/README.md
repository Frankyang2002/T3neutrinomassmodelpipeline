# Scientific utilities and repository tools

`tools/` is distinct from the physics implementation. All 18 original Python files here were inspected read-only at GitHub `main` commit `806cdf08608dbe3f5e0dcd8af785cfad16f26af6` (8 November 2026). They are *not* all obsolete, and none is a substitute for native Matchete/RGBeta matching.

## Saved-data reporting: retain

| File | Actual role in inspected source |
| --- | --- |
| `RebuildCurrentReports.py` | Rebuild the supported current figures and HTML from saved numerical data; imports figure redirect and label helpers. |
| `RebuildReports.py` | Rebuild selected fermion-threshold and Weinberg figures from saved results. |
| `RebuildSavedStudyPlots.py` | Selectively rebuild saved-data comparison plots. |
| `RedirectStudyFigures.py` | Context manager for routing Matplotlib files into a presentation directory; used by reporting utilities. |
| `LabelSensitivityDirections.py` | Adds scientific parameter-deformation definitions to saved sensitivity displays; imported by `RebuildCurrentReports.py`. |
| `PlotDisplayWrapper.py` | Applies consistent model labels to an existing plotting module; invoked by report rebuild workflows. |
| `RestoreClassicRGEReport.py` | Reconstructs an older, coloured RGE dashboard from saved data; historical presentation support, not physical computation. |
| `RestoreSharedA2LegacyReport.py` | Reconstructs the saved shared-A2 RGE dashboard from an existing template; specialised historical reporting. |

## Audit and scientific-source support: retain for now

| File | Role |
| --- | --- |
| `AuditScriptDependencies.py` | First-pass read-only inventory of `scripts/` imports, textual references and saved provenance. |
| `RefineScriptDependencyAudit.py` | AST-based second-pass classification of imports and literal invocations; imports the first-pass auditor. |
| `ConsolidateFittedBenchmarks.py` | Selects and validates canonical saved fit candidates; one-time-style data operation, but potentially important to benchmark reproducibility. |
| `ThesisRepositoryAudit.py` | New local-only per-file inventory, safe review and optional hash-verified external archive. |

## One-time migration tools: review individually for external archival

| File | Inspected purpose | Archival caveat |
| --- | --- | --- |
| `AuditDocumentationRetirement.py` | Historical note retirement with reference/manifest checks | May still be cited by prior audit instructions |
| `AuditLegacyNotes.py` | Old root README/study-note archival | Preserve prior audit provenance |
| `AuditOneTimeArtifacts.py` | One-time root log and terminal-capture archival | Its tests or docs may still refer to it |
| `ConsolidateHistoricalDocs.py` | Historical research-document consolidation | Old tests and reproducibility notes may depend on the tool |
| `FinishHistoricalDocsConsolidation.py` | Final migration of retained research notes into thematic volumes | Old tests and historical record references may still depend on the tool |
| `RenameNumberedResearchFiles.py` | Completed migration to descriptive script names | Old provenance may depend on the mapping and migration manifest |
| `MoveStudyFigures.py` | One-time PNG migration for full-study figures | Existing saved displays and documented instructions may depend on old paths |

**These are audit candidates, not approved deletions.** The local audit checks all tracked text, tests, code and saved outputs. Exact incoming references, local edits, and incomplete coverage prevent an automatic recommendation. Referencing docs alone requires a separate deliberate override, so migration history is not silently broken.

`tools/RebuildCurrentReports.py` is a good supported reporting entry point. Avoid using file names with `Legacy` or `Restore` to decide that a program is valueless: these can be the only sources that regenerate a figure used in the thesis.

See [Thesis file audit](../docs/THESIS_FILE_AUDIT.md) for preview and restore commands.
