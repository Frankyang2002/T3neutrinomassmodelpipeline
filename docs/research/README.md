# Technical research records

Supplementary T3 calculations, investigations and implementation notes are
organised here by scientific topic rather than by the sequence of software
revisions. The [research document index](INDEX.md) lists relocated documents
and their current relative paths.

| Folder | Scientific content |
| --- | --- |
| `model-and-symmetry/` | Gauge representations, group-theory invariants, Yukawa interactions and symmetry completeness |
| `eft-matching/` | Integrating out heavy fermions, scalar mixing and thresholds, operator matching and Weinberg coefficients |
| `numerical-studies/` | UV/intermediate running, parameter scans, benchmark comparisons and sensitivity analyses |
| `scientific-reports/` | Figures, tables and thesis-oriented reporting investigations |
| `research-notes/` | Additional model-specific technical studies |
| `project-records/` | Historical test, build and repository-management records |

The principal physics reference remains [Model and calculation contract](../MODEL_AND_PIPELINE.md),
with [Physics](../PHYSICS.md), [Scalar sector](../SCALAR_SECTOR.md),
[Numerical methods](../NUMERICAL.md) and [Results](../RESULTS.md).
The consolidated [research history](../RESEARCH_HISTORY_INDEX.md) is preserved.

The intended EFT sequence is **UV → integrate out all heavy fermions →
scalar-containing intermediate EFT → integrate out all heavy scalars → SMEFT**.
Archived experiments or conditional numerical outputs are not promoted to
certified general-model predictions by being moved into these folders.

Only documents with no detected incoming source references and no relative
file links are moved automatically; all other files remain at their original
paths. File contents and names are preserved for scientific traceability.

## Moving files in a local checkout

From the repository root, create readable indexes, then preview the proposed
relocations. The plan is stored under `output/maintenance/` and must be
reviewed before any move.

```powershell
python -X utf8 scripts/maintenance/OrganizeResearchFiles.py --index
python -X utf8 scripts/maintenance/OrganizeResearchFiles.py
```

After reviewing the approved source and destination paths:

```powershell
python -X utf8 scripts/maintenance/OrganizeResearchFiles.py `
    --apply --confirm ORGANIZE_RESEARCH_FILES
.\scripts\RunCoreTests.ps1
```

For an explicit, optional review of the isolated experimental script files,
use `--include-experiments` on the **preview**. The approved script paths, if
any, are then included in the same reviewed apply operation. Do not use this
option for the ordinary documentation-only reorganisation.

Restoration is reversible while the relocated files are unchanged:

```powershell
python -X utf8 scripts/maintenance/OrganizeResearchFiles.py --restore
```

The local manifest records the original and new path and SHA-256 of each
moved file. The tool never commits, pushes, deletes generated physics outputs,
or changes source file contents.
