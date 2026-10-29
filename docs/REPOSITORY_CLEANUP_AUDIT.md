# T3 repository cleanup: historical-document reference audit

**Phase:** documentation audit and optional external retirement.  **Source baseline:**
`f21d15d2d38f4f313a0f2214bc52f90e1469dce2` (verified read-only on GitHub).
The installed local checkout may contain the uncommitted authoritative documents
`docs/SCALAR_SECTOR.md`, `docs/NUMERICAL.md`, and `docs/RESULTS.md`; this utility
uses the **actual local files**, not the remote snapshot, for its decisions.

## Scope and safety boundary

The audited group consists of 14 older explanatory documents from the scalar
relevance batches, zero-boundary probes, and reduced-scalar final reporting.
The intended successors are:

- `docs/SCALAR_SECTOR.md`: physical scalar potential, 50 additional quartic and
  eight quadratic operator slots, 46/4 fixed-order vertex selection, radiative
  generation, quadratic mixing and scope limitations.
- `docs/NUMERICAL.md`: numerical evolution, Yukawa fitting, report generation,
  interactive RGE and sensitivity diagnostics, convergence/provenance rules.
- `docs/RESULTS.md`: final frozen reduced-scalar benchmark outcomes, comparison
  with common A(-2) Yukawas, and historical UV-running diagnostic limitations.

**The archive is outside the repository**, by default at:

`C:\thesistemp\T3_cleanup_archive\documentation_v1\`

Nothing in `Lagrangian/`, `RGE/`, `Numerical/`, `physics/`, `scripts/`,
`tests/`, `configs/`, `output/`, `Reports/` or the fitted result directories is
retired by this tool. In particular, production and historical scripts are
not classified by filename. Some apparently historical scripts remain imports
of current reporting code; deleting them now would break reproducibility.

The audit searches current local text files, both tracked and untracked, for
historical Markdown filenames and stems. A candidate is blocked by references
from a **retained** file. References between documents archived together are
allowed. A candidate which becomes blocked can in turn block another candidate:
this transitive dependency is checked before any move.

The three authoritative master documents often list old filenames solely as
**historical provenance**. Plain mentions in those three are recorded as
non-blocking; active Markdown links or navigation instructions remain blockers.
This special rule does not apply to any other document, Python script, test,
configuration, or report. The utility excludes itself, this audit guide and
its generated report files from reference counting, and never follows symlinked
source trees. Dynamic path construction or references hidden inside binary
artifacts cannot be proved absent by text search. Review the report before
archiving even files marked safe.

## Dry-run commands (PowerShell; no unsigned .ps1 execution)

From the project root with `.venv` active:

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
python -X utf8 -m tools.AuditDocumentationRetirement
Get-Content -Encoding UTF8 .\output\cleanup_audit\documentation_v1\summary.md
```

The audit creates only:

```text
output/cleanup_audit/documentation_v1/
    audit.json          # candidate hashes, source commit, all reference occurrences
    references.csv      # references, source filenames, line numbers, type
    summary.md          # safe/blocked/missing list and example blockers
```

The generated files are audit evidence, **not** physics inputs. A dry run moves
nothing. It does not affect the previously verified Python test baseline
(2,144 passed / 10 skipped / 1 warning in the user's checkout).

## After reviewing the report

If every safe file has been reviewed, the optional explicit command is:

```powershell
python -X utf8 -m tools.AuditDocumentationRetirement --apply --confirm ARCHIVE_SAFE_ONLY
```

This archives **only** the currently safe subset. It refuses to run if any
master document is absent or fails the SHA-256 fingerprint of the delivered
consolidated versions (with LF/CRLF normalised), if a source file could not be scanned,
if the archive destination is inside the checkout or already exists, or if a
candidate changed after the audit. For each archive member, its original SHA-256
is recorded; the copy's SHA-256 is checked before its original is removed.
The external archive writes `manifest.json` and keeps the original project
relative path (e.g. `docs/T3_SCALAR_VERTEX_BATCH2.md`).

Do **not** run `--apply` merely because the scanner reports some files safe:
first check `summary.md`, particularly whether a planned thesis source link
should remain accessible inside the repository. Files flagged blocked stay
in place. This tool does not automatically rewrite their references.

## Restoring an archived document

The archive does not modify a historical file's bytes. For a single retired
file, restore it manually from the archived path, for example:

```powershell
Copy-Item -LiteralPath 'C:\thesistemp\T3_cleanup_archive\documentation_v1\docs\T3_SCALAR_VERTEX_BATCH2.md' `
  -Destination '.\docs\T3_SCALAR_VERTEX_BATCH2.md'
```

Before restoring, ensure that destination is not already present and compare
the restored file's `Get-FileHash -Algorithm SHA256` against the manifest.
The original is retained until its verified archive copy has been written.

## Required post-archival validation

```powershell
git --no-pager status --short
python -X utf8 -m pytest -q --tb=line
```

The expected prior baseline is **2,144 passed, 10 skipped, one warning**.
Any new failure or changed skip count requires investigation. Following this
phase, use the reference audit results to update the remaining historical
links with **complete replacement files**, then audit again; do not delete
historical batch scripts or source numerical data as part of a documentation
cleanup.

**Installation prerequisite:** Both replacements from
`T3_Numerical_Results_Consolidation_f21d15d.zip` and the earlier
`docs/SCALAR_SECTOR.md` must already be installed. A dry run still reports
references if they are not, but `--apply` is disabled until all three master
documents match their delivered versions. This prevents retiring their
original sources while an old master document is still installed.
