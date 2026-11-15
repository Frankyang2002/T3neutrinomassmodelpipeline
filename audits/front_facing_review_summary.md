# Repository presentation and consolidation review

**GitHub source checked:** `Frankyang2002/T3neutrinomassmodelpipeline`, `main` at `8719090f19e9f6a224b0e4eb1de1c642056277f0`.

**Inventory:** the 313 tracked filenames uploaded from the local Windows checkout on 15 November 2026. The original list excludes files beneath `output/` and `Reports/output/`. This is **not a full repository-wide dependency audit**. The local uncommitted files were not accessible.

## Classification (313/313 files)

| Proposed action | Files | Meaning |
| --- | ---: | --- |
| CONSOLIDATE_REVIEW | 77 | Review document content; extract scientific conclusions into canonical pages, then consider archiving redundant source notes. |
| REVIEW_CODE_DEPENDENCIES | 154 | Preserve Python/Wolfram/PowerShell implementation pending callers, import and output dependency checks. |
| PRESERVE_TEST | 72 | Retain regression tests until mapped to verified replacement paths. |
| PRESERVE_PROVENANCE_REVIEW | 8 | Retain JSON data/schema evidence pending source provenance review. |
| REVIEW_ARCHIVE_BACKUP | 2 | Two tracked `.before_v2` copies: candidate backups, not automatic deletions. |
| **Total** | **313** | All source entries classified. |

The largest concentration of development terminology is in executable sources and their tests, not the public-facing scientific workflow. Renaming modules directly would break imports or commands. The first cleanup should therefore target *reader-facing documents*, not scientific code.

### Priority 1: Consolidate documentation

Map the 77 historical documents to the existing permanent pages identified in `front_facing_review.csv`. Do **not** claim the contents have all been merged: only selected GitHub source documents were read in detail. Existing detailed derivations and checks should stay available until extracted and verified.

Examples checked against the source:

- `docs/WORKER4_OPERATOR_CLOSURE.md` contains information about the seven-term Wilson evolution and potentially omitted operator classes; carry it into `docs/IMPLEMENTATION_STATUS.md` before archival.
- `docs/WorkstreamC_FinalConsolidated_v23_README.md` explicitly describes a restricted approximate benchmark, not a verified physical matching coefficient.
- `docs/T3_THESIS_FULL_GENERAL_BATCH1_HANDOVER.md` covers UV and scalar-EFT running but does not establish general one-loop scalar/fermion matching.
- `Lagrangian/RunT3Batch2LaurentAudit.wl` is a Laurent audit, not a renormalised matching calculation.

### Priority 2: Keep scientific implementations stable

There are 154 code-related entries and 72 test entries. In particular, `Numerical/experimental/Batch14EFT1Running.py` imports `Batch13GeneralUVRunning.py`; renaming just one would break the code. For public usage, leave the implementation under the hood and surface a small number of scientifically named commands.

### Priority 3: Check exact dependencies locally

Run the included `scripts/maintenance/AuditFrontFacingReferences.py` from the repository root. It creates a read-only reference report under `audits/`. A zero-reference result is *not* permission to delete scientific modules: dynamic loads, Wolfram script calls, configuration references and saved result provenance may not be visible to a literal scan.

**Nothing in this package automatically archives or removes the 313 files.** The CSV is a review manifest and its `archive_approved` column is `NO` throughout.

## Scientific public interface

Preserve the real EFT sequence:

```text
UV theory (SM + fermions + scalars)
    -> integrate out all heavy fermions
Scalar EFT (SM + remaining heavy scalars + effective operators)
    -> integrate out all heavy scalars
SMEFT (SM + relevant Wilson operators)
    -> evolve and extract neutrino observables
```

The above is the intended physical hierarchy; present source limitations mean it should not be advertised as a universally verified, complete, fully general one-loop calculation for all 16 models.

## Scope of included replacements

- `README.md` and `docs/README.md` are **full front-facing replacements** based on inspected GitHub versions. They do not rename or alter calculation paths.
- `audits/front_facing_review.csv` classifies all 313 filenames without inventing usage relationships.
- `scripts/maintenance/AuditFrontFacingReferences.py` audits references locally and writes a CSV. It performs no deletions.

Review any local modifications before installing the replacement documentation.
