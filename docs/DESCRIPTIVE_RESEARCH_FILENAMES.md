# Descriptive filenames for numbered T3 research stages

## Scope and source

Read-only inspected baseline: Git commit `fe0f511104cce8d7d914d114e3fd884eb47c6a95` (29 October 2026). This migration executes **only on the user's local checkout**. It does not access or modify GitHub. The supplied migration is a **full standalone Python file**, not a patch. It generates renamed copies of the full original files from the user's own checkout.

The existing numbered research names refer to genuinely different calculations. The rename changes **source filenames and their exact references**, not the calculations themselves.

| Original source | Descriptive filename |
|---|---|
| `scripts/Stage30CoupledA0UVRunning.py` | `scripts/CoupledA0UVRunning.py` |
| `scripts/Stage31CoupledUVRunning.py` | `scripts/CoupledUVRunning.py` |
| `scripts/Stage33FermionThresholdHandoff.py` | `scripts/FermionThresholdHandoff.py` |
| `scripts/Stage34IntermediateEFTRunning.py` | `scripts/IntermediateEFTRunning.py` |
| `scripts/Stage35ScalarThresholdHandoff.py` | `scripts/ScalarThresholdHandoff.py` |
| `scripts/Stage36FinalWeinbergRunning.py` | `scripts/FinalWeinbergRunning.py` |
| `scripts/RunGeneralT3Stage16.py` | `scripts/RunGeneralT3FermionFirstThresholds.py` |
| `scripts/RunGeneralT3Stage17.py` | `scripts/RunGeneralT3EFT1ScalarPairMixing.py` |
| `scripts/RunGeneralT3Stage18.py` | `scripts/RunGeneralT3EFT1WilsonSourceTransport.py` |
| `scripts/RunGeneralT3Stage19.py` | `scripts/RunGeneralT3WeinbergSourceTransport.py` |
| `scripts/Stage22TensorAlgebra.py` | `scripts/RGBetaTensorAlgebra.py` |
| `docs/GENERAL_T3_STAGE17_STATUS.md` | `docs/GENERAL_T3_EFT1_SCALAR_PAIR_MIXING.md` |
| `tests/python/test_general_t3_stage36_final_weinberg.py` | `tests/python/test_general_t3_final_weinberg.py` |

The full mapping of all eligible `StageNN` filenames (source scripts, PowerShell/Wolfram launchers, research notes, tests) is generated from actual **tracked files** in the local checkout. The tool refuses naming collisions or unexpected target paths. Supporting source references to the exact module, script and document filenames are rewritten locally.

### Intentionally excluded

- `output/**` and `Reports/output/**`, especially **saved scientific results**, historical stage-named RGE/matching exports and original numerical provenance.
- `docs/archive/**` and test fixtures, whose content is historical evidence.
- Batch-numbered filenames, since `Batch` labels are not part of this request.
- Narrative historical labels (`Stage 16`, `Stage 17`), physics variable names, and old *saved-data* paths (`output/general_t3_stage17/...`): the chronology and data contracts must remain reconstructible.
- JSON inputs and exported data, which must retain their original source fingerprints.

**Important:** Because some historical JSON files record script or documentation paths, renaming the active scripts cannot retroactively make every old provenance path valid. Original source files are preserved in the external backup and in Git history; do not rewrite saved numerical provenance to hide those historical references. This migration therefore does not certify old artifact replay. Regenerate reports from current paths only after running tests and checking source hashes.

## Local PowerShell procedure

Start from the repository root and use the existing Python virtual environment:

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline

# Optional: verify audit tool tests first.
python -X utf8 -m pytest -q tests/python/test_descriptive_research_file_names.py --tb=line

# Read-only plan, with a machine-readable mapping saved under output/.
python -X utf8 -m tools.RenameNumberedResearchFiles `
    --report output/cleanup_audit/descriptive_research_filenames/plan.json

# Review this JSON before applying.
Get-Content -Encoding UTF8 output/cleanup_audit/descriptive_research_filenames/plan.json |
    Select-Object -First 75

# Apply (modifies local checkout; NO GitHub writes).
python -X utf8 -m tools.RenameNumberedResearchFiles `
    --apply --confirm RENAME_SOURCE_FILES
if ($LASTEXITCODE -ne 0) { throw "Rename failed; inspect backup and error" }

# All existing tests; no symbolic recalculation or optimization.
python -X utf8 -m pytest -q --tb=line
if ($LASTEXITCODE -ne 0) { throw "Tests failed; do not commit rename" }

git status --short
```

The command rejects a changed HEAD, source files with uncommitted local modifications, existing rename destinations and a pre-existing backup location. It permits unrelated uncommitted deletions from the previous documentation cleanup. **Do not run `git reset --hard`** to make the rename work.

## Backup and restoration

Before editing anything, the tool copies every affected original file to the directory adjacent to the project:

`C:\thesistemp\T3_stage_filename_backup_fe0f511\originals\...`

The external backup's `manifest.json` records the SHA-256 of each original, new destination and rewritten bytes. The copy is verified before source renaming begins. On a failed application, the tool attempts to restore the original files automatically and leaves the backup intact.

To restore manually, first inspect `manifest.json`, remove only the renamed destinations listed under `renames`, and then copy the contents of `originals/` back into the project root while preserving paths. Do not delete unrelated user changes. If tests fail, retain the backup and do not commit until the failure is diagnosed.

## Testing scope

The package's seven tests cover naming maps, Python import/CLI path rewrites, preservation of saved-data paths, document titles, and full backup/application in a synthetic Git checkout. They are **not** a replacement for the Windows repository's full Python tests or external Wolfram/Matchete/RGBeta regression. Paths constructed dynamically from stage numbers may require manual follow-up if the full tests reveal them.

### Files added by this package

- `tools/RenameNumberedResearchFiles.py` — standalone complete local migration tool.
- `tests/python/test_descriptive_research_file_names.py` — seven offline tests.
- `docs/DESCRIPTIVE_RESEARCH_FILENAMES.md` — this guide.

No previously archived documents, physics scripts, saved data, numerical configurations or GitHub files are included in the ZIP.
