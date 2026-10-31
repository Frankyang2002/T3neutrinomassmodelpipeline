# Final repository entry-point consolidation

Install all files into `C:\thesistemp\T3neutrinomassmodelpipeline\` preserving paths.
`README.md` in the ZIP is **a proposed full replacement**; the safer installation
method is the script below, which archives the current README before changing it.
Do not directly extract the ZIP over your existing README unless you have already
saved the original. For the safe method, extract **only**
`scripts/maintenance/FinishRepositoryConsolidation.py` and
`README_REPLACEMENT.md` and this guide from the ZIP. (PowerShell commands below
extract them automatically.)

```powershell
$zip = 'C:\path\to\T3_Final_Source_Consolidation.zip'
$repo = 'C:\thesistemp\T3neutrinomassmodelpipeline'
Expand-Archive $zip -DestinationPath "$repo\_consolidation_install" -Force
New-Item -ItemType Directory -Force "$repo\scripts\maintenance" | Out-Null
Copy-Item "$repo\_consolidation_install\scripts\maintenance\FinishRepositoryConsolidation.py" "$repo\scripts\maintenance\"
Copy-Item "$repo\_consolidation_install\README_REPLACEMENT.md" "$repo\README_REPLACEMENT.md"
python "$repo\scripts\maintenance\FinishRepositoryConsolidation.py" --root $repo
Get-Content "$repo\consolidation_preview.json"
python "$repo\scripts\maintenance\FinishRepositoryConsolidation.py" --root $repo --apply
```

The script leaves old code, generated matching files, tests and configurations
untouched. It inventories residual temporary candidates rather than deleting
files that might be runtime dependencies.

After successful regeneration of the saved numerical inputs (or restoration
from archive), verify:

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
$env:PYTHONUTF8='1'
python scripts/maintenance/SmokeTestFullPipeline.py
git status --short
```

**Important:** archiving `output/` and `Reports/output/` can make the smoke test
fail because it consumes saved generated inputs. This package is a documentation
and repository-layout cleanup, not a full clean-room Mathematica/RGBeta rebuild.
Do not claim full regenerability until the required source generation commands
have been executed in the local Wolfram environment.
