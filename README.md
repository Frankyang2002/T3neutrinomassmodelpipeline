# Batch 12 — consolidate completed exploratory all-16 benchmark

This is the final **data aggregation** step. It reads (does not rerun or modify) the previously completed output directories. It verifies `Computed` for 16 models and each model's 16 U/E/W/S cases. It produces a compact CSV, manifest and Markdown overview.

## Installed paths

- `scripts/ConsolidateT3Exploratory16.py` — new
- `scripts/TestT3Exploratory16Consolidation.py` — new

## PowerShell

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -Path "$HOME\Downloads\T3ExploratoryAll16ConsolidationBatch12.zip" -DestinationPath . -Force
python -X utf8 -m scripts.TestT3Exploratory16Consolidation
python -X utf8 -m scripts.ConsolidateT3Exploratory16
```

Output directory: `output/studies/general_scalar_extended/exploratory_full_quartic_all16_consolidated_v1/` with `all16_cases.csv`, `all16_manifest.json`, `all16_overview.md`.

If an output directory exists, supply a fresh `--out` path; files are not overwritten.

## Physics caveats

- Matchete-to-RGBeta quartic mapping is provisional by invariant order.
- The six recovered ordinary models use frozen unsupported complex quadratic couplings.
- Nonmixed new extra quartics were set to zero at the F threshold rather than evolved from the UV.
- Simultaneous EFT1 Wilson running and matching therefore remain an **exploratory approximation**, not a complete UV scalar RGE treatment or normalized physical prediction.

The EFT sequence remains UV → integrate out all fermions → EFT1 → integrate out both scalars → SMEFT.
