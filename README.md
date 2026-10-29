# T3 Batch 1 — bulk continuation (existing repository)

## Exact destination

Copy `scripts/RunT3FullScalarBatch1.ps1` into the **existing** repository at `scripts/RunT3FullScalarBatch1.ps1`. No existing files need to be replaced. This deliberately does **not** override the established T3 production Lagrangian, Matchete or RGBeta code.

## What it does

Runs existing `scripts.AuditGeneralT3QuadraticBasis`, the existing 16-model general quartic RGBeta generator, and `scripts.RunGeneralT3QuadraticAudit` in one batch. When passed `-ExportBetas`, it asks locally installed Wolfram/RGBeta to regenerate UV and EFT1 beta expressions for **all 16 models**. When passed `-RunTests`, it runs all Python tests. Results are written under `output/audits/full_scalar_bulk_batch1/`, and the generator writes its established outputs under `output/general_t3_stage20/`. Audit scripts may update their existing audit artifact paths. None is a physical full-sequential-matching implementation.

## PowerShell from repository root

```powershell
powershell -ExecutionPolicy Bypass -File .\scripts\RunT3FullScalarBatch1.ps1 -ExportBetas -RunTests
```

If Wolfram is unavailable, run the Python structural stages only:

```powershell
powershell -ExecutionPolicy Bypass -File .\scripts\RunT3FullScalarBatch1.ps1 -RunTests
```

Inspect `output/audits/full_scalar_bulk_batch1/batch1_status.json` afterward.

## Known scientific limits

- The additional 50 quartic entries are summed over 16 **different models**, not 50 couplings per model.
- Eight additional bilinears remain outside the complete coupled scalar RGE/matching pipeline.
- Matchete/RGBeta coupling normalisation is not established for the complete basis.
- Native RGBeta export is not a physical one-loop C5-matching certificate.
- No existing neutrino-fit results are relabelled as fully matched.

Read-only source snapshot used for the batch runner: `119f5d46f89159470af59bff3189ab65424ef280`. Run against the user's latest local checkout; confirm differences if any.
