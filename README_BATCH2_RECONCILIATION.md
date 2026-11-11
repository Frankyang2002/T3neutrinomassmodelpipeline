# Batch 2: A(-4)/D(-2) source inventory reconciliation

**Scope:** Read-only reconciliation of the uploaded Matchete source evidence, not a renormalized one-loop matching implementation.

## Evidence checked

| Model | Full ledger terms | One-loop hbar terms | Candidate LLSS | Candidate LLHH | Other EFT1 | Laurent pole-residue terms |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| A(-4) | 135 | 131 | 48 | 6 | 81 | 45 |
| D(-2) | 1200 | 1194 | 352 | 12 | 836 | 403 |

The one-loop term counts and Laurent-source counts agree. Source sector assignments are those of the saved ledgers and should not be assumed to form a complete independent physical EFT basis.

## Verified limitations (saved evidence flags)

- The Laurent audit has not performed UV-hard minus EFT1-soft subtraction or added complete local MSbar counterterms.
- The EFT1 scalar bubble's pole is cancelled by its own counterterm. This does not establish full fermion matching.
- Direct comparison of the native pole and the isolated EFT1 scalar counterterm has relative discrepancy approximately 0.463 for A(-4) and 1.681 for D(-2). This is **not** a valid complete renormalization consistency test: relevant UV/EFT counterterms and maps remain absent.
- The saved pole-closure report has six failed physical matching gates for both benchmark models: verified UV/EFT1 hard-soft matching, complete scalar threshold, complete hard topology inventory, native UV counterterms, identification of the native expression as bare UV 1PI, and the signed GammaCC-Weyl map.

## Next actual physics implementation

1. Fix the exact spinor operator and scalar component normalisation by comparing a chosen LLSS amplitude in the Matchete and EFT1 conventions, including signs and factors of two.
2. Identify the amplitude/scheme represented by the saved Matchete `MatchedEFT` expression; do not assume it is a bare UV amplitude.
3. Compute UV and EFT1 amplitudes using the same external states and subtraction convention, project onto the retained EFT1 basis, account for other relevant dim-5 classes, and derive a finite Wilson boundary.
4. Only after independent checks, inject that boundary into the existing DOP853 solve. Do not refit or rerun all 16 models before A(-4) and D(-2) are certified.

## Install and rerun read-only validation (PowerShell)

```powershell
$Project = 'C:\thesistemp\T3neutrinomassmodelpipeline'
$Downloads = (New-Object -ComObject Shell.Application).Namespace('shell:Downloads').Self.Path
$Zip = Join-Path $Downloads 'T3_Batch2_Reconciliation.zip'
$Evidence = Join-Path $Downloads 'T3_Batch2_Matchete_Evidence.zip'
Expand-Archive -LiteralPath $Zip -DestinationPath $Project -Force
Set-Location -LiteralPath $Project
$Out = Join-Path $Project ('output\audits\batch2_source_reconciliation_' + (Get-Date -Format 'yyyyMMdd_HHmmss') + '.json')
python -X utf8 scripts\AuditT3MatchingEvidenceArchive.py --archive "$Evidence" --out "$Out"
Get-Content -LiteralPath $Out
```

This archive adds only `scripts/AuditT3MatchingEvidenceArchive.py`, `README_BATCH2_RECONCILIATION.md` and `BATCH2_EVIDENCE_AUDIT.json`. No production script is replaced. No numerical outputs are modified.
