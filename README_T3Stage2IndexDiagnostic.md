# T3 hierarchical stage-2 index diagnosis and controlled common-mass experiment

## Cause identified from the actual `T3_Ap0_stage2_failure.zip`

The fresh-kernel input `input_tree.wxf` and `input_transition_tree.wxf` are the same
43-term tree-level EFT1 expression. Four monomials (1-based positions 33, 34,
42, 43) carry `MF[r]^(-1) yA[p,r] y2[q,r]` with **three occurrences of the
same typed `NFlavor` index `r`**, causing the confirmed
`FindDummyIndices::trippleindex` diagnostic at stage-2 tree matching.
Two of the terms involve `y1Conj`; two have baseline `y1`. Therefore setting
`y1Conj=0` does not by itself repair the representation of indexed masses.

The physical expression is a finite generation sum, for example
`W[p,q]=Sum[y1Conj[p,r] y2[q,r]/MF[r],{r,1,3}]`. Replacing it by a properly
defined two-Flavor-index Wilson tensor is the exact long-term repair. It
must also be carried consistently through inherited one-loop EFT matching and
EFT1 Wilson running; this package **does not invent that implementation**.

## What this package adds

- `scripts/DiagnoseT3Stage2Index.py` — safely parses uncompressed WXF; audits
  the actual 43-term tree EFT with no Wolfram installation or dependence on
  the model's numeric benchmark. Reports the four contractions and source
  couplings; refuses unsupported WXF tokens, truncated input, and changed
  audited term structure.
- `Lagrangian/RunT3ExtendedHierarchyNative.wl` — full-file replacement of the
  prior experimental driver with optional `common` mass mode. `indexed`
  remains the default; no change to production matching files.
- `scripts/RunT3ExtendedHierarchyBatch.py` — full-file replacement of prior
  batch runner, adds `--common-mass-diagnostic`, uses a separate output
  directory, and refuses common-mass reports unless they explicitly deny
  certification of independent heavy masses.
- `scripts/TestT3Stage2IndexDiagnostic.py` — synthetic and source-contract
  regression, plus the real WXF check when present.

## Commands (PowerShell, repository root)

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "$HOME\Downloads\T3Stage2IndexDiagnostic_Local.zip" -DestinationPath . -Force

python scripts\DiagnoseT3Stage2Index.py --self-test
python scripts\TestT3Stage2IndexDiagnostic.py
python scripts\DiagnoseT3Stage2Index.py --assert-a0 --output-json Reports\output\workstream_a\extended_c5_hierarchy\stage2_index_diagnosis.json

# CONTROLLED EXPERIMENT: three degenerate heavy flavours (one common MF).
# Do not confuse it with the previous three *independent* MF[r] benchmark.
python scripts\RunT3ExtendedHierarchyBatch.py --pilot --stage1 --common-mass-diagnostic
if ($LASTEXITCODE -eq 0) {
    python scripts\RunT3ExtendedHierarchyBatch.py --pilot --hierarchy --common-mass-diagnostic
}
$LASTEXITCODE
```

The common-mass experiment starts from the UV action and regenerates stage 1:
`T3_EXPERIMENTAL_FLAVOR_MASSES=0` but `NFlavor=3` is retained. Therefore it
removes indexed `MF[r]` through a consistent Lagrangian rebuild rather than
silently discarding it from an already-matched EFT. In that mode
`MF_1=MF_2=MF_3`, and `1/MF` is a scalar; the four identified triple-index
monomials are avoided. The Matchete stage-2 outcome **is not known until
this local command runs**. Other independent stage-2 issues may remain.

Common-mass outputs are isolated under:
`Reports/output/workstream_a/extended_c5_hierarchy/common_mass_diagnostic/`.
The original indexed-mass stage-1 results and blocked stage-2 outputs are
left intact. The new runner still accepts the original indexed-mass commands,
though stage 2 should not be rerun in that mode until the Wilson tensor issue
is addressed.

## Status and scientific scope

- **Confirmed from user WXF:** exact 43-term tree, four triple heavy-generation
  contractions, one loop log's explicit `FindDummyIndices::trippleindex`.
- **Tested in this environment:** restricted WXF parser, real uploaded WXF
  audit, Python compile and synthetic provenance checks.
- **Not locally verified:** new native common-mass Matchete run (requires
  Windows Wolfram/Matchete installation).
- **Not certified:** arbitrary three independent mass eigenvalues, complete
  intermediate-EFT running, UV-pole cancellation, absolute Matchete–RGBeta
  coupling conversion, finite physical hierarchical C5 or low-energy masses.

No GitHub changes. No existing production source files are replaced.
