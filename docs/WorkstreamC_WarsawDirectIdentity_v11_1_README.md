# Workstream C: Native Warsaw direct-identity probe v11.1

## Verified failure in v11

The installed Matchete package loads, but v11 reported `SMModelLoadFailed` before building either operator. The original v11 code called `LoadModel["SM"]` without the `ResetAll[]` initialization that the project's working `Lagrangian/model/LagrangianBuilder.wl` performs first. It also wrapped the model load in `Quiet@Check`, which suppressed useful diagnostic messages and could interpret a nonfatal message as a load failure. This is a **probe failure**, not evidence of incorrect B/C matching.

## Replacement

- `Lagrangian/ProbeT3BCWarsawDirectIdentity.experimental.wl` — **full replacement**, v11.1. It calls `ResetAll[]` before `LoadModel["SM"]` and, only if that fails, retries with `UsingFrontEnd` after another reset. It retains model-load messages. It compares the direct Warsaw source structure with `NCM[Transp[l], DiracProduct[GammaCC,Proj[-1]], l]` using the Matchete simplifier. Its JSON output retains the v11 schema accepted by the existing Python consumer.
- `docs/WorkstreamC_WarsawDirectIdentity_v11_1_README.md` — these notes.

Existing `scripts/ResolveT3BCWeinbergOperator.experimental.py` from v11 is unchanged, as are all production matching files.

## PowerShell

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "$HOME\Downloads\T3_WorkstreamC_WarsawDirectIdentity_v11_1.zip" -DestinationPath . -Force

$dir = "output\audits\t3_bc_y_warsaw_direct_identity_v11"
New-Item -ItemType Directory -Force $dir | Out-Null
wolframscript -file Lagrangian\ProbeT3BCWarsawDirectIdentity.experimental.wl `
    "$dir\native_identity.json"


$probeExit = $LASTEXITCODE
Get-Content "$dir\native_identity.json"
if ($probeExit -eq 0) {
    python scripts\ResolveT3BCWeinbergOperator.experimental.py `
        --matchete-identity "$dir\native_identity.json" `
        --output "$dir\summary_live.json"
}
```

`ExactPlusEquality` is the only accepted live certification for the direct operator identity. If `NotProvenByNativeSimplifier`, `ExactMinusEquality`, or a load failure is reported, **do not** force certification; the `native_identity.json` report contains the compared expressions and/or loader diagnostics. Upload that report for diagnosis. Capture `$LASTEXITCODE` immediately after `wolframscript` before running `Get-Content`.

## Remaining boundaries

An exact native operator identity fixes only the relative native/Warsaw operator convention. It does not prove original EOMEFT versus independent canonical EFT1 equivalence, full EFT1/EFT2 counterterms, the subtracted finite coefficient, X-threshold matching or physical C5. The source-level Weyl/Grassmann identity already passed under v11's Python self-test.
