# T3 native extended C5 — Matchete flavor-index fix

## Native failure established on Windows

After correcting Matchete's load-before-parse issue, the A(0) pilot reached
`build-full-general-fermion-quartic-Lagrangian` and emitted:

`DefineFlavorIndex::IndexDimension: The flavor index dimension '1' is not a positive integer larger than 1 or an undefined Symbol.`

The installed source `Lagrangian/model/T3Fields.wl` calls
`DefineFlavorIndex[NFlavor, Lookup[f, "Multiplicity", 3], ...]`. The prior
native pilot specified heavy flavor multiplicity 1, which the installed
Matchete flavor-index API disallows. The ordinary numerical UV state uses
three heavy generations.

## Implemented

- The complete native WL driver accepts only 2 or 3 heavy flavors and retains
  the preceding separate-top-level Matchete `Needs` load fix.
- The complete Python batch runner defaults to **3 heavy flavors**; 2 is an
  optional reduced-cost native test. It now refuses 1 **before** subprocess
  execution. A future API change can be diagnosed by the new regression.
- The complete synthetic integration test now uses the supported
  three-heavy-flavor case; a dedicated regression mocks all six native
  model invocations with both supported dimensions and rejects 1.
- **No** production RGBeta file, benchmark config, matching formula, original
  UV trajectory or physical Matchete↔RGBeta coupling dictionary is changed.

## Files in ZIP (relative to repository root)

1. `Lagrangian/RunT3ExtendedC5Native.wl` — FULL replacement.
2. `scripts/RunT3ExtendedMatchingBatch.py` — FULL replacement.
3. `scripts/TestT3ExtendedPhysicalStackSynthetic.py` — FULL replacement of earlier synthetic test.
4. `scripts/TestT3ExtendedFlavorIndexFix.py` — NEW native argument regression.
5. `README_T3ExtendedFlavorIndexFix.md` — this note.

## PowerShell

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
$zip = "$HOME\Downloads\T3ExtendedFlavorIndexFix_Local.zip"
Expand-Archive -LiteralPath $zip -DestinationPath . -Force

python scripts\RunT3ExtendedMatchingBatch.py --self-test
python scripts\TestT3ExtendedFlavorIndexFix.py
python scripts\TestT3ExtendedMatcheteLoadBoundary.py
python scripts\TestT3ExtendedPhysicalStackSynthetic.py

# Small native validation: two flavors satisfy the installed Matchete API.
python scripts\RunT3ExtendedMatchingBatch.py --pilot --heavy-flavours 2
$pilotExit = $LASTEXITCODE

# All six with physically relevant three-heavy-generation input, only
# after the small native pilot succeeds.
if ($pilotExit -eq 0) {
    python scripts\RunT3ExtendedMatchingBatch.py --all --heavy-flavours 3
}
$LASTEXITCODE
```

If the pilot fails at a **later stage**, consult
`Reports/output/workstream_a/extended_c5_native/logs/A_p0_h2.log` and
`Reports/output/workstream_a/extended_c5_native/native/A_p0_h2.json`.
When successful, the three-flavor runs have suffix `_h3` rather than `_h1`.

The existing numerical flavor proxy uses three heavy generations independently
of the optional two-flavor Matchete smoke test; it is not a physical C5
comparison or a coupling-normalization certificate.

## Actual test scope

The complete files passed Python compilation, guarded source checks, WL
bracket check, both unsupported/accepted argument checks, mocked six-model
native invocation (2 and 3), and the prior load-order regression. **No
live Matchete one-loop calculation is available in this container.** The local
pilot may reveal a later independent matching problem.
