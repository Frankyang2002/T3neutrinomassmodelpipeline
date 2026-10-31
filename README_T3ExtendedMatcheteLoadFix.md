# Extended Matchete C5 native load-order correction

## What failed

Six A/D/E native attempts stopped at `load-Matchete` with `Bar::shdw`,
`Coupling::shdw`, and `Loading failed!`. The numerical sensitivity proxies and
16-model RGBeta UV integration were not implicated.

## Root cause and fix

Wolfram Language resolves symbol contexts while reading an entire expression.
The previous `RunT3ExtendedC5Native.wl` loaded `Matchete`` from inside the same
large `Module` that contained `Bar[Coupling[...]]`. Those symbol names were
read before `Needs["Matchete`"]` executed. This caused premature `Global``
symbol creation and shadowing warnings.

The full replacement loads Matchete in its own **completed top-level
expression** before the large matching `Module` is read. If loading fails,
it saves a blocked native JSON result and exits nonzero. All matching and C5
extraction operations remain unchanged.

## Files and project paths

- `Lagrangian/RunT3ExtendedC5Native.wl`: **complete replacement** of the
  previously supplied native matching driver.
- `scripts/TestT3ExtendedMatcheteLoadBoundary.py`: new source-level regression
  test; checks bracket structure, ordering, negative control, and unchanged
  matching/extraction interface.

## PowerShell

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
$zip = "$HOME\Downloads\T3ExtendedMatcheteLoadFix_Local.zip"
Expand-Archive -LiteralPath $zip -DestinationPath . -Force

# Source-level test only (does not execute Matchete):
python scripts\TestT3ExtendedMatcheteLoadBoundary.py
python scripts\RunT3ExtendedMatchingBatch.py --self-test

# Real live Matchete load + A(0) pilot:
python scripts\RunT3ExtendedMatchingBatch.py --pilot

# Only if pilot returns 0:
if ($LASTEXITCODE -eq 0) {
    python scripts\RunT3ExtendedMatchingBatch.py --all
}
$LASTEXITCODE
```

This fix has **not** been executed against the user's Windows Matchete
installation. The first pilot might reveal a subsequent matching-stage failure.
The one-step Matchete result, even if successful, is **not** a certified
hierarchical threshold calculation or Matchete↔RGBeta coupling dictionary.
