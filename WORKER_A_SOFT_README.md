# Worker A scalar-soft Laurent audit

Source: local uploaded `worker_a_scalar_soft_20261112_004007.zip`.

## Installed file

- `scripts/WorkerAScalarSoftLaurent.py`: read-only source parser and symbolic Laurent verifier. Requires `sympy`.

`WORKER_A_SOFT_RESULTS.json` is a run artifact, not installed source.

## Findings

Let `W = conj(T3EFT1W12)`, `WE = conj(T3EFT1WExtra)`, `L = lambdaT3`, `LE = conj(lambdaT3ExtraQ020101I1)`, `A = MS1`, `B = T32MinusMass`, `C = T32PlusMass`, `Q = mubar^2`. Laurent coefficients below set Matchete's formal `hbar` to one. Physical loop normalization not independently certified.

- A(-4) soft pole is `-2/sqrt(3) * L*W`, and its finite function is the source's two-mass scalar bubble.
- D(-2) soft pole is `-2/sqrt(3) * (L*W + LE*WE)`.
- D(-2) soft finite expression factorizes as `f_diag*(L*W+LE*WE) + f_mix*(LE*W+L*WE)`, with **nonzero f_mix**. Hence setting `WE=0` still leaves a `LE*W` contribution when the extra quartic is nonzero.
- The D(-2) continuation and mass-repair source files are byte-identical.

**Important:** These are scalar-soft Matchete expressions from a joint scalar threshold calculation at the later scale; they do not by themselves establish the fermion-threshold UV-hard minus EFT1-soft subtraction. Do not combine pole norms computed in different normalizations/scales, or use this soft finite term as a boundary at `mu_F` without an independently checked EFT1 operator mapping. No counterterm or full operator-basis cancellation is certified. No production Wilson arrays are generated.

## Commands (PowerShell)

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
$zip = Get-ChildItem "$env:USERPROFILE\Downloads\worker_a_scalar_soft_verified*.zip" | Sort-Object LastWriteTime -Descending | Select-Object -First 1
if (-not $zip) { throw "Verification ZIP not found" }
Expand-Archive -LiteralPath $zip.FullName -DestinationPath . -Force
$source = 'output\audits\worker_a_scalar_soft_20261112_004007.zip'
$runId = Get-Date -Format 'yyyyMMdd_HHmmss'
python scripts\WorkerAScalarSoftLaurent.py --source-zip $source --output "output\audits\worker_a_soft_laurent_$runId.json"
```

If source ZIP has a different name/location, change `$source` to that path.
