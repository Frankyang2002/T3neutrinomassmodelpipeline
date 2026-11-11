# Worker A: source-exact scalar-soft numerical operator kernel

## New code path

- `scripts/WorkerASoftMixingKernel.py` (full standalone new Python source).

This script reads the original native scalar-soft expressions in the user's
`worker_a_scalar_soft_20261112_004007.zip` and the scalar endpoints in
`worker_a_sources_20261112_003127.zip`. It recalculates A(-4) and D(-2)
Laurent pole and finite structures symbolically, verifies source reconstruction,
checks diagonal/off-diagonal coefficient symmetry, checks the logarithmic scale
derivatives, and evaluates the **source-basis scalar-soft** coefficient matrices.

## Benchmarks

At the scalar matching endpoint muS = 1e3 GeV, with mu_bar = muS, in
Matchete source's hbar=1 convention:

| Model | f_diag | f_mix | coefficient of conjugated W12 at WExtra=0 |
|---|---:|---:|---:|
| A_-4 | 0.469517224388 | 0 | 0.0484097955966 |
| D_-2 | 0.386651177977 | -0.00890256316054 | 0.0368955355874 |

For D(-2), the source contains
`f_diag * (lambdaT3 * W12bar + lambda_extra_bar * WExtrabar)
 + f_mix * (lambda_extra_bar * W12bar + lambdaT3 * WExtrabar)`.
The evaluated coefficient of W12bar includes a nonzero mixed-extra-quartic term,
even after setting WExtra=0.

### Strict limits

- This is a **joint scalar-threshold** source calculation at muS, not the
  EFT1-soft subtraction required for fermion matching at muF=1e5 GeV.
- Scalar mass eigenvalue and quartic mappings to native Matchete sources are
  still provisional. The code reports both as uncertified.
- Matchete hbar bookkeeping to canonical MSbar `1/(16*pi*pi)` normalization
  and GammaCC-to-Weyl normalization are not certified.
- No finite physical C5, no renormalised fermion matching and **no production
  boundary** produced. The existing solver and all history remain unchanged.
- The source ZIP archives are required inputs; no need to put them in Git.

## PowerShell

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
$zip = Get-ChildItem "$env:USERPROFILE\Downloads\worker_a_soft_mixing_kernel*.zip" | Sort-Object LastWriteTime -Descending | Select-Object -First 1
if (-not $zip) { throw 'Kernel ZIP missing in Downloads' }
Expand-Archive -LiteralPath $zip.FullName -DestinationPath . -Force
python -m py_compile scripts\WorkerASoftMixingKernel.py
$runId = Get-Date -Format 'yyyyMMdd_HHmmss'
python scripts\WorkerASoftMixingKernel.py `
    --soft-zip 'output\audits\worker_a_scalar_soft_20261112_004007.zip' `
    --sources-zip 'output\audits\worker_a_sources_20261112_003127.zip' `
    --output "output\audits\worker_a_mixing_kernel_$runId.json"
```

The `--sources-zip` archive might instead reside in Downloads depending on
where you saved the exporter result. Adjust only this path if necessary.
