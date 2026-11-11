# Worker A — A(-4) candidate EFT1 scalar bubble at the fermion scale

This is a **new calculation**, distinct from the earlier scalar-threshold
numerical kernel: it substitutes **mu_F = 100,000 GeV** and **the scalar
quadratic masses from the UV benchmark at mu_F**, then evaluates the scalar
bubble multiplying the tree `W12` operator. It also evaluates the six direct
`LLHH` UV Matchete Laurent terms with all three heavy flavours.

## Source and implementation

- `scripts/WorkerAFermionEFT1SoftA4.py`: complete standalone calculation;
  reads the two input ZIPs without altering them; verifies SHA256 for the UV
  raw expression against the full 135-term ledger; validates both conjugate
  hard Laurent triplets; parses and symbolically reconstructs the original
  A(-4) native scalar-soft formula; computes and exports complex 3x3 ordered
  candidate tensors and their norms.
- `scripts/TestWorkerAFermionEFT1SoftA4.py`: regression tests against the
  actual uploaded source bundles, including degenerate scalar limit and
  explicit fixed-coupling scale cancellation.
- `output/parallel_t3/worker_a_fermion_soft_candidate_v1/A_-4.json`:
  benchmark output calculated directly from the two uploaded source ZIPs.

## Interpretation

All amplitudes use the **ordered GammaCC convention** and an explicitly
assumed `hbar = 1/(16*pi^2)` for reporting GeV^-1 values. The soft expression
comes from Matchete's **later** joint-scalar threshold. The same topology and
tree-W12 insertion are evaluated at `mu_F`; this is useful evidence but is
**not proof** of the EFT1-soft subtraction in full UV→EFT1 matching.

The scalar mass-squared coordinate mapping is taken from the UV input
`mS1Sq`, `mS2Sq` without a certified scalar Hessian / Matchete mass dictionary.
The complex heavy-mass matrix is almost diagonal at the given benchmark;
this script rejects non-diagonal inputs rather than assuming an SVD rotation.

**Not done:** full EFT1 operator basis/loops, MSbar counterterms, matching
scheme, GammaCC-to-Weyl map, finite renormalised one-loop Wilson arrays.
The diagnostic `hard + soft` is **not** promoted to a production coefficient.

## Run (PowerShell, project root)

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
$zip = Get-ChildItem "$env:USERPROFILE\Downloads\worker_a_fermion_soft_at_F_A4*.zip" | Sort-Object LastWriteTime -Descending | Select-Object -First 1
if (-not $zip) { throw 'ZIP not found in Downloads' }
Expand-Archive -LiteralPath $zip.FullName -DestinationPath . -Force
python -m py_compile scripts\WorkerAFermionEFT1SoftA4.py scripts\TestWorkerAFermionEFT1SoftA4.py
$runId = Get-Date -Format 'yyyyMMdd_HHmmss'
python scripts\WorkerAFermionEFT1SoftA4.py --sources-zip 'output\audits\worker_a_sources_20261112_003127.zip' --soft-zip 'output\audits\worker_a_scalar_soft_20261112_004007.zip' --output "output\audits\worker_a_F_soft_A4_$runId.json"
python scripts\TestWorkerAFermionEFT1SoftA4.py --sources-zip 'output\audits\worker_a_sources_20261112_003127.zip' --soft-zip 'output\audits\worker_a_scalar_soft_20261112_004007.zip'
```

No production solver files or other workstreams are changed.
