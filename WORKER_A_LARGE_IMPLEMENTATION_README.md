# Worker A: combined direct-channel calculation at the fermion threshold

## Scope

At mu_F = 100,000 GeV, all 3 heavy fermions are integrated simultaneously while **both S1 and S2 are retained in EFT1**. Neither scalar is integrated at mu_F. This package provides **source-verified direct LLHH and candidate EFT1 scalar-loop diagnostics**, **not** the full renormalised one-loop Wilson boundary.

## Files (relative to project root)

- `scripts/WorkerAFermionDirectThresholdBoth.py` — new full combined A(-4) and D(-2) calculation, JSON output with source hashes, ordered 3x3 complex arrays, pole/finite separation and explicit blockers.
- `scripts/TestWorkerAFermionDirectThresholdBoth.py` — new numerical, symbolic, and negative regression tests.
- `scripts/WorkerAFermionEFT1SoftA4.py` — complete existing A(-4) source-checked implementation included for reproducible imports (same implementation as earlier package).
- `scripts/TestWorkerAFermionEFT1SoftA4.py` — existing A(-4) test suite included for reproducible imports.
- `scripts/WorkerASoftMixingKernel.py` — complete existing scalar-soft symbolic source parser included for reproducible imports.
- `output/parallel_t3/worker_a_direct_threshold_v1/benchmarks.json` — tested A(-4)/D(-2) numerical output (diagnostic only).

## Execution (PowerShell)

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
$zip = Get-ChildItem "$env:USERPROFILE\Downloads\worker_a_fermion_direct_both_large*.zip" | Sort-Object LastWriteTime -Descending | Select-Object -First 1
if (-not $zip) { throw 'ZIP not found in Downloads' }
Expand-Archive -LiteralPath $zip.FullName -DestinationPath . -Force
python -m py_compile scripts\WorkerAFermionDirectThresholdBoth.py scripts\TestWorkerAFermionDirectThresholdBoth.py scripts\WorkerAFermionEFT1SoftA4.py scripts\WorkerASoftMixingKernel.py
$runId = Get-Date -Format 'yyyyMMdd_HHmmss'
python scripts\WorkerAFermionDirectThresholdBoth.py --sources-zip 'output\audits\worker_a_sources_20261112_003127.zip' --soft-zip 'output\audits\worker_a_scalar_soft_20261112_004007.zip' --output "output\audits\worker_a_F_direct_both_$runId.json"
python scripts\TestWorkerAFermionDirectThresholdBoth.py --sources-zip 'output\audits\worker_a_sources_20261112_003127.zip' --soft-zip 'output\audits\worker_a_scalar_soft_20261112_004007.zip'
```

## Source-backed derivations

- A(-4): six hard LLHH terms; repeats existing A(-4) benchmark to protect established results.
- D(-2): twelve direct LLHH terms, **six survive** restriction of the independent `y2Conj` coupling. Complex conjugation of the surviving `y1`, `y2` matrices is preserved.
- Both: direct UV hard Laurent pole and finite logarithmic terms, plus the formal EFT1-scalar-bubble topology evaluated at mu_F from original Matchete scalar-soft expressions; all three heavy masses retained.
- D(-2): scalar split masses taken as `mS2Sq +/- 2*abs(muT3GeneralS2Sq)` (candidate Matchete-RGBeta dictionary), with source-derived diagonal and offdiagonal operator-space mixing. Mixed quartic survives with `WExtra=0`.
- Both: within the ordered direct channel, opposite pole residues and cancellation of explicit d/d(log mu^2) of the paired diagnostic. This is **not** an all-operator UV/EFT1 matching proof.

## Outstanding physics for full Worker A completion

1. Calculate/identify and project the **full EFT1-soft** one-loop amplitude at mu_F in the same operator basis as the UV one-loop calculation. The scalar-soft source reused here was emitted for the later joint-scalar threshold, so topology transfer is **not certified**.
2. Include the renormalisation-counterterm/evanescent-operator prescription and prove UV-EFT1 finite matching with pole cancellation for the **full restricted operator inventory** (A: 135 total raw operators; D: 359 retained after Yukawa restriction), not just direct LLHH. Verify whether non-LLHH sectors feed LLSS mixing.
3. Certify `GammaCC` -> two-component Weyl coefficient conventions, SU(2) CG contractions and scalar-mass/coupling dictionary; then construct and test the finite renormalised EFT1 Wilson arrays. No such arrays are exported by this package.

The project sequence remains `UV -> integrate out ALL fermions -> EFT1 with BOTH scalars -> integrate out BOTH scalars -> SMEFT`. No production solver, RGE, or physical C5 files are modified.
