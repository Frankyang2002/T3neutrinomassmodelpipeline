# Worker C — restricted D(-2) numerical finite-part diagnostic (v2)

This replaces the first overlap runner, which incorrectly required the isolated
one-loop hard Matchete coefficient to be pole-free. A hard-only term generally
has a nonzero 1/epsilon pole; the saved source's **combined hard + inherited +
running** result has independently verified exact NF=1 pole cancellation.

## Contents and paths

- `scripts/RunT3WorkerCNumericOverlap.py` — full replacement runner.
- `tests/test_t3_worker_c_numeric_overlap.py` — full replacement test module.
- `WORKER_C_NUMERIC_OVERLAP_README.md` — this file.

## Behaviour

- Reuses saved restricted `hard` and `total` WXF candidates and the unmodified
  `Lagrangian/RunGeneralDMinus2C5RoughBenchmark.wl` evaluator. No new Match[].
- Creates a fresh isolated benchmark with `y2Conj=0` and
  `PoleHandling="FinitePartDiagnostic"` for both candidate evaluations.
- Requires both Wolfram outputs to declare `FinitePartOnly=True`,
  `RenormalizationCertified=False`, and valid 3x3 finite matrices.
- Records hard and total numerical pole residues separately. A nonzero hard
  pole is expected. The combined result must have numerical residue <=
  `max(1e-24, 1e-10 * hard_pole_residue)` GeV^-1, justified as a diagnostic
  floating-point check by the saved exact NF=1 total pole cancellation.
- The extracted hard finite part is NOT a renormalized standalone Wilson
  coefficient. It is NOT safe to add to the existing scalar bubble.
- Optional `--bubble-input` remains conditional: must match benchmark model,
  scales, and scalar eigenmasses. This alone does not certify the independent
  Matchete-to-real-tensor coupling normalization. Never add its difference
  to the production calculation.

## PowerShell

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "$HOME\Downloads\T3_worker_C_numeric_overlap_POLE_FIX_20261112.zip" -DestinationPath . -Force
python -m unittest discover -s tests -p "test_t3_worker_c_numeric_overlap.py" -v
$restricted = "output\audits\worker_c_restricted_D_m2_20261112_052721_aaa5ec32"
python -m scripts.RunT3WorkerCNumericOverlap --restricted-dir $restricted --preflight
python -m scripts.RunT3WorkerCNumericOverlap --restricted-dir $restricted
```

A new directory under `output/audits/worker_c_numeric_overlap_D_m2_*`
contains `result.json`, `hard_wolfram.log`, `total_wolfram.log`,
`hard/.../profile_c5_rough_benchmark_result.json` and corresponding total
result and CSVs. Previously blocked output directories are untouched.

## Verification

8 offline Python tests passed in the generation environment. The Wolfram
numeric evaluation has **not** been executed here and must be checked locally.
No complete one-loop scalar threshold or physical operator normalization
certification is claimed.

## NF=1 zero-Takagi-mode diagnostic (2026-11-12)

The saved hard-only D(-2) finite part numerically evaluates all nine C5
entries, but the original `RunGeneralDMinus2C5RoughBenchmark.wl` then exited
14 at `SingularValueList`, before exporting the CSVs. A rank-two symmetric
Weinberg matrix of the form `y1 y2^T + y2 y1^T` has one zero singular value;
the exact cause of Wolfram exit 14 has not been independently verified.

This runner now reads the original evaluator (SHA256 recorded) and writes a
separate `worker_c_benchmark_evaluator.wl` in its new audit directory. It
replaces **only** the single legacy mass-extraction line with the three-value
mathematically equivalent construction

```wolfram
massesEV = Sort[N[Sqrt[Clip[Re[Eigenvalues[
    ConjugateTranspose[mnuEV].mnuEV]], {0, Infinity}]], 12]];
```

The `M^dag M` eigenvalues are the squared Takagi singular values, including
any zero mode. The input C5, flavour treatment, finite-part extraction,
benchmark parameters, and all source Matchete files are untouched. The
runner rejects unexpected changes to the original evaluator source instead
of applying a blind replacement. The original `Lagrangian/` evaluator file
is **not** replaced. This is diagnostic mass postprocessing, not physics
renormalization or a certified nonbubble correction.
