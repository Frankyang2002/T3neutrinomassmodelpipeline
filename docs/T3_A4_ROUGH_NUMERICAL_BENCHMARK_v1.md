# A(-4) rough numerical neutrino benchmark (not a completed matching calculation)

This patch deliberately prioritises numerical outputs over independent source certification.

## Physics workflow

UV -> integrate all three heavy fermions together at 100,000 GeV -> EFT1 with
both scalar multiplets -> integrate both scalar multiplets together at
1,000 GeV -> SM + Weinberg -> run to 100 GeV -> Takagi / PMNS / NuFIT.

The tentative Weinberg coefficient is constructed from the ALREADY COMPUTED
finite ordered direct hard/soft channel:

  A_pq = (2 lambdaT3)/(sqrt(3) 16 pi^2) sum_r (y1_pr* y2_qr*/M_r)
         [mS1^2 log(mS1^2/M_r^2) - mS2^2 log(mS2^2/M_r^2)]/(mS1^2-mS2^2)

  C5_approx,pq = f_GammaCC (A_pq+A_qp)

The default **assumption** is `f_GammaCC=1`. This value is not source-derived.
All three flavours and the Batch 1 running scalar masses are used. The direct
channel retains F Yukawas/lambda at F. The distinct formal EFT1 Wilson RG
term is **not** added to avoid unverified double counting of the same
F/S logarithmic channel. The EFT1 RG remains a separate useful diagnostic.

The SM + Weinberg solver uses `atol=1e-24`, because K is around 1e-15 GeV^-1 and the repository default 1e-11 atol does not control its relative integration error.

This is a rough estimate, not the full renormalised one-loop result. In
particular, other F-hard dimension-five operators, scalar finite diagrams,
MSbar counterterms, running of all operators and the independent operator
normalisation/dictionary are not included. All result files flag this.

## Files (full files, not patches)

- `Numerical/parallel_t3/complete/a4_rough_benchmark.py` — new
- `scripts/RunT3A4ApproxBenchmark.py` — new
- `tests/test_t3_a4_rough_benchmark.py` — new
- `docs/T3_A4_ROUGH_NUMERICAL_BENCHMARK_v1.md` — new

## Execute locally in PowerShell

```
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "$HOME\Downloads\T3_A4_Rough_Numerical_Benchmark_v1.zip" -DestinationPath . -Force
python -X utf8 -m unittest discover -s tests -p test_t3_a4_rough_benchmark.py -v
python -X utf8 -m scripts.RunT3A4ApproxBenchmark
Get-Content output\parallel_t3\a4_rough_neutrino_benchmark_v1\A_-4\benchmark.json
```

Prerequisite: previously installed `a4_direct_hard_soft.py` from Patch 2A
and current `SMWeinbergEvolution.py`, `rough_numerics.py`, the Batch 1
numerical endpoints, and NuFIT target already in the repository.

Optional sensitivity check (separate output):

```
python -X utf8 -m scripts.RunT3A4ApproxBenchmark --operator-factor -1 --out output\parallel_t3\a4_rough_neutrino_benchmark_v1_minus\A_-4
```

A factor -1 changes the overall neutrino Majorana phase but leaves masses
and most observables unchanged; for a scale sensitivity use factor 2.

Outputs: `benchmark.json` and `numerical_matrices.npz`.
No outputs are written to the strict full-physical-C5 input folder.
