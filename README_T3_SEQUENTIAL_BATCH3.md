# T3 sequential matching — Batch 3 (ordinary three-flavour hard scalar template)

## What is implemented

For **D(+2)** and **A(−4)** only, the existing Wolfram Match stage from Batch 2 has already produced an independent scalar-only one-loop hard C5 *per heavy eigenvalue* template. This Batch 3 code consumes those saved expressions and the corresponding `output/general_observable/all16_v1/<tag>/running.json` values, and performs:

1. Input/provenance checks (`HardScalarCandidateExported`, `CheckLagrangianPassed`, the exact six-term Matchete expression fingerprint). An optional stricter flag verifies the original first-threshold source file hashes as well.
2. The explicit sum over the three **nondegenerate**, positive, diagonal heavy-fermion masses at muF, with the two three-by-three flavour Yukawa matrices at muF.
3. Complex-conjugate Yukawa convention `Bar[y1] Bar[y2]` with `lambdaT3` unbarred. The symmetric flavour projection is explicitly documented; its off-diagonal physical normalization remains provisional.
4. The exact finite scalar one-loop kernel at muS, with numerically stable equal-mass limit.
5. The direct leading EFT1 LLSS log and its fixed-order matching-scale cancellation diagnostic. The raw 1/epsilon residue is retained separately.
6. An **optional, independent numerical cross-check** against the repository's archived `lo_joint_scalar_matching/exploratory_all16_v6` matrices. The archived LLSS term is **never added** to the new result, because that would double-count it.
7. Individual per-heavy-eigenstate source matrices, total C5 matrices, conditional Takagi masses, CSV and JSON summary. No Wolfram is run.

## Formula derived from the uploaded Matchete template

Let

```
P_ij = -2/sqrt(3) * (lambdaT3/(16*pi^2))
       * 1/2 * sum_{r=1}^3 [conj(y1_{i r})conj(y2_{j r})
                           +conj(y1_{j r})conj(y2_{i r})] / MF_r .
```

For `a=mS1^2`, `b=mS2^2`, define

```
F(mu) = 1 + [a log(mu^2/a) - b log(mu^2/b)]/(a-b),
F(mu; a=b) = log(mu^2/a).
```

Then the exact exported source template, *with its explicit UV pole*, becomes

```
C5_hard,unsubtracted(muS) = P / epsilon + P * F(muS),
C5_EFT1,direct-log(muS)   = -2 P * log(muS/muF),
C5_LLSS,finite-subset     = P * [F(muS) - 2 log(muS/muF)].
```

`hbar` has been set to `1/(16*pi^2)`, exactly as in the existing D(-2) numerical benchmark. The direct-log coefficient and sign agree with the repository's saved D(+2) and A(-4) EFT1 running expressions. At fixed scalar masses/couplings, `d/dlog(muS) [hard_finite + direct_log] = 0`.

**Pole warning:** The raw hard coefficient contains `P/epsilon`. The script records the residue, but does **not** claim an independently proved MSbar counterterm or renormalized physical observable. A formal subtraction of `-P/epsilon` is not the same as showing a source-derived complete counterterm cancellation.

**Completeness warning:** The inherited O(hbar) Wilson sector from the first threshold and other full EFT1 contributions are not added. The output is the **LLSS hard-finite + direct-log subtotal**, which is **not complete fixed-order SMEFT C5** and is not an independently normalized physical neutrino mass prediction. The archived provisional all-16 calculation contains this same LLSS subtotal; do not add the two results together.

## Project file locations

All ZIP entries install directly in the repository root:

- `Numerical/running/T3OrdinaryHardFlavor.py` — source-locked numerical reduction and checks.
- `scripts/RunT3OrdinaryHardFlavorBatch.py` — command-line runner, resume, outputs.
- `tests/python/test_t3_ordinary_hard_flavor.py` — offline unit tests.
- `tests/fixtures/t3_ordinary_hard_c5_candidate.txt` — exact source expression test fixture (from the successful local hard matching).
- `README_T3_SEQUENTIAL_BATCH3.md` — this file.

## PowerShell

From `C:\thesistemp\T3neutrinomassmodelpipeline`:

```powershell
Expand-Archive -LiteralPath "$HOME\Downloads\T3_Sequential_Batch3.zip" -DestinationPath . -Force
python -X utf8 -m pytest -q tests/python/test_t3_ordinary_hard_flavor.py
python -X utf8 -m scripts.RunT3OrdinaryHardFlavorBatch --check-only
python -X utf8 -m scripts.RunT3OrdinaryHardFlavorBatch --run-id hard_flavor_01 --require-reference-parity
```

No previous hard matching is rerun. If archived reference parity files are absent, omit `--require-reference-parity` to calculate without reference comparisons. `--verify-fermion-sources` additionally checks source hashes from the saved provenance when the F-only release is present locally.

For an interrupted run, resume with the *same* options:

```powershell
python -X utf8 -m scripts.RunT3OrdinaryHardFlavorBatch --run-id hard_flavor_01 --require-reference-parity --resume
```

Output:

- `output/general_sequential_pilot/ordinary_hard_flavor/hard_flavor_01/D_p2/result.json`
- `output/general_sequential_pilot/ordinary_hard_flavor/hard_flavor_01/A_m4/result.json`
- `output/general_sequential_pilot/ordinary_hard_flavor/hard_flavor_01/summary.json`
- `output/general_sequential_pilot/ordinary_hard_flavor/hard_flavor_01/comparison.csv`

Input paths are configurable with `--hard-root`, `--running-root`, `--reference-root`, and `--output-root`.

## Validation performed here

- 14/14 offline tests passed.
- The complete runner was executed against the **two uploaded successful Matchete expressions** and their **two uploaded RG snapshots**; it performed all three eigenstate sums, direct logs, and fixed-order scale checks without rerunning Wolfram.
- Against separately retrieved existing archived LLSS component numbers for the same two benchmark points, all three matrix components (hard finite, EFT1 logarithm and sum) agreed to better than `5e-14` relative Frobenius precision.
- The runner's `--resume` path was tested and reused both results.

Expected conditional (not certified physical) masses (ascending, eV):

| Model | m1 | m2 | m3 |
|---|---:|---:|---:|
| D(+2) | 0.0196944948422731 | 0.0218828367230930 | 0.0262594043481166 |
| A(−4) | 0.0194265340013435 | 0.0215851015754182 | 0.0259021221670917 |

## Next physics work

To obtain a genuinely complete fixed-order C5, the **first-threshold one-loop inherited Wilson sector**, full EFT1 Wilson running and **independent counterterm/normalization audit** must be obtained and matched in the same conventions. Only then use the resulting total C5 for physical neutrino-mass claims. The mixed scalar B(-1) model is not covered by this ordinary mass-basis template.
