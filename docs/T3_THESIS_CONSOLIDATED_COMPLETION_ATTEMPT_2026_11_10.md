# T3 A(-4)/D(-2) consolidated completion — honest execution contract

**Repository inspected read-only:** `Frankyang2002/T3neutrinomassmodelpipeline`, GitHub `main` commit `e0a667bc0023c76c45f17332c871aceb5ee773cf` on 10 November 2026. The user's local directory is `C:\thesistemp\T3neutrinomassmodelpipeline`. Nothing has been changed on GitHub.

## Actual physics result and outstanding work

This package **does not finish full-general two-threshold matching, numerical Yukawa optimisation, or physical neutrino predictions**. The mathematical gap is *substantial*, not a missed program flag. The earlier implementation plan was too optimistic, and splitting Laurent analysis into an extra named batch was unnecessary. This consolidated runner therefore refuses to label a restricted calculation complete.

- **Executed on user's Windows machine** (per their reported logs; results inspected on GitHub): full native general UV/EFT1 *renormalisable* one-loop RG for both A(-4) and D(-2); complete raw all-three-fermion Matchete `F`-only matching with both scalar multiplets retained; exact loop-order and explicit Laurent split. The numerical RG trajectories and symbolic WXF files are retained.
- **Implemented and offline-tested in this package:** manifest/hash and source-inventory validation, D(-2) generated five extra quartics and scalar self-bilinear at the S endpoint, numerical *raw unrenormalised* LLHH source channel from the exact saved Matchete operator ledger (A: six terms; D: twelve), a source-contract input interface for genuinely finite `C5(MS)`, and actual SMEFT solver/Takagi/PMNS/NuFIT adapters that execute **only** after a valid source-derived boundary exists.
- **Not implemented or executed** for these B1/B2 values: physical Matchete-to-canonical Wilson basis and signed spinor normalization; complete UV-hard minus EFT1-soft subtraction with local MSbar counterterms and pole cancellation; numerically evolved **complete** EFT1 Wilson tensor; *joint* full one-loop integration of both scalars with all general quartics and D(-2) self-bilinear; 18-real-Yukawa **full-general** refit. The earlier D(-2) one-loop scalar matching prototype supports a provisional one-heavy-flavour contraction and must not be substituted for three-heavy-flavour results.

**Formal direct hard contribution for thesis input**, computed from source `one_loop_operator_ledger.json` and the true Batch 1 `MF` endpoint by interpreting the *native ordered* LLHH coefficients and neglecting the RG-generated MF off-diagonals only when their relative norm is below `1e-8`:

| Model | 3-mass diagonal [GeV] | Raw ordered LLHH regular Frobenius [GeV^-1] | Raw ordered LLHH pole Frobenius [GeV^-1] |
|---|---|---:|---:|
| A(-4) | 110354.78, 121390.26, 132425.74 | ~1.45e-16 | ~2.45e-16 |
| D(-2) | 105426.09, 115968.70, 126511.31 | ~1.93e-16 | ~2.55e-16 |

These numbers are **not renormalised Wilson coefficients**, **not scalar-threshold Weinberg coefficients**, and **not neutrino masses**. They omit EFT1 soft/counterterm cancellation and must not be inserted into SMEFT as `C5`. D(-2) permits `y2Conj`; at this original UV-zero benchmark it remains zero, so its direct extra-LLHH source channel is zero. The D(-2) five generated quartics and nonzero self-bilinear are nevertheless retained in the scalar endpoint evidence and require full scalar matching.

## Install in local project

Every ZIP member is a **complete file** and the original relative project path is retained:

- `Numerical/parallel_t3/complete/__init__.py`
- `Numerical/parallel_t3/complete/full_general_completion.py`
- `Numerical/parallel_t3/complete/hard_llhh_source.py`
- `scripts/RunT3ThesisFullGeneralConsolidated.py`
- `tests/test_t3_thesis_full_general_complete.py`
- `docs/T3_THESIS_CONSOLIDATED_COMPLETION_ATTEMPT_2026_11_10.md`

Prerequisites: previous completed Batch1 and Batch2 raw/Laurent v2 outputs in their existing project paths. No rerunning symbolic generation or Matchete is needed. Do not overwrite prior outputs.

```powershell
$repo = 'C:\thesistemp\T3neutrinomassmodelpipeline'
Set-Location $repo
$zip = 'C:\path\to\T3_FullGeneral_Consolidated_2026_11_10.zip'
Expand-Archive -LiteralPath $zip -DestinationPath $repo -Force
python -X utf8 -m unittest discover -s tests -p test_t3_thesis_full_general_complete.py -v
python -X utf8 -m scripts.RunT3ThesisFullGeneralConsolidated --all --run-id thesis_full_general_consolidated_v1
```

The final invocation writes `output/parallel_t3/full_general_consolidated/thesis_full_general_consolidated_v1/{A_-4,D_-2}/formal_F_LLHH_not_a_Wilson.json` and `completion_report.json`, plus the common `summary.json`. It returns exit code **2** if finite matching is missing. This is an **intentional, scientifically necessary** failure, not a Matchete crash. No neutrino masses or fit files are created in this condition.

## Full-general matching completion handoff

If in future the genuinely complete, source-derived and numerically evaluated MSbar Weinberg coefficient becomes available, place **independent** records at:

`output/parallel_t3/full_general_complete_input/A_-4/renormalized_c5_at_MS.json`

`output/parallel_t3/full_general_complete_input/D_-2/renormalized_c5_at_MS.json`

The `C5_FORMAT` in `Numerical/parallel_t3/complete/full_general_completion.py` defines the strict input contract, including the F/Scalar integration order, source hashes, full-EFT1-Wilson RG provenance, and a scalar matching WXF source hash. A true matched `3x3 C5_MS_GeV_inverse` plus the existing numerical SMEFT state triggers the actual `SMWeinbergEvolution.evolve_sm_weinberg` solver from `1e3` to `100` GeV and the 5-observable NuFIT Gaussian chi-square. The 18-variable *refit* still needs a per-trial full general RG+matching evaluator; a fixed benchmark Wilson matrix cannot be used to fit varying Yukawas. The script reports this without inventing a result.

**No additional named sub-batch is recommended** merely to regenerate already working symbolic audits. The only useful completion work is implementing the actual matching and corresponding parameter-to-observable evaluator.
