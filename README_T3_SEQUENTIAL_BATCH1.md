# T3 first-batch numerical handoff (D(+2), A(-4), B(-1))

## Purpose and scope

This is a **small, source-backed implementation batch**, extracted after
inspecting the supplied **local D(-2) consolidated runner** and the current
GitHub `main` source (read-only). Existing files are not replaced. No GitHub
change is made. It deliberately **does not claim complete D(+2) sequential
one-loop matching**: the difficult model-dependent hard one-loop threshold
has not yet been ported from the D(-2) Wolfram runner.

The numerical pipeline reused here is:

`UV -> integrate ALL F -> EFT1(S1,S2) -> direct LLSS->C5 one-loop log
-> simultaneously integrate ALL scalars (finite bubble) -> conditional SMEFT C5`

The **dynamic** backend uses the repository's existing
`scripts/RunT3GeneralSequentialLeadingLog.py:process` and native source
boundary, source-resolved HHSS tensor, general quartic registry, and RG
snapshots. It recomputes conditional LLSS contributions at the saved
reference state; no Wolfram matching or fitting is repeated.

The **archived-v6** backend independently reconstructs/checks the saved
conditional log, bubble, and Takagi magnitudes; it does not recalculate
matching diagrams. This is useful for a quick preflight and comparison.

Both backends explicitly preserve:

- the unknown signed Matchete `GammaCC -> Weyl` normalization (`--spinor-factor
  1` is a convention, **not certification**);
- conditional source-to-RGBeta quartic aliases and special scalar bases;
- no independently computed hard one-loop fermion/scalar contributions other
  than the LLSS finite scalar bubble;
- no full dimension-five operator-completeness certification;
- no claim that the generated masses are a neutrino-data fit;
- the complete independent D(-2) source and result files, unchanged.

The numerical consistency checks include 3x3 complex flavour symmetry,
`C5_total = C5_log + C5_bubble`, finite nonnegative Takagi masses, SU(2)
projection residual, and fixed-order matching-scale cancellation. These are
necessary consistency checks, **not** the D(-2) exact `1/epsilon` pole proof.

## New files and intended paths

- `Numerical/running/T3SequentialBatchNumerics.py`
- `scripts/RunT3SequentialNumericalBatch.py`
- `tests/python/test_t3_sequential_batch_numerics.py`
- `README_T3_SEQUENTIAL_BATCH1.md`

## Install locally (PowerShell)

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "$HOME\Downloads\T3_Sequential_Batch1.zip" -DestinationPath . -Force
python -X utf8 -m pytest -q tests/python/test_t3_sequential_batch_numerics.py
```

## Preflight (no writes / no matching)

```powershell
python -X utf8 -m scripts.RunT3SequentialNumericalBatch --check-only
```

The default models are D(+2), A(-4), B(-1). This reports source paths and
whether the existing all-16 boundary, general quartic seeds, RG snapshots,
and registry are available. The read-only D(-2) manifest check looks at
`output/general_sequential_pilot/D_m2_stage1/data/general_d_m2_scalar_threshold/sequential_run/sequential_matching_manifest.json`.

## Fresh conditional LLSS C5 calculation — first three

```powershell
python -X utf8 -m scripts.RunT3SequentialNumericalBatch --run-id first3_v1
```

This uses the dynamic source-backed calculator. It does **not** run or overwrite
expensive Matchete/Wolfram matching. If interrupted or blocked by one model:

```powershell
python -X utf8 -m scripts.RunT3SequentialNumericalBatch --run-id first3_v1 --resume
```

To include the existing *conditional* SMEFT low-scale evolution in a fresh run:

```powershell
python -X utf8 -m scripts.RunT3SequentialNumericalBatch --run-id first3_smeft_v1 --with-smeft
```

## Quick archived baseline (no matching)

```powershell
python -X utf8 -m scripts.RunT3SequentialNumericalBatch --backend archived-v6 --run-id first3_archived_v1
```

The repo's archived v6 numerical results (not new results from this package)
include these **eV per |N_GammaCC|, conditional only**, ascending:

| Model | m1 | m2 | m3 |
|--|--:|--:|--:|
| D(+2) | 0.01969449 | 0.02188284 | 0.02625940 |
| A(-4) | 0.01942653 | 0.02158510 | 0.02590212 |
| B(-1) | 0.01776825 | 0.01974256 | 0.02369107 |

These are *not* independently hard-matched sequential physical predictions.
The dynamic calculation can differ because it includes general HHSS terms
not present in the older v6 illustrative background.

## All 16 — conditional numerical baseline

```powershell
python -X utf8 -m scripts.RunT3SequentialNumericalBatch --all16 --run-id all16_conditional_v1
```

No fit, new Mathematica matching, or D(-2) rerun is triggered. Missing input
files or failed projections produce a `Blocked` entry, not a fabricated C5.

## Outputs

For run `first3_v1`, outputs are in:

```
output/general_sequential_pilot/all16_numerical_batch/first3_v1/
  settings.json
  summary.json
  comparison.csv
  D_p2/result.json
  A_m4/result.json
  B_m1/result.json
  <model>/error.json  (only on failure)
```

All `result.json` files contain full complex 3x3 LLSS-induced C5 matrices,
separate direct-log and finite-bubble pieces, scalar-mass class, Takagi
singular-value magnitudes and consistency residuals. Source hashes protect
`--resume` against inadvertent source changes.

## What is still required for full D(+2) matching

The D(-2) local runner contains the separate **physical-basis hard Matchete
one-loop scalar Match**, EOM reduction, inherited C5 sector, fixed-order
assembly, Warsaw reference and `1/epsilon` pole audit. It is tightly bound to
the D(-2) `NewScalar1,T3S2R,T3S2I` physical basis and hard-coded term counts;
renaming D(-2) to D(+2) would be mathematically incorrect. D(+2) has ordinary
complex scalar multiplets with a different SU(2) CG and threshold registry.

A complete D(+2) generalization requires a new physical EFT1 WXF continuation,
source-exact all-F boundary including full quartic inventory, matching of both
scalar multiplets to SMEFT and normalization of that hard result against the
LLSS direct running and inherited terms. This package leaves those terms
**uncomputed**, rather than silently identifying a one-scale C5 with an
independently derived hard sequential contribution.
