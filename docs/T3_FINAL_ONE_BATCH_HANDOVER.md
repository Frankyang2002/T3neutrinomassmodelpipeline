# A(-4) and D(-2): single-package full-general completion attempt

Checked GitHub `main` at `e0a667bc0023c76c45f17332c871aceb5ee773cf` on 2026-11-10. No remote changes.

## Exact EFT sequence

UV at 1e7 GeV -> integrate **all three heavy fermion flavours together** at one F threshold (1e5 GeV) -> EFT1 with **both** S1 and S2 -> integrate **both scalar multiplets together** at one S threshold (1e3 GeV) -> SMEFT -> Weinberg coefficient at 100 GeV -> Takagi/PMNS -> fit 18 UV real Yukawa parameters.

## Source status

- The local Batch 1 general UV/EFT1 **renormalisable** RG is complete for the thesis initial benchmark.
- Batch 2's first F-threshold full general Matchete symbolic source and Laurent split are complete.
- The independent EFT1 *Wilson* anomalous-dimension construction, finite subtraction and mass/coupling mapping have not been physically assembled. The original raw Laurent output cannot simply be relabeled a finite Wilson coefficient.
- The new Wolfram script attempts to run the **existing native sequential Matchete** two-threshold matching, not analytic triangle/bubble formulas; for D(-2), it labels the result as unrotated holomorphic S2 basis. That does not establish an exact physical neutral-scalar mass eigenbasis when the running-generated bilinear is nonzero.
- The native C5 evaluator attempts three independent heavy-flavour sums, retains complex running Yukawas, and refuses an ambiguous flavour map, uncancelled regulator pole, or unmapped Matchete coupling. Its potential numeric result is an explicitly unverified candidate. It does not insert the missing EFT1 Wilson RG. This code cannot be executed in the current container because Matchete is not installed here.
- The optional `--fit` flag hooks the 18-real Yukawa optimizer to the native candidate: **each trial reruns full general UV/EFT1 renormalisable RG and re-evaluates the native C5 expression**, with no triangle/bubble substitution. The optimizer executes only if the native C5 has passed the explicit numeric/pole check. Such fits remain **provisional**, because the complete EFT1 Wilson RGE, independent finite subtraction and D(-2) mass-basis map are not yet established.

## Project files

- `Lagrangian/RunT3FullGeneralTwoThresholdNative.wl` — full native sequential Matchete two-stage candidate; reruns symbolic first-stage F matching as part of one generic sequential call
- `Lagrangian/RunT3FullGeneralNativeC5Numerics.wl` — exact-index guarded numerical candidate conversion
- `scripts/RunT3ThesisFinalOneBatch.py` — one CLI; native numerical continuation, SMEFT solver and optional provisional fit
- `Numerical/parallel_t3/complete/final_one_batch_fit.py` — 18-real-yukawa objective/least-squares module
- `tests/test_t3_thesis_final_one_batch.py` — offline source and numerical interface tests
- `docs/T3_FINAL_ONE_BATCH_HANDOVER.md` — this file

## Local PowerShell

```
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath 'C:\path\to\T3_Final_One_Batch_2026_11_10.zip' -DestinationPath (Get-Location) -Force
python -X utf8 -m unittest discover -s tests -p test_t3_thesis_final_one_batch.py -v
python -X utf8 -m scripts.RunT3ThesisFinalOneBatch --all --run-id thesis_final_one_batch_v1 --fit --fit-maxeval 40
```

Output: `output/parallel_t3/full_general_final_one_batch/thesis_final_one_batch_v1/<MODEL>/`.

**The exit code remains 2 even if a provisional candidate and fit are produced; the complete finite Wilson matching, EFT1 Wilson RG and model-specific scalar mass basis are not established. This is NOT a finished numerical pipeline.** Do not treat any formal native Matchete output or trial fit as completed model predictions.
