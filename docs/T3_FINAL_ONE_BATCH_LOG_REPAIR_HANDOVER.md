# Existing final-one-batch pipeline: verified log repair

Source: uploaded November 10 run logs. This is NOT a new physics batch.

## Cause, not speculation

1. A(-4) F matching succeeds (three heavy fermion flavours together). The *joint* S1/S2 match fails with `FindDummyIndices::trippleindex` on `y1[p,r] y2[q,r] / MF[r]`. The raw flavoured mass inverse cannot be passed uncontracted through Matchete's next threshold. The existing repository has `RunT3EFT1WilsonStageLauncher.wl` + `RunThresholdStageEFT1Wilson.wl` for independent `T3EFT1W12[p,q]` coefficients. The replacement native script uses this exactly at TREE LEVEL to attempt the scalar-soft part. It saves the full F hard O(hbar) source separately. It DOES NOT pretend the missing one-loop hard Wilson reparameterisation or EFT1 Wilson RGE has been performed.
2. D(-2) F generic matching fails with an unrotated holomorphic self-bilinear. The previous successful Batch2 first-threshold Matchete run **did** use `T3Scalar2Minus/Plus` and split masses, with three heavy generations. The replacement checks that source and refuses the generic unrotated match; the physical scalar-mass-basis second threshold is still unimplemented in this path. Do not treat an unrotated run as physics.
3. The old Python wrapper incorrectly reported `Sequential Matchete returned failure: Success` when stage 1 had status Success but stage 2 failed, because the sequential runner's status carried stage 1 in the report. The repair instead identifies the exact stage and preserves diagnostics.

## Files in ZIP (all complete replacements)

- `Lagrangian/RunT3FullGeneralTwoThresholdNative.wl`
- `scripts/RunT3ThesisFinalOneBatch.py`
- `tests/test_t3_thesis_final_one_batch.py`
- `docs/T3_FINAL_ONE_BATCH_LOG_REPAIR_HANDOVER.md` (new)

All written paths are relative to `C:\thesistemp\T3neutrinomassmodelpipeline`. The prior `thesis_final_one_batch_v1` output is reused READ-ONLY for A F WXF; no UV/F RGBeta rerun and no GitHub writes. For D, reuse the saved physical Batch2 F source.

## Commands

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath 'C:\path\to\T3_Final_One_Batch_Log_Repair.zip' -DestinationPath (Get-Location) -Force
python -X utf8 -m unittest discover -s tests -p test_t3_thesis_final_one_batch.py -v
python -X utf8 -m scripts.RunT3ThesisFinalOneBatch --all --run-id thesis_final_one_batch_logrepair_v2 --debug
```

Do not use `--fit` for this diagnostic run. Preserve the successful F and Laurent sources. Expect **exit code 2** until the hard matching source, EFT1 Wilson RG, and D physical scalar continuation are implemented. If A scalar soft succeeds, it exports `A_scalar_soft_only_unrenormalized.wxf` and keeps `A_full_F_hard_one_loop_preserved.wxf`, with `ScalarSoftOnly_HardInheritanceAndWilsonRGPending`, but not physical C5.

**Validation:** Python parser/tests passed in the container; Wolfram/Matchete not installed in the local execution environment and not run here. Numerical physical results are not yet available.
