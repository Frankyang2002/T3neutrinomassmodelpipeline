# Batch 2 — conditional numerical result runner

Place the archive contents in the repository root. No existing files are overwritten.

Run in PowerShell from `C:\thesistemp\T3neutrinomassmodelpipeline`:

```powershell
$env:PYTHONUTF8 = "1"
python scripts/maintenance/RunBatch2Numerics.py --model A_p0
```

If automatic configuration discovery reports missing scenario/study files, supply the **existing** JSON files with `--scenario` and `--study`. This runner refuses to invent RGBeta-to-Matchete convention multipliers.

The results go to `output/audit/batch2/conditional_observables.json`; detailed findings and failures go to `output/audit/batch2/batch2_result.json`, with per-test logs in the same directory.

This uses the existing six-model UV -> direct one-scale C5 -> SM Weinberg low-energy calculation. It is **conditional** on explicitly unverified coupling maps and sudden decoupling, and **does not** certify sequential EFT1 matching. It is not a complete physical neutrino-mass prediction.

The four previously failing selected test commands reduce to two issues:
1. Eight failures in the hard-pole tests use Windows cp932 as the implicit encoding for Unicode `ϵ`; `PYTHONUTF8=1` addresses this test-environment issue without modifying scientific code.
2. The three quartic/triplet integration-test failures descend from `verify_matchete_source` in `scripts/QuarticActionIdentity.py`. The source guard expects an exact older Lagrangian assembly text; the current builder includes `assembled`, `Expand` and `RelabelIndices`. Do **not** bypass this gate without verifying the physical convention.

This archive deliberately does not change or certify the scientific source. The user's repository and Wolfram output are unavailable for full execution inside this environment; the local command is the authoritative numerical run.
