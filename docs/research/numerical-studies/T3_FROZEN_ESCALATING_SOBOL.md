# Frozen UV Sobol escalation

Full replacement: `scripts/RunT3FrozenEstablishedSobol.py`.

For every frozen model/scenario, preserve accepted source fits. For unresolved cases, sample independent scrambled Sobol batches, doubling the number of points per round (up to the cap), applying the established multi-start 5/8/12 active-parameter local refinement to retained candidates. Retain the lowest chi2 in checkpoint `round_NNN.json`. A case advances only if chi2 <= target; if all rounds fail, exit code 2 and do not advance to the next case. `--continue-on-failure` optionally overrides. A physically blocked case, including the legacy C5 parity failure, terminates after the bounded budget; more Sobol cannot fix the kernel.

Default escalation: 128, 256, 512, 1024, 2048, 4096, 4096, 4096 (8 rounds). The existing fitted benchmark and normal scalar physical checks are unchanged.

Use a **new run ID** if previous `result.json` is already in the directory. Existing `result.json` is read without repeating calculations. Resume interrupted rounds using the same run ID. The source campaign is still `overnight_six_v1_bulk48`.

```powershell
python -u -X utf8 -m scripts.RunT3FrozenEstablishedSobol --source-campaign overnight_six_v1_bulk48 --run-id frozen_escalating_v1 --sobol-points 128 --max-sobol-points 4096 --rounds 8 --seed-count 5
```

Status `ExistingSpecialFit` refers to B(-1)/C(-1) independently validated cases using a separate matching evaluator, not a new search. Results use conditional general-quartic candidate RGEs and matching, not complete general scalar EFT.
