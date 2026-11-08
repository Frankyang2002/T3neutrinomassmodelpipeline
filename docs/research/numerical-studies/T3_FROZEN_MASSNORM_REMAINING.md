# Remaining frozen models, mass-normalized search

Install `scripts/RunT3FrozenMassNormalizedRemaining.py` in the repository; retain the already installed `scripts/SearchT3FrozenPersistent.py` from the successful A(+0) mass-normalized package.

The campaign uses `output/studies/general_scalar_extended/bulk_48/overnight_six_v1_bulk48/campaign.json` to select unresolved cases. It excludes UV-running, the B(-1)/C(-1) special mixed-scalar evaluator, and four legacy C5 parity-unverified models A(-4), B(-3), C(-3), E(-4). It skips already accepted baseline cases, and recognizes accepted checkpoint states in `--previous-run`. It preserves all model-level checkpoints and fail-stops on nonacceptance. A future parity correction would be required before the four excluded models can be fitted reliably.

Each selected model runs the **same mass-normalized** Sobol + all-18 least-squares fitting method that accepted A(+0). No UV or matching physics changed. Set `--max-cycles 0` for unlimited cycles, interrupted with Ctrl+C. When restarting, pass `--resume` and the identical inputs. A different run-prefix starts a separate campaign.

Run tests:

```powershell
python -X utf8 -m pytest -q tests\test_t3_frozen_massnorm_remaining.py tests\test_t3_frozen_persistent.py
```

Inspect selected models without starting fits:

```powershell
python -u -X utf8 -m scripts.RunT3FrozenMassNormalizedRemaining --run-prefix frozen_massnorm_remaining_v1 --previous-run frozen_massnorm_a0_v1 --dry-run
```

Run:

```powershell
python -u -X utf8 -m scripts.RunT3FrozenMassNormalizedRemaining --run-prefix frozen_massnorm_remaining_v1 --previous-run frozen_massnorm_a0_v1 --from-run frozen_massnorm_a0_v1 --points 256 --seeds 3 --nfev 300 --radii 0.02,0.05,0.10 *> .\frozen_massnorm_remaining_v1.log
```

For automatic continuation after an interruption or an unfinished model, run the **same** command with `--resume`. Do not run two instances with the same run prefix concurrently.
