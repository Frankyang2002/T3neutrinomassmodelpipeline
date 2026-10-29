# Frozen-UV Yukawa fit recovery

Purpose: improve currently unaccepted, mass-safe frozen-UV T3 candidates via fresh, independent scrambled Sobol runs and local least-squares Yukawa optimisation. Full UV running is excluded with `--cases`.

## Project files
- `scripts/ContinueT3GeneralScalarSobol.py`: complete replacement (adds `--cases` scenario filtering to the existing runner; original physics calls retained).
- `tests/test_t3_frozen_continuation.py`: tests for filtering and retention.
- `docs/T3_FROZEN_UV_FIT_RECOVERY.md`: this file.

The source `overnight_six_v1_bulk48` must exist. Existing accepted cases are copied and preserved; poor-fit results are replaced only if the new finite chi-square is lower. Neither the scalar UV input nor the C5 normalization is changed. The search varies the 18 real Yukawa entries within the existing `search_box`, using mass/perturbativity screens and the existing matching and fitting implementations.

**Important restriction:** This campaign cannot solve extreme-hypercharge cases A(-4), B(-3), C(-3), E(-4) whose legacy fast C5 parity test fails, and cannot guarantee a good fit in every model. B(-1)/C(-1) special cases are reused from separately validated matches. The full 48-case CSV is emitted for compatibility; selected scenarios only are newly evaluated.

## Windows PowerShell

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
python -X utf8 -m pytest -q tests\test_t3_frozen_continuation.py
python -u -X utf8 -m scripts.ContinueT3GeneralScalarSobol `
    --source-campaign overnight_six_v1_bulk48 `
    --run-id frozen_recovery_v1 `
    --cases small_uv_frozen_5e-9,large_uv_frozen_0p5 `
    --models A_+0,A_-2,B_+1,C_+1,D_+0,D_+2,D_-2,E_+0,E_+2,E_-2 `
    --sobol-points 256 --max-nfev 200 --max-rounds 12 `
    *> .\frozen_recovery_v1.log
```

Resume by appending `--resume` to exactly the same arguments. The run stores `continuation_progress.json`, `campaign.csv`, `campaign.json`, `summary.md`, and individual best `fit.json`/`fitted_config.json` under `output/studies/general_scalar_extended/bulk_48/frozen_recovery_v1/`. Intermediate rounds are saved in `attempts/` after each round. `--max-rounds` is the total budget for a run and cannot be changed on resume; use a new run ID for another budget. Running full frozen recovery may take considerable time.

## Status
Conditional general-quartic EFT1 only. UV RGEs are frozen for both targeted scenarios. No general holomorphic quadratic extensions or full sequential Wilson EFT matching. Successful neutrino fit scores are not validation of complete model physics.
