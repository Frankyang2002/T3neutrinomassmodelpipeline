# T3 general-scalar targeted screen/parity diagnosis

This is a read-only numerical diagnostic for the existing candidate 16x3 general-quartic campaign. It does not change RGE equations, screening bounds, matching, or any saved fit.

New project files:
- `scripts/DiagnoseT3GeneralScalarScreens.py`
- `tests/test_t3_general_scalar_diagnostics.py`
- `docs/T3_GENERAL_SCALAR_DIAGNOSTICS.md`

Run from project root:

```powershell
python -X utf8 -m pytest -q tests/test_t3_general_scalar_diagnostics.py
python -u -X utf8 -m scripts.DiagnoseT3GeneralScalarScreens --run-id screen_parity_v1
```

Default evaluates historical fitted seed for A(+0), A(-2), B(+1) low/full UV, plus A(-4) high/frozen UV for exact-vs-fast legacy C5 matrix parity. The full UV sample may be slow; no Sobol search is launched by default.

Optional 4 fresh Sobol diagnostic points for A(-2):

```powershell
python -u -X utf8 -m scripts.DiagnoseT3GeneralScalarScreens --run-id a2_screen_4_v1 --cases A_-2:small_uv_running_5e-9 --sobol-points 4
```

Outputs: `output/studies/general_scalar_extended/diagnostics/<run-id>/diagnostics.json` and `diagnostics.csv`. The JSON includes C5 matrices for any mass-safe legacy parity test. Requires saved candidate beta files in `bulk_48/overnight_six_v1_bulk48/betas` or `bulk_48/bulk48_v1/betas`.

**Interpretation**: `ParityFailed` means the saved fast legacy C5 kernel and the authoritative symbolic evaluator disagree. Do not relax tolerance or use those models for fits until the discrepancy is resolved. `ScreenRejected` specifies individual failed inequalities. This is candidate general-quartic work; it is not certified full EFT matching.
