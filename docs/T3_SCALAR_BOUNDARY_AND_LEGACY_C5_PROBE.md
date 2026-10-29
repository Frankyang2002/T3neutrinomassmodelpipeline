# Diagnostic probe: UV scalar boundary and legacy C5 matrix parity

New files only; no existing physics files or fitted benchmarks overwritten.

- `scripts/ProbeT3ScalarBoundaryAndLegacyC5.py`
- `tests/test_t3_scalar_boundary_legacy_probe.py`
- `docs/T3_SCALAR_BOUNDARY_AND_LEGACY_C5_PROBE.md`

## Run (PowerShell from repository root)

```powershell
python -X utf8 -m pytest -q tests/test_t3_scalar_boundary_legacy_probe.py
python -u -X utf8 -m scripts.ProbeT3ScalarBoundaryAndLegacyC5 --run-id a2_scalar_probe_v1 --beta-campaign overnight_six_v1_bulk48
```

Optional other full-UV models:

```powershell
python -u -X utf8 -m scripts.ProbeT3ScalarBoundaryAndLegacyC5 --run-id fulluv_scalar_probe_v1 --models A_-2 A_+0 B_+1 --factors 0.25 0.5 0.75 1.0 1.25
```

## What it does

It multiplies each *UV scalar mass* by factors 0.25, 0.5, 0.75, 1.0, 1.25; internally the mass-squared inputs are multiplied by the factor squared. The 18 historical fitted real Yukawas, topology quartic, extra quartic boundaries, RG beta functions, matching sources, and all screening thresholds are unchanged. It integrates candidate UV 1e7->1e5 GeV and EFT1 1e5->1e3 GeV for each trial. It reports scalar masses and each screening inequality. **This is a diagnostic boundary scan, not a valid fit or automatic replacement for the original benchmark boundary.** An endpoint that passes the screens is not a full EFT certification.

Separately for A(-4) high-scalar/frozen UV, it reruns the saved EFT1 candidate background and compares fast legacy C5 and authoritative `evaluate_symbolic_c5`. It saves both complex matrices and tests possible global sign, conjugation, transpose and scalar-normalisation relationships *only as diagnostics*. No matching tolerances or formulas change.

Output:

`output/studies/general_scalar_extended/diagnostics/<run-id>/probe.json`

`output/studies/general_scalar_extended/diagnostics/<run-id>/scalar_boundary.csv`

The saved `probe.json` includes `legacy_parity_probe.analysis`. Inspect the numeric values before deciding whether a mass boundary or C5 kernel correction is justified. The full calculation requires the existing saved beta exports, numerical Python environment and Mathematica-generated artifacts; no Wolfram call is performed by this script.
