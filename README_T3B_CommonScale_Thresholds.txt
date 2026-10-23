T3-B common-scale, mass-resolved threshold-log boundary
======================================================

Location: C:\thesistemp\T3neutrinomassmodelpipeline
Extract this ZIP into the project root, preserving paths.

Files supplied:
  Numerical/diagnostics/T3BExploratoryCoupledBenchmark.py
      FULL REPLACEMENT of previous exploratory benchmark.
  tests/test_t3b_common_scale_threshold_logs.py
      New numerical regression tests.
  README_T3B_CommonScale_Thresholds.txt
      This documentation.

Requires existing files from previous installations (not replaced):
  Numerical/diagnostics/T3BHeavyLogAlignment.py
  RGE/running/intermediate/T3BC12FlavorRGE.py

PowerShell tests:
  python -m pytest -q tests/test_t3b_common_scale_threshold_logs.py

PowerShell fitted benchmark using saved fixed-order and normalization inputs:
  python -m Numerical.diagnostics.T3BExploratoryCoupledBenchmark `
    --config "Numerical\benchmarks\fitted_reference\T3_dS1_2_dS2_2_dF_1_alpha_m1.json" `
    --uv-rgbeta "output\single\T3_B_alpha_m1\data\uv_rgbeta_rge.json" `
    --control "output\t3b_fitted_fixed_order_control.json" `
    --normalization "output\t3b_fitted_unfitted_normalization.json" `
    --threshold-scheme common_scale_log_boundary `
    --compare-legacy-threshold `
    --output "output\t3b_fitted_common_scale_threshold.json"

The default scheme is common_scale_log_boundary.
To reproduce the previous approximation explicitly, use:
    --threshold-scheme endpoint_extrapolation

Physics and numerical conventions:
  * The established production flavour boundary is sum_r C12_r,
    where C12_r[p,q] = (y1*[p,r] y2*[q,r] + y2*[p,r] y1*[q,r])/(2 Mr).
  * At common MF, seed ordered C5 with
      C5(MF) = (4 sqrt(3)/3) lambdaT3(MS)/(16 pi^2)
               * sum_r C12_r ln(MF/Mr).
  * Integrate the coupled (C12, ordered-C5) system strictly MF -> MS
    using the saved scalar-only renormalisable background.
  * Convert ordered C5 to the physical Majorana coefficient by factor 2,
    consistent with the earlier absolute normalization audit.
  * Fixed-order sum of boundary log and ordinary EFT running log:
      ln(MF/Mr) + ln(MS/MF) = ln(MS/Mr).
    This identity is enforced against the production C5 direct matrix.
  * Subtract the production direct C5 before adding the coupled contribution:
      C5_exploratory = C5_production - C5_production_direct + 2*C5_coupled_ordered.

Meaning of new output fields:
  threshold_scheme
  threshold_scheme_details.mass_resolved_log_mf_over_mr
  threshold_scheme_details.heavy_boundary_c5_ordered
  common_boundary_fixed_order_direct_residual
  legacy_threshold_comparison (when --compare-legacy-threshold is used)

LIMITATIONS:
  This is NOT sequential matching at each heavy mass. The reference
  repository has only all-heavy UV RGEs and a scalar-only EFT below a single
  matching scale. The scheme instead encodes the known hierarchical one-loop
  logarithms as a Wilson boundary at that common scale. Above-MF scalar-only
  running is no longer extrapolated. Threshold-scheme difference is a
  sensitivity, not a controlled uncertainty or a physical correction.
  Generated eC-sector operators, full 3-generation EFT closure and higher-order
  threshold matching still require treatment. The reported masses are
  exploratory and at MS, not validated low-energy predictions.

No GitHub files are edited or pushed.
