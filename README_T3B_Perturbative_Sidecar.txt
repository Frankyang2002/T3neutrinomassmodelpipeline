T3-B COMMON-SCALE COUPLED RGE: PERTURBATIVE ORDER BOOKKEEPING
=============================================================

Model: ordinary T3-B, alpha=-1, three heavy Majorana fermions integrated
TOGETHER at common MF = 100 TeV. This bundle does not implement sequential
decoupling or change the standard numerical pipeline's physical predictions.

Files installed (NEW: neither overwrites your existing code):
  Numerical/diagnostics/T3BPerturbativeAccuracy.py
  Numerical/orchestration/T3BCoupledLLSidecar.py
  tests/test_t3b_perturbative_accuracy.py
  README_T3B_Perturbative_Sidecar.txt

Dependencies: the versions of these previously installed local files used
in the project already:
  Numerical/diagnostics/T3BExploratoryCoupledBenchmark.py
  Numerical/diagnostics/T3BHeavyLogAlignment.py
  RGE/running/intermediate/T3BC12FlavorRGE.py

Run in PowerShell from C:\thesistemp\T3neutrinomassmodelpipeline:
  python -m pytest -q tests/test_t3b_perturbative_accuracy.py

  python -m Numerical.orchestration.T3BCoupledLLSidecar `
      --config "Numerical\benchmarks\fitted_reference\T3_dS1_2_dS2_2_dF_1_alpha_m1.json" `
      --uv-rgbeta "output\single\T3_B_alpha_m1\data\uv_rgbeta_rge.json" `
      --control "output\t3b_fitted_fixed_order_control.json" `
      --normalization "output\t3b_fitted_unfitted_normalization.json" `
      --coupled "output\t3b_fitted_common_scale_threshold.json" `
      --low-energy "output\t3b_fitted_low_energy_comparison.json" `
      --output "output\t3b_fitted_coupled_ll_sidecar.json"

This is an OPT-IN sidecar adjacent to the existing standard pipeline outputs.
The production pipeline, standard C5 and fitted predictions remain untouched.
There is no option to promote its C5 to production automatically.

Perturbative expansion of the implemented one-loop Wilson RGE:
   X = X[0] + X[1] + X[2] + X[3] + ...
   X[0](MF) = (sum_r C12_r, 0)
   C5[1](MF) = [4 sqrt(3)/3 * lambdaT3(MS)/(16 pi^2)]
                  * sum_r C12_r * ln(MF/Mr)
   d X[n]/d ln mu = hbar * Gamma_1(mu) X[n-1]  (n>=1)

Both Wilson matrices are 3x3, and C5 uses the physical Majorana factor of 2
in the audit. The renormalisable background is kept at its already integrated
one-loop values; its own perturbative expansion is not recomputed.

Wilson-loop order:
   C5[1]: direct C12 -> C5 plus the heavy-mass log boundary.
   C5[2]: first effects of C12 and C5 diagonal Wilson running.
   C5[3]: first possible C5 -> C12 -> C5 reverse-feedback effect.

The standalone audit tests one-loop frozen-coupling direct production agreement,
series reconstruction of the FULL COUPLED numerical C5, and reverse order.
The O(hbar) integrated result differs from the production frozen-coupling term
by higher-order renormalisable-running contributions; do NOT calibrate it away.

Crucial limitation: this is leading-log RESUMMATION of a subset of one-loop
EFT RGEs, NOT a complete O(hbar^2) Wilson coefficient. At the same formal
second order there can also be one-loop C12 threshold matching, RGE of the
one-loop hard C5 piece, genuine two-loop mixing terms and higher-order scalar
threshold matching. Full 3-generation operator closure also remains unproven.
The stage explicitly marks physical_prediction=false and will not replace
the established one-loop fixed-order C5. No GitHub changes.
