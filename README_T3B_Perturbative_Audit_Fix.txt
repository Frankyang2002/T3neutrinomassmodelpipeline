T3-B PERTURBATIVE AUDIT: FALSE REVERSE-MIXING FAILURE FIX
===========================================================

Problem
-------
The previous perturbative_accuracy() computed the reverse-on and reverse-off
Taylor towers in two independent scipy solve_ivp integrations. The higher
Wilson orders change the adaptive integration mesh; comparing tiny O(hbar)
and O(hbar^2) outputs from separate meshes can give a false nonzero residual
in those orders. The exact order counting, with C5[0]=0, requires their
reverse-mixing contributions to be identically zero. Reverse feedback first
reaches C5 at O(hbar^3).

Fix
---
Integrate the reverse-on and reverse-off Taylor towers together in ONE
solve_ivp call with a shared adaptive mesh. Preserve the existing output
schema and the single-branch _solve_expansion() function for callers.
The sidecar command and matching scheme are unchanged.

Full replacement project paths:
  Numerical/diagnostics/T3BPerturbativeAccuracy.py
  tests/test_t3b_perturbative_accuracy.py
New instruction file:
  README_T3B_Perturbative_Audit_Fix.txt

Extract ZIP into:
  C:\thesistemp\T3neutrinomassmodelpipeline

PowerShell from repository root:
  python -m pytest -q tests/test_t3b_perturbative_accuracy.py

  python -m Numerical.orchestration.T3BCoupledLLSidecar `
      --config "Numerical\benchmarks\fitted_reference\T3_dS1_2_dS2_2_dF_1_alpha_m1.json" `
      --uv-rgbeta "output\single\T3_B_alpha_m1\data\uv_rgbeta_rge.json" `
      --control "output\t3b_fitted_fixed_order_control.json" `
      --normalization "output\t3b_fitted_unfitted_normalization.json" `
      --coupled "output\t3b_fitted_common_scale_threshold.json" `
      --low-energy "output\t3b_fitted_low_energy_comparison.json" `
      --output "output\t3b_fitted_coupled_ll_sidecar.json"

This is a numerical-audit fix, not a change to the C12/C5 coupled physical
model or production matching. The sidecar remains explicitly exploratory.
All three heavy fermions continue to decouple simultaneously at MF.

Tested here: 7 passed (including a fitted-scale weak complex-Yukawa regression)
in the available local environment. Full benchmark can run with the user's
saved fitted control JSON and RGBeta output.
