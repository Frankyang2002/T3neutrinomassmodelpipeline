T3-B SOURCE-DERIVED C12 FLAVOUR RGE -- LOCAL REPLACEMENT
======================================================

Extract the ZIP into:
  C:\thesistemp\T3neutrinomassmodelpipeline

Files (relative to the project root):
  RGE/running/intermediate/T3BC12FlavorRGE.py        NEW
  Numerical/diagnostics/T3BExploratoryCoupledBenchmark.py  COMPLETE REPLACEMENT
  tests/test_t3b_c12_source_derived.py              NEW

Source-derived projected beta:
  16*pi^2 * d C12/d ln(mu) =
     (-3*g2^2 + 2*lambda12 - lambda12Adj/2)*C12
   + (Ye*Ye^dagger*C12 + C12*(Ye*Ye^dagger)^T)/2
   + (8*sqrt(3)/3)*conjugate(lambdaT3)*C5_ordered

The charged-lepton flavour term is inherited from the fermion collinear
anomalous dimension in the generic master tensor RGE. The gauge/quartic and
reverse-mixing factors are the already-audited T3-B one-generation projection.

This code does NOT prove full operator-space closure or sequential heavy
threshold decoupling. The fit-point neutrino masses remain exploratory.

PowerShell (run from repository root):

  python -m pytest -q tests/test_t3b_c12_source_derived.py

  python -m RGE.running.intermediate.T3BC12FlavorRGE `
      --verify-master `
      --output "output\t3b_c12_flavour_master_verification.json"

  python -m Numerical.diagnostics.T3BExploratoryCoupledBenchmark `
      --config "Numerical\benchmarks\fitted_reference\T3_dS1_2_dS2_2_dF_1_alpha_m1.json" `
      --uv-rgbeta "output\single\T3_B_alpha_m1\data\uv_rgbeta_rge.json" `
      --control "output\t3b_fitted_fixed_order_control.json" `
      --normalization "output\t3b_fitted_unfitted_normalization.json" `
      --output "output\t3b_fitted_c12_source_derived.json"

The scalar-threshold result is expected to remain numerically very close to
(or identical to) the previous fitted exploratory output, since the previous
projected C12 matrix action was already algebraically correct. This upgrade
supplies a direct source derivation and independent complex Yukawa regression;
it does NOT claim a new numerical effect.
