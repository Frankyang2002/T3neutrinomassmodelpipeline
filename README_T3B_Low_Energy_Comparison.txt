T3-B fixed-order versus exploratory coupled neutrino observables at low energy
=============================================================================

Extract the ZIP into your existing local project:
  C:\thesistemp\T3neutrinomassmodelpipeline
Do not modify GitHub. All paths below are relative to this project root.

Files supplied:
  Numerical/diagnostics/T3BLowEnergyComparison.py   NEW COMPLETE MODULE
  tests/test_t3b_low_energy_comparison.py            NEW TESTS
  README_T3B_Low_Energy_Comparison.txt                THIS FILE

The existing T3BExploratoryCoupledBenchmark.py is NOT replaced. Input artifacts
must already exist from the fitted T3-B fixed-order control and coupled run.

PowerShell:
  python -m pytest -q tests/test_t3b_low_energy_comparison.py

  python -m Numerical.diagnostics.T3BLowEnergyComparison `
    --config "Numerical\benchmarks\fitted_reference\T3_dS1_2_dS2_2_dF_1_alpha_m1.json" `
    --uv-rgbeta "output\single\T3_B_alpha_m1\data\uv_rgbeta_rge.json" `
    --eft1-rgbeta "output\single\T3_B_alpha_m1\data\eft1_rgbeta_rge.json" `
    --control "output\t3b_fitted_fixed_order_control.json" `
    --coupled "output\t3b_fitted_common_scale_threshold.json" `
    --ordering NO `
    --output "output\t3b_fitted_low_energy_comparison.json"

What the code does:
  1. Reads fixed-order and exploratory C5(MS) from your saved JSON outputs.
  2. Checks the existing production decomposition:
       fixed = hard + direct; exploratory = hard + coupled.
  3. Re-runs the SAME renormalisable UV and scalar-only EFT RGEs from the
     provided fitted config and RGBeta inputs, then verifies its MS boundary
     agrees with the saved fixed-order control's intermediate samples.
     This step is required because the control JSON did not save gY or g3.
  4. Builds matched SM+Weinberg states using the established
       Numerical.running.ScalarThresholdBoundary.build_sm_weinberg_initial_conditions
  5. Evolves both full 3x3 symmetric complex matrices via
       Numerical.running.SMWeinbergEvolution.evolve_sm_weinberg
     from MS=1000 GeV to mu_low=100 GeV with identical SM inputs.
  6. Applies the repository's m_nu = -(v^2/2) C5 convention, rotates into the
     charged-lepton mass basis and calls calculate_neutrino_observables.
  7. Saves masses, Delta m^2, |U_PMNS|, sin^2(theta_ij), angles in degrees,
     Takagi residuals, C5 low matrices, and coupled-minus-fixed differences.

Use --mu-low VALUE to override config scales.mu_low_gev; --vev-gev VALUE
for an explicit electroweak vev, and --ordering NO, IO or AUTO. The fitted
reference is expected to have normal ordering, so NO is the default.

Important scope:
  * The low-energy evolution itself uses the current production SM+Weinberg
    beta functions and is NOT a new or approximate RGE implementation.
  * The upstream coupled C5(MS) remains EXPLORATORY (not completely matched
    sequential heavy EFT; 3-generation operator closure not proven).
  * The repository uses an unbroken SM+Weinberg evolution to 100 GeV; this
    script does not implement electroweak-threshold matching or decoupling.
  * Therefore the results are exploratory low-scale numerical observables,
    NOT validated physical predictions.
  * Production source files and the existing control/coupled JSON are never
    modified. The script writes only the --output JSON.

Local authoring tests: 7 passed, 1 real-repository integration test skipped
because the full repository is not present in the authoring environment.
Please run the test and full benchmark in your local repository.
