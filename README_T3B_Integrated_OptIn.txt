T3-B: unified fixed-order pipeline with opt-in exploratory C12<->C5 LL
==================================================================

Project root: C:\thesistemp\T3neutrinomassmodelpipeline

Files (complete new files; production code untouched):
  Numerical/orchestration/T3BIntegratedCoupledPipeline.py
  tests/test_t3b_integrated_coupled_pipeline.py
  README_T3B_Integrated_OptIn.txt

Dependencies: This addon calls the repository's standard
Numerical.orchestration.PipelineNumericalResults.run_pipeline_numerical_results.
For --enable-coupled it also needs the modules installed in earlier work:
  Numerical/diagnostics/T3BHeavyLogAlignment.py
  Numerical/diagnostics/T3BExploratoryCoupledBenchmark.py
    (the COMMON-SCALE LOG BOUNDARY version)
  Numerical/diagnostics/T3BPerturbativeAccuracy.py
    (the REVERSE ORDER FIX version)
  RGE/running/intermediate/T3BC12FlavorRGE.py

The run's data/ directory must contain previously generated:
  uv_rgbeta_rge.json
  eft1_rgbeta_rge.json
  final_weinberg_coefficient.json

Test:
  python -m pytest -q tests/test_t3b_integrated_coupled_pipeline.py

One command, standard fixed-order only:
  python -m Numerical.orchestration.T3BIntegratedCoupledPipeline --config "Numerical\benchmarks\fitted_reference\T3_dS1_2_dS2_2_dF_1_alpha_m1.json" --output-dir "output\single\T3_B_alpha_m1"

One command, standard fixed-order + opt-in exploratory LL:
  python -m Numerical.orchestration.T3BIntegratedCoupledPipeline --config "Numerical\benchmarks\fitted_reference\T3_dS1_2_dS2_2_dF_1_alpha_m1.json" --output-dir "output\single\T3_B_alpha_m1" --enable-coupled --max-order 4

Outputs:
  output/single/T3_B_alpha_m1/data/running_diagnostics.json
    (normal, authoritative fixed-order production output)
  output/single/T3_B_alpha_m1/data/t3b_coupled_ll_exploratory.json
    (only with --enable-coupled, new separate artifact)

Interpretation:
  - All 3 heavy fermions integrated out together at common mu_F.
  - Keep individual heavy masses in the one-loop matching logarithms.
  - Evolve C12 and C5 simultaneously over mu_F to mu_S.
  - Subtract existing direct one-loop term to avoid double counting.
  - Use the standard SM+Weinberg RGE to run fixed and exploratory C5 to mu_low.
  - Include the perturbative Wilson-order audit with max_order=4 by default.
  - Production matching and final C5 are NOT replaced by exploratory values.
  - Complete higher-order matching and full 3-generation closure remain open.
  - The CLI is an ADDITIVE FRONT END, not a change to the existing pipeline's
    normal call sites. It invokes the standard numerical stage unchanged.
  - The standard numerical stage itself regenerates its usual diagnostics
    and figures, just as it would without this wrapper.

Test coverage: unit tests use stubs for standard pipeline dispatch and
check default-off behavior, separate-artifact writing and failure handling.
The full fitted benchmark must be exercised in the user's local repository.
