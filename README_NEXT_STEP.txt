Standalone T3-B coupling-trajectory bridge (no production modifications).

REPLACEMENT FILE (replaces previous standalone bundle version):
  Numerical/running/CoupledT3BWilson.py

NEW FILES:
  Numerical/running/T3BWilsonBackground.py
  tests/test_t3b_wilson_background.py

Run from project root:
  python -m pytest -q tests/test_t3b_wilson_background.py tests/test_coupled_t3b_mixing.py

The existing Numerical.running.IntermediateScalarRunner.run_intermediate_scalar_segment(...)
returns an IntermediateScalarRunningResult object. In the Python code where
that result is available, add (for a local diagnostic only):

  from Numerical.running.T3BWilsonBackground import export_background
  export_background(intermediate_result, 'output/t3b_scalar_background.json')

Do not insert this call until the current code path providing
'intermediate_result' has been identified in your local project. This module
intentionally refuses 3x3 Yukawa matrices, because projecting 3 flavours onto
one diagonal element would be unphysical. It is therefore compatible only
with a genuine one-generation intermediate-scalar state.

If you have a genuine one-generation background JSON, run:

  python -m Numerical.running.CoupledT3BWilson \
    --validation output/c12_c5_validation.json \
    --mapping output/c12_c5_exact_quartic_map.json \
    --background output/t3b_scalar_background.json \
    --s0 '[1e-5,0]' --kappa0 '[0,0]' \
    --output output/t3b_coupled_diagnostic.json

The example s0 is not your matched physical Wilson coefficient.
Use your actual matched one-generation coefficient for quantitative comparisons.
Nothing in this bundle connects to threshold matching or 3-flavour numerics.
