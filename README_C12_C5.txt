T3-B alpha=-1 independent combined-workflow prototype.

Paths:
  validation/DeriveT3BCrossQuarticMap.py
  Numerical/running/CoupledT3BWilson.py
  tests/test_coupled_t3b_mixing.py

Step 1: python -m pytest -q tests/test_coupled_t3b_mixing.py
Step 2: python -m validation.DeriveT3BCrossQuarticMap --data-dir "output\single\T3_B_alpha_m1\data" --output "output\c12_c5_exact_quartic_map.json"
Step 3: provide real scalar running background JSON with samples [{mu_gev,g2,lambdaH,lambda12,lambda12Adj,lambdaT3,yu,yd,ye},...]. Yukawas here are single-generation complex numbers, NOT full matrices.
Step 4: python -m Numerical.running.CoupledT3BWilson --validation output\c12_c5_validation.json --mapping output\c12_c5_exact_quartic_map.json --background PATH_TO_BACKGROUND.json --s0 "[1e-5,0]" --kappa0 "[0,0]" --output output\c12_c5_coupled_trial.json

This module is a standalone one-generation test; does not alter the production result, three-generation flavour matching, or scalar threshold.
