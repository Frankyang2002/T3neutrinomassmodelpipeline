T3-B operator-closure scan -- one-generation production tensor RGE
=================================================================

New project files:
  RGE/running/intermediate/T3BOperatorClosure.py
  tests/test_t3b_operator_closure.py

Required existing (unchanged) project files:
  RGE/running/intermediate/AuditC12C5Mixing.py
  RGE/running/intermediate/ScalarOnlyWilsonTensorRGE.py
  RGE/general/WilsonTensorRGE.py
  RGE/general/AnomalousDimensions.py

Extract the ZIP at:
  C:\thesistemp\T3neutrinomassmodelpipeline

Unit tests:
  python -m pytest -q tests/test_t3b_operator_closure.py

Full one-generation L/eC tensor scan, source-resolved:
  python -m RGE.running.intermediate.T3BOperatorClosure `
    --data-dir "output\single\T3_B_alpha_m1\data" `
    --scope llec `
    --output "output\t3b_operator_closure.json"

If the first scan finds additional components, explicitly test feedback from
that residual tensor into the retained LLSS/LLHH channels:
  python -m RGE.running.intermediate.T3BOperatorClosure `
    --data-dir "output\single\T3_B_alpha_m1\data" `
    --scope llec `
    --feedback `
    --output "output\t3b_operator_closure_feedback.json"

Both commands reconstruct the matched one-generation Wilson seed and run the
production master one-loop tensor RGE, including its scalar/fermion anomalous
dimensions, Yukawa terms, gauge terms and scalar-quartic terms.

IMPORTANT ABOUT INTERPRETATION
------------------------------
- The one-generation L/eC scan tests only the implemented psi^2 phi^2 tensor
  sector. It is not a complete EFT classification or a three-generation proof.
- The 'pivot_coefficients' are algebraic coordinates. If the residual is
  nonzero, they do not define physical anomalous dimensions for a closed 2x2
  operator basis.
- A nonzero residual component marks failure of exact two-template closure.
- --feedback applies the full RGE a second time to the unretained part and
  checks its LLSS/LLHH contribution. These are beta-numerator iterations, not
  physical Wilson coefficients; no logarithm or 16pi^2 factors have been
  inserted into a resummed result.
- The calculation never modifies the production numerical pipeline.
- When results are available, inspect 'results.C12' and 'results.C5' separately.
- No neutrino-mass prediction or revised physical 2.51% claim is made.

The seven pure symbolic projection tests pass in this environment. The full
production tensor scan requires the user's local repository and JSON exports.
