# Stage 10: numerical **pilot**, not completed general T3 model

## Confirmed from user's Stage 9 run

- 44 Python tests passed.
- Matchete tensor/canonical SU(2) factors equal **1** for both C5-active quartic contractions in ordinary B(-1) and C(-1).
- Earlier checks: 16/16 quartic enumeration, formal zero-extra C5 regression, linear response extraction and additive RGBeta source terms succeeded.
- **NOT** verified: RGBeta's `delS2` tensor, overall coupling-potential sign/normalisation, or a Matchete↔RGBeta map.

## Implemented in this package

- `Numerical/experimental/GeneralT3C5Pilot.py`: reads the existing formal `delta_c5_atoms.json` for B/C(-1), parses the *actual Matchete derivative expressions* through `RGE.matching.MatcheteC5Parsing`, requires a single `conjugate(y1)^2` or `conjugate(y2)^2` flavour structure, substitutes the hard matching scale `μbar2=MF²`, refuses ultraviolet regulator poles, and evaluates a symmetric flavour response matrix for diagonal, positive heavy masses. Input extra couplings are **Matchete** coefficients, not RGBeta couplings. This is the **hard-boundary response only**, not the final `C5` after sequential running.
- `Numerical/experimental/GeneralT3RGSourcePilot.py`: reads the independently verified RGBeta source expressions at vanishing six special quartics, evaluates one-loop `β` via the **existing** `RGBetaEvaluator.evaluate_inputform_expression`, and reports a frozen-background leading-log step, `δλ = β^(1)/(16π²) log(μIR/μUV)`. The input `rgbeta_lambdaT3` is an RGBeta coefficient.
- Two CLI wrappers, one PowerShell multi-model runner, 10 new Python tests and a clearly marked **toy numerical input**.
- No production matching, RGE or fitted-parameter code was modified by Stage 10. Existing experiment files are retained in the consolidated archive.

## Commands (PowerShell; repository root)

Place `T3_General_Scalar_Batch10.zip` in the repository root:

```powershell
Expand-Archive -LiteralPath .\T3_General_Scalar_Batch10.zip -DestinationPath . -Force
python -m pytest -q tests\python\test_general_scalar_trial.py tests\python\test_scalar_quadratic_prototype.py tests\python\test_general_scalar_c5_proof.py tests\python\test_general_rgbeta_coverage.py tests\python\test_general_t3_stage6.py tests\python\test_general_t3_stage10.py
powershell -ExecutionPolicy Bypass -File .\scripts\RunGeneralT3Stage10.ps1
```

Inputs are the saved Batch 8 `output/general_t3_trial/{B_m1,C_m1}/general_c5_algebraic/delta_c5_atoms.json` and `output/general_rgbeta_source_proof.json`. The PowerShell command uses the synthetic `docs/general_c5_toy_benchmark.json` and writes:

```
output/general_t3_numerical_pilot/B_m1_matchete_hard_c5.json
output/general_t3_numerical_pilot/C_m1_matchete_hard_c5.json
output/general_t3_numerical_pilot/B_m1_rgbeta_source.json
output/general_t3_numerical_pilot/C_m1_rgbeta_source.json
```

For real inputs, copy the toy JSON, change the values and pass `-Config path\to\your.json`. **Do not** equate `matchete_extra_couplings` with RGBeta's special couplings until the tensor and Lagrangian normalisations have been checked.

### Interpretation and limitations

1. The hard C5 response is evaluated at `μ=MF` in the matching calculation's conventions. It is **not** the fully renormalised hierarchical final coefficient or a neutrino-mass prediction; the existing complete pipeline handles the original contribution separately.
2. The formal Matchete derivative must be free of UV regulator poles after the hard-scale substitution. If a `1/ϵ` term appears, **the pilot refuses to evaluate it**. Proper pole/RGE subtraction must be validated first.
3. The flavour lift assumes a *single symmetric same-species Yukawa bilinear* and diagonal heavy fermion mass eigenvalues (matching the existing numerical heavy-mass-basis handling). Any other flavour dependence is rejected.
4. The frozen-background RGBeta source is a first-order leading-log estimate, not integrated running. Its numerical value is independent of the Matchete coupling input.
5. Zero new quartics at the UV input are not an RG invariant restriction in B/C(-1): the measured source terms are nonzero in general. Inserting their loop-generated values in a one-loop mass is formally of higher order.
6. The other 14 ordinary models have no new quartic dependence in the *observed one-loop C5 artifacts*; they are unchanged by this numerical pilot. Quadratic scalar mass terms and general mass mixing remain unimplemented.

## Next full-theory steps

- Extract and calibrate RGBeta's symmetric `delS2` invariant coefficients, including both field orientations and overall factors in `AddQuartic`. Matchete/canonical = 1 **alone does not provide this map**.
- Complete the missing RGBeta quartic basis and expose its corresponding UV/intermediate numerical states and solvers with tested transformations.
- Extend full hierarchical matching and numerical final-C5 construction to allow new independent hard contributions and RG effects; validate the zero-new-coupling benchmark and scale consistency.
- Treat the symmetry-allowed quadratic mass matrices (and special real/complex scalar choices) separately; this changes propagators, thresholds and mass bases.

The broader objective remains a general T3 framework, but the current status should be stated precisely as **general scalar-quartic symbolic matching plus an experimental numerical hard-response/leading-log study**.
