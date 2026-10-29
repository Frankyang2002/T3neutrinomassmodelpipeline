# General T3 extension: validated physics and next numerical gate

## Verified on the user's 16-model Windows Matchete run (Batch 5)

- 16/16 expected extra quartic bases were constructed and validated.
- All baseline and extended one-loop symbolic matching jobs completed.
- Formal zero-extra-coupling limit was proved for all 16 `C5` outputs.
- The only extra quartics observed in the formal one-loop Weinberg coefficient
  were `lambdaT3ExtraQ020020I1` and `lambdaT3ExtraQ020200I1` in B/C(-1).
- The other 14 models had text-identical one-loop `C5` output in the comparison.

These claims concern additional **scalar quartics** in the separately defined
ordinary two-complex-scalar models. They do not establish the complete
renormalisable UV theory including independent quadratic mixings and all
possible special-case fermionic couplings.

## RGBeta structural status

- Of 50 additional complex scalar quartics summed across the 16 separate
  model points, 12 have structurally analogous vertices already in RGBeta.
- 38 are still missing from the current RGBeta interaction set.
- No Matchete-to-RGBeta numerical coupling normalisation is verified yet.
- For the leading-order C5 effect, only two B/C(-1) structures require an
  immediate matching-to-RGBeta map. All other omitted structures may affect
  the renormalisable RG running even when absent from the fixed-order C5.

## Stage 6 (new)

1. Confirm the linear response of the C5 coefficient to the new couplings
   and identify which Yukawa structures multiply them.
2. Confirm additive sources in their existing RGBeta one-loop beta functions.
3. Export the exact Matchete tensor convention relative to an explicit
   SU(2) doublet tensor. This determines a Matchete-to-canonical coefficient,
   not yet the RGBeta-to-canonical coefficient.

Successful Stage 6 results do **not** certify the enlarged numerical pipeline.
The remaining numerical implementation should not equate differently
normalised symbolic couplings without verifying RGBeta's `delS2` conventions.

## PowerShell

```powershell
Expand-Archive -LiteralPath .\T3_General_Scalar_Batch6.zip -DestinationPath . -Force
python -m pytest -q tests\python\test_general_scalar_trial.py tests\python\test_scalar_quadratic_prototype.py tests\python\test_general_scalar_c5_proof.py tests\python\test_general_rgbeta_coverage.py tests\python\test_general_t3_scalar_quartic_response_audit.py
powershell -ExecutionPolicy Bypass -File .\scripts\RunGeneralT3ScalarQuarticResponseAudit.ps1
```

To revert the experimental version entirely, switch back to the pre-trial Git
branch or revert the experiment commit. Avoid using `git restore` on files
with independent local changes. No GitHub repository changes are made by
this package.
