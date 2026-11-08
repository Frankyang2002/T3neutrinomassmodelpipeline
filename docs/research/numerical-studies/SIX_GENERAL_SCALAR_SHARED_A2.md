# Six-model general-scalar shared-A(-2) comparison — preflight

**Source-only audit; not new physics results.**

Selected models: `A_-2`, `D_-2`, `B_-1`, `C_-1`, `B_-3`, `C_-3`.

Preferred starting scenario: `large_uv_frozen_0p5`, because the saved general-scalar `A_-2` reference is labelled `CandidateConstrainedAccepted` there. `small_uv_frozen_5e-9` has a `CandidateConstrainedPoorFit` A(-2) source, so do not use it as the main *optimised acceptable* reference.

From the repository root in PowerShell:

```powershell
python -X utf8 -m scripts.AuditSixGeneralScalarSharedA2
Import-Csv .\output\studies\general_scalar_extended\six_model_shared_a2_audit\preflight.csv | Format-Table -AutoSize
Get-Content .\output\studies\general_scalar_extended\six_model_shared_a2_audit\reference_A_m2_yukawas.json
```

This checks that the general-scalar campaign contains 48 distinct cases, verifies the A(-2) representation, reads its actual saved 3x3 `y1`/`y2`, records the source hashes, and prepares a six-model audit report. It does **not** refit or perform physical evolution. It explicitly labels missing shared-input results as `NOT_CALCULATED`.

**To perform the actual study:** use each target model's existing *general-scalar* matching/RGE evaluator, replace **only** its Yukawa input arrays with these reference arrays, preserve its other target-model inputs, and compute the full flavour-resolved `C5`, Takagi masses and PMNS observables at a common low-energy scale. Validate the matching scheme and scalar-mass threshold before calling outputs physical. Do not substitute the reduced-scalar evaluator. In the saved 48-case study, B/C(-3) had no mass-safe fit, and B/C(-1) use a conditional one-scale mixed-scalar matching approximation.

For an analytical comparison, check matrix proportionality: `k=Tr(C5_A^† C5_model)/Tr(C5_A^† C5_A)` and `epsilon=||C5_model-k*C5_A||_F/||C5_model||_F`. Only compute those from **actual newly evaluated** common-input C5 matrices.
