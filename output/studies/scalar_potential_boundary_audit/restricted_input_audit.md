# Three historical T3 studies — scalar boundary audit

Read-only inspection of saved fitted configurations and fit records. No fits, RGBeta integrations, Matchete matching or historical outputs were rerun.

| Study | lambdaT3 A–D (E) | Input boundary | Models | Chi2 <= 5.86 and fit Success | Explicit reused quartics zero | New quartics absent | Quadratics absent |
|---|---:|---|---:|---:|---:|---:|---:|
| `small_uv_running_5e-9` | 5e-09 (2.5e-09) | mu_UV=1e7 GeV (UV on) | 16 | 16 | 12 / 12 | 38 / 38 | 8 / 8 |
| `small_uv_frozen_5e-9` | 5e-09 (2.5e-09) | mu_F=1e5 GeV (UV frozen) | 16 | 15 | 12 / 12 | 38 / 38 | 8 / 8 |
| `large_uv_frozen_0p5` | 0.5 (0.25) | mu_F=1e5 GeV (UV frozen) | 16 | 13 | 12 / 12 | 38 / 38 | 8 / 8 |

## Meaning of the checks

- **Explicitly zero:** one of the 12 model-slots of additional quartics already represented in the old fitted input state, recorded as exactly zero real and imaginary parts.
- **Absent / unrepresented:** one of the 38 newer quartic slots or eight additional bilinears is not a parameter in the saved fitted input. It is **not** an independently checked zero-valued coefficient in a general-potential model.
- UV-frozen fitted inputs retain `base_state.mu_gev=1e7` metadata, but their numerical values are **imposed at mu_F=1e5**; UV equations were not integrated.
- Failed / poor-fitting cases remain in the output, identified by the fit record and chi2 threshold. Matching/fit success and operator completeness are different tests.

## Scope of the scalar-potential conclusion

The original restricted scalar potential remains a well-defined fixed-order leading-one-loop benchmark when absent additional terms are **assumed zero by model restriction**. The audit does not reconstruct a fully general scalar potential or certify quantum corrections, nonzero coupling maps or scalar vacuum stability.

## Separate future general-scalar numerical studies

Store new results only under `output/studies/general_scalar_extended/` with distinct run identifiers, input hashes, RG/matching scope and original-comparison references. **Never write to** `output/studies/lambda_5e-9_refit/`, `output/studies/fermion_start_lambda_comparison/` or the saved Stage 36 / Batch 7–8 output trees. See `docs/T3_GENERAL_SCALAR_STUDY_SEPARATION.md`.

## Reproducibility

- Source JSON files fingerprinted: **99** (SHA-256 keys use portable `/` separators).
- Audit status: `ExplicitReusedQuarticsZeroAcrossAllCases`.
- A separate CSV contains the 48 model-study records, including every poor-fit status and the distinction between explicit zero and absent parameters.
