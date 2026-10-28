# Scalar-potential boundary audit and separation of historical versus extended T3 studies

**Purpose:** preserve the original fitted benchmarks, check the saved restricted scalar inputs, and reserve a separate namespace for any future calculations with additional general scalar operators. No historical fit, C5 artifact or numerical pipeline is modified by this package.

## Historical studies audited

The audit reads all 16 fitted configurations and their fit records in each of the following folders:

| Case key | Existing study input folder | `lambdaT3` classes A–D | Class E | Actual numerical boundary |
|---|---|---:|---:|---|
| `small_uv_running_5e-9` | `output/studies/lambda_5e-9_refit/` | `5e-9` | `2.5e-9` | UV at `1e7` GeV; normal UV RGEs enabled |
| `small_uv_frozen_5e-9` | `output/studies/fermion_start_lambda_comparison/5e-9/` | `5e-9` | `2.5e-9` | Fermion threshold `1e5` GeV; UV RGEs disabled |
| `large_uv_frozen_0p5` | `output/studies/fermion_start_lambda_comparison/0p5/` | `0.5` | `0.25` | Fermion threshold `1e5` GeV; UV RGEs disabled |

The original `base_state.mu_gev` field in the UV-frozen configurations may still read `1e7` GeV; the code uses `frozen_result(initial, mu_f, uv=True)`, imposing the values at `1e5` GeV. It must not be called UV running.

A fourth available study (`5e-10`, E `2.5e-10`) is **not included** in this three-study audit. It can be added in a separate follow-up without changing any existing results.

## What a zero-boundary audit can and cannot prove

The candidate scalar registry (`RGE/general/generated/general_rgbeta_quartic_registry.json`) has **50 additional quartic model-slots**: 12 pre-existing structural quartic slots in the ordinary numerical states and 38 newly registered quartic slots absent from them. There are eight additional quadratic model-slots not implemented by the original state.

The audit checks all stored reused quartics as **exactly zero** (`real=imag=0`) and reports model/study exceptions. It records missing new quartic and quadratic parameters as **unrepresented**, not as independently certified zero-valued inputs. Original restricted-potential calculations effectively set these absent terms to zero **by model definition**, not by evolving or checking a general model.

The original baseline scalar self-, Higgs-portal and scalar–scalar-portal couplings remain present. Do not confuse them with the additional quartics catalogued by the general registry.

Refer to `docs/T3_REDUCED_SCALAR_POTENTIAL_CONCLUSION.md` (prior Batch 8 package) and `output/scalar_relevance/final_leading_scalar_potential.md` for the fixed-order relevance derivation. Relevant findings:

- Four extra quartic slots (two each in B(-1)/C(-1)) can enter the one-loop Weinberg matching **if independently nonzero**.
- The first RG-induced extra quartics from an all-extra-zero boundary are beyond strict fixed-order leading-one-loop `C5` when inserted in that graph.
- Allowed additional quadratic mixing can affect the one-loop amplitude at leading order if nonzero; it is **not** generically forbidden by the gauge and Z2 symmetries.
- The B/C numerical magnitude from Batch 7 is an **illustrative fitted-reference point**, **not** a validated universal bound on the three historical input studies (particularly `lambdaT3=0.5`).

## Saved fitting quality

The source fit records have both `status` and `chi2_fit`. The audit keeps every model, including `PoorFit`, and reports the number satisfying **both** `status=Success` and `chi2_fit <= 5.86`. A comparison summary may use a different status convention (`Success` for numerical evaluation) even if the resulting chi2 exceeds the fitting screen: these are not synonymous.

Saved comparisons based on independently fitted Yukawa matrices are not same-seed perturbations of a scalar coupling. Refit only when the requested new physics calculation actually needs a new fitted input; retain original fits unchanged for direct provenance.

## Separate output trees (mandatory convention)

Historical directories are read-only and **never overwritten**:

- `output/studies/lambda_5e-9_refit/`
- `output/studies/fermion_start_lambda_comparison/`
- `output/general_t3_stage22/` (Stage 30–36)
- `output/scalar_relevance/` (Batches 1–8)

This **input-only** audit writes only:

- `output/studies/scalar_potential_boundary_audit/restricted_input_audit.json`
- `output/studies/scalar_potential_boundary_audit/restricted_input_audit.csv`
- `output/studies/scalar_potential_boundary_audit/restricted_input_audit.md`

**All future computations including the additional general scalar operators** must be written into their own run folders, e.g.:

```text
output/studies/general_scalar_extended/
    small_uv_running_5e-9/<run_id>/
    small_uv_frozen_5e-9/<run_id>/
    large_uv_frozen_0p5/<run_id>/
```

`<run_id>` must be unique for every new numerical attempt and identify its source code/input snapshot. Keep original restricted benchmarks as cross-references only. Do not run a general-scalar calculation into one of the old fitted, comparison, or Stage 36 output folders. Future code should default to **no overwrite**, accepting existing run directories only after an explicit user-controlled choice.

Each future run must record:

1. Study key, model, UV-on or UV-frozen boundary scale, full fitted-input config pathname and SHA-256.
2. All 50 relevant quartic **model-slots** and eight bilinears as applicable to the chosen model, identifying explicit zero, nonzero, and unsupported/unrepresented coefficients.
3. Actual beta inventory and UV/EFT1 field content; whether fermion and scalar thresholds are treated at tree or loop level.
4. Exact `C5` source, map/basis convention, and whether the Matchete/RGBeta identification is *certified* or *candidate*.
5. Fit/source provenance, `chi2`, numerical tolerances, successful solver status and limitations.
6. A link/path to the **unchanged** restricted baseline used for comparison.

The existence of this reserved folder is **not** a claim that the fully general scalar sector has been simulated. A numerical calculation using additional quartics is not physical `C5` matching until its coupling normalization and threshold operators are validated.

## PowerShell

From the local repository root after extracting the ZIP:

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
python -m pytest -q tests\python\test_t3_saved_scalar_boundaries.py
python -m scripts.AuditT3SavedStudyScalarBoundaries
Get-Content -Encoding UTF8 output\studies\scalar_potential_boundary_audit\restricted_input_audit.md
```

No RG, optimisations, or matching calculations run. Inputs are SHA-256 fingerprinted with platform-independent `/` keys. The audit's own outputs may be regenerated in place; historical study inputs/outputs are never modified.

**Inspected repository snapshot:** GitHub `main` commit `a02e04c6def5cf8522b2f612ccd6bff7d9963b91` (read-only); local source files take precedence if they differ from this snapshot.
