# Frozen fitted-reference benchmarks (16 T3 models)

These files are **NuFIT-calibrated reference points**, selected from the first-pass feasibility search and its retries. All 16 candidates were independently reevaluated by the existing numerical pipeline, and their reported chi-squared values agreed. They are intended as fixed starting points for within-model and between-model T1/T2/T3 sensitivity comparisons.

## What was fitted
Five normal-ordering oscillation observables: delta_m21_sq_ev2, delta_m3l_sq_ev2, sin2_theta12, sin2_theta13, sin2_theta23. These are fitted targets, **not independent predictions**. The solutions are not unique, and other phenomenological constraints have not been established.

## Coupling-regime warning
The input configurations originated in `smallY_smallL`, but the feasibility optimisation adjusted Yukawa entries and held lambdaT3 fixed. The labels `smallY` and `smallL` have **not been independently validated against scenario boundaries**. Treat this set as `fitted_reference` rather than proof of a coupling regime. Manifest records each lambda magnitude and maximum Yukawa entry for checking against the intended definitions.

## Class E
The class-E `(3,3,2)` configurations retain their original representation metadata. No new multiplicity correction was inserted. The numerical machinery used to evaluate them must continue to apply its own representation-specific treatment.

## Usage
Extract ZIP at the repository root. These are immutable reference input configurations. Do not overwrite original `configs/generated_model_sets/comparison` files or previous fit/sensitivity diagnostics. `benchmark_search.enabled` and `sensitivity.enabled` are false in the frozen files to avoid accidental refitting or scans. The existing fitted sensitivity output remains the verification record.

## Source provenance
`benchmark_manifest.json` records each chosen feasibility stage, chi-squared, independently verified chi-squared, predictions, coupling norms and SHA-256 of the exact installed config bytes. All 16 models have fit chi-squared below 0.01.
