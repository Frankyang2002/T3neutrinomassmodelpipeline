# Continuous T3 scalar running → model-specific C5 → 100 GeV neutrino results

## Installation

Extract the ZIP into the root of `T3neutrinomassmodelpipeline`, preserving directories. **Do not copy the ZIP's documentation to `README.md`.** No old physics sources, fitted reference files or historical output directories are overwritten.

New files:

- `scripts/RunT3FullScalarNeutrino.py`
- `tests/test_t3_full_scalar_neutrino.py`
- `docs/T3_FULL_SCALAR_NEUTRINO_BATCH.md`

This requires the previous successful 16-model run from `scripts/RunT3FullScalarContinuous.py` with `--run-id continuous_v1` (and its per-model files in `output/studies/general_scalar_extended/continuous_full_scalar/continuous_v1/`).

## Run (PowerShell at project root)

```powershell
python -X utf8 -m pytest -q tests/test_t3_full_scalar_neutrino.py
python -m scripts.RunT3FullScalarNeutrino --all --run-id continuous_v1
```

To resume interrupted work **without rerunning already computed models**, use:

```powershell
python -m scripts.RunT3FullScalarNeutrino --all --run-id continuous_v1 --resume
```

If changing the script, matching source or upstream trajectory, choose a new `--output-id`, e.g.:

```powershell
python -m scripts.RunT3FullScalarNeutrino --all --run-id continuous_v1 --output-id neutrino_v2
```

Outputs:

```
output/studies/general_scalar_extended/continuous_full_scalar/continuous_v1/neutrino_v1/
    summary.json
    summary.csv
    A_p0_neutrino.json
    ... 16 model JSON files where successful ...
```

The script checkpoints after each model. Failures are recorded individually and the command returns exit code 2 if any model is blocked. Saved results are only reusable if the source trajectory, endpoint, matching source and runner code hashes have not changed.

## Numerical calculation actually performed

- **All 16 cases:** Read continuous coupled RGBeta UV → 100 TeV fermion threshold → 1 TeV scalar EFT1 background. Reconstruct the *evolved* UV Yukawa and heavy mass matrices, as well as the actual EFT1 scalar endpoint. This means the new stage no longer uses frozen reference Yukawas for the B/C calculation.
- **B(-1), C(-1):** Diagonalise the evolved scalar mass matrix with a nonzero complex off-diagonal bilinear. Evaluate the project's saved mass-eigenstate **direct one-scale** Matchete ordered coefficient `A5`, using actual evolved Yukawas and heavy-fermion masses. Symmetrise `C5=A5+A5^T`. Include two direct additional quartics using the existing limited source-level `lambda_Matchete=-lambda_RGBeta` identification. Other quartics affect the coupled RG evolution, but no generic 17-coupling physical normalization certificate is claimed.
- **Eight cases without an additional quadratic:** Evaluate the repository's existing model-specific hierarchical hard+running kernel (or designated four legacy common-scale kernels) on the evolved threshold states. Additional quartics enter their RG running, but are excluded from the selected direct single-quartic Weinberg topology by the established diagonal-propagator field-flow selection.
- **Six A/D/E self-bilinear cases:** **Nonzero self-bilinear matching is not implemented.** The numerical output is the original no-mixing `C5` formula on the evolved scalar-running background **solely as a restricted control**, not the general-scalar neutrino mass. Each output marks `additional_quadratic_c5_omitted=true`. The scalar spectrum saved upstream *does* include the self-bilinear; that spectrum is not spuriously inserted into a diagonal propagator formula.
- **All successful C5 evaluations:** Evolve the resulting `C5` and SM couplings from 1 TeV to 100 GeV using the repository's `Numerical.running.SMWeinbergEvolution`, convert `mnu=-v² C5/2`, diagonalize the charged-lepton basis, and output three neutrino singular masses and five normal-ordering observables without refitting Yukawas.

### What this does **not** do

It does **not** calculate the missing A/D/E scalar mass-basis one-loop matching; it does **not** construct genuine fermion-first `C12` and `C5` Wilson operators from the entire general scalar Lagrangian or a full two-threshold renormalized matching; and the unconstrained general-scalar Matchete/RGBeta CG and coupling conventions remain uncertified. For B/C, interpreting the original **simultaneous** all-heavy one-scale coefficient at the scalar endpoint is conditional. The 100 GeV RG evolution is real, but cannot cure missing higher-threshold matching. `physical_full_general_scalar_matching_complete` is always `false`.

**Conclusion:** This is an actual model-specific numerical continuation, not another inventory audit. It delivers candidate low-energy neutrino quantities for all supported kernel types, with the six self-mixed models visibly classified as restricted controls. It must not be presented as 16 fully matched general-scalar neutrino predictions.

## Validation provenance

- Source interfaces inspected from GitHub `Frankyang2002/T3neutrinomassmodelpipeline` at tree SHA `357e82cf2bacda023cc933e26c882895fd60e918` (29 Oct 2026).
- Locally, this batch's new Python files compile. Nine isolated contract tests pass against explicitly stubbed prior-stage module interfaces. These are **not** equivalent to running the full new numerical pipeline in your checkout.
