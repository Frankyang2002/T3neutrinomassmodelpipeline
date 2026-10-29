# General T3: existing model-specific C5 and final SM+Weinberg running

## Scope

Connect the successful Stage 33 / 35 *numerical states* to existing model-specific one-loop matching artifacts. No new symbolic matching is derived. The authoritative code is reused:

- `Numerical/running/FinalC5TrajectoryAdapter.py`: heavy mass-basis conversion and construction of C5 inputs from the heavy threshold at 100 TeV and scalar masses and `lambdaT3` at 1 TeV;
- `Numerical/running/SMWeinbergStage.py`: evaluation of an **existing** symbolic C5 artifact;
- `Numerical/running/ScalarThresholdBoundary.py`: initial SM+Weinberg state;
- `Numerical/running/SMWeinbergEvolution.py`: one-loop running 1 TeV -> 100 GeV;
- `physics/NeutrinoMass.py`, `physics/NeutrinoTrajectory.py`, `physics/NeutrinoObservables.py`: low-energy neutrino mass and charged-lepton basis observables.

Full prior UV and intermediate RGEs are **not** rerun. SHA-256 validations tie Stage 33/34/35 outputs to the inputs. For each model the coefficient file is resolved from its *own* `(class,alpha)` path, never by recycling another model's matching artifact.

## Existing C5 sources

For twelve model/alpha values, the successful hierarchical fixed-order source is:

`output/full/hypercharge/T3_<family>_alpha_<pN_or_mN>/data/final_weinberg_coefficient.json`

The source must be `status=Success`, `scheme=fixed_order_one_loop_MSbar_hierarchical`, have a successful renormalized hard threshold and a validated manifestly symmetric direct-running contribution. `FinalWeinbergAdapter` combines the hard and direct terms.

Four extreme-hypercharge cases (`A_-4`, `B_-3`, `C_-3`, `E_-4`) currently do not have this hierarchical artifact in the inspected repository. They have the *different* legacy one-loop output:

`output/hypercharge/T3_<family>_alpha_<mN>/data/c5_coefficient.txt`

This uses the existing `evaluate_symbolic_c5` parser and *does not* reproduce the same hierarchical hard-plus-direct scheme. To include these four rough benchmarks the caller must explicitly supply `--allow-legacy`. The combined 16-model result is labelled **mixed matching sources**, and source kind is recorded per model. Do not treat it as a uniform precision comparison.

## Missing physics and convention caveats

- No re-evaluation of the full Matchete matching with the 38 new Stage 29 scalar quartic vertices. Those new terms have no certified EFT1 basis transfer and are absent below 100 TeV.
- The 12 reused quartic slots were evolved in Stage 34; they enter the scalar threshold state, but there is no new Wilson-operator matching.
- Existing matched C5 coefficients are evaluated on the Stage 33/35 RGE backgrounds. That is a conditional approximate prediction, **not** a new matched EFT calculation.
- Model-specific mass eigenvalues differ from the grouped matching scales (100 TeV and 1 TeV); this is an explicit common-scale approximation.
- The matching convention / coupling map between Matchete and RGBeta extra quartics remains uncertified.
- The `NO` mass eigenstate labelling follows the fitted-reference benchmark convention. No refit to oscillation data is performed.

## Files / project paths

- `scripts/FinalWeinbergRunning.py` — new complete runner.
- `tests/python/test_general_t3_final_weinberg.py` — new offline regression tests.
- `docs/GENERAL_T3_FINAL_C5_AND_WEINBERG_RUNNING.md` — this document.

All generated JSON outputs are placed under `output/general_t3_stage22`.

## Commands — PowerShell from repository root

Install ZIP (place it in repository root first):

```powershell
Expand-Archive -LiteralPath .\T3_General_Scalar_Stage36_ExistingC5_WeinbergRunning.zip -DestinationPath . -Force
```

Tests:

```powershell
python -m pytest -q tests\python\test_general_t3_final_weinberg.py
```

First model (hierarchical C5):

```powershell
python -m scripts.FinalWeinbergRunning --model A_+0
```

Run all 16 **with explicit opt-in for the four legacy approximations**:

```powershell
python -m scripts.FinalWeinbergRunning --all --allow-legacy
```

Regenerate/check the combined summary without rerunning:

```powershell
python -m scripts.FinalWeinbergRunning --summary
```

## Output

Individual model:

`output/general_t3_stage22/stage36_A_p0_low_energy_neutrinos.json`

All-model summary:

`output/general_t3_stage22/stage36_16_model_low_energy_neutrino_summary.json`

The per-model record carries provenance hashes, matching scheme, heavy/scalar numerical matching inputs, full complex C5 at 1 TeV and 100 GeV, low-energy neutrino mass matrix, `NO` neutrino mass eigenvalues, squared splittings and absolute PMNS entries. The combined summary differentiates the twelve hierarchical and four legacy artifacts. All files use the existing benchmark without refitting.

## Offline verification

The Stage 36 unit tests mock the external numerical services but validate all source routing, missing-artifact failures, alpha isolation, opt-in policy, and handoff contract. Existing Stage 30/31/33/34/35 Python regressions also pass in a local isolated staging tree. *The actual Stage 36 matching evaluations and numerical integration must run in the user's full Windows repository.*
