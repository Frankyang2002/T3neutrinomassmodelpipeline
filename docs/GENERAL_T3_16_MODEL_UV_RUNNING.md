# General T3: 16-model convention-labelled one-loop UV running

## Scope

Stage 30 completed `A_+0` with a live, class-level projected RGBeta quartic
beta system and coupled DOP853 integration. Stage 31 extends that **isolated**
method to the other 15 ordinary T3 models, not into production.

- 33 newly registered quartic slots among the remaining 15 models are represented
  as complex sidecar ODE coordinates.
- 12 *reused* extra-quartic slots (six in `B_-1`, six in `C_-1`) are **not**
  double-counted: their Stage 29 boundary values replace the matching existing
  complex scalar entries in a *copy* of the fitted benchmark state.
- Four models (`A_-4`, `D_+2`, `E_+2`, `E_-4`) have zero extra quartic slots.
  Their standard UV beta functions are still exported and integrated.
- Together with Stage 30's five new `A_+0` slots, the inventory is **16 models,
  38 new slots + 12 reused slots = 50 additional quartic slots**.

## Equations and implementation

RGBeta's **native** projected class-level `BetaTerm[Quartic,1]` is used for
all quartics, not unprojected individual `BetaTerm[symbol,1]`. The complete beta
system is regenerated after registering the appropriate extension:

\[
  \frac{dX}{d\ln\mu} = \frac{\beta_X^{(1)}}{16\pi^2},\quad
  \mu:10^7\to10^5\ \mathrm{GeV}.
\]

All base and new beta expressions use the production numerical `RGBetaEvaluator`
and `BetaVector` with a temporary extra-complex sidecar, integrated by SciPy
`solve_ivp(method='DOP853',rtol=1e-8,atol=1e-11)`.

The output uses frozen fitted UV input matrices *without refitting* and Stage 29
illustrative quartic boundary candidates. The extra quartic signs and relative
magnitudes are conditional on the Stage 27 provisional triplet convention.
**The Matchete/RGBeta physical component identification is not certified.**

No threshold matching, EFT running, neutrino-mass calculation, stability survey,
or physical viability claim is included.

## Files and installation (repository root)

- `scripts/ExportGeneralT3AllModelsUVOneLoop.wl`: generalised **per-model** RGBeta
  export with a fail-closed projector and full inventory checks.
- `scripts/CoupledUVRunning.py`: single-model coupled UV integrator and
  16-model summary validator.
- `scripts/RunGeneralT3CoupledUVRunning.ps1`: sequential local driver, stops on failures
  instead of reporting false completion.
- `tests/python/test_general_t3_coupled_uv_running.py`: synthetic
  integration/inventory/provenance checks, not a native RGBeta substitute.

Run the local Stage 31 tests, then the entire 15-model batch:

```powershell
python -m pytest -q tests\python\test_general_t3_coupled_uv_running.py
.\scripts\RunGeneralT3CoupledUVRunning.ps1
```

To run only one model initially (recommended to check the native exporter):

```powershell
.\scripts\RunGeneralT3CoupledUVRunning.ps1 -Models 'A_-2'
```

Re-running all models overwrites only the isolated Stage 31 output JSON files.
The Stage 30 `A_+0` numerical result must remain available at
`output/general_t3_stage22/stage30_A_p0_coupled_uv_running.json` for final summary.

For a given model such as `A_-2`, the generated files are:

- `output/general_t3_stage22/stage31_A_m2_full_uv_betas.json`
- `output/general_t3_stage22/stage31_A_m2_coupled_uv_running.json`

After all 15 have succeeded:

- `output/general_t3_stage22/stage31_16_model_uv_summary.json`

Stage 31 does not modify the production RGBeta model, UV ODE state, numerical
configurations, or the remote GitHub repository.

## Verification scope

Python tests use controlled synthetic fixtures to verify conservation of
coupling inventory and the numerical handoff contracts. The actual 15-model
RGBeta projector and numerical integrations can only be verified by running
`wolframscript` and the Python pipeline locally with the user's native package.
Failure, including a genuinely singular quartic projection, must not be
bypassed. Retain the raw native failure output if any model fails.
