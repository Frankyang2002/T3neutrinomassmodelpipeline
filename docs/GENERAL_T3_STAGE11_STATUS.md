# Stage 11: finite mass-ratio logarithms in the experimental C5 pilot

## Source of the Stage 10 failure

In the user's Stage 10 run, **54 Python tests passed** but the first real
Matchete derivative artifact was rejected with `Unresolved matching logarithm;
hard-scale convention unknown`. The Stage 10 implementation rejected every
remaining `Log` after setting `μbar2 = MF^2`, even though terms such as
`Log[MF^2/MS1^2]` are finite **physical mass-ratio logarithms** and must
not be discarded. This was an overly restrictive validation rule. The exact
remaining logarithm(s) in the user's real artifact are not available in this
package, so the failure's precise expression is not yet established.

## Implemented change

- `Numerical/experimental/GeneralT3C5Pilot.py`: keeps finite logarithms of
  **known, positive physical scalar masses**, and substitutes the renormalisation
  squared scale `μbar2 = MF_r^2` separately for each positive fermion mass.
  It requires dimensionless log arguments, refuses unresolved unknown symbols,
  and continues to **reject UV regulator poles**. It will not arbitrarily
  zero `Log` terms. `MS1`, `MS2`, `MS` are masses in GeV; `mS1Sq`, `mS2Sq`
  are mass-squared inputs in GeV^2. The squared values are derived from the
  masses when omitted, with a consistency check if both are provided.
- `docs/general_c5_toy_benchmark.json`: includes toy scalar masses
  `MS1=1350 GeV`, `MS2=1780 GeV`. These are *not* fitted model parameters.
- `scripts/RunGeneralT3NumericalPilot.py`: writes a machine-readable
  `status=Blocked` diagnostic containing the actual response atoms if the
  real Matchete expression cannot be interpreted safely.
- `scripts/RunGeneralT3Stage10.ps1`: executes both B/C hard-C5 pilots and
  both independent RGBeta source pilots, aggregating errors at the end so
  one failure cannot hide the other three reports.
- `tests/python/test_general_t3_stage10.py`: tests mass-ratio logs, missing
  scalar masses, inconsistent masses, unsafe dimensional logs and unknown
  symbols in addition to the earlier Stage 10 regressions.

**Scope:** This remains an experimental *one-loop Matchete hard-boundary*
response and a **separate** one-step RGBeta beta-source computation. No
Matchete-to-RGBeta normalization, hierarchical final C5, complete quadratic
mass mixing or RG-improved neutrino-mass calculation is claimed.

## Run (PowerShell, from repository root)

```powershell
Expand-Archive -LiteralPath .\T3_General_Scalar_Batch11.zip -DestinationPath . -Force
python -m pytest -q tests\python\test_general_scalar_trial.py tests\python\test_scalar_quadratic_prototype.py tests\python\test_general_scalar_c5_proof.py tests\python\test_general_rgbeta_coverage.py tests\python\test_general_t3_stage6.py tests\python\test_general_t3_stage10.py
powershell -ExecutionPolicy Bypass -File .\scripts\RunGeneralT3Stage10.ps1
```

Outputs are unchanged:

```
output/general_t3_numerical_pilot/B_m1_matchete_hard_c5.json
output/general_t3_numerical_pilot/B_m1_rgbeta_source.json
output/general_t3_numerical_pilot/C_m1_matchete_hard_c5.json
output/general_t3_numerical_pilot/C_m1_rgbeta_source.json
```

The hard-C5 JSON records the chosen matching-scale convention, finite
mass-ratio log handling, and supplied scalar masses. If another symbolic
structure occurs, the Matchete hard-C5 JSON instead contains the exact
`reason` and `raw_response_atoms`. This is intentional: no extrapolated
numerical result is returned for an unsupported symbolic expression.

## Verification performed here

- 59 pytest cases passed using local isolated copies of the previously
  examined parser/evaluator as test fixtures. The production repository
  itself was not modified.
- A four-command CLI integration check succeeded with **synthetic** B/C
  response artifacts containing both a matching-scale log and a finite
  mass-ratio log, plus synthetic RGBeta source artifacts. This does not
  prove the as-yet-unseen real Matchete derivative expression will parse.
- The ZIP contains full replacement files and all earlier experimental files.
