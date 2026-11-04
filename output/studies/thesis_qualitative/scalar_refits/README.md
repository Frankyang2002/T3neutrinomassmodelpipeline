# Qualitative general-T3 scalar refits

**Implementation baseline (verified from GitHub):** `main` at
`6ebe9937a72880c6f180bcd026261779800ac71f` (unchanged from supplied SHA).
The accepted `continue_general16_v1/merged_all16.json` selects 11 original and 5
continuation fits. The original fitted `lambdaT3` is 0.1 for A/B/C/D and 0.05
for E. All 16 source-specific fitted configs were inspected via the GitHub
connector. **No numerical refit in this package is claimed to have run.**

## Inputs and physics

The new script calls `RunT3GeneralZeeSensitivity.fit_selections` and `.sources`
to validate the manifest, source result/config/Yukawas and, for nonmixed
models, the source-specific `general_uv_overrides.json`, canonical
`full16_general_v1` UV/EFT1 RGBeta data, and full-general Matchete C5.
It runs `RunT3GeneralAll16FromReference.NonMixedEvaluator` for the 14 nonmixed
models and `RunT3GeneralScalarBCExploratoryFit.CandidateEvaluator` for B(-1),
C(-1). The old restricted evaluator and `RunT3GeneralScalar48Campaign.py`
are **not** used. The historical B/C tracked reference is consulted by the
existing helper **only for registry and representation inventory**, not for
physical Yukawas or scalar boundaries; the source is replaced by the accepted
merged fitted config.

Both backends are **conditional one-scale** C5 calculations; they do not
perform certified UV -> integrate out all fermions -> integrate out all scalars
-> SMEFT sequential matching. B/C additionally uses the conditional
`lambda_M = -lambda_RGBeta` map for two direct scalar quartics. Extra allowed
conjugate Yukawas use the original saved boundary values and are not fitted
independently.

## Scalar policies (never conflated)

- `--mode lambda-only` (default): change *only* UV `lambdaT3`, retaining
  physical values of independent quartics and scalar bilinears. For B/C,
  because each of seven extra quartics is parameterised by
  `quartic_fraction * abs(lambdaT3)`, compensate with
  `q_new = q_old * abs(lambda_old)/abs(lambda_new)`. Reject if `q_new > 1`.
- `--mode joint-quartics`: additionally scale **extra general UV quartics**
  with `abs(lambdaT3)` while keeping ordinary quartics and mixed scalar
  bilinears fixed. For B/C this means holding the fraction `q_old` constant.

For B(-1) and C(-1), the saved quartic fraction is 0.1 at `lambdaT3=0.1`.
Therefore `lambda-only` would require fraction 20 at `0.0005`, fraction 2
at `0.005`, and fraction 0.02 at `0.5`. The first two cases are explicitly
unsupported, not silently clamped or converted to joint scaling.

## Fit method and interpretation

At each supported benchmark `0.0005`, `0.005`, `0.5`, perform a new bounded
SciPy `least_squares` optimisation over all 18 **real** entries of `y1/y2`.
The starting vector is `y_original * sqrt(lambda_original/lambda_new)`, a
heuristic seed only. There are 5 fitted oscillation observables for 18
parameters; a small optional seed-local regulariser (default 0.001) chooses
one underdetermined branch. These Yukawa matrices are not unique.

Fit acceptance uses the source's approximate 5-observable chi-square cutoff
`chi2 <= 5.86`. A supported fit may be `PoorFit` or `Failed`. Missing inputs,
failed baseline reproduction, incompatible scalar configurations and missing
physical backends are explicitly reported, never imputed. The code also
records optimizer outcome separately from physical fit acceptance.

Successful per-point output appears at
`<run-id>/<model-tag>/lambda_<value>/result.json`, including matrices, norms,
5 predictions, neutrino masses, diagnostic C5, provenance hashes and status.
Saved `fitted_config.json` and (where applicable) `general_uv_overrides.json`
are included alongside. The run root contains `comparison.csv` and
`summary.json` with accepted-point log-log norm slopes. There are no bundled
synthetic physics results.

## PowerShell (from repository root)

```powershell
Expand-Archive -LiteralPath .\t3_qualitative_scalar_refits.zip -DestinationPath . -Force
python -X utf8 -m unittest discover -s tests/python -p 'test_t3_qualitative_scalar_refits.py' -v
python -X utf8 -m scripts.RunT3QualitativeScalarRefits --models pilot --mode lambda-only --preflight
python -X utf8 -m scripts.RunT3QualitativeScalarRefits --models pilot --mode lambda-only --run-id pilot_lambda_only_v1 --max-nfev 30
# Only after checking pilot results:
python -X utf8 -m scripts.RunT3QualitativeScalarRefits --models all --mode lambda-only --run-id all16_lambda_only_v1 --max-nfev 30
# Distinct scalar hypothesis; especially useful for unsupported B/C small lambdas:
python -X utf8 -m scripts.RunT3QualitativeScalarRefits --models B_-1,C_-1 --mode joint-quartics --run-id bc_joint_quartics_v1 --max-nfev 30
```

Use a **new** `--run-id` if an output folder already exists. The script never
overwrites an existing run directory. Nonzero exit code 2 means an input/model
was blocked or at least one physically supported benchmark was not accepted.
Explicitly unsupported B/C lambda-only points are recorded, not counted as
numerical fit failures.

**Validation available in the package:** Python syntax and 8 isolated unit tests
passed locally. This does **not** establish that actual A(-4)/B(-1) UV/EFT1/
Matchete pilot runs succeeded; those need the user's complete local checkout,
canonical generated artifacts, and existing physics dependencies.
