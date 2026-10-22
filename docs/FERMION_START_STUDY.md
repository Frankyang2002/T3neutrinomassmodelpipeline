# Fermion-threshold-start T3 benchmark study

**Project file:** `Numerical/orchestration/RunFermionStartLambdaStudy.py`

This is a **new** boundary-value study: the parameters taken from the 5e-9 fit (or rescaled for 5e-10) are **imposed at the fermion matching scale**. The UV RGEs are disabled; the fermion threshold is matched, the scalar-only EFT1 is run to the scalar threshold, the scalar threshold is matched, and the SM+Weinberg EFT is run to low energy. This must **not** be described as a full-UV RGE result.

Two benchmark quartics per study:
- `--case 5e-9`: A–D 5e-9; E 2.5e-9
- `--case 5e-10`: A–D 5e-10; E 2.5e-10

## PowerShell commands (project root)

Extract the ZIP at the repository root. Verify the installed module:

```powershell
python -m Numerical.orchestration.RunFermionStartLambdaStudy --help
```

Prepare both input sets (does not integrate):

```powershell
python -m Numerical.orchestration.RunFermionStartLambdaStudy prepare --case 5e-9
python -m Numerical.orchestration.RunFermionStartLambdaStudy prepare --case 5e-10
```

Pilot using A(-2) with 64 Sobol points per case:

```powershell
python -m Numerical.orchestration.RunFermionStartLambdaStudy fit --case 5e-9 --model T3-A-m2 --sobol-points 64 --resume
python -m Numerical.orchestration.RunFermionStartLambdaStudy fit --case 5e-10 --model T3-A-m2 --sobol-points 64 --resume
```

Inspect `output/studies/fermion_start_lambda_comparison/<case>/fits/T3_dS1_1_dS2_3_dF_2_alpha_m2.json` for `status` and `chi2_fit`. The shared donor benchmark requires status `Success` (chi2 <= 5.86).

Fit all 16 models, **after pilot**, with 256 Sobol samples per model:

```powershell
python -m Numerical.orchestration.RunFermionStartLambdaStudy fit --case 5e-9 --sobol-points 256
python -m Numerical.orchestration.RunFermionStartLambdaStudy fit --case 5e-10 --sobol-points 256
```

**Important:** Changing `--sobol-points` changes checkpoint identity. With `--resume`, incompatible checkpoints are refused. A 64-point pilot checkpoint for A(-2) must be backed up/moved before that model is rerun with 256. Other already successful models may be reused with `--resume` without using their Sobol checkpoint.

For completed fitted sets:

```powershell
python -m Numerical.orchestration.RunFermionStartLambdaStudy compare --case 5e-9 --mode optimised
python -m Numerical.orchestration.RunFermionStartLambdaStudy compare --case 5e-9 --mode shared
python -m Numerical.orchestration.RunFermionStartLambdaStudy compare --case 5e-10 --mode optimised
python -m Numerical.orchestration.RunFermionStartLambdaStudy compare --case 5e-10 --mode shared
python -m Numerical.orchestration.RunFermionStartLambdaStudy display --case 5e-9
python -m Numerical.orchestration.RunFermionStartLambdaStudy display --case 5e-10
Start-Process 'Reports/output/display_fermion_start/5e-9/index.html'
Start-Process 'Reports/output/display_fermion_start/5e-10/index.html'
```

## Output paths

- Inputs: `output/studies/fermion_start_lambda_comparison/<case>/prepared/`
- Fit records: `.../<case>/fits/`, `.../<case>/fitted/`
- Sobol checkpoints: `.../<case>/sobol/`
- Numerical comparisons: `.../<case>/comparison/optimised/` and `.../comparison/shared_a_m2/`
- Endpoint HTML comparisons: `Reports/output/display_fermion_start/<case>/index.html`

## Scope and limits

- Reuses prior fitted 5e-9 results and saved symbolic artifacts, not the symbolic pipeline.
- The same model's A(-2) **new fermion-start fit** supplies the shared Yukawa matrices at each quartic. Not the old UV-start A(-2) fit.
- Uses the repository's real Sobol sampler and staged 5/8/12 active-parameter least squares.
- Screens with chi2 <= 5.86, not a calibrated statistical confidence region.
- **This package creates fitted configurations, comparison CSV/JSON, and endpoint comparison HTML. It does not yet create the former energy-dependent interactive RGE-running dashboard or the old pipeline's `running_diagnostics.json`.** Do not apply `InteractiveModelComparison.py` to it without a new, fermion-start-aware trajectory serializer.
- Physics has not been executed in the author's local environment; syntax compilation is the current validation.
