# Workstream D — Stage 2 packed-bilinear correction

## Verified failure

The previous `AuditGeneralT3FermionCrossings` run completed all three full-UV
fermion mass crossings for 8/16 models and reported `Invalid packed RGBeta
MF/layout` for the eight models with an additional holomorphic scalar bilinear.
This failure is a Stage 2 audit state-size assumption, not evidence of an
unphysical heavy mass matrix or a source RGE failure.

The existing full-scalar runner carries an ODE vector of length:

- `layout.size` for `quadratic_kind == 'none'`;
- `layout.size + 2` otherwise, because Re(mu2) and Im(mu2) are appended **after**
  the baseline and extra-quartic blocks.

`StateVectorLayout.block('MF')` continues to cover the unchanged complex
fermion mass matrix. The Stage 2 extractor previously rejected the appended
bilinear because it required the entire ODE vector to have length `layout.size`.

## Changes

1. `Numerical/experimental/GeneralT3RunningFermionThresholds.py` now accepts an
   **explicit** `trailing_components=2` in `mass_matrix_from_vector`, validating
   the whole ODE vector and extracting `MF` using the original named block.
   It retains zero as the default, rejecting undeclared trailing coordinates.
2. `scripts/AuditGeneralT3FermionCrossings.py` supplies exactly 2 trailing
   coordinates for `self_s1`, `self_s2`, and `mixed` models, and 0 for `none`.
3. `tests/test_general_t3_running_fermion_thresholds.py` adds two tests for
   correct real/imaginary `MF` extraction with an appended complex bilinear,
   and rejection of wrong/nonfinite appended coordinates.

No modifications to UV beta functions, initial parameters, scalar potentials,
threshold equations, existing scalar data or legacy B/C benchmark runners.
Neither physical EFT matching nor per-fermion decoupling is performed.

## Apply on Windows (PowerShell)

Download the fix ZIP, then from the **existing** local worktree run:

```powershell
$repo = 'C:\thesistemp\T3neutrinomassmodelpipeline'
$zip = "$env:USERPROFILE\Downloads\T3_WorkstreamD_Stage2_PackedBilinear_Fix.zip"
Expand-Archive -LiteralPath $zip -DestinationPath $repo -Force
Set-Location $repo
python -m unittest discover -s tests -p 'test_general_t3_running_fermion_thresholds.py' -v
python -m scripts.AuditGeneralT3FermionCrossings --all --output 'output/studies/general_scalar_extended/continuous_full_scalar/continuous_v1/workstream_d_fermion_audit_v2/audit_summary.json'
```

The previous Stage 2 audit is deliberately left untouched. The default audit
command refuses to overwrite existing results; use the new `--output` path.

Expected unit tests: **13 passing** (11 original + 2 new) when applied to the
unchanged project dependency versions. Model-level UV replay of all sixteen
cases must be checked on the user's actual Windows repository; no such local
replay was performed in the isolated test harness.

## Outputs and interpretation

The corrected audit retains the original JSON schema and `candidate` flags.
The 8 previously successful models should preserve their reported roots;
the 8 previously blocked models should now pass the packed-layout check.
Additional independent checks can still flag genuine numerical failures.

The computed fermion roots continue to use the **undecoupled** full-UV RG
trajectory, and the scalar roots continue to use the original grouped-F EFT1
trajectory. Their combined order is a diagnostic, not a physically consistent
sequence of decoupled EFT evolutions. Finite matching remains missing.
