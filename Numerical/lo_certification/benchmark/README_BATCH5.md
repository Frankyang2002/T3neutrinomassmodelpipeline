# Batch 5 — A_-4 evolved-input direct UV loop versus conditional EFT1+joint scalar matching

**Scope:** numerical reconciliation, not a physical one-loop UV/EFT matching certificate.

This batch uses the already-completed **Batch 4** dense fitted-input results. It does not rerun the RGEs or alter any existing releases. It **reuses** unchanged repository functions:

- `scripts/BuildT3DirectOneLoopInputs.py::construct` — native SU(2) CG, charge/fermion-mass, Higgs-background and source Weyl-vertex translation, still **provisional**;
- `Numerical/lo_certification/direct_uv/native_uv_loop_reference.py::independent_reference` — full spectral, finite, zero-external-momentum UV loop;
- `Numerical/direct_one_loop/core.py::evaluate` — independently implemented production UV loop;
- `Numerical/matching/GeneralT3FlavourC5.py::prepare_heavy_flavour` — heavy fermion mass-basis checks;
- `Numerical/lo_certification/benchmark/common_input_batch4.py` — validated, dense, 20-coupling frozen input and common-input EFT summary.

## Output files: new only, no source replacement

Repository base: `C:\thesistemp\T3neutrinomassmodelpipeline`

Files in the ZIP, each **complete**, located at the paths shown relative to that root:

```
Numerical/lo_certification/benchmark/common_point_batch5.py
Numerical/lo_certification/benchmark/README_BATCH5.md
scripts/RunLOCommonPointBatch5.py
tests/test_lo_common_point_batch5.py
```

Run from PowerShell:

```powershell
$repo = 'C:\thesistemp\T3neutrinomassmodelpipeline'
$zip = "$env:USERPROFILE\Downloads\T3_Batch5_Common_Point_UV_EFT.zip"
if (-not (Test-Path $zip)) { throw "ZIP not downloaded: $zip" }
Expand-Archive -LiteralPath $zip -DestinationPath $repo -Force
Set-Location $repo
$env:PYTHONPATH = (Get-Location).Path
python -X utf8 -m unittest discover -s tests -p 'test_lo_common_point_batch5.py' -v
python -X utf8 scripts\RunLOCommonPointBatch5.py --check-only --batch4-run-id source01
python -X utf8 scripts\RunLOCommonPointBatch5.py --batch4-run-id source01 --run-id source01
```

If output exists, use `--run-id source02` (never overwrite an old report). `--batch4-run-id` must refer to an actually completed Batch 4 folder.

Successful output location:

```
output/general_observable/lo_certification_benchmark/batch5_uv_eft_source01/
    summary.json
    report.md
    original_uv_recomputed_reference.json
    UV_direct_input.json
    UV_direct_independent.json
    UV_direct_production.json
    F_direct_input.json
    F_direct_independent.json
    F_direct_production.json
    S_direct_input.json
    S_direct_independent.json
    S_direct_production.json
```

The runner evaluates in memory and publishes atomically into a new directory. Original fitted benchmarks, RGBeta exports, Wilson release, direct UV reference, and Batch 4 files are **read-only**. The GitHub repository is not modified.

## Exact numerical stages

Common fitted benchmark at 10^7 GeV; heavy fermions at 10^5 GeV; EFT1 scalar endpoint at 10^3 GeV.

1. Validate all 20 original UV couplings against the `fitted_running.json` initial state, source SHA256, RG success, `F`/`S` scales and the dense-flavour `LLSS` factor in Batch 4. Independently validate that the Batch 4 scalar response reproduces its `C5` and its `mν = v² C5/2` matrix, with `v=246.22` GeV.
2. Load the original provisional direct-UV component vertices and the native CG component, generator and covariance audits. Require original frozen SHA256 and source-audit alignment. **Reconstruct the original UV inputs using the existing constructor**, then require agreement with the published independent reference matrix, and agreement of independent spectral UV and production UV one-loop calculations. If this baseline fails, do **not** accept new outputs.
3. Recompute the full exact **direct UV** one-loop Majorana mass matrix from three explicitly different parameter assignments, all with the Batch 4 fermion-sector `y1`, `y2`, `MF` evolved to `μF`:

| Stage | Scalar quadratic and Higgs couplings | Interpretation |
|---|---|---|
| `UV` | Original frozen values at `μUV=10^7 GeV` | Heavy F evolved, frozen scalar sector: diagnostic mixed-scale input |
| `F` | UV-running scalar values at `μF=10^5 GeV` | Numerically common running scale for available couplings; unimplemented extra scalar operators projected out |
| `S` | EFT1 scalar values at `μS=10^3 GeV` | Heavy F at `μF`, scalars at `μS`: deliberately **hybrid** matching-scale sensitivity check |

4. For every stage, check finite UV pole cancellation, positivity and hierarchy, symmetric full matrix, and agreement between both established UV evaluators. Compare the complete complex 3×3 matrix against the existing **conditional** Batch 4 EFT matrix (not only singular values), fitting just one descriptive complex multiplier:

\[
  n_\star=\frac{\operatorname{Tr}[M_{\rm EFT}^\dagger M_{\rm UV}]}{\operatorname{Tr}[M_{\rm EFT}^\dagger M_{\rm EFT}]},
  \qquad
  r=\frac{\|M_{\rm UV}-n_\star M_{\rm EFT}\|_F}{\|M_{\rm UV}\|_F}.
\]

`n_star` is **not** identified with physical `N_GammaCC`. The EFT matrix is per unit unknown `N_GammaCC`; the direct-UV normalisation is also source-uncertified. The two calculations do **not** contain the same perturbative terms.

## Why the F/S direct-UV scenarios are *projected*

The existing direct-UV input constructor explicitly rejects nonzero adjoint Higgs portals and various extra scalar operators. RGEs can generate some of those couplings despite them being zero in the frozen UV benchmark. Batch 5 does **not** invent their source-to-component mass/interaction map. Where the requested scalar stage has such couplings, Batch 5 records the **actual evolved numbers** in `unsupported_evolved_scalar_couplings_dropped` and projects those couplings out **only for that UV diagnostic scenario**. These are not complete evolved UV Lagrangians. `lambdaS2Adj` is a scalar self quartic, not a Higgs portal; every omitted source coupling is listed by name. No omission is made to the existing Batch 4 EFT calculation.

This release does **not** establish equality of running full UV and EFT at a common threshold. A true scheme- and order-matched comparison requires the native scalar/Higgs Hessian dictionary, the signed `GammaCC`→master Weyl conversion and correct hard-one-loop fermion plus EFT soft/renormalisation matching. Results here only locate numerical scale and flavour sensitivity, respecting the required sequence:

```
UV -> integrate out ALL fermions together -> EFT1 with BOTH scalars
   -> integrate out BOTH scalars together -> SMEFT.
```

## Validation boundary

On the authoring environment, **19 tests passed and one full-checkout integration test skipped** (no copy of the user's Windows project/Batch 4 result). Neither the new full UV calculation nor the new numerical residuals can be claimed as executed against the user directory until the PowerShell commands above complete there. All physical certification flags stay false.
