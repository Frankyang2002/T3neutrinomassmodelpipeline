# T3 Batch 12 — end-to-end *conditional* numerical neutrino output and full-UV fit

Models **A(−4), D(−2)**. Installs as additions to the local project. Does not write to GitHub or overwrite Batches 1–11.

## What is implemented

- Read source-verified **Batch 11** finite native LLHH hard/soft *sector*, check its status and source SHA-256, and independently reconstruct its 3×3 coefficient directly from the archived one-loop operator ledger plus the original Matchete 10/68-term scalar-soft expressions. The compiled SymPy regular/pole functions evaluate all nine complex flavour entries without repeated symbolic substitutions.
- Verify native pole cancellation and numerical agreement with the saved Batch 11 result to relative `2e-7`, so neither Batch 9's scalar bubble nor Batch 10's extra reference tensors are appended a second time.
- Apply a **user-specified signed factor** `--signed-factor 1` to translate the *source-convention sector* to a **hypothetical** SMEFT coefficient. **1 is only an explicit numerical convention assumption**; it is not a derived GammaCC/Weyl factor.
- Run the repository's existing one-loop **SM + Weinberg** RG from 1000 to **100 GeV**, calculate `mnu=-v² C5/2`, Takagi masses, charged-lepton-corrected PMNS angles, and five NuFIT 6.1 residuals (uncorrelated Gaussian approximation). Default is **baseline only**: no UV rerunning or fit.
- With `--fit`, fit **all 18 real UV Yukawa entries** in the two 3×3 matrices defined at **10⁷ GeV**, not 18 arbitrary scalar-threshold effective parameters. Each trial calls the repository's *actual* full-inventory RGBeta **UV** and **EFT1** ODE solvers, computes the F mass and Yukawa endpoint at 10⁵ GeV, and recomputes both scalar masses and quartic at 10³ GeV. The source hard and scalar contributions are recalculated from these values; then the SM+Weinberg RGE and observable residuals are evaluated. Source coefficients are never fitted independently of the Yukawa matrices.
- The UV fit requires an unchanged canonical symbolic beta export and zero unprovided renormalizable F matching shifts. Trials outside the supported positive, approximately diagonal heavy-fermion basis, or with nonzero D `y2Conj`, are **rejected**, not approximated as valid. A newly integrated initial UV point must reproduce the saved Batch 11 benchmark within relative `5e-3` before fitting.

## Important limitation

This is an **end-to-end numerical demonstration for the *included source sector***, **not** a completed physical T3 neutrino mass prediction. Still unresolved: complete finite fermion matching and subtraction of other operator channels; complete EFT1 dimension-5 anomalous dimensions and operator closure (the Batch 11 pole-derived log is retained, not duplicated); native Matchete/RGBeta scalar invariant and Wilson basis conversion for D; and the signed GammaCC-to-left-Weyl spinor normalization. `--signed-factor 1` does not remove that uncertainty. Fit chi-squared values must not be quoted as physical NuFIT fits of the complete T3 model.

The calculation always preserves **UV → integrate out all three F at one μF → EFT1(S1,S2) → integrate out both scalars at one μS → SMEFT**, never reverses or separates these thresholds. For D the split neutral singlets remain both heavy fields at the common scalar threshold.

## Files to install

| Project-relative path | Role |
|---|---|
| `Numerical/parallel_t3/complete/two_model_conditional_uvfit_batch12.py` | Cached symbolic Laurent kernel, pole-closed three-flavour sector, complete UV/EFT1 RG trial runner, SMEFT observables and optimizer |
| `scripts/RunT3ConditionalUVFitBatch12.py` | Preflight, execution, separate output owners and per-model results |
| `tests/test_t3_conditional_uvfit_batch12.py` | Fast source-convention and vectorisation checks; real Matchete source comparisons when archived outputs exist |
| `README_T3_CONDITIONAL_UV_FIT_BATCH12.md` | This note |

Output (default): `output/parallel_t3/thesis_conditional_uvfit_batch12_v1/A_-4/result.json`, similarly `D_-2/result.json`, plus `summary.json`.

## PowerShell

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline

Expand-Archive `
  -LiteralPath "$HOME\Downloads\t3_conditional_uvfit_batch12.zip" `
  -DestinationPath . -Force

.\.venv\Scripts\python.exe -m pytest -q `
  tests\test_t3_conditional_uvfit_batch12.py

.\.venv\Scripts\python.exe -X utf8 `
  -m scripts.RunT3ConditionalUVFitBatch12 --all --preflight

# Baseline: C5(100 GeV), masses, PMNS and chi2 for BOTH models
.\.venv\Scripts\python.exe -X utf8 `
  -m scripts.RunT3ConditionalUVFitBatch12 --all --signed-factor 1 --debug

# Optional: independent output folder, actual UV and EFT1 RG on every trial
.\.venv\Scripts\python.exe -X utf8 `
  -m scripts.RunT3ConditionalUVFitBatch12 --all --fit `
  --signed-factor 1 --max-nfev 18 --yukawa-bound 0.005 `
  --output output/parallel_t3/thesis_conditional_uvfit_batch12_fit_v1 `
  --debug
```

You can use `--model A_-4` to fit A first, then `--model D_-2` with another output folder. For a different signed factor use a new output folder. **Never overwrite the prior report**, as it records which hypothetical normalization was chosen.

## Verification and provenance

Standalone construction environment: 6 isolated tests passed; 2 tests requiring archived native files were skipped. The user's successful Batch 11 output has been reported, but integration against their Windows repository is not executable in this construction environment. This release validates actual source status and hashes at runtime. It does **not** certify matching physics by a test pass.
