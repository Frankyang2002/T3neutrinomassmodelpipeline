# Batch 11 — exploratory B(-1)/C(-1) general-scalar Yukawa refit

## Scope

**New independent numerical studies only.** `scripts/RunT3GeneralScalarBCExploratoryFit.py` reuses the verified historical fitted inputs and the **existing** Sobol / sensitivity-ranked (5/8/12 active-parameter) Yukawa optimiser. The scalar extra-coupling boundary is **fixed**, not fitted. There are 18 real Yukawa fit variables (`y1` and `y2`, three by three each).

Each **trial** (including its starting historical Yukawa point) executes:

1. Deep-copy the historically fitted model configuration and replace only its 18 real Yukawa entries in memory.
2. Set each of seven candidate general-scalar quartics to `--quartic-fraction * abs(lambdaT3)` at the appropriate boundary. Set the additional holomorphic quadratic mixing to `--mu2-fraction * min(mS1Sq,mS2Sq)` GeV². Defaults are `0.1` and `0`, respectively. These are **illustrative choices**.
3. For `small_uv_running_5e-9`, integrate candidate **UV** RGBeta from 10⁷ to 10⁵ GeV, retaining Yukawa and heavy-fermion mass running. For `small_uv_frozen_5e-9` and `large_uv_frozen_0p5`, **do not** integrate UV: apply the original parameter values directly at 10⁵ GeV.
4. Diagonalize the fermion mass matrix at 10⁵ GeV with the repository's existing Majorana Takagi / vectorlike SVD routines, rotating `y1,y2` consistently. Project to the scalar-only EFT1 and evolve the candidate RGBeta quartics and mass mixing 10⁵ → 10³ GeV.
5. Diagonalize the actual 2×2 mixed-doublet quadratic mass matrix, and evaluate the **saved** Matchete **one-scale ordered** coefficient using the RG endpoint. Insert the **conditional** source-supported mapping `lambda_M=-lambda_RGBeta` for the two additional direct-C5 quartics, without claiming certified physical field identification.
6. Reuse the existing SM+Weinberg integrator for 10³ → 10² GeV, charged-lepton basis change, neutrino mass conversion, NuFIT 6.1 (five-observable) approximate χ².

This is a *conditional one-scale* benchmark, **not** newly matched general-potential sequential EFT1 C5, physically certified cross-package coupling normalisation, or an independent neutrino prediction. No production numerical module, Matchete source, or RGE source is modified.

## Output segregation

Run IDs are mandatory and must be new. Each run writes only to:

`output/studies/general_scalar_extended/<case>/<run-id>/`

- `B_m1_exploratory_fit.json` or `C_m1_exploratory_fit.json`: numerical fit diagnostics and source hashes.
- `B_m1_exploratory_fitted_config.json` or `C_m1_exploratory_fitted_config.json`: complete **new** fitted input copy with updated Yukawas and explicit exploratory provenance. Inherited old study metadata refers only to historical provenance.
- `summary.json`, `summary.md`: run results and scope limitations.

For project completeness there is a single full new script, a new test module and this document. Previous `output/studies/lambda_5e-9_refit/`, `output/studies/fermion_start_lambda_comparison/`, Stage 36, Batch 7/8, Batch 9, and Batch 10 files are never opened for writing.

## Local PowerShell

Extract the ZIP at `C:\thesistemp\T3neutrinomassmodelpipeline`. Run from that directory in the project `.venv`:

```powershell
python -m pytest -q tests\python\test_t3_general_scalar_bc_exploratory_fit.py

# Verify exact historical seeds before optimisation (two independent runs).
python -m scripts.RunT3GeneralScalarBCExploratoryFit `
  --case small_uv_frozen_5e-9 --model B_-1 --run-id b_small_seed_v1
python -m scripts.RunT3GeneralScalarBCExploratoryFit `
  --case large_uv_frozen_0p5 --model C_-1 --run-id c_large_seed_v1

# After both probes succeed: separate local refits. Default: 8 Sobol points
# then the existing staged 5/8/12-parameter refinement if necessary.
python -m scripts.RunT3GeneralScalarBCExploratoryFit `
  --case small_uv_frozen_5e-9 --model B_-1 --fit `
  --sobol-points 8 --run-id b_small_fit_v1
python -m scripts.RunT3GeneralScalarBCExploratoryFit `
  --case large_uv_frozen_0p5 --model C_-1 --fit `
  --sobol-points 8 --run-id c_large_fit_v1

Get-Content -Encoding UTF8 output\studies\general_scalar_extended\small_uv_frozen_5e-9\b_small_fit_v1\summary.md
Get-Content -Encoding UTF8 output\studies\general_scalar_extended\large_uv_frozen_0p5\c_large_fit_v1\summary.md
```

`--fit` is opt-in; without it, the command evaluates the historical Yukawa seed in the *new* scalar scenario, still performing RG and low-energy evolution. Increase `--sobol-points` to 32 or 64 for a broader, more expensive search after the pilot. `--bound` (default `3.0`) is an absolute per-entry Yukawa bound; the B/C small-scalar historical Yukawa entries are O(1), so do **not** use the high-scalar fitted-study default `0.01` for the small-scalar cases.

Optional explicit scalar-mixing input: `--mu2-fraction 0.01` **and a new run ID**. The one-scale matrix handles its nonzero eigenbasis but does not certify a physical matching result.

## Validation status

Offline tests exercise the immutable output path, parameter shape, copied trial state, original seeded inputs, stage-interface calls, and a complex mixed-mass one-scale coefficient under the negative signed map. They do **not** substitute for full numerical integration against the user's local UV/EFT1 RGBeta beta files, source coefficient, or Matchete output. The first local two-seed probes are therefore a required numerical gate before interpreting any refit.
