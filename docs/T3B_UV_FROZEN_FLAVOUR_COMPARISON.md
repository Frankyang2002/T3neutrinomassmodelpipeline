# T3-B-m1: UV-frozen fitted flavour comparison

## What this runs

`Numerical/analysis/CompareFermionStartT3BFlavour.py` is an **independent, read-only** diagnostic of two existing, individually Sobol-optimised configurations:

- `output/studies/fermion_start_lambda_comparison/5e-9/fitted/T3_dS1_2_dS2_2_dF_1_alpha_m1.json`
- `output/studies/fermion_start_lambda_comparison/0p5/fitted/T3_dS1_2_dS2_2_dF_1_alpha_m1.json`

It also requires the corresponding `fits/<model>.json` records with `status=Success` and `chi2_fit<5.86`. **It will not substitute a prepared/unfitted 0.5 seed, silently refit parameters, or guess missing matrices.** The case label `0p5` denotes `lambdaT3=0.5` for class B.

The workflow is:

1. Load the existing fitted config and model-specific symbolic JSON artifacts. Require matching model ID, case, coupling, accepted fit and explicit fermion-start marker.
2. Impose saved initial couplings directly at `mu_F` via the repository's `frozen_result(..., uv=True)`; assert `uv.nfev=0`. Run actual scalar-only EFT1 renormalisable beta functions from `mu_F` to `mu_S` and collect sampled background.
3. Evaluate fixed-order scalar-threshold `C5 = C5_hard + C5_direct` using the existing symbolic hierarchical coefficient; cross-check with `FinalC5TrajectoryEvaluator`.
4. Check the existing *unfitted* heavy-log normalization of the direct contribution; refuse unexplained discrepancies.
5. Independently call the repository's **opt-in exploratory** full-flavour two-operator `C12`–`C5` leading-log system with its common-scale matching-log boundary. It does not overwrite the fixed-order `C5`.
6. Run the same final SM+Weinberg evolution from `mu_S` down to `mu_low` for fixed-order and exploratory matrices.
7. Quantify `C12`/`C5` alignment, the correction perpendicular to fixed `C5`, Yukawa SVD misalignment, individual heavy-column `C12` contributions and the transformed neutrino mass correction `U.T @ delta(mnu) @ U` in the charged-lepton mass basis.

All matrix projections use complex Frobenius products. Takagi diagonalisation is performed using the repository's full 3-generation complex-matrix routine, **not** a real-only eigenvalue shortcut. Ten-percent Yukawa-column variation tests are explicitly *tree `C12` texture-only*, not new low-energy neutrino predictions.

## Local execution — PowerShell

Unzip the package **at the repository root**:

```powershell
Expand-Archive -LiteralPath .\T3B_UV_FROZEN_FLAVOUR_COMPARISON.zip -DestinationPath . -Force
```

Run the already-fitted `5e-9` reference only (does not require a high-scalar fit):

```powershell
python -m Numerical.analysis.CompareFermionStartT3BFlavour --cases 5e-9
Start-Process "output\studies\t3b_uv_frozen_flavour_alignment\index.html"
```

Then, when `0p5/fitted/` and `0p5/fits/` exist from the previous high-scalar Sobol workflow, compare both:

```powershell
python -m Numerical.analysis.CompareFermionStartT3BFlavour --cases 5e-9,0p5
Start-Process "output\studies\t3b_uv_frozen_flavour_alignment\index.html"
```

If your model-specific artifact locations differ, pass the same `--artifact-map PATH.json` used by the existing fermion-start study. Adjust EFT1 sample count with `--samples 121` (default 121).

If the 0.5 point has **not** been fitted, first run the previously provided high-scalar workflow (do **not** use an unfitted rescaled seed as a completed benchmark):

```powershell
python -m Numerical.orchestration.RunHighScalarHalfFermionStart preflight --check-artifacts
python -m Numerical.orchestration.RunHighScalarHalfFermionStart prepare
python -m Numerical.orchestration.RunHighScalarHalfFermionStart fit --resume --model T3-B-m1
```

`fit --model T3-B-m1` uses the existing local workflow; it may require its Sobol search and EFT runs and is separate from this **no-refit diagnostic**.

## Output

Files are created locally under `output/studies/t3b_uv_frozen_flavour_alignment/`:

- `index.html`, `comparison.json`, `comparison.csv`: status and side-by-side metrics.
- `5e-9/summary.json`, `5e-9/index.html`, plus diagnostic PNGs.
- `0p5/summary.json`, `0p5/index.html`, plus diagnostic PNGs **only if the fitted input exists and the evaluation passes**.
- `saved_5e-9_yukawa_only/`: bundled independently verified *algebraic* input benchmark/figures, **no RGEs and no coupled `C5`**.

### Key metrics

For each case, the complex projection is

`alpha = vdot(C5_fixed, C5_coupled-C5_fixed)/vdot(C5_fixed,C5_fixed)`;

`eps_perp = ||deltaC5-alpha*C5_fixed||_F/||C5_fixed||_F`.

The neutrino mass-basis check is

`eps_mass = ||offdiag(U_fixed.T deltaMnu U_fixed)||_F/||Mnu_fixed||_F`.

The physical charged-lepton rotation and Takagi phase conventions are applied **first**. For complex symmetric Majorana matrices the transpose in the last expression, *not* Hermitian conjugation, is intentional.

### Interpretation boundaries

- **Verified from repository source:** `RunFermionStartLambdaStudy` freezes UV with `nfev=0` while running EFT1 and Weinberg; `T3BExploratoryCoupledBenchmark` explicitly has incomplete full operator closure; `FinalC5ContributionDiagnostics` splits fixed-order matching into hard and direct.
- **Generated by the new script:** optional exploratory coupled correction, its mass eigenbasis and plots, generated locally *only after the above inputs are found*. These are fresh diagnostics, not historical output files.
- **Not established:** a validated next-to-leading-order physical neutrino prediction or a perturbativity proof for `lambdaT3=0.5`.
- The GitHub-saved `5e-9` fit is an *available* UV-frozen baseline; it should not be mislabeled as the user's otherwise unavailable historical "2.29%" run. To analyze that different run, the corresponding local config, full scalar-threshold matrices and coupled matching scheme are needed.
- Missing `0p5` fitted files appear as `MissingFittedInput` in the comparison rather than being silently substituted.

Source inspected on GitHub `main` at commit `9abe36f583e3cb2d7a3d1ce709032af7bb357a7f` (24 October 2026). No source changes were sent to GitHub.
