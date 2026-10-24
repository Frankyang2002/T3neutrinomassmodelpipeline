# T3-B-m1 Yukawa flavour basis: read-only next diagnostic

Inspected GitHub `Frankyang2002/T3neutrinomassmodelpipeline`, `main`, tree SHA `9abe36f583e3cb2d7a3d1ce709032af7bb357a7f` (24 October 2026). Install **locally only**. The script does not change source fits or save into the production result folders.

## Intended project paths

- `Numerical/analysis/AnalyseT3BYukawaBasis.py` — complete standalone read-only diagnostic.
- `docs/T3B_YUKAWA_BASIS.md` — instructions.
- `output/studies/t3b_yukawa_basis/example_fitted_reference/` — precomputed example from an *excerpt* of the saved fitted-reference config, not full original config or a validated low-energy physics run. Rebuild locally using full config.

## Run (PowerShell, repository root)

```powershell
# Check for local file conflicts before extracting.
git status --short
Expand-Archive -LiteralPath .\T3B_YUKAWA_BASIS_NEXT.zip -DestinationPath . -Force
python -m py_compile Numerical\analysis\AnalyseT3BYukawaBasis.py

# Pilot on the repository's original T3-B-m1 fitted-reference input.
# The values are interpreted at mu_F, NOT run from mu_UV.
python -m Numerical.analysis.AnalyseT3BYukawaBasis `
    --config Numerical\benchmarks\fitted_reference\T3_dS1_2_dS2_2_dF_1_alpha_m1.json `
    --output output\studies\t3b_yukawa_basis\fitted_reference_frozen_boundary `
    --fraction 0.1
Start-Process 'output\studies\t3b_yukawa_basis\fitted_reference_frozen_boundary\index.html'
```

To repeat the same algebraic check on a **genuine previously saved UV-frozen optimised** configuration, supply its actual config path with `--config`. It must contain the Yukawa pair *imposed at* `mu_F`; do not substitute a different full-UV-running benchmark just because its couplings have a similar label.

Once the separate high-scalar fermion-start Sobol study has actually completed, the expected new T3-B-m1 fit is:

```powershell
python -m Numerical.analysis.AnalyseT3BYukawaBasis `
    --config output\studies\fermion_start_lambda_comparison\0p5\fitted\T3_dS1_2_dS2_2_dF_1_alpha_m1.json `
    --output output\studies\t3b_yukawa_basis\lambda_0p5_optimised
Start-Process 'output\studies\t3b_yukawa_basis\lambda_0p5_optimised\index.html'
```

This script **does not** perform the 0.5 fit or the coupled-C12/C5 integration. A high-scalar fit is not available just because the command is documented.

## Numerical definitions

In the heavy fermion mass basis, assuming positive diagonal `MF`:

\[
C_{12}^{(r)}=\frac{Y_{1,r}^*Y_{2,r}^{*T}+Y_{2,r}^*Y_{1,r}^{*T}}{2M_r},\qquad C_{12}=\sum_rC_{12}^{(r)}.
\]

Use the exact real/imag matrix storage format of the existing numerical configurations. The example has `MF=diag(100000,110000,120000) GeV`; the routine supports complex diagonal MF by the physical phase redefinitions needed for positive masses. It deliberately refuses non-diagonal MF rather than introducing an unverified heavy-basis rotation.

The SVD comparison is

\[
Y_1=U_1D_1V_1^\dagger,\qquad \tilde Y_2=U_1^\dagger Y_2V_1.
\]

It measures `||offdiag(tildeY2)||F/||Y2||F` and the best-complex-scalar shape residual between `Y2` and `Y1`. This SVD is only a *mathematical* comparison; rotating the lepton and Majorana indices this way generally changes the form of `Ye` and `MF`. The matched C12 calculation uses the original physically specified charged-lepton/heavy-mass basis, **not** the arbitrary SVD frame.

For a controlled test, separately multiply each `Y2` heavy column by 1.1, and calculate the associated new **tree-level C12 texture**. Its mixing *proxy* comes from the eigenvectors of `M^dagger M` in the charged-lepton mass basis, where `M` is proportional to the C12 texture. Charged-lepton rotation uses the same left SVD and sorted singular-value convention as `physics/NeutrinoTrajectory.py`. The normal-ordering index convention is used. Near-degenerate texture singular values are rejected.

The +10% **global** Y2 rescaling must leave these angles unchanged, since C12 simply rescales. Column-specific shifts change the *relative* heavy-flavour contributions and can rotate the texture. These are **not final PMNS predictions** or the effect of real coupled RGEs. They are controlled interventions showing how to break the leading flavour-scaling symmetry.

## Example from the saved fitted reference (historical λT3=0.1 input)

| Diagnostic | Numeric result |
|---|---:|
| `||Y2||F/||Y1||F` | `0.777484` |
| Y2 off-diagonal fraction in Y1-SVD basis | `0.935077` |
| Sum of three heavy-state C12 norms divided by norm of their sum | `2.627523` |
| Leading C12 texture θ12, θ13, θ23 | `(34.1872°, 8.7574°, 43.2323°)` |
| Global `Y2 *= 1.1` | Angles unchanged (roundoff) |
| Change only heavy column 1 of Y2 by +10% | `(Δθ12, Δθ13, Δθ23) = (+5.0697°, +0.3961°, +0.6241°)` |
| Change only heavy column 2 | `(-9.1780°, -0.3577°, -1.6986°)` |
| Change only heavy column 3 | `(+3.8845°, -0.1074°, +1.0017°)` |

These numbers derive from the fitted-reference **input matrices**, algebraically treated as a fermion-threshold boundary. The example does not reproduce the original historical UV-running trajectory; the angles are C12-only texture diagnostics and must **not** be equated to the previously calculated 2.29%, 2.56%, or 2.4925% coupled low-energy observables.

## Output files

The selected `--output` directory contains `summary.json`, `index.html`, `yukawa_svd.png`, `heavy_c12_terms.png`, and `column_texture_sensitivity.png`. Nothing is modified elsewhere. Verify `status=ReadOnlyT3BYukawaTextureComplete` and the input path/hash in `summary.json`.

## Interpretation

1. A small second Yukawa is **not** supported by this fitted-reference pair: its Frobenius norm is 78% of the first, and most of it remains off diagonal in the Y1-SVD frame.
2. There is substantial cancellation among different heavy-state contributions to C12; suppressing or perturbing individual columns can change its flavour eigenvectors.
3. Global Yukawa rescalings that preserve the relative bilinear structure are invisible to the mixing angles at this leading texture level. This provides a natural reason why the repository's T1/T2/T3 sensitivity directions show weak angular response, without claiming the full EFT Wilson correction is exactly aligned.

For the **exact** local 2.29% UV-frozen coupled benchmark, also rerun the earlier `AnalyseT3BFlavourAlignment.py` using the three saved *matching* `--coupled`, `--control`, and `--low-energy` JSON files. Those precise files were not present among the provided attachments or publicly saved GitHub outputs, so no exact 2.29% alignment numbers have been invented.
