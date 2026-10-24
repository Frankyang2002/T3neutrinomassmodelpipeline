# T3-B-m1 flavour alignment: read-only saved-matrix diagnostic

Checked against GitHub `main` tree `9abe36f583e3cb2d7a3d1ce709032af7bb357a7f` (24 October 2026). This is a **new independent analysis script**, not a replacement for existing physics code. It reads JSON already produced by the repository and **never executes UV running, intermediate running, Weinberg running, matching, fitting, or symbolic programs**.

## Install / execute (PowerShell; project root)

Review any local changes before extracting. Extract `T3B_FLAVOUR_ALIGNMENT.zip` into your local repository root.

```powershell
Expand-Archive -LiteralPath .\T3B_FLAVOUR_ALIGNMENT.zip -DestinationPath . -Force
python -m py_compile Numerical\analysis\AnalyseT3BFlavourAlignment.py

# Recalculate alignment from the complete, *local* previously saved JSON:
python -m Numerical.analysis.AnalyseT3BFlavourAlignment `
  --coupled output\t3b_fitted_common_scale_threshold.json `
  --control output\t3b_fitted_fixed_order_control.json `
  --low-energy output\t3b_fitted_low_energy_comparison.json `
  --output output\studies\t3b_flavour_alignment\fitted_reference_local

Start-Process output\studies\t3b_flavour_alignment\fitted_reference_local\index.html
```

Uses `numpy` and `matplotlib` already used by the project. It outputs `summary.json`, `index.html`, `c12_c5_alignment.png`, and `mass_basis_delta.png` in the selected output directory.

The ZIP also contains a precomputed example under `output/studies/t3b_flavour_alignment/fitted_reference_saved/` based on **selected numeric fields transcribed from saved GitHub files**. The three `source_excerpts/*.json` files are explicitly excerpts, not copies of the complete original diagnostics. The included HTML/PNGs and `summary.json` were generated with the same script from these excerpts; the normal command above recomputes from the complete local originals.

## Benchmark identity / important limitation

The example uses `Numerical/benchmarks/fitted_reference/T3_dS1_2_dS2_2_dF_1_alpha_m1.json`, where `lambdaT3=0.1`, and the saved `output/t3b_fitted_common_scale_threshold.json` coupled study. Its result is `2.4925034125%` at the scalar threshold. This is **not** the later local optimised benchmark reported as `2.29%`, nor the shared comparison reported as `2.56%`. The original saved example had UV running in the upstream calculation. The diagnostic itself performs no running, but this does not mean the input data came from a frozen-UV pipeline. Do not conflate these.

For your **exact** optimised 2.29% result, pass the corresponding `--coupled`, `--control`, and `--low-energy` files produced from the same benchmark and matching scheme. It will refuse missing or inconsistent files. For a UV-frozen study, all three source files must come from that UV-frozen evaluation; changing a CLI label is not sufficient.

The coupled treatment remains exploratory and does not establish complete three-generation operator-space closure or complete higher-order matching.

## Calculations

At the same physical flavour basis and scalar-threshold scale `mu_S=1000 GeV`:

- `fixed = hard + direct`, checked against the separate saved fixed-order control.
- `delta = coupled - fixed` and `||delta||F / ||fixed||F`.
- For reference `A` and candidate `B`, best complex scalar `alpha=Tr(A^dag B)/Tr(A^dag A)`, then shape residual `||B-alpha A||F/||B||F`.
- `||delta-alpha fixed||F/||fixed||F` distinguishes aligned versus flavour-changing correction.

At `mu_low=100 GeV`:

- Uses independently saved, evolved low-energy `C5` matrices; **does not** evolve them.
- For the real-symmetric saved benchmark, the eigenvectors of fixed `C5` are valid for the CP-conserving Takagi mixing directions up to diagonal phases. It checks `max||U|-saved |PMNS||` and refuses a mismatch larger than `1e-6`.
- Computes `U.T @ (C5_coupled-C5_fixed) @ U`, its off-diagonal relative norm, and signed-eigenvalue first-order rotation estimates for diagnostics. Always compare with the saved exact mixing angle differences.
- Refuses significantly complex cases instead of silently applying a real-symmetric diagonalisation to the wrong problem. For general complex Yukawa/charged-lepton mixing, use the project's full Takagi and charged-lepton-basis routines.

## Example result (from saved fitted reference, not the new UV-frozen 0.5 study)

| Quantity | Calculated result |
| --- | ---: |
| `||delta C5(MS)|| / ||C5_fixed(MS)||` | 0.0249250341 |
| shape residual `C12(MS)` vs `C5_fixed(MS)` | 0.0093153660 |
| shape residual `C12(MS)` vs `C5_direct(MS)` | 0.0082999120 |
| parallel factor of `delta` relative to fixed `C5` | -0.0249178355 |
| perpendicular part `||delta_perp|| / ||C5_fixed||` | 0.0005990010 |
| mass basis `||offdiag(U.T delta C5 U)|| / ||C5_fixed||` | 0.0003092792 |
| `|PMNS|` reconstruction max deviation | 3.34e-10 |

Saved exact low-energy angle differences: `theta12 = +0.0278269 deg`, `theta13 = +0.0089040 deg`, `theta23 = -0.0031850 deg`.

The data support a dominant rescaling and a much smaller rotation-generating matrix component. This is an observed alignment in the saved exploratory example, **not** a proof of exact proportionality or a validated high-order prediction.
