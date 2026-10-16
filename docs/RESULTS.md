# T3 numerical inputs, benchmark results and interpretation

A compact record of the inputs, verified observations and limitations documented in `INFO_T3_Numerical_Inputs_Sources_and_Justification.md`, `INFO_T3_NUMERICAL_RESULTS.md`, `INFO_T3_BENCHMARK_SELECTION_AND_NUFIT_FEASIBILITY.md`, `INFO_T3_MULTIPLICITY.md` and `Reports/notes/T3_Numerical_Findings.md`.

**Scope:** reported benchmark and numerical outputs, not a claim of completed phenomenological validation. Use the saved JSON/CSV and exact frozen input manifest for quantitative reproduction.

## 1. Input provenance

Use distinct categories: **literature input (L)**, **chosen benchmark (B)**, **fitted quantity (F)**, **derived quantity (D)** and **numerical convention (N)**.

| Quantity | Documented input | Category and interpretation |
|---|---|---|
| UV scale | \(10^7\) GeV | B: chosen UV boundary |
| Gauge couplings at UV | \(g_Y=0.38265,\ g_2=0.59493,\ g_3=0.7649\) | L: rounded SM-running reference values |
| Higgs VEV | \(v=246.22\) GeV | conventional input; project uses \(m_\nu=-v^2C_5/2\) |
| Thresholds | \(\mu_F=10^5\), \(\mu_S=10^3\), \(\mu_{\rm low}=10^2\) GeV | B: fermion-first hierarchy |
| Heavy-flavour masses | \(\mathrm{diag}(10^5,1.1\times10^5,1.2\times10^5)\) GeV | B: chosen masses, not measured/pole masses |
| Higgs quartic | \(\lambda_H=0.25\) | B: **not** taken from the SM-running table |
| Scalar self-quartics | \(\lambda_{S_1}=0.10,\lambda_{S_2}=0.12\) | B |
| Higgs portals | \(\lambda_{H1}=0.02,\lambda_{H2}=0.03\) | B |
| Scalar–scalar portal | \(\lambda_{12}=0.01\) | B |
| Reference \(\lambda_{T3}\) | 0.1 or 1.0 fixed grid; E factor 1/2 | B; multiplicity-only normalisation |
| New Yukawas \(y_1,y_2\) | per-model fitted matrices in `Numerical/benchmarks/fitted_reference/` | F |
| ODE tolerances | `rtol=1e-8`, `atol=1e-11` | N; not by themselves a convergence validation |

The source for rounded SM gauge and charged-fermion Yukawa parameters is Antusch, Hinze & Saad (2026), Table II at \(10^7\) GeV. That table's \(g_1\) is GUT-normalised; the project uses \(g_Y=\sqrt{3/5}\,g_1\). Diagonal benchmark quark Yukawa matrices omit a realistic CKM texture. Additional T3 couplings and threshold masses are chosen, not determined by that SM reference.

## 2. Fitting target and caveats

The documented NuFIT 6.1 normal-ordering target input file is `configs/nufit_6_1_ic24_no.json`:

| Observable | Target | Approximate 1σ |
|---|---:|---:|
| \(\Delta m_{21}^2\) (eV²) | \(7.537\times10^{-5}\) | \(9.7\times10^{-7}\) |
| \(\Delta m_{31}^2\) (eV²) | \(2.511\times10^{-3}\) | \(2.05\times10^{-5}\) |
| \(\sin^2\theta_{12}\) | 0.3088 | 0.00665 |
| \(\sin^2\theta_{13}\) | 0.02248 | 0.00057 |
| \(\sin^2\theta_{23}\) | 0.470 | 0.0155 |

The fit objective uses an approximate independent-error diagonal Gaussian form,

\[
\chi^2_{\rm approx}=\sum_{i=1}^5\left[\frac{O_i^{\rm model}-O_i^{\rm target}}{\sigma_i}\right]^2.
\]

This is not the official NuFIT multidimensional likelihood. The documented fitted reference manifest contains independently re-evaluated candidates for 16 models with very small fitted \(\chi^2\), which establishes numerical reproducibility of calibrated inputs, **not predictive validation**. Check the NuFIT dataset/variant label before citing it in the thesis.

## 3. Fixed comparison completion versus model viability

The reported October 2026 numerical-only resume completed 63 successes plus one previously completed case, with no failed numerical jobs across the 64-case fixed grid. It reused generated matching/RGE artifacts rather than rerunning Matchete. A subsequent saved-results display rebuild succeeded for all 64 model rows (384 PNGs and four HTML dashboards). Neither result implies that all points agree with oscillation data. Some large-coupling reference cases yielded neutrino masses far outside the measured scale.

`Numerical/diagnostics/ExtractFullComparisonResults.py` extracts saved diagnostics into `Reports/output/full/t3_full_comparison_observables.csv` and `.json`; it does not refit or rerun the symbolic pipeline.

## 4. Recorded numerical sensitivity findings

The saved `smallY_smallL` staged-running scan covers **12 fitted models × five heavy-mass factors = 60 successful points**; four additional models are marked `MissingArtifacts`. The full-running stage at baseline factor 1 reproduced accepted fitted observables (maximum documented chi-square 0.0037006394). These are **data coverage limitations**, not theoretical exclusions.

The cumulative stage toggles are:

| Stage | UV beta evolution | Intermediate EFT beta evolution | Final Weinberg evolution |
|---|---|---|---|
| `000` | off | off | off |
| `100` | on | off | off |
| `110` | on | on | off |
| `111` | on | on | on |

An off stage freezes beta evolution but **retains matching**. The final comparison uses the fully evolved charged-lepton basis, so these are conditional staged-running diagnostics, not scheme-independent decompositions of quantum contributions.

**Median incremental percentages at the fitted reference \(r=1\)** over the 12 covered models:

| Stage | Solar mass-squared splitting | Atmospheric mass-squared splitting |
|---|---:|---:|
| UV `000 → 100` | +19.951% | +19.892% |
| EFT1 `100 → 110` | +0.549% | +0.549% |
| Weinberg `110 → 111` | -12.336% | -12.336% |
| Total `000 → 111` | +6.268% | +6.079% |

The medians in each column are computed independently and **must not be summed**. Individual stage ratios multiply; their logarithms add. The reported additive-log consistency residual is below \(3\times10^{-14}\) percentage points.

EFT1's median is small but its *model-level* baseline solar contribution spans about -4.438% to +23.420%. The final Weinberg solar correction stays around -12.35% to -12.32% in the investigated heavy-mass-factor range, reflecting its fixed scalar-to-low-energy scale interval. This is **not** a universal correction factor.

At baseline, the median UV effect on \(\sin^2\theta_{13}\) was approximately +0.1974%, versus about +0.000214% in the final Weinberg stage; EFT1 mixing changes were at the numerical-noise level in this data set. These are contingent on the fitted models, input textures and basis prescription.

For nonbaseline heavy-mass factors the couplings were **not refitted**, so shifted points are sensitivity probes rather than alternative good-fit benchmarks. Four missing-artifact cases must not be interpreted as phenomenological failure. A complete independent row-by-row comparison to an earlier running CSV was explicitly left pending in the original numerical-findings note.

## 5. Interpretation and unresolved validations

1. Check exact NuFIT release/subset label and whether the 5-observable approximation is sufficient for the thesis claim.
2. Justify mass-threshold choices; distinguish running Lagrangian masses, matching scales and physical pole masses.
3. Verify Higgs quartic convention and higher-scale matching; \(\lambda_H=0.25\) is a benchmark assumption.
4. Document diagonal quark Yukawa/CKM approximation.
5. Do not infer vacuum stability, perturbativity across all scales or full phenomenological viability solely from small initial quartic values and successful integration.
6. Save exact benchmark config paths and SHA-256 provenance whenever generating final comparison figures.
7. In staged-RG plots distinguish the LLSS fixed-order diagnostic from the actual direct \(LLSS\to O_5\) contribution.
8. Preserve scenario-specific input configurations; a numerical resume can accidentally reuse previously resolved fit parameters and wash out nominal small/large distinctions.

## 6. Primary references and reproducibility pointers

- Antusch, Hinze & Saad, *Updated running quark and lepton parameters at various scales* (2026), https://arxiv.org/abs/2510.01312 .
- NuFIT 6.1, https://www.nu-fit.org/?q=node/309 ; target PDF: https://www.nu-fit.org/sites/default/files/v61.tbl-parameters.pdf .
- Restrepo, Zapata & Yaguna (2013), https://arxiv.org/abs/1308.3655 .
- Ma (2006), https://arxiv.org/abs/hep-ph/0601225 .
- Frozen fitted inputs and manifest: `Numerical/benchmarks/fitted_reference/benchmark_manifest.json`.
- Authoritative per-study outputs: `Reports/output/full/{staged_rg_comparison,weinberg_running_comparison,fermion_scale}/smallY_smallL/scan.json`, where available locally.
- Detailed findings source: `Reports/notes/T3_Numerical_Findings.md`.

This summary reproduces selected documented findings. For complete scan rows, use the saved numerical output, not this prose.
