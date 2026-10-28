# T3 numerical results: accepted benchmarks, comparisons and limits

**Authoritative results-level summary for the honours project.** Read with [PHYSICS.md](PHYSICS.md), [SCALAR_SECTOR.md](SCALAR_SECTOR.md) and [NUMERICAL.md](NUMERICAL.md). This document separates **completed numerical evaluations**, **calibrated neutrino fits**, **conditional exploratory scalar models**, and **physical conclusions not yet established**.

**Evidence inspected:** read-only repository commit `f21d15d2d38f4f313a0f2214bc52f90e1469dce2`, particularly `output/studies/reduced_scalar_final_reports/final_reduced_v2/{summary.md,manifest.json}`, the frozen Sobol source summary, original `docs/RESULTS.md`, the reduced reporting docs, and the frozen general-scalar 48-case CSV/report provenance. Figures and row-wise numerical results should always be cited against their **saved data and source hashes**, not this narrative alone. This documentation update does not recompute, modify or certify any physics output.

## 1. Summary of what the project demonstrates

| Dataset / question | Recorded result | What it demonstrates |
|---|---|---|
| Original 16 × 4 fixed benchmark grid | **64/64 numerical cases completed** (one earlier completed and 63 successful resumption cases) | The configured numerical jobs ran; **not** 64 accepted neutrino fits |
| Primary zero-additional-scalar frozen study | **32/32 individually refitted cases accepted**, over 16 models and two lambda scenarios | The existing fixed-order T3 C5 machinery can be calibrated to the five NuFIT targets within the declared benchmark, not an independent prediction |
| Transfer of fixed A(−2) Yukawas to other representations | **10/32 shared-reference evaluations accepted**, five in each scenario | Model-dependent Yukawa adjustment generally matters; failure of the shared input does **not** exclude a representation |
| Historical running-UV scalar boundary | **0/16** unchanged historical seeds satisfy the retained 0.5–5 TeV scalar-window screen | The **chosen high-scale boundary** does not remain inside the assumed benchmark hierarchy; not a theorem against all UV boundaries |
| General-scalar 16 × 3 frozen 48-case atlas | Frozen campaign statuses **20 conditional accepted / 6 constrained poor-fit / 22 no mass-safe fit**; **two** accepted running-UV B/C cases used a relaxed topology boundary | Exploratory general-quartic numerics under incomplete physical matching, **not** 20 experimentally viable complete models |
| Original `smallY_smallL` staged-running scan | **12 models × 5 factors = 60 successful points**; four models had `MissingArtifacts` | Conditional running sensitivity within the available saved source coverage |

The primary honours result is the **restricted leading-one-loop, zero-additional-scalar-boundary comparison**. Neither the older general-scalar 48 cases nor the nominal 64-grid completions should be silently aggregated into its 32 individual acceptances.

## 2. Benchmark inputs and classification

Maintain the distinction between literature-motivated inputs (L), chosen benchmark parameters (B), fitted quantities (F), derived observables (D) and numerical conventions (N):

| Quantity | Recorded value | Role |
|---|---|---|
| Gauge couplings at `10^7 GeV` | `gY=0.38265`, `g2=0.59493`, `g3=0.7649` | L: rounded SM-running reference values, not a new project fit |
| Higgs VEV | `v=246.22 GeV` | N: used with `m_nu=-v^2 C5/2` |
| Standard scales | `muUV=1e7`, `muF=1e5`, `muS=1e3`, `muLow=1e2 GeV` | B: selected fermion-first threshold hierarchy |
| Original illustrative heavy masses | `diag(1e5,1.1e5,1.2e5) GeV` | B: not measured pole masses or every final-fit mass spectrum |
| Scalar ordinary inputs (original reference) | `lambdaH=0.25`, `lambdaS1=0.10`, `lambdaS2=0.12`, `lambdaH1=0.02`, `lambdaH2=0.03`, `lambda12=0.01` | B: explicit model choices, **not** inferred from the SM reference table |
| Original fixed grid | `Y=0.005/0.5`, `lambdaT3(ref)=0.1/1.0` | B: nominal control combinations |
| Reduced frozen scenarios | `lambdaT3=5e-9` and `0.5` at `muF` (class E half each) | B: fixed input quartic boundary |
| Two 3×3 real Yukawa matrices | 18 real fitted entries per model/scenario | F: numerical optimiser result, not independent predicted mixing |
| Saved ODE control example | `rtol=1e-8`, `atol=1e-11` | N: precision parameters, not proven universal convergence |

The older gauge/charged-fermion starting values cite Antusch, Hinze & Saad's 2026 running-parameter table at `10^7 GeV`. Their `g1` is GUT-normalised while the project tabulates `gY=sqrt(3/5)*g1`. The original diagonal quark Yukawa background is not a realistic CKM fit. The Higgs quartic `0.25` is a separate benchmark choice. Do not generalise these selected parameters to every final-fit scalar state: the **per-run fitted configuration and manifest** define the actual numerical input.

### Five fitted observables

The file `configs/nufit_6_1_ic24_no.json` supplies the normal-ordering NuFIT 6.1 reference used in the saved study:

| Observable | Target | Approximate 1σ |
|---|---:|---:|
| `Delta m21^2` [eV²] | `7.537e-5` | `9.7e-7` |
| `Delta m31^2` [eV²] | `2.511e-3` | `2.05e-5` |
| `sin²(theta12)` | `0.3088` | `0.00665` |
| `sin²(theta13)` | `0.02248` | `0.00057` |
| `sin²(theta23)` | `0.470` | `0.0155` |

The saved objective uses the sum of five squared individually normalised residuals and accepts `chi2_approx <= 5.86`. **That is a project screening cut, not an official NuFIT likelihood or confidence level.** The fit has 18 real Yukawa search parameters for five observables, so neither a tiny chi-square nor a successful optimization implies a unique physical solution, a prediction of those fitted observables, or a phenomenological exclusion of other cases.

## 3. Final reduced-scalar frozen comparison — primary thesis result

**Source files:**

```text
output/studies/reduced_scalar_sobol/reduced_frozen_overnight_v2/comparison.json
output/studies/reduced_scalar_sobol/reduced_frozen_overnight_v2/comparison.csv
output/studies/reduced_scalar_final_reports/final_reduced_v2/summary.md
output/studies/reduced_scalar_final_reports/final_reduced_v2/frozen_comparison.csv
output/studies/reduced_scalar_final_reports/final_reduced_v2/oscillation_observable_pulls.csv
output/studies/reduced_scalar_final_reports/final_reduced_v2/manifest.json
```

The two frozen-UV scenarios impose the reduced scalar boundary at `muF = 1e5 GeV`, so **no UV running is performed in this primary frozen comparison**. One-loop candidate scalar-only EFT1 evolution is evaluated as needed. Both `S1` and `S2`, their ordinary self/portal interactions, scalar masses, and `lambdaT3` remain in the model; the *additional* gauge-allowed quartics/bilinears are zero at the input boundary. The scalar reduction is explained in [SCALAR_SECTOR.md](SCALAR_SECTOR.md).

| Frozen input scenario | Individual accepted / tested | Shared A(−2) accepted / evaluated |
|---|---:|---:|
| `small_uv_frozen_5e-9` | **16/16** | **5/16** |
| `large_uv_frozen_0p5` | **16/16** | **5/16** |
| **Both** | **32/32** | **10/32** |

For direct transparency, the **saved `chi2_approx` results** are reproduced below. A(-2) is the newly optimised **scenario-matched** reference, not the original historic reference and not a refit in the target models. All other model inputs and group factors are those of the respective target case.

| Model | Small: individual | Small: shared A(−2) | Large: individual | Large: shared A(−2) |
|---|---:|---:|---:|---:|
| `A_+0` | `1.121489e-11` | `0.00817089` | `0.005162948` | `0.007585785` |
| `A_-2` | `5.301052e-17` | `5.301052e-17` | `3.539857e-11` | `3.539857e-11` |
| `A_-4` | `2.081946e-16` | `11527.68` | `2.860277e-14` | `11464.09` |
| `B_+1` | `1.198898e-21` | `1.473904` | `3.446356e-11` | `2.187401` |
| `B_-1` | `8.937932e-18` | `37.39544` | `1.067709e-08` | `35.45188` |
| `B_-3` | `8.297949e-15` | `12177.75` | `1.482725e-09` | `12124.08` |
| `C_+1` | `6.602229e-20` | `9273.324` | `3.484557e-20` | `9256.323` |
| `C_-1` | `2.286049e-10` | `9749.82` | `2.453295e-18` | `9739.222` |
| `C_-3` | `8.040784e-11` | `17818.81` | `8.57308e-15` | `17797.39` |
| `D_+0` | `2.575589e-12` | `0.0002023746` | `1.046566e-08` | `0.0002090227` |
| `D_+2` | `4.125362e-12` | `232.8243` | `8.591204e-11` | `238.6801` |
| `D_-2` | `1.203500e-10` | `0.01288426` | `3.455100e-11` | `0.01225028` |
| `E_+0` | `1.368353e-09` | `3922.104` | `4.127976e-11` | `4106.154` |
| `E_+2` | `2.435921e-17` | `2932.508` | `5.586120e-11` | `3106.756` |
| `E_-2` | `1.185518e-09` | `3915.549` | `3.026984e-15` | `4099.545` |
| `E_-4` | `8.270065e-12` | `15394.62` | `8.651976e-11` | `15437.31` |

**Physics interpretation:** The tested input Yukawa matrices can be individually adjusted to reproduce the selected five neutrino oscillation inputs for every representation in these two benchmark scenarios. A shared Yukawa structure has much weaker transferability across representations. This supports the usefulness of model-specific numerical fitting and the representation dependence of the amplitudes, **not** independent prediction, a 32-model exclusion test or universal full-general-T3 validity.

The scalar mass screen in the frozen numerical study was **0.5–5 TeV**, and the heavy-fermion hierarchy screen was `rhoF <= 0.5`; neither bound was relaxed to produce the final accepted cases. The final C5 sources are **12 hierarchical hard+direct** and **four legacy one-scale** schemes. They are documented as **nonuniform**: the comparison is meaningful only within this stated numerical approximation. Radiatively generated candidate B/C additional quartics are omitted from the fixed-order C5; generic nonzero extra quartics and scalar mixing are **not** solved.

## 4. UV-running scalar-mass diagnostic — not a failed model fit

**Source:** `final_reduced_v2/running_uv_formal_endpoints.csv` and the original historical `reduced_48_v1` source inputs. The existing historical scalar and Yukawa boundary is held at `10^7 GeV` and evolved to the fermion threshold (`10^5 GeV`) without refitting or scalar boundary retuning. The final report also replayed **16 candidate UV renormalisable trajectories** between those two scales. `running_uv_trajectories.csv` and `running_uv_scalar_trajectories.png` are replay diagnostics, not accepted physical scalar-threshold continuations.

The saved nominal-EFT1 scalar endpoints at `muS=1e3 GeV` are **formal diagnostics**. Their scalar masses span approximately **30–169 TeV**, so **all 16 unchanged historical seeds fail** the original **0.5–5 TeV scalar window**. Several points also challenge the imposed heavy-fermion/scalar decoupling hierarchy. When scalars run above the assumed nominal matching threshold, the 1 TeV scalar-EFT endpoint **cannot** be interpreted as an independently valid physical low-energy neutrino prediction without new threshold identification and rematching.

This is a result **about the chosen high-scale input boundary**, not a proof that every UV-running T3 theory or every boundary choice is impossible. In particular, these are *not* 16 additional failed oscillation fits; the study did not perform new Yukawa fitting under physically revised mass thresholds.

## 5. Earlier fixed grid and staged-running sensitivity

The original `16×4` fixed-grid numerical-only resumption reached **64 complete numerical cases**, sourced from the stored matching/RGE artifacts. Presentation extraction yielded complete display rows and historical dashboards, but numerical completion should never be called neutrino-fit acceptance; some large-coupling nominal references give masses far outside the oscillation scale.

The `smallY_smallL` staged-running scan covered **12 models × five heavy-mass factors = 60 successful points**, with four `MissingArtifacts` model rows. The stages cumulatively switch on UV, scalar EFT1 and SM+Weinberg **beta evolution**, **without switching off matching**:

| Stage | UV beta | EFT1 beta | Final Weinberg beta |
|---|---|---|---|
| `000` | frozen | frozen | frozen |
| `100` | on | frozen | frozen |
| `110` | on | on | frozen |
| `111` | on | on | on |

Reported **median incremental percentage changes at the fitted heavy-mass reference factor `r=1`**, over the 12 covered models:

| Comparison | Solar `Delta m21²` | Atmospheric `Delta m31²` |
|---|---:|---:|
| `000 -> 100` (UV) | `+19.951%` | `+19.892%` |
| `100 -> 110` (EFT1) | `+0.549%` | `+0.549%` |
| `110 -> 111` (Weinberg) | `-12.336%` | `-12.336%` |
| `000 -> 111` (total) | `+6.268%` | `+6.079%` |

These separately taken medians **do not add**; exact stage ratios multiply and their logarithms add. The archived numerical findings reported an additive-log identity residual below `3e-14` percentage points. The EFT1 median is small, but the baseline model-level solar EFT1 increment spanned about `-4.438%` to `+23.420%`. The final Weinberg solar increment was approximately `-12.35%` to `-12.32%` in that restricted heavy-mass scan; this is **not** a general correction factor. For `sin²(theta13)` the documented median UV change at reference was about `+0.1974%`, compared with about `+0.000214%` for final Weinberg evolution, while EFT1 mixing effects were near the numerical-noise level for these inputs.

At nonbaseline heavy-mass factors the fitted Yukawas were **not refitted**. Hence these points are controlled sensitivities, not alternative accepted fits. Four missing source sets are *coverage gaps*, not failed representations. Results depend on the final charged-lepton basis and numerical inputs; this is not a scheme-independent partition of quantum corrections.

**Saved numerical sources:**

```text
Reports/output/full/staged_rg_comparison/smallY_smallL/scan.json
Reports/output/full/weinberg_running_comparison/smallY_smallL/scan.json
Reports/output/full/fermion_scale/smallY_smallL/scan.json
output/full/comparison/...
Numerical/benchmarks/fitted_reference/benchmark_manifest.json
```

Some historical numerical source files live under `Reports/output/full/` even though their path contains `Reports`: they are **data**, not expendable duplicate figures. Preserve them and their hashes. The original `Numerical/diagnostics/ExtractFullComparisonResults.py` reads saved diagnostics; it does not perform a new physical fit.

## 6. Conditional general-quartic campaign (separate exploratory study)

The general-scalar 16×3 campaign is **not** the final reduced benchmark. Its frozen source is:

```text
output/studies/general_scalar_extended/bulk_48/bulk48_v1/campaign_frozen.csv
Reports/output/general_scalar_48_final/provenance.json
Reports/output/general_scalar_48_final/tables/benchmark_comparison.csv
```

The reporting script `scripts/BuildGeneralScalar48FinalReports.py` validates the **exact** frozen campaign CSV hash, its 48 unique model/scenario entries and **20/6/22** result classification **before writing**. Two accepted B(-1)/C(-1) UV-running fits had explicitly **relaxed `lambdaT3` UV boundaries**; those must be kept separate from original-boundary acceptances. Some accepted B/C cases were imported from independently validated exploratory fits, not newly refitted by the 48-case bulk runner. Class E input `lambdaT3` retains the half-normalisation.

The candidate general-quartic integrations include additional renormalisable quartic running where exported, but **do not** provide the full extra holomorphic quadratic sector, certified Matchete/RGBeta scalar map, complete coupled intermediate Wilson system or finite general-scalar threshold matching. Four extreme-hypercharge models use legacy one-scale matching. Therefore 20 conditional acceptance labels do **not** establish 20 physically validated general-scalar models, and `NoMassSafeNeutrinoFit` is not the same outcome as `PoorFit` (the former never reached a valid objective evaluation).

The frozen campaign CSV stores **endpoint records**, not complete scale series or sensitivity covariance. Its PDF/HTML atlas is a reporting layer, not evidence for newly run sensitivity ovals. Use saved supplementary B/C six-fit trajectory samples for their own labelled source family; do not relabel an earlier reduced-scalar or fixed-grid RGE dashboard as a general-scalar trajectory.

## 7. Figures, interactive reports and what they actually depict

| Artifact | Source / meaning |
|---|---|
| `final_reduced_v2/frozen_chi2_comparison.png` | 32 individual and 32 shared-reference saved chi-square values |
| `final_reduced_v2/running_uv_formal_mass_endpoints.png` | Historical high-scale scalar boundary failure against the original 0.5–5 TeV window |
| `final_reduced_v2/running_uv_scalar_trajectories.png` | **UV-only** candidate renormalisable replay, `1e7 -> 1e5 GeV`, not a physical scalar-threshold plot |
| `Reports/output/reduced_scalar_lambda_style_v1/` | Saved-result HTML/PDF/tables using original lambda-study reporting style; not new fits or RG computation |
| `Reports/output/display_reduced_scalar/reduced_frozen_overnight_v2/` | Reduced-scalar SM+Weinberg RGE and T1/T2/T3 sensitivity ovals, with explicit failed gamma-point accounting |
| `Reports/output/general_scalar_48_final/` | Conditional general-quartic 48-case status matrix, saved observables and frozen-input atlas |

The reduced interactive gamma scans **re-evaluate forward EFT1 and low-energy diagnostics** rather than invent them from frozen campaign endpoints. High-gamma solver failures are retained in source JSON as `NumericalFailure`; the oval renderer omits models lacking sufficient valid T1/T2/T3 points and never interpolates failed ones. **Ovals are empirical visualization envelopes, not confidence contours.** Where a dashboard contains fewer than 16 models, this means diagnostic coverage was insufficient, **not** physical model exclusion.

Do not delete `output/studies/reduced_scalar_interactive_diagnostics/reduced_frozen_overnight_v2/` merely because an HTML file already exists. Preserve the source `running_diagnostics.json`, `sensitivity_diagnostics.json` and `comparison_summary.json` to rebuild and audit the figures.

## 8. Reproducibility, scientific status and limitations

**Primary reproducibility anchors** are:

- The original matched hard/direct C5 artifacts and the model-specific legacy cases; relevant native RGBeta beta JSON exports and scalar-coupling registry.
- The frozen campaign `comparison.json`/`.csv` and individual case `fit.json` files, including the 18 fitted Yukawa parameters and statuses.
- `final_reduced_v2/manifest.json`, which stores exact source SHA-256 values, accepted counts and UV replay solver costs. Existing data must not be LF/CRLF re-encoded without updating/verifying its provenance.
- The saved `configs/nufit_6_1_ic24_no.json` and frozen input benchmark manifests.
- The original full/staged numerical `scan.json`/`.csv` data, not only their PNG/HTML displays.

**What is verified:** saved-file identities where manifests check hashes; the documented numerical fit and report contracts; the recorded 32/10 acceptance counts; successful candidate UV replay as reported by `final_reduced_v2`; existing diagnostic-source handling and the user's most recent complete Python regression (**2144 passed, 10 skipped, 1 warning** with `python -X utf8 -m pytest -q`). Tests were run after the earlier test fixes and installation of scalar documentation; this numerical/RESULTS package still requires local regression after installation.

**What remains conditional or unverified:** full physical Matchete/RGBeta scalar-coupling and component normalisation, arbitrary nonzero extra quadratic mass mixing and its C5 effect, genuinely complete general-scalar fermion-first Wilson running and finite scalar threshold matching, a consistent full two-loop truncation, convergence of all unexplored regimes, vacuum stability, electroweak-broken scalar spectra, and phenomenological confidence regions. A successful numerical solver, approximate five-observable fit or source-level tensor audit does not by itself prove those items.

**Defensible thesis statement:** *For a deliberately restricted zero-additional-scalar boundary and the existing fixed-order matching prescriptions, all 32 frozen model/scenario cases admit individually fitted Yukawa inputs satisfying the chosen five-observable screen, while only 10 of 32 evaluations succeed when the scenario-specific A(-2) Yukawas are held fixed across representations. Separately, the unchanged chosen UV-running boundary drives the scalar masses outside the assumed TeV window; this diagnoses the input scenario, not the entire class of T3 theories.*

### Quick source inspection (Windows PowerShell)

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Get-Content -Encoding UTF8 .\output\studies\reduced_scalar_final_reports\final_reduced_v2\summary.md
Import-Csv .\output\studies\reduced_scalar_final_reports\final_reduced_v2\frozen_comparison.csv |
    Format-Table model,case,chi2_individual,chi2_shared_A_m2 -AutoSize
Get-Content -Encoding UTF8 .\output\studies\reduced_scalar_final_reports\final_reduced_v2\manifest.json |
    Select-Object -First 24
python -X utf8 -m pytest -q --tb=line
```

No command above refits Yukawas or regenerates Matchete/RGBeta. See [NUMERICAL.md](NUMERICAL.md) for optional report-only and forward-diagnostic commands, which have separate costs and output directories.

## 9. References and historical documentation ownership

- NuFIT 6.1 normal ordering: https://www.nu-fit.org/?q=node/309 ; target table https://www.nu-fit.org/sites/default/files/v61.tbl-parameters.pdf . Dataset variant is fixed by the saved `configs/nufit_6_1_ic24_no.json` input.
- Antusch, Hinze & Saad, *Updated running quark and lepton parameters at various scales* (2026), https://arxiv.org/abs/2510.01312 . Rounded SM reference values are starting inputs, not a T3 prediction.
- Restrepo, Zapata & Yaguna (2013), https://arxiv.org/abs/1308.3655 .
- Ma (2006), https://arxiv.org/abs/hep-ph/0601225 .
- Saved original detailed findings: `Reports/notes/T3_Numerical_Findings.md`, older `INFO_T3_NUMERICAL_RESULTS.md`, `INFO_T3_Numerical_Inputs_Sources_and_Justification.md`, `docs/T3_REDUCED_SCALAR_FINAL_RESULTS.md`, `docs/T3_REDUCED_SCALAR_SOBOL_COMPARISON.md` and historical 48-case notes. Their critical results, limits and source pointers are consolidated here. **Do not retire those files until references have been audited and an external hash-preserving archive exists.**
