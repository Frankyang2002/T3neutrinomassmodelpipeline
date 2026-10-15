# T3 neutrino mass model — numerical findings log

Updated: 2026-10-15. Scenario: `smallY_smallL`. Local pipeline, **fermion-first EFT**.

## Scope and provenance

- Fitted 16-model sensitivity comparison: user-provided comparison summary and embedded saved trajectories in `fitted_model_envelopes_all.html`.
- Fermion threshold scan: reviewed previous user-uploaded `smallY_smallL.zip` (60 successful model/factor points, 4 unavailable fitted cases), as described in the previous numerical analysis.
- Final SM+Weinberg comparison: **rechecked directly** against the most recent user-uploaded `smallY_smallL.zip`, containing `scan.csv`, `scan.json`, and five `weinberg_running_*.png` figures.
- The *staged* UV/EFT1/Weinberg decomposition below is **implemented but not run on the user's local numerical artifacts yet**; do not quote stage contributions until its output is checked.

## Model and EFT structure

The numerically implemented trajectory is

\[
\mathrm{UV\ T3}\xrightarrow[\mu_F]{\mathrm{integrate\ F}}\mathrm{SM}+S_1+S_2
\xrightarrow[\mu_S]{\mathrm{integrate\ scalars}}\mathrm{SM}+C_5
\xrightarrow[\mu_{\rm low}]{\mathrm{RGE}}m_\nu.
\]

The corresponding scale ordering is `mu_UV > mu_F > mu_S > mu_low`.
The project uses existing RGBeta UV/EFT1 beta functions, the existing fermion/scalar boundary projections, and the existing `FinalC5TrajectoryEvaluator` for scalar-threshold matching. No new symbolic matching expressions are generated in these studies.

## Completed numerical results

### Coupling-rescaling sensitivity

Each model is independently fitted at `gamma = 1`. Three compensated UV deformations preserve the nominal leading product `lambdaT3 * (y1*y2 + y2*y1)` at the input scale:

- T1: `lambdaT3 -> gamma lambdaT3`, `y1,y2 -> gamma^(-1/2) y1,y2`.
- T2: `lambdaT3 -> gamma lambdaT3`, `y2 -> gamma^(-1) y2`.
- T3: `y1 -> gamma y1`, `y2 -> gamma^(-1) y2`.

The 16 fitted baseline points were verified against the saved fitted predictions in the user's comparison summary, with zero reported maximum prediction difference. The broad scan covered `10^-3 <= gamma <= 10^3`; some non-baseline T1/T2 points failed, so figures must **not interpolate across failures**. `0.95 <= gamma <= 1.05` was too narrow to show useful excursions relative to NuFIT ovals in the Zee-style observable-pair dashboard.

**Interpretation:** these are coupling sensitivities **with the full RG pipeline enabled**, not a direct measurement of an RG correction.

### Fermion mass / threshold variation

The physical heavy-fermion mass matrix and integration scale were both multiplied by a factor `r`, retaining fixed scalar threshold and other fitted input couplings:

\[
M_F\mapsto rM_F,\qquad\mu_F\mapsto r\mu_F,
\quad r\in\{0.5,0.75,1,1.5,2\}.
\]

- 12 matching-artifact-supported models produced 60/60 successful numerical points.
- Four additional fitted models lacked matching artifacts in the available calculation sets: `(1,3,2), alpha=-4`; `(2,2,1), alpha=-3`; `(2,2,3), alpha=-3`; `(3,3,2), alpha=-4`.
- `r=1` reproduced the fitted predictions for all 12 included models.
- The resulting mass-squared splittings changed strongly as `r` varied, while mixing angles changed comparatively little.

**Interpretation:** this combines mass dependence in matching and RG evolution. It does **not** isolate running contributions. Comparisons across independently fitted representations are descriptive rather than controlled single-parameter representation effects.

### Isolated final SM+Weinberg EFT running

For the same matched `C5(mu_S)` and a *common fully evolved low-scale charged-lepton basis*, compare the fully evolved Weinberg coefficient to the coefficient held fixed at the scalar threshold:

\[
\delta_O^{\mathrm{Weinberg}}=100\left[\frac{O(K(\mu_{\rm low}))}
{O(K(\mu_S)\ \mathrm{frozen})}-1\right]\%.
\]

**Verified from the 60 successful rows in the uploaded Weinberg comparison CSV:**

| Observable | Median relative correction | Minimum | Maximum |
|:--|--:|--:|--:|
| `delta_m21_sq_ev2` | -12.341907% | -12.563515% | -12.200077% |
| `delta_m3l_sq_ev2` | -12.341989% | -12.563428% | -12.200163% |
| `sin2_theta12` | -0.00006337% | -0.00010841% | +0.00001473% |
| `sin2_theta13` | +0.00020859% | -0.00317567% | +0.00024962% |
| `sin2_theta23` | -0.00019233% | -0.00039893% | -0.00017557% |

The ~-12.34% mass-squared correction corresponds to an approximately common neutrino-mass factor `sqrt(1 - 0.1234) ~ 0.936`, i.e. roughly a 6.4% reduction in masses, if common mass-matrix rescaling is an adequate approximation.

**Interpretation:** this isolates *running of Weinberg K after the scalar threshold*, while holding the final charged-lepton basis fixed. It does not isolate UV or intermediate-EFT beta-function contributions, nor all quantum corrections to the UV theory.

## Completed staged-RG ablation (verified 2026-10-15)

**Source:** user-uploaded staged scan archive `smallY_smallL.zip`, files `scan.csv`, `scan.json`, five `staged_rg_*.png` figures. The raw archive is not included in this note.

**Coverage:** 60/60 evaluated points succeeded (12 distinct fitted models, five mass factors each); four additional models are explicitly `MissingArtifacts`. The 12 full-running `111` evaluations at factor `r=1` reproduce the fitted benchmarks with maximum chi-squared `0.0037006394235763714` (< 0.01 acceptance threshold). Independent row-by-row comparison against the previous full-running CSV has **not** been performed here because that CSV was not supplied alongside this upload.

The four cumulative configurations, identical original input and existing matching functions, are:

| Stage | UV RGE | EFT1 RGE | Final Weinberg RGE |
|:--|:--:|:--:|:--:|
| `000` | off | off | off |
| `100` | on | off | off |
| `110` | on | on | off |
| `111` | on | on | on |

An off UV/EFT1 stage is a **zero-beta state transport**, not a deletion of its threshold matching. The final-EFT comparison freezes `K` but uses the same fully evolved low-energy charged-lepton basis for every configuration. Therefore this is a *fixed-basis, conditional staged-running diagnostic*, not a completely independent no-running SM calculation or an unambiguous decomposition into individual loop corrections. Activating an upstream stage changes the quantities subsequently used in matching.

For each positive observable `O`, the incremental stage percentage is `100*(O_after/O_before - 1)`; the total compares `111/000`. The individual percentages multiply rather than add. Alternatively, `100*ln(O_after/O_before)` yields an exactly additive decomposition. The maximum absolute row-wise discrepancy between the sum of all three log contributions and the total log ratio in the uploaded results is below `3e-14` percentage points.

### Mass-squared differences: median incremental percentages

**At the fitted benchmark `r=1` (12 models):**

| Stage | Solar `delta_m21_sq_ev2` | Atmospheric `delta_m3l_sq_ev2` |
|:--|--:|--:|
| UV `000 -> 100` | +19.951% | +19.892% |
| Intermediate EFT1 `100 -> 110` | +0.549% | +0.549% |
| Final Weinberg `110 -> 111` | -12.336% | -12.336% |
| Total `000 -> 111` | +6.268% | +6.079% |

These medians were taken independently per column and **must not be directly summed**. The main pattern is a positive UV contribution partially cancelled by negative final-EFT evolution; EFT1 has a smaller typical contribution but its sign and size vary significantly by model.

**Across all 60 successful points:**

| Stage | Median solar correction | Solar range |
|:--|--:|--:|
| UV | +21.361% | +5.792% to +62.58% |
| EFT1 | +0.937% | -4.442% to +28.52% |
| Weinberg | -12.342% | -12.564% to -12.200% |
| Total | +6.774% | -11.38% to +69.21% |

These are results for particular independently fitted models and mass ratios, not a universal representation-only scaling law.

### Dependence on the fermion threshold

All figures below are **medians over the 12 models** and use the solar mass-squared splitting. The common scalar threshold is `mu_S=1000 GeV`; the baseline fermion threshold is `mu_F=100000 GeV`.

| Mass factor `r` | `mu_F/mu_S` | UV | EFT1 | Weinberg | Total |
|--:|--:|--:|--:|--:|--:|
| 0.50 | 50 | +23.568% | -0.004% | -12.350% | +9.092% |
| 0.75 | 75 | +21.405% | +0.303% | -12.342% | +7.481% |
| 1.00 | 100 | +19.951% | +0.549% | -12.336% | +6.268% |
| 1.50 | 150 | +18.012% | +0.937% | -12.328% | +4.740% |
| 2.00 | 200 | +16.678% | +1.240% | -12.322% | +3.760% |

Weinberg running varies little with `mu_F` because its integration limits `mu_S -> mu_low` are fixed. UV and EFT1 effects are more sensitive to the fermion matching scale, as expected from their changing running intervals and threshold inputs.

### Flavour-mixing diagnostics

At `r=1`, median percentage changes were:

| Stage | `sin2_theta12` | `sin2_theta13` | `sin2_theta23` |
|:--|--:|--:|--:|
| UV | -0.0163% | +0.1974% | -0.00758% |
| EFT1 | approximately 0 | approximately 0 | approximately 0 |
| Weinberg | -0.0000650% | +0.000214% | -0.000192% |
| Total | -0.0163% | +0.1976% | -0.00778% |

The EFT1 variations in mixing angles are at floating-point noise level in this dataset, consistent with an intermediate-stage change predominantly to the overall neutrino mass scale in this setup. This is an observation about the saved calculations; it is not a proof that EFT1 can never change mixing in other scenarios.

### Representation-level differences at `r=1`

The solar EFT1 increment spans **-4.438% to +23.420%** over the 12 fitted models, much broader than the final Weinberg increment (**-12.519% to -12.201%**). For example, `(2,2,3), alpha=-1` gives EFT1 **-4.438%**, whereas `(3,3,2), alpha=+2` gives **+23.420%**. Because each representation uses an independent fit, the difference combines representation and parameter dependence.

### Validation and limitations

1. Confirmed 60 successful numerical rows plus 4 explicitly missing-artifact model rows.
2. Confirmed `111` baseline chi-squared values agree with the fitted benchmarks to the expected small level (maximum < 0.01).
3. Confirmed exact log-additivity numerically (< `3e-14` percentage points residual).
4. The `111` vs `110` median Weinberg correction is consistent with the previous final-EFT calculation (~ -12.34%), but a full **row-by-row independent comparison against earlier CSV data remains pending**.
5. Perturbativity and the relationship between running masses and fixed matching scales were not independently revalidated from the archived output; the scan enforces explicit scale ordering.
6. `000` is **not** a fully consistent all-betas-zero renormalisation scheme: the charged-lepton basis is held fixed to the `111` endpoint. Treat the stage contributions as comparative diagnostics.

## Figure and project paths

- Numerical script: `Numerical/orchestration/RunFittedStagedRGComparison.py`.
- Numerical outputs: `Reports/output/full/staged_rg_comparison/smallY_smallL/scan.csv`, `scan.json`, and five `staged_rg_*.png` plots.
- This note: `Reports/notes/T3_Numerical_Findings.md`.

### Reproduce (PowerShell from project root)

```powershell
python -m Numerical.orchestration.RunFittedStagedRGComparison
python -m Numerical.orchestration.RunFittedStagedRGComparison --run --limit 2
python -m Numerical.orchestration.RunFittedStagedRGComparison --run
```

## Next numerical step

No additional RGE computation is required for rough results. The next useful analysis is to select representative rows/plots and compare the staged corrections by representation, with controlled benchmark caveats. If a truly *all-running-off* reference is later required, define its charged-lepton and matching-scale conventions separately rather than relabelling this fixed-basis diagnostic.
