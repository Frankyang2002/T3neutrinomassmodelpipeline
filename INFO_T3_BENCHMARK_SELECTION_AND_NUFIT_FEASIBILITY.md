# INFO — T3 benchmark selection, NuFIT feasibility and sensitivity comparisons

**Project:** T3 neutrino mass generation  
**Status:** Methodology and preliminary numerical result; the 16-model feasibility results are not yet available here.  
**Scope:** Ordinary T3 model classes A–E, normal ordering (NO) in the completed pilot.

## 1. Motivation

The project compares neutrino observables **within** a T3 model under controlled parameter variations and **between** different T3 representations. The original numerical comparison used four nominal UV coupling regimes:

| Scenario | Nominal Yukawa regime | Nominal scalar-coupling regime |
|---|---|---|
| `smallY_smallL` | small | small |
| `smallY_largeL` | small | large |
| `largeY_smallL` | large | small |
| `largeY_largeL` | large | large |

Here `Y` denotes the Yukawa couplings and `L` the T3 scalar coupling, conventionally `lambdaT3`. These labels describe **input benchmark categories**, not evidence that the points reproduce neutrino data. The exact numeric small/large boundaries must be taken from the relevant generated configuration files; they have not been independently established in this note.

The original 16 representations × 4 scenarios = **64 fixed-benchmark cases** were used for numerical evolution and sensitivity plots. Successful numerical evolution does **not** imply a successful fit to oscillation measurements. Some fixed-benchmark predictions were far outside the NuFIT region.

## 2. Why introduce fitted benchmarks?

For a radiative T3 neutrino mass, the leading dependence is schematically

\[
m_\nu \sim \frac{v^2}{16\pi^2 M}\,\lambda_{T3}\,(y_1\,\mathcal F\,y_2^T + y_2\,\mathcal F^T\,y_1^T),
\]

where the actual loop, threshold, representation and RGE dependence is evaluated by the numerical pipeline rather than this schematic expression. Thus the overall mass scale is sensitive to the coupling product, while the mixing angles and relative mass splittings depend on the **flavour textures**.

Comparing sensitivity curves based on arbitrary points far from the observed region mainly reveals properties of those arbitrary inputs. We therefore introduced a **NuFIT-targeted feasibility fit** to find physically relevant starting points before constructing the final within-model and between-model sensitivity comparisons.

This is a deliberately **data-calibrated benchmark selection**. The close agreement of the fitted central point with NuFIT is a consequence of fitting to NuFIT, **not an independent prediction or validation**. The subsequent sensitivity trajectories, model-dependent differences, parameter constraints and any outliers are the objects of comparison.

## 3. Feasibility procedure actually implemented

The standalone pilot is `Numerical/orchestration/RunFitFeasibility.py`; the batch wrapper is `Numerical/orchestration/RunBatchFitFeasibility.py`. They reuse the existing T3 point evaluator, the model-specific saved `uv_rgbeta_rge.json`, `eft1_rgbeta_rge.json` and `final_weinberg_coefficient.json` artifacts, and the NuFIT 6.1 target JSON.

For each selected model:

1. **Mass-scale calibration:** multiply both initial Yukawa matrices by the same positive factor, testing 15 logarithmically spaced factors between \(10^{-7}\) and \(1\). Retain the factor giving the smallest evaluated oscillation \(\chi^2\).
2. **Flavour exploration:** independently vary the **18 real entries** of the two \(3\times3\) Yukawa matrices using 256 scrambled Sobol samples. Each entry is bounded within three times the largest absolute entry of its calibrated matrix, with either sign allowed. The calibrated baseline is also retained as a candidate.
3. **Limited local refinement:** probe coordinate sensitivities and locally optimise six selected real Yukawa entries, with a default maximum of 30 least-squares function evaluations.
4. **Record:** save the best-fit \(\chi^2\), the five fitted oscillation observables, the fitted Yukawa parameters, calibration history and any failures. The batch wrapper records all models, including outliers, in CSV and JSON summaries.

The pilot **holds `lambdaT3` fixed** at the value in each input configuration. It does **not** enforce strict `smallY`/`largeY` bounds after calibration and flavour fitting. Accordingly, these are **feasibility fits**, not yet a controlled comparison of four disjoint coupling regimes. A subsequent constrained fit must define and enforce each regime's numeric bounds before interpreting fitted results as belonging to that regime.

The full repository also contains `Numerical/fitting/SobolBenchmarkSearch.py`, which supports a broader Sobol search and multistage local least-squares refinement. The feasibility runner is a cheaper preliminary test; its \(\chi^2\) minimum is not guaranteed to be global.

## 4. Verified pilot result: one model only

The uploaded `fit_feasibility.zip` contains a completed pilot for

- **Class A**, \(\alpha=-2\), representation \((d_{S_1},d_{S_2},d_F)=(1,3,2)\);
- model key `T3_dS1_1_dS2_3_dF_2_alpha_m2`;
- scenario label `smallY_smallL`;
- normal mass ordering.

| Quantity | Result from uploaded feasibility summary |
|---|---:|
| Best \(\chi^2\) | \(3.7042522209895953\times10^{-9}\) |
| \(\Delta m^2_{21}\) | \(7.5369991070\times10^{-5}\ \mathrm{eV}^2\) |
| \(\Delta m^2_{3\ell}\) | \(2.5109999837\times10^{-3}\ \mathrm{eV}^2\) |
| \(\sin^2\theta_{12}\) | 0.3088000221 |
| \(\sin^2\theta_{13}\) | 0.0224800020 |
| \(\sin^2\theta_{23}\) | 0.4700009295 |
| Best uniform Yukawa calibration factor | 0.01 |
| \(\chi^2\) after calibration, before flavour optimisation | 166174.9261 |
| Fixed \(\lambda_{T3}\) in fitted config | 0.1 (real) |
| Local-refinement evaluations reported | 30 |

The calibrated Yukawa magnitude **alone** did not give a good fit; allowing the flavour structure to vary was essential. This establishes feasibility for **one representation and one configuration** only. It does not show that the other 15 models or the other three coupling scenarios have small \(\chi^2\).

The reported fit contains five oscillation observables. In particular, \(\delta_{\rm CP}\) is **not** included in the saved best-fit prediction dictionary, so this result should not be described as a fit to every measured neutrino observable. The extremely small \(\chi^2\) indicates close numerical agreement with the configured five target values; it does not by itself establish model predictivity or experimental viability beyond them.

## 5. Planned 16-model feasibility comparison and outliers

The first batch covers the 16 representations under the `smallY_smallL` input scenario. All models should use the same NuFIT ordering/target, the same fitting algorithm and numerical conventions, but retain their **own representation-specific RGE and matching artifacts**. Each model's Yukawa flavour texture is fitted independently.

The comparison should report, for **every** model, the best \(\chi^2\), five predicted observables, best Yukawa scale, fit status and fitted parameters. Cases that fail to reach the target region must remain visible as **outliers**, rather than being discarded or forced to agree.

A low \(\chi^2\) is evidence of fit feasibility **within the explored parameter space**. A high \(\chi^2\) may indicate limited search coverage, restrictive bounds, numerical issues or genuinely poor agreement; it is not proof of exclusion without further checks. Screening values such as \(\chi^2<10\) are practical flags, not formal confidence statements.

## 6. Class E multiplicity

Class E has a distinct multiplicity/charge-channel structure compared with the other T3 classes. Earlier project benchmark conventions include an effective normalisation for Class E (for example `lambdaT3_effective = lambdaT3_reference / 2` in the separate full-study benchmark convention). **That convention must not be blindly multiplied into the numerical fit a second time.**

The feasibility runners reuse each representation's existing RGE and final Weinberg-coefficient artifacts through the standard evaluator. They do **not** introduce a new universal multiplicity factor for E. This is the correct approach for avoiding an *additional* ad hoc factor, but the feasibility scan **does not independently verify** that the upstream Class E matching and multiplicity implementation is correct.

When comparing fitted coupling values between Class E and other classes, explicitly state whether a reported coupling is the input `lambdaT3`, an effective normalised coupling, or another convention. The representation-specific matching should be checked independently before attributing a sensitivity difference solely to multiplicity.

## 7. Within-model and between-model comparisons

**Within a model:** choose its fitted benchmark and apply the same T1/T2/T3 transformations over the same \(\gamma\) range (initially 0.95–1.05). Plot observable–observable trajectories and compare their size, direction and curvature against the NuFIT uncertainty regions. These transformations approximately preserve leading coupling products, so small shifts can be physically meaningful.

**Between models:** fit each representation independently to the same oscillation target, then compare trajectories around each model's own fitted point. Do **not** impose identical Yukawa matrices across inequivalent representations, and do **not** artificially force all \(\chi^2\) values to be equal. Model-specific fitting success, outliers and sensitivity remain results.

Because the fitted central values are intentionally near NuFIT, proximity of the central points is **not** an independent prediction. The interesting results are the allowed parameter textures, residual fit quality, representation dependence and behaviour under controlled perturbations.

The existing 64 fixed-benchmark sensitivity outputs should remain separate from the fitted-benchmark sensitivity outputs. No fitted sensitivity scan is claimed to have been completed here.

## 8. Reproducibility and outputs

From the repository root in PowerShell:

```powershell
python -m Numerical.orchestration.RunBatchFitFeasibility
python -m Numerical.orchestration.RunBatchFitFeasibility --limit 2 --run
python -m Numerical.orchestration.RunBatchFitFeasibility --run
```

The default scenario is `smallY_smallL`. The batch output is

`Reports/output/full/fit_feasibility/smallY_smallL/batch_summary.csv`

and `batch_summary.json`, with each model's `feasibility_summary.json` and `fitted_candidate_config.json` in a model-specific subdirectory. Completed summaries are skipped unless `--rerun` is requested.

These files are **candidate fits** and should not silently overwrite the original `configs/generated_model_sets/comparison/.../input` benchmarks or `running_diagnostics.json` files. A later fitted-sensitivity runner must explicitly load the saved fitted candidate for the relevant model.

## 9. Interpretation and next decision

1. Run the 16-model feasibility batch and inspect the \(\chi^2\) distribution, including Class E and any outliers.
2. Independently re-evaluate promising saved fitted configurations and check units, mass ordering and observable extraction.
3. If most models fit, run within-/between-model sensitivity around those fitted points.
4. Only when comparing the four `small/large Y` and `small/large lambda` scenarios as **distinct physical regimes**, enforce explicit bounds on both Yukawa magnitudes and `lambdaT3` and repeat the fitting.

**Current conclusion:** one Class A model has demonstrated a very close five-observable fit. The feasibility of the other representations and the scientific comparison of constrained coupling regimes remain to be established numerically.
