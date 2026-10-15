---
title: T3 Numerical Inputs — Sources, Assumptions and Justification
tags: [T3, numerics, references, benchmark, provenance]
status: audit-of-current-main-2026-10-14
---

# T3 numerical inputs — sources, assumptions and justification

> [!important] Scope
> This is a provenance register for the *current published `main` branch* of [T3neutrinomassmodelpipeline](https://github.com/Frankyang2002/T3neutrinomassmodelpipeline), checked 2026-10-14. It distinguishes literature inputs from chosen benchmark parameters, fitted quantities, calculations, and numerical tolerances. This document does **not** mean that all chosen T3 parameters have been experimentally constrained.

## 1. Classification scheme

- **L — literature-derived input**: value taken from an identified external analysis, with stated scale/convention.
- **B — benchmark assumption**: an illustrative, controlled choice; not a measurement.
- **F — fitted value**: output of a specified optimisation against experimental targets; not an independent prediction.
- **D — derived value**: calculated from source inputs or the model, using an explicit equation.
- **N — numerical setting**: implementation/integration choices, not physical observables.

## 2. Input provenance table

| Quantity | Used value / definition | Type | Source in repository | Justification and caveat |
|---|---|---|---|---|
| UV scale | $\mu_{\mathrm{UV}}=10^7\,\mathrm{GeV}$ | B | `Numerical/benchmarks/fitted_reference/*.json`, `base_state.mu_gev` | Chosen high-energy initial-condition scale; selected to match an available literature SM running-coupling table. Not a measured T3 scale. |
| SM UV gauge couplings | $g_Y=0.38265$, $g_2=0.59493$, $g_3=0.7649$ | L (rounded) | Fitted-reference JSON `base_state.sm` | Antusch–Hinze–Saad (2026), Table II, $10^7$ GeV, SM $\overline{\mathrm{MS}}$ running; hypercharge converted to the project's convention. |
| UV fermion Yukawas | $y_u=4.39\times10^{-6}$, $y_c=0.00222$, $y_t=0.6462$; $y_d=9.7\times10^{-6}$, $y_s=1.93\times10^{-4}$, $y_b=0.00961$; $y_e=2.8681\times10^{-6}$, $y_\mu=6.0421\times10^{-4}$, $y_\tau=0.010264$ | L (rounded) | Fitted-reference JSON `base_state.sm.yu/yd/ye` | Central values from Antusch–Hinze–Saad Table II, $10^7$ GeV. **Diagonal real matrices are an additional texture/basis simplification** and do not encode CKM mixing. |
| SM trace | $T=\mathrm{Tr}(3Y_u^\dagger Y_u+3Y_d^\dagger Y_d+Y_e^\dagger Y_e)$ | D | Matrices above | Approximately $1.253$ from the numerical matrices. Recompute from exact config when needed. |
| Higgs VEV | $v=246.22\,\mathrm{GeV}$ | L / conventional derived | Fitted-reference JSON `numerical.vev_gev`; `physics/NeutrinoMass.py` | Standard electroweak convention $v=(\sqrt{2}G_F)^{-1/2}$. Distinguish from $v/\sqrt2\simeq174$ GeV and the repository convention $m_\nu=-v^2 C_5/2$. |
| Fermion thresholds | $\mu_F=10^5\,\mathrm{GeV}$ | B | Fitted-reference JSON `scales.mu_fermion_threshold_gev` | Benchmark EFT threshold; no direct experimental determination or T3 fit to this scale. |
| Scalar thresholds | $\mu_S=10^3\,\mathrm{GeV}$ | B | Fitted-reference JSON `scales.mu_scalar_threshold_gev` | Benchmark hierarchy $\mu_F>\mu_S$ compatible with implemented fermion-first EFT path. |
| Low scale | $\mu_{\mathrm{low}}=100\,\mathrm{GeV}$ | B/N | Fitted-reference JSON `scales.mu_low_gev` | Reference endpoint near electroweak scale, not the actual neutrino-oscillation experiment scale. |
| Heavy flavour masses | $M_F=\operatorname{diag}(10^5,1.1\times10^5,1.2\times10^5)\,\mathrm{GeV}$ | B | Fitted-reference JSON `base_state.ordinary.MF` | Chosen slightly nondegenerate flavour benchmark; do not call these measured/pole masses. Matrix structure and allowed Majorana/vectorlike treatment depend on model. |
| Scalar mass parameters | $m_{S_1}^2=1\,306\,122.448979592\,\mathrm{GeV}^2$, $m_{S_2}^2=10^6\,\mathrm{GeV}^2$ | B | Fitted-reference JSON `base_state.ordinary` | Chosen Lagrangian mass-squared parameters, not independently established physical pole masses. |
| Higgs quartic | $\lambda_H=0.25$ | B | Fitted-reference JSON `base_state.sm.lambdaH` | Explicit benchmark assumption, **not** obtained from Table II of the SM-running paper. Do not claim measured Higgs-mass matching at $10^7$ GeV. |
| Scalar self-quartics | $\lambda_{S_1}=0.10$, $\lambda_{S_2}=0.12$ | B | Fitted-reference JSON | Selected illustrative perturbative values; no cited experimental determination. |
| Higgs portals | $\lambda_{H1}=0.02$, $\lambda_{H2}=0.03$ | B | Fitted-reference JSON | Selected small nonzero portal couplings; detailed vacuum stability and direct-detection phenomenology are not established by this choice. |
| Scalar–scalar portal | $\lambda_{12}=0.01$ | B | Fitted-reference JSON | Small benchmark coupling. |
| Extra quartics | Optional adjoint/cross quartics set to 0 in baseline | B | Fitted-reference JSON `base_state.ordinary` | **Initial conditions**, not generally symmetry-enforced or necessarily RG-invariant. |
| Mixing quartic | $|\lambda_{T3}|=0.1$ (A–D), $0.05$ (E) | B | Fitted-reference JSON; `INFO_T3_MULTIPLICITY.md` | A–D nominal small-L reference; E halved for two equal-weight charge channels in benchmark. This is *only* a benchmark normalization, not a replacement for CG/group factors. |
| New Yukawas $y_1,y_2$ | Distinct $3\times3$ real fitted matrices for each of 16 models | F | `Numerical/benchmarks/fitted_reference/<model>.json` | Calibrated to NuFIT oscillation targets; not predicted uniquely and not justified by a single universal $Y=0.005$ or $0.5$. |
| Target oscillation data | Five normal-ordering parameters listed below | L, then Gaussian approximation | `configs/nufit_6_1_ic24_no.json` | NuFIT 6.1 data source, with arithmetic-average asymmetric 1-sigma uncertainties; correlation matrix / published non-Gaussian likelihood not used. |
| ODE tolerances | `rtol=1e-8`, `atol=1e-11`, 120 stored scale points | N | Fitted-reference JSON `numerical` | Solver accuracy/sampling choices. Numerical convergence under tightening tolerances must be tested rather than inferred from settings. |
| One-loop normalization | $16\pi^2\,dX/d\ln\mu=\beta_X^{(1)}$ | D / convention | `Numerical/core/BetaVector.py`; `INFO_NUMERICAL.md` | Evaluator produces one-loop numerator; factor of $1/(16\pi^2)$ applied once in vector integration. |

## 3. SM source: exact citation and hypercharge convention

**Primary source:** Stefan Antusch, Kevin Hinze and Shaikh Saad, *Updated running quark and lepton parameters at various scales*, **Physical Review D 113, 095011 (2026)**, [arXiv:2510.01312](https://arxiv.org/abs/2510.01312), [published PDF](https://dirros.openscience.si/Dokument.php?id=44286&lang=eng). **Table II**, column $\mu=10^7$ GeV, uses 2024 PDG inputs and SM $\overline{\mathrm{MS}}$ running.

The paper's Table II gives $g_1=0.493977$ in GUT normalization. The numerical config uses

$$
g_Y=\sqrt{\frac35}\,g_1\simeq0.38265.
$$

Consequently

$$
g_Y^2\simeq0.14642,\qquad g_2^2\simeq0.35394,\qquad g_3^2\simeq0.58507.
$$

The UV SM inputs are **SM-only extrapolations**, not a fully self-consistent determination including the extra T3 fields above their thresholds. These are reasonable reference boundary conditions for a demonstrative comparison, provided this distinction is stated.

The input SM matrices are diagonal in the benchmark. The paper also supplies quark mixing parameters; their omission means this benchmark does not implement a realistic CKM matrix in $Y_u,Y_d$ simultaneously.

## 4. Oscillation target and fitting procedure

**Primary target:** [NuFIT 6.1 parameter table](https://www.nu-fit.org/sites/default/files/v61.tbl-parameters.pdf); [NuFIT release page](https://www.nu-fit.org/?q=node/309). The *repository's* target file is [`configs/nufit_6_1_ic24_no.json`](https://github.com/Frankyang2002/T3neutrinomassmodelpipeline/blob/main/configs/nufit_6_1_ic24_no.json). Dataset label in this file reads “IC24 with Super-Kamiokande atmospheric data.” **Recheck the selected table subsection/data-release label before final thesis citation**, since the externally indexed NuFIT PDF also lists an IC23 variant; do not conflate alternatives.

| Observable | Target central value | Approximate $1\sigma$ used |
|---|---:|---:|
| $\Delta m_{21}^2$ | $7.537\times10^{-5}\,\mathrm{eV}^2$ | $9.7\times10^{-7}\,\mathrm{eV}^2$ |
| $\Delta m_{3\ell}^2$ (NO: $\Delta m_{31}^2$) | $2.511\times10^{-3}\,\mathrm{eV}^2$ | $2.05\times10^{-5}\,\mathrm{eV}^2$ |
| $\sin^2\theta_{12}$ | $0.3088$ | $0.00665$ |
| $\sin^2\theta_{13}$ | $0.02248$ | $0.00057$ |
| $\sin^2\theta_{23}$ | $0.470$ | $0.0155$ |

The code [`Numerical/fitting/OscillationFit.py`](https://github.com/Frankyang2002/T3neutrinomassmodelpipeline/blob/main/Numerical/fitting/OscillationFit.py) builds a diagonal covariance matrix from supplied independent errors when the JSON contains no covariance. Therefore its fitting objective is approximately

$$
\chi^2_{\rm approximate}=\sum_{i=1}^{5}\left(\frac{O_i^{\rm model}-O_i^{\rm NuFIT}}{\sigma_i}\right)^2.
$$

**Not the official NuFIT $\Delta\chi^2$ surface.** An excellent fitted $\chi^2$ is expected by construction when the fitting parameters are sufficiently flexible. It is *not* a test of predictive power. The repository's frozen-reference README says all 16 selected candidates were independently re-evaluated with agreement, and all report $\chi^2<0.01$; this checks reproducibility of the implemented evaluator, not broad phenomenological viability.

**Optimisation provenance:** [`Numerical/orchestration/RunFitFeasibility.py`](https://github.com/Frankyang2002/T3neutrinomassmodelpipeline/blob/main/Numerical/orchestration/RunFitFeasibility.py) initially scales Yukawas logarithmically, then explores 18 real $y_1/y_2$ matrix entries using seeded Sobol sampling and limited local least-squares refinement. The chosen final candidate for each model is recorded in [`benchmark_manifest.json`](https://github.com/Frankyang2002/T3neutrinomassmodelpipeline/blob/main/Numerical/benchmarks/fitted_reference/benchmark_manifest.json), with stage, predictions, norms, $\chi^2$, independent confirmation and SHA-256 hashes. Some models use retry/focused stages, not just the first pilot. The 16 frozen reference files are the definitive fitted input matrices.

## 5. Benchmark spectrum and coupling justification

The mass hierarchy in the numerical setup is

$$
\mu_{\rm UV}=10^7\,\mathrm{GeV}>\mu_F=10^5\,\mathrm{GeV}>\mu_S=10^3\,\mathrm{GeV}>\mu_{\rm low}=10^2\,\mathrm{GeV}.
$$

This permits the **fermion-first EFT sequence actually implemented** in the repository. It does **not** establish that T3 particles physically live at these masses. The choice allows study of a hierarchy with separate threshold stages. The heavy-flavour mass entries $100,110,120$ TeV and scalar mass-squared values are benchmark inputs; distinguish Lagrangian running masses from pole masses and from nominal matching scales. For detailed scale sensitivity, threshold variations would be an additional study, not something already proven by the baseline.

For the scalar couplings, the defensible statement is: *we select perturbative UV test values, with optional interactions initially zero to reduce parameter dimensionality.* It would be incorrect to say that $\lambda_H=0.25$, $\lambda_{S_i}$ or portals are experimentally fixed by the cited SM-running table. The Higgs quartic in particular needs an explicit potential-normalization check and Higgs-mass/RGE matching study before calling it a physical SM determination.

Class-E mixing-quartic initialization follows the charge-channel audit in [`INFO_T3_MULTIPLICITY.md`](https://github.com/Frankyang2002/T3neutrinomassmodelpipeline/blob/main/INFO_T3_MULTIPLICITY.md): A–D have one compatible channel; E has two equal-weight channels under the unbroken-SU(2) multiplet-mass benchmark. Hence $\lambda_{T3}^{E}=\lambda_{T3}^{\rm ref}/2$. The true representation-specific Clebsch–Gordan factors remain in model calculations. Do **not** multiply the physical E amplitude by an extra factor of two on top of what the pipeline calculates.

## 6. Model and physics sources (background, not numeric parameter determinations)

- [Restrepo, Zapata & Yaguna, *Models with radiative neutrino masses and viable dark matter candidates* (2013)](https://arxiv.org/abs/1308.3655): one-loop topology and model-class context. It does **not** determine this project's benchmark mass spectrum or fitted Yukawa matrices.
- [Ma, *Verifiable radiative seesaw mechanism of neutrino mass and dark matter* (2006)](https://arxiv.org/abs/hep-ph/0601225): illustrative one-loop radiative mass mechanism; not a substitute for the representation-specific T3 matching.
- [Aebischer, Bresciani & Selimović, *Anomalous Dimension of a General Effective Gauge Theory II: Fermionic Sector*](https://arxiv.org/abs/2512.16890): context for general EFT anomalous-dimension calculations; the project's actual exported RGBeta equations and group-factor validations remain authoritative for its implemented T3 RGE.
- Repository implementation: [`INFO_LAGRANGIAN.md`](https://github.com/Frankyang2002/T3neutrinomassmodelpipeline/blob/main/INFO_LAGRANGIAN.md), [`INFO_RGE.md`](https://github.com/Frankyang2002/T3neutrinomassmodelpipeline/blob/main/INFO_RGE.md), [`INFO_NUMERICAL.md`](https://github.com/Frankyang2002/T3neutrinomassmodelpipeline/blob/main/INFO_NUMERICAL.md).

## 7. Reproducing the documented input points

From the repository root in PowerShell:

```powershell
# List the frozen fitted reference benchmarks
Get-ChildItem Numerical\benchmarks\fitted_reference\T3_*.json

# Show the manifest and its provenance metadata
python -c "import json; p=json.load(open('Numerical/benchmarks/fitted_reference/benchmark_manifest.json')); print('models:',len(p['models'])); [(print(x['model'], x['fit_chi2'], x['y1']['frobenius_norm'], x['y2']['frobenius_norm'])) for x in p['models']]"

# Inspect exact UV inputs for a selected model
python -c "import json; p=json.load(open('Numerical/benchmarks/fitted_reference/T3_dS1_2_dS2_2_dF_1_alpha_m1.json')); print(json.dumps(p['base_state'],indent=2))"
```

To reproduce *computed results*, not just inspect inputs, reuse the existing model-specific UV/intermediate RGBeta and matched-$C_5$ artifacts plus the numerical runner. Do not assume that the complete artifact set is contained in the public GitHub tree. Reproducing a NuFIT fit requires the target file and the fit-stage configurations. Avoid overwriting frozen reference configurations.

## 8. Open justification points to resolve before thesis submission

1. **Verify NuFIT dataset variant and release-date text.** Match the exact intended table subsection to `ic24_no.json` and explain why five observables are used (CP phase and absolute neutrino mass are not constrained by this objective).
2. **Document why these threshold masses were chosen**, beyond computational convenience. They are scenario assumptions, not fitted mass estimates.
3. **Check Higgs quartic convention and matching:** $\lambda_H=0.25$ should not be presented as SM $\lambda_H(10^7\,\mathrm{GeV})$ from the running-coupling paper.
4. **Record the CKM approximation:** diagonal $Y_u,Y_d$ omits realistic quark-flavour mixing. Decide whether it is numerically acceptable at the desired precision.
5. **Check vacuum stability/perturbativity if making phenomenological claims.** Small quartics and vanishing optional quartics do not, by themselves, prove stability or RG-invariance.
6. **Assess threshold/pole-mass distinction and dependence on chosen scales** before treating the numerical outputs as physical constraints.
7. **Preserve exact fitted-reference SHA-256/config paths** with all figures and tables, so the presented benchmark can be rerun without silent parameter drift.

## 9. Thesis-ready interpretation

> The Standard Model gauge and charged-fermion Yukawa couplings at the ultraviolet scale $\mu_{\rm UV}=10^7\,\mathrm{GeV}$ were initialized using the $\overline{\rm MS}$ running values reported by Antusch, Hinze and Saad (2026), Table II. The additional T3 scalar couplings, heavy-field masses and threshold scales were chosen as benchmark parameters to define a reproducible fermion-first multi-scale scenario. The new Yukawa matrices $y_1$ and $y_2$ were fitted independently for each representation to five normal-ordering NuFIT 6.1 oscillation observables using an approximate diagonal Gaussian likelihood. The resulting agreement with those fitted observables constitutes numerical calibration, not an independent model prediction. Model-dependent anomalous dimensions, group factors and threshold matching are supplied by the pipeline's generated RGE and matching expressions. The benchmark's interpretation is consequently restricted to the selected parameter values, approximations and validated stages of the numerical calculation.

## 10. Primary links

- Antusch, Hinze & Saad (2026): https://arxiv.org/abs/2510.01312 ; https://dirros.openscience.si/Dokument.php?id=44286&lang=eng
- NuFIT 6.1: https://www.nu-fit.org/?q=node/309 ; https://www.nu-fit.org/sites/default/files/v61.tbl-parameters.pdf
- PDG: https://pdg.lbl.gov/2025/listings/particle_properties.html
- T3 repository: https://github.com/Frankyang2002/T3neutrinomassmodelpipeline
- Fitted benchmark manifest: https://github.com/Frankyang2002/T3neutrinomassmodelpipeline/blob/main/Numerical/benchmarks/fitted_reference/benchmark_manifest.json
