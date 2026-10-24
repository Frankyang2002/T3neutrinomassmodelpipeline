# T3 radiative-neutrino-mass thesis: proposed narrative and chapter plan

**Status:** planning document, 24 October 2026. The packaged plots use the **saved optimised, UV-frozen** `lambdaT3 = 5e-9` benchmarks. The high-scalar (`0.5`) and controlled full-UV versus UV-frozen comparisons remain to be verified. Do not present proposed conclusions as established results.

## Central research question

**How do threshold matching and renormalisation-group (RG) evolution modify the magnitude and flavour structure of radiatively generated neutrino masses across T3 representations, and how do these quantum effects constrain the underlying Yukawa and scalar couplings?**

The thesis should separate **the flavour orientation** of the neutrino mass matrix (mixing angles), **its eigenvalue scale and splittings**, and **the degeneracies of the underlying parameters**. These three questions are related but not identical.

## Logical spine

1. **The T3 mechanism supplies a Majorana mass.** Integrating out fermions creates the tree-level LLSS coefficient `C12`, whose flavour tensor is the mass-weighted symmetrised Yukawa product. Subsequent scalar-threshold matching and running yield the Weinberg coefficient `C5`, which becomes `m_nu = -v^2 C5 / 2` after EWSB.
2. **At leading order, matching depends on combinations of parameters.** A simplified dependence is `C5 ∝ lambdaT3 × Sym(Y1* W Y2†)`, where `W` contains the heavy-state mass/loop weights. Changing `lambdaT3` while compensating Yukawa magnitudes can preserve this combination. This motivates, but does not by itself prove, equivalent physical predictions at the complete matched order.
3. **Running may distinguish matching-degenerate points.** RGEs for Yukawas, scalar couplings, and gauge interactions are not generally invariant under those compensating rescalings. Determine whether departures arise primarily in the full UV interval, the scalar-only EFT, or the final SM+Weinberg EFT.
4. **Mass corrections need not produce large angle corrections.** The production one-loop `C12 -> C5` source is a scalar flavour-blind prefactor times the *same* `C12` tensor. Such a term transfers the flavour texture without rotating its eigenvectors. Hard matching and mass-resolved logarithmic weights may misalign the full `C5` relative to `C12`.
5. **Representations change RG and matching coefficients.** Compare the 16 T3 variants under both independently optimised and shared inputs. A claim that representation dominates over coupling strength requires a controlled cross-representation numerical comparison, not only similarly small fit chi-squares.

## Suggested chapter structure

### Chapter 1 — Introduction

- Motivation: neutrino masses and mixing beyond the minimal SM, radiative mass generation and EFT techniques.
- Research question above, clearly limited to the T3 family, chosen benchmarks, fixed-order calculation and available numerical stages.
- Three concrete objectives: quantify running effects; test matching degeneracies; explain mass-versus-mixing sensitivity.
- Scope: common fermion threshold; **do not** claim sequential decoupling or physical full coupled-Wilson NLO predictions.

### Chapter 2 — Physical and model background

- Three-flavour mixing, charged-lepton basis, PMNS parameterisation, Majorana versus Dirac masses, `Delta m^2` versus absolute mass scale.
- Weinberg operator, `C5 = C5^T`, Takagi diagonalisation for a complex symmetric matrix. Distinguish symmetric matrices from diagonal matrices.
- T3 field content and allowed SU(2)/hypercharge representation choices; introduce 16 model variants and their group-theory consequences.
- Loop mass generation and threshold hierarchy `mu_F > mu_S > mu_low`. Discuss literature and related radiative mechanisms.

### Chapter 3 — Analytical matching and RGE structure

**3.1 T3 Lagrangian and one-loop mass generation.** State definitions, index conventions, sign conventions and thresholds; write the full implemented matching expression where available (avoid replacing loop kernels by `1/M` except when labelled schematic).

**3.2 Fermion-threshold matching.** Define the symmetric LLSS flavour tensor:

`(C12)_{ij} = (1/2) sum_r [Y1_{ir}* Y2_{jr}* + Y2_{ir}* Y1_{jr}*]/M_r`.

Derive why the *symmetrised* form is required for a Majorana mass; stress that symmetry does not make the matrix diagonal. Show the heavy-state decomposition `C12 = sum_r F_r`.

**3.3 Scalar-only EFT and Weinberg generation.** Show the production one-loop source `16 pi^2 dC5/dln(mu) ⊃ (4 sqrt(3)/3) lambdaT3 C12` for T3-B in the project's operator conventions; separate **hard scalar matching** and **direct one-loop C12 -> C5** transport. Explain why the current production code contains the direct log in the matched C5 and should not double count it with an added numerical solver.

**3.4 Final SM+Weinberg evolution.** Explain beta-function terms proportional to C5 and charged-lepton flavour-dependent terms. State which observables can consistently be evaluated at which stage.

**3.5 Scaling directions and predictions.** Derive `lambda -> gamma lambda`, `Y1,Y2 -> gamma^{-1/2}Y1,Y2` and the two other T1/T2/T3 directions in the existing `BenchmarkSensitivity.py`. Explicitly state what is exact for the *leading bilinear* and what is merely a test for the complete matched coefficient with running.

**3.6 Flavour orientation versus eigenvalue rescaling.** Introduce `F = Sym(Y1* W Y2†)`, `U^T C5 U = D5`, and a nonalignment measure `||offdiag(U^T C12 U)||_F/||C12||_F`. Describe charged-lepton basis conventions. Note that `U` is fixed by `C5`; diagonal C5 in that basis is tautological, whereas near-diagonal C12 is not.

### Chapter 4 — Numerical framework and reproducibility

**4.1 Model inputs:** representations, heavy/scalar masses, boundary scales, Yukawa matrices, scalar couplings, whether values are real/complex.

**4.2 Pipeline:** UV beta functions (where enabled), common fermion matching, numerical renormalisable scalar-only EFT solution, symbolic fixed-order Wilson transport/hard matching, numerical Weinberg EFT evolution.

**4.3 Fits:** define target observables, objective `chi2`, mass ordering, charged-lepton basis, constraint conventions. Distinguish a small `chi2` from evidence that a model is preferred; the fitted Yukawa space is underdetermined.

**4.4 Controlled studies:** (i) UV on/off, (ii) different `lambdaT3` seeds, with or without refitting as separately labelled experiments, (iii) optimised/shared Yukawas across representations, (iv) matching-preserving rescaling directions.

**4.5 Diagnostics:** Frobenius norms and projection residuals, mass-basis off-diagonal fractions, heavy-state cancellation ratio, mass-squared drift and mixing-angle drift. Report input/output JSON paths and code versions.

### Chapter 5 — Results and physical interpretation

**5.1 Fit viability:** baseline 16-model optimised fit table and `fig03`. Note that *two* saved C-sector points have appreciable chi-square (~0.814 and ~10.87), so a `Success` flag does not imply an excellent fit.

**5.2 Coupling rescaling / matching degeneracy:** analytical `fig02`, then controlled numerical sensitivity scan(s). Test scaling at fixed matching and after running. Do not interpret plot `fig02` as a numerical result.

**5.3 Where quantum corrections matter:** inspect full UV-on vs UV-off on matched inputs, and decompose hard/direct/EFT running. `fig07` quantifies the UV-frozen T3-B scalar-threshold decomposition, **not** a general UV-vs-EFT result. Provide matched observables at each actual stage, not invented intermediate PMNS values.

**5.4 Scalar-coupling strength:** test high vs small `lambdaT3` for fixed Yukawa inputs **and**, separately, after re-optimising Yukawas. Larger `lambdaT3` alone increases the direct source, but compensating Yukawas may cancel that increase. High-scalar `0.5` conclusions are **pending** until fit outputs exist.

**5.5 Representation dependence:** compare the model-family variation at shared inputs against variation along the rescaling directions. Control for refitting to targets. Put claims about which is dominant **after** this comparison.

**5.6 Why mixing angles barely move:** begin with `fig04` showing symmetrisation plus fitted cancellation, `fig05` showing matrices, `fig06` showing heavy-state cancellation, `fig07` showing matching contributions. Then contrast `fig08` (~8.08% change in each splitting over 1 TeV -> 100 GeV) with `fig09` (largest angle change ~8.3e-5 degrees). Keep these final two figures explicitly restricted to the **SM+Weinberg** interval. The scalar-EFT direct term preserves flavour structure at leading fixed-one-loop order, but the full matched C5 need not be exactly proportional to C12.

### Chapter 6 — Conclusions and outlook

- What the benchmark *demonstrates*: matching-scale flavour tensor from Yukawas, coupling degeneracies at leading order, direct Wilson mixing structure, measured Weinberg-EFT differences between masses and mixing.
- What still depends on study completion: strong/weak scalar comparison, full-UV effects, representation-dominance claim.
- Limitations: real-Yukawa reference fits, common fermion threshold rather than sequential thresholds, perturbative order, operator-basis closure of exploratory coupled solver, fit degeneracy and parameter selection.

### Appendices

A. Representation tables and group-theoretical coefficients. B. Wilson-coefficient conventions and Matchete consistency checks. C. Full fit summary (16 models). D. Detailed numerical/solver settings, source hashes, generated tables. E. Optional exploratory C12/C5-coupled investigation clearly not labelled production NLO.

## Interpretation: the newly computed flavour result

For the UV-frozen T3-B(-1) fitted benchmark, in the real orthogonal eigenbasis of the full matched `C5`:

- `||offdiag(U^T A U)|| / ||A|| = 78.61%`, where `A = Y1* diag(1/M) Y2†` is **unsymmetrised**.
- `||offdiag(U^T C12 U)|| / ||C12|| = 1.943%`, for `C12 = (A+A^T)/2`.
- `sum_r ||F_r||/||sum_r F_r|| = 4.598` shows sizeable heavy-state cancellations.
- Across one saved model per class, the symmetrised fractions range from **1.28% to 4.11%**; this is not a theorem of Majorana neutrinos, and fits were constrained to reproduce oscillation data.
- In the 1 TeV -> 100 GeV **Weinberg EFT**, `Delta m^2` values change by **-8.08%** and the largest `theta_ij` shifts by about `8.3e-5` degrees. This numerical comparison does **not** establish a physical intermediate-scale PMNS mixing-angle curve in the scalar-only EFT.

Interpretation: the *fitted* symmetric Yukawa flavour tensor is already well aligned with the complete scalar-threshold mass matrix. The direct `C12 -> C5` source transfers this structure with a scalar prefactor. The observed small angle drift is therefore consistent with approximate eigenvector alignment, **not** a consequence of symmetry alone.

## Minimal remaining numerical work (in priority order)

1. **Complete the `lambdaT3=0.5` UV-frozen fitted comparison** for representative models and collect actual saved outputs. This tests whether rescaling equivalence survives the change in scalar coupling.
2. **Analyse existing controlled UV-on/off sensitivity runs** at matching-degenerate inputs: tabulate relative changes in mass splittings and angles; only then finalise the statement about UV breaking of degeneracy.
3. **Evaluate a shared-input representation comparison**, alongside optimised fits, to quantify whether representation-dependent deviations dominate over the tested coupling variations.

Do **not** expand to sequential fermion decoupling, speculative neutrino-sector NLO running, or full complex-Yukawa reconstruction without a separate project objective.
