# Implementation status — current code audit (10 November 2026)

**Current operational status:** The repository implements substantial general T3 Lagrangian, one-loop renormalisable RGE, Matchete matching, Wilson-evolution, two-scalar threshold, and Yukawa fitting machinery. **It does not currently provide one executable, fully general-Lagrangian, complete one-loop, fermion-first -> scalar-second numerical prediction-and-fit pipeline for A(-4) and D(-2).** This is an **implementation/connectivity** finding, not a demand for experimental or independent physics verification.

**Source provenance:** Read-only GitHub `main`, commit `480262541959ff2c8ed2f3392bebfaf9c971ff64`, inspected 10 November 2026. The separately supplied local pilot v3 ZIP (`scripts/RunT3ThesisTwoModelPilot.py`) was also inspected. The user's uncommitted Windows checkout and generated Wolfram outputs were not accessible for this audit. No full Wolfram matching, EFT rerun or Yukawa optimisation was executed for these notes. Older study results and test counts below remain historical, not new runs.

## Target and conventions

```text
UV (SM + F1,F2,F3 + S1,S2; general permitted interactions)
   -> integrate out ALL heavy fermions at MF = 1e5 GeV
EFT1 (SM + BOTH scalar multiplets + induced Wilson operators)
   -> run renormalisable couplings and appropriate EFT1 Wilson coefficients
   -> integrate out BOTH scalar multiplets together at MS = 1e3 GeV
SMEFT (SM + Weinberg operator)
   -> run to mu_low = 100 GeV
   -> C5, m_nu = -v^2 C5/2, Takagi masses, PMNS mixing
```

For the thesis pilot, `muUV=1e7 GeV`, `lambdaT3(UV)=0.1`, ordinary source scalar-coupling benchmark `0.01`, and the thesis Appendix B real Yukawas are initial inputs. The old class-E factor `1/2` is **historical to certain frozen-scalar study inputs**; it is not to be applied to the new thesis benchmark. All scales and starting inputs must be read from the individual run's config/manifest; older completed Batch 8 data use a different benchmark and (by default) `mu_low=173 GeV`.

## Actual implementations and what they do

| Source | Implemented operation | Boundary / caveat |
|---|---|---|
| `Lagrangian/model/LagrangianBuilder.wl`, `Lagrangian/interactions/GeneralScalarOperators.wl`, `Lagrangian/RunFullGeneralScalarBuild.wl` | Opt-in construction and source validation of symmetry-allowed general fermion/scalar interactions | A symbolic Lagrangian is not itself matching, running, or a numerical C5 |
| `RGE/running/rgbeta/T3RGBetaModel.wl`, `RGE/running/rgbeta/RunGeneralT3QuadraticRGBeta.wl`, `scripts/RunT3CompleteGeneralScalarSymbolic.py` | One-loop UV/EFT1 renormalisable beta exports, including extra scalar quartics and applicable quadratic masses | Extended beta export is experimental; does not establish a Matchete/RGBeta mass-basis dictionary or Wilson evolution |
| `scripts/RunT3FullScalarContinuous.py` | Candidate continuous UV -> fermion-threshold -> EFT1 running of general scalar quartics and quadratics | Its output sets `fermion_threshold_wilson_matching_performed=False`, `scalar_threshold_wilson_matching_performed=False` and `full_sequential_neutrino_mass_prediction=False` |
| `Lagrangian/RunMatching.wl`, `Lagrangian/RunGeneralPhysicalMatching.wl`, `scripts/RunT3GeneralMatcheteMatching.py` | Matchete matching at `LoopOrder -> 1`; general one-scale branches and selected self-bilinear mass bases | This does NOT hand off a full general fermion-first, run-EFT1, scalar-second numerical Wilson trajectory |
| `Lagrangian/RunGeneralT3FermionFirstEFT1OneLoop.wl`, `scripts/RunGeneralT3FermionFirstEFT1OneLoop.py` | Genuine F-only tree and one-loop Matchete matched-EFT-difference export with both scalars light | The Python driver restricts this experiment to B(-1)/C(-1); raw LLHH is unrenormalised and not a physical C5 |
| `Numerical/parallel_t3/background/runner.py` | Existing one-loop candidate UV/EFT1 RGBeta numerical running, with Matchete-sourced F-tree LLSS | Benchmark override is normally development-specific; the thesis pilot injects new values at runtime |
| `Numerical/parallel_t3/wilson/fermion_hard_triangle_batch1.py`, `scripts/RunT3FHardTriangleBatch1.py` | Source-resolved, all-heavy-fermion analytic one-loop triangle contribution | Only an included hard topology; Stage 1 calls `builder.Options(..., reduced_potential=True)` |
| `Numerical/parallel_t3/rough_inputs/builder.py` | Native scalar-coupling to canonical quartic-tensor conversion | In reduced mode, unmapped nonzero scalar quartics are omitted; this is not the full general potential |
| `Numerical/parallel_t3/wilson/strict_one_loop_batch2.py`, `scripts/RunT3StrictOneLoopBatch2.py` | First-order EFT1 full `psi2phi2` Wilson tensor using a candidate seven-channel ADM | Candidate ADM; full one-loop UV Wilson boundary and all operators are not established |
| `Numerical/parallel_t3/threshold/matching.py`, `Numerical/parallel_t3/threshold/pipeline.py` | One common threshold for both scalars, computed LLSS-quartic bubble, then SM+Weinberg running and observables | Does not include all possible one-loop threshold topologies, renormalisation/EOM terms and general-coupling mass-basis effects |
| `scripts/RunT3PoleClosureBatch7.py`, `Numerical/parallel_t3/wilson/native_pole_closure_batch7.py`, `Numerical/parallel_t3/wilson/final_repro_batch8.py` | Source-preserving pole/counterterm diagnostic and consolidated three-prescription numerical output across 16 models | Explicit source flags reject complete physical matching; scenarios are conditional |
| `scripts/RunT3ReducedScalarSobolComparison.py` | Historical 18-real-dimensional Sobol + bounded SciPy least-squares fitting | Uses historical restricted matching/evaluator, not the new full-general EFT |
| `scripts/RunT3ThesisTwoModelPilot.py` (local v3) | Sobol 18D initial sampling, limited-parameter least-squares, re-evaluates Workstream B -> analytic hard triangle -> strict EFT1 -> scalar bubble -> SMEFT for A(-4)/D(-2) | **Restricted conditional numerical evaluator**; not a full general one-loop match |

## A(-4) versus D(-2) general interaction content

From `RGE/general/generated/general_rgbeta_quartic_registry.json` at the audited commit:

| Model | Additional scalar quartic slots | Additional scalar quadratic/bilinear slots | Extra general Yukawa interaction |
|---|---:|---:|---|
| `A_-4` | **0** | **0** | None from the opt-in `y1Conj/y2Conj` criterion |
| `D_-2` | **5** | **1** (neutral S2 self-bilinear) | `y2Conj` is allowed in the opt-in general Yukawa sector |

Thus A(-4)'s renormalisable **additional-scalar** inventory is simpler, but that alone does not complete its one-loop F threshold. D(-2)'s full general calculation must transport the extra scalar couplings and potentially nonzero bilinear, and handle the mass-basis/tensor mapping at the scalar threshold. Setting a permitted coupling to zero is a legitimate **benchmark boundary choice**, not evidence it has been included in matching for general nonzero values.

## Execution evidence: do not conflate completed runs with complete physics

- Earlier reduced-potential frozen benchmark results: the documented **32/32** individually fitted cases are completed **historical restricted** calculations; they are not general one-loop fits. The class-E half applies only to their saved input convention.
- Workstreams Batches 1–8 have previously reported **16/16** completion and conditional low-energy results with three distinct prescriptions. `analytic_triangle`, `native_formal`, and `native_as_UV_minus_soft_hypothesis` are **not** equivalent physical MSbar determinations. Batch 8 is consolidation/audit of saved sources, not source regeneration of the full general Lagrangian.
- The locally reported D(-2) thesis baseline produced conditional masses approximately `[0.00080619, 0.01508806, 0.08696690] eV` and an approximate five-observable `chi2=85161.4`. It **did not** use the full-general matching pipeline and is not an optimised result. A(-4)'s successful current baseline/fit has not been evidenced in this audit.
- Source-level code audit is **not** a fresh numerical execution. Do not claim a general numerical run has succeeded without its matching/trajectory/output files and model-specific manifests.

## Feasibility for the current two-model request

The user wants fitted Yukawas for `A_-4` and `D_-2` **from the full symmetry-allowed Lagrangian, all relevant one-loop F and S matching, candidate one-loop RG evolution across EFT1 and SMEFT, and final mass/mixing observables**. Independent verification is not a prerequisite for preliminary results, but the computational chain itself must execute the requested terms.

**Decision (10 November 2026): NOT READY AS A SINGLE FULL-GENERAL END-TO-END PIPELINE.** Existing full-general beta/RGE sources and Matchete matching experiments do not feed the same full-general two-threshold numerical Wilson solver. Running the present v3 Sobol fit will optimise the restricted analytic-triangle approximation only. Calling it a full one-loop general result would be incorrect, regardless of chi-square.

The shortest physics-development route is: (1) establish the general first-threshold matched Wilson source including its finite/pole treatment for A(-4), then D(-2), (2) hand off to the existing Wilson RGE without discarding relevant general scalar couplings, (3) account for D(-2) extra Yukawa and scalar bilinear mass basis as applicable, (4) implement/assemble the required second-threshold one-loop terms, and (5) wire that evaluator to the existing Sobol/least-squares fitter. Reusing the existing solvers is desirable; the missing matching/handoff is **not** resolved by rerunning tests or increasing optimisation iterations.

See [two-model code-only feasibility audit](T3_A4_D2_GENERAL_PIPELINE_FEASIBILITY_2026_11_10.md) for the precise go/no-go decision. Older [numerical methods](NUMERICAL.md), [results](RESULTS.md), and report READMEs document **their dated historical studies**, and should not be interpreted as describing a working full-general version of the new two-model pilot.
