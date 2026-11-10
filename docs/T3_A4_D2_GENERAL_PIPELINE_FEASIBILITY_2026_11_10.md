# A(-4) / D(-2) general-Lagrangian one-loop pipeline: code audit

**Audit date:** 10 November 2026. **Repository:** `Frankyang2002/T3neutrinomassmodelpipeline`, GitHub `main` commit `478c658c70a397cff7f04d397e6f7dd448179d01` (current `main` at follow-up audit); additionally inspected the locally supplied v3 thesis-pilot script. **No Wolfram executable, locally generated Matchete artifact, or complete user Windows checkout was available to rerun; claims below describe source implementation and the user's reported D(-2) baseline only.** No GitHub modifications.

## Requested deliverables and present feasibility

1. Two sets of **optimised UV Yukawas** for `A_-4` and `D_-2` using the **full general symmetry-allowed Lagrangian**: **not available today** from a completed full-general sequential EFT evaluator. The existing Sobol optimiser is not the blocker.
2. For each fitted point, **one-loop F matching, EFT1 running, simultaneous one-loop two-scalar matching, SMEFT running to 100 GeV, and neutrino masses/mixing**: **not available today** with all requested general interactions. The existing pipeline performs a **restricted** conditional subset and can produce approximate preliminary fits if explicitly labelled as such.
3. **Verification or certification**: not required by the present request. Removing the verification requirement does not create missing Wilson matching terms or connect currently separate symbolic and numerical branches.

## Exact code relationships

```text
GENERAL LAGRANGIAN / RGBeta side                     EXISTING FITTED NUMERICAL side
-------------------------------------------------   --------------------------------
Lagrangian/RunFullGeneralScalarBuild.wl               local scripts/RunT3ThesisTwoModelPilot.py
scripts/RunT3CompleteGeneralScalarSymbolic.py           -> Numerical/parallel_t3/background/runner.py
  -> UV/EFT1 general quartic & quadratic betas         -> scripts/RunT3FHardTriangleBatch1.py
scripts/RunT3FullScalarContinuous.py                     -> builder.Options(reduced_potential=True)
  -> candidate coupled renormalisable trajectories      -> analytic hard fermion triangle ONLY
  -> NO F or S Wilson-threshold matching                -> scripts/RunT3StrictOneLoopBatch2.py
                                                          -> strict candidate EFT1 Wilson
Separate full-general Matchete branches                  -> joint S scalar bubble ONLY
Lagrangian/RunGeneralPhysicalMatching.wl                 -> SMEFT RGE -> masses/PMNS -> chi2
  -> one-scale general C5 (or special self-bilinear)
Lagrangian/RunGeneralT3FermionFirstEFT1OneLoop.wl
  -> raw one-loop first-threshold Matchete output
  -> Python caller currently B_-1/C_-1 ONLY
  -> NO numerical full-general A_-4 / D_-2 bridge
```

**Important existing implementation:** `Lagrangian/RunMatching.wl` already defines `RunSequentialT3Matching`, with separate tree and one-loop stage matching, inherited one-loop terms propagated through tree matching at the next stage, and fresh-kernel threshold handling. Therefore symbolic two-stage matching need not be reinvented. The gap is full-general A(-4)/D(-2) invocation, running-scale insertion, operator coefficient extraction, general mass/tensor mapping, and numerical handoff.

**No source entry point was found that connects these general-Lagrangian branches to full general fermion-first matching, general EFT1 coefficient evolution, general two-scalar threshold, and the two-model optimizer in one numerical evaluation.** This is a claim about inspected code paths, not proof no unshared local work exists.

## Detailed blockers, by model

### `A_-4`

- From `RGE/general/generated/general_rgbeta_quartic_registry.json`: **0 additional quartics and 0 additional scalar bilinears**. The opt-in extra-fermion Yukawa criterion in `T3RGBetaModel.wl` selects only `yF=+/-1/2`, so A(-4) (`yF=-3/2`) has no `y1Conj/y2Conj` extension.
- This makes the **renormalisable general scalar inventory relatively simple** and may permit reuse of much of the existing tensor builder at the thesis benchmark.
- **Still missing** for the requested full-one-loop claim: a connected, all-relevant-terms F-only Matchete one-loop hard Wilson boundary, proper associated UV/EFT1 soft subtraction and counterterms, complete finite second-threshold matching, and handoff to the existing fitter. A calculated analytic triangle is not the sum over all one-loop F matching topologies.

### `D_-2`

- From the same registry: **5 additional quartic slots** and **1 allowed extra quadratic self-bilinear** for the neutral scalar S2. `T3RGBetaModel.wl` also permits **`y2Conj`** in the opt-in general Yukawa basis (`yF=-1/2`). The v3 thesis pilot sets extra conjugate Yukawa coefficients and permitted bilinears to **zero at its UV input**, and the matching tensor builder runs in reduced-potential mode.
- The full general version needs a candidate scalar tensor/bilinear mass-basis handoff that does not silently discard allowed nonzero terms, plus the first- and second-threshold one-loop contributions. Zero input is permitted, but running-induced terms and any deliberate nonzero general benchmark must be treated consistently.
- `scripts/RunT3GeneralMatcheteMatching.py` supports a **one-scale** self-bilinear mass-basis candidate for D(-2), but not the requested two-threshold Wilson trajectory.

## What is executable now without claiming full general one-loop

- **Restricted thesis-point evaluator:** the user's reported `D_-2` baseline completed the existing UV + EFT1 + analytic hard triangle + scalar bubble + SMEFT calculation. Reported conditional masses: `0.0008061911`, `0.0150880588`, `0.0869668958` eV. Reported approximate five-observable chi-square: `85161.4086`. This is not a fitted or full-general result.
- **Restricted fitting:** `scripts/RunT3ThesisTwoModelPilot.py --fit ...` does Sobol sampling of 18 real Yukawas and limited-parameter bounded SciPy least-squares on that same restricted evaluator. A low chi-square would not change the physical scope.
- **General scalar renormalisable evolution in isolation:** `scripts/RunT3FullScalarContinuous.py` has a `--model` path for the candidate all-coupling scalar RG calculation when the required generated beta artifacts are present. It explicitly does not compute F/S Wilson threshold matching or the full sequential neutrino mass.
- **Symbolic one-scale matching in isolation:** `scripts/RunT3GeneralMatcheteMatching.py --model A_-4` and `--model D_-2` name the two models; those one-scale Matchete outputs are **not** a drop-in replacement for fermion-first + scalar-second matching.
- **16-model Batch 8:** final deterministic *conditional* report across analytic and native/hypothesis prescriptions, using existing earlier saved sources; not full-general input regeneration.

Do not infer new A(-4) fitted Yukawas, fitted D(-2) Yukawas, chi2, PMNS angles, or neutrino masses for a full-general one-loop evaluator: **none were generated by this audit**.

## Minimum next implementation checkpoint

For implementation scope, see [four large batches](T3_A4_D2_FULL_GENERAL_IMPLEMENTATION_BATCHES.md). First choose **A(-4)** because its additional scalar-interaction inventory is empty. Build a *single-point, source-derived* full-general two-threshold `C5` evaluator that consumes the general matching/RGE outputs and visibly includes the required first- and second-threshold contributions. It need not be externally verified for this numerical pilot, but its operator inventory, input definitions, and pole/finite prescriptions must be explicit rather than assuming away missing terms.

Only after that single point runs through F matching -> EFT1 coefficient evolution -> joint scalar matching -> SMEFT should the existing 18-real-Yukawa Sobol/least-squares search be attached. Extend the same contract to D(-2), including its extra scalar interactions and `y2Conj`. This is an actual implementation step, not a documentation toggle, and likely requires new matching/adapter code.

## Direct code references (GitHub, read-only)

- `Lagrangian/model/LagrangianBuilder.wl`, `Lagrangian/interactions/GeneralScalarOperators.wl`, `Lagrangian/RunFullGeneralScalarBuild.wl`
- `RGE/running/rgbeta/T3RGBetaModel.wl`, `RGE/running/rgbeta/RunGeneralT3QuadraticRGBeta.wl`, `RGE/general/generated/general_rgbeta_quartic_registry.json`
- `scripts/RunT3CompleteGeneralScalarSymbolic.py`, `scripts/RunT3FullScalarContinuous.py`
- `Lagrangian/RunMatching.wl`, `Lagrangian/RunGeneralPhysicalMatching.wl`, `scripts/RunT3GeneralMatcheteMatching.py`
- `Lagrangian/RunGeneralT3FermionFirstEFT1OneLoop.wl`, `scripts/RunGeneralT3FermionFirstEFT1OneLoop.py`
- `scripts/RunT3FHardTriangleBatch1.py`, `Numerical/parallel_t3/rough_inputs/builder.py`, `Numerical/parallel_t3/wilson/strict_one_loop_batch2.py`, `Numerical/parallel_t3/threshold/matching.py`
- `Numerical/parallel_t3/wilson/native_pole_closure_batch7.py`, `Numerical/parallel_t3/wilson/final_repro_batch8.py`
- Local v3 `scripts/RunT3ThesisTwoModelPilot.py` and `README_T3_THESIS_TWO_MODEL_PILOT.md` (separately supplied, not asserted to exist on `main`)
