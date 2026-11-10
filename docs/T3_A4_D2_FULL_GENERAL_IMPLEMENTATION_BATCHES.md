# A(-4) and D(-2): full-general one-loop implementation batches

**Status: proposed work, not completed work.** Updated 10 November 2026 following read-only code inspection of GitHub `main` at `478c658c70a397cff7f04d397e6f7dd448179d01`. The user's Windows checkout and generated Wolfram artifacts were not executed for this note.

## Objective and acceptance criterion

Produce preliminary, numerically usable results for **A(-4)** and **D(-2)** using the **general symmetry-allowed T3 renormalisable Lagrangian** and a one-loop, fermion-first two-threshold chain. No independent physical certification is requested. However, a numerical computation cannot be called full-general if the evaluator silently discards general interaction coefficients or replaces an unimplemented one-loop contribution with a known restricted formula.

```text
UV (muUV = 1e7 GeV): SM + F(three flavours) + S1 + S2
  -> run full applicable renormalisable UV RGBeta one-loop system
  -> integrate out ALL fermions together at MF = 1e5 GeV
EFT1: SM + BOTH scalar multiplets + induced Wilson operators
  -> run full applicable renormalisable EFT1 couplings and the relevant Wilson coefficients
  -> integrate out BOTH scalar multiplets at ONE MS = 1e3 GeV
SMEFT: run C5 and SM couplings to muLow = 100 GeV
  -> mnu = -v^2 C5/2, Takagi masses, PMNS mixing and five-observable chi2
```

Use the thesis benchmark for fixed inputs (in particular `lambdaT3=0.1`, no historical class-E adjustment) and vary the two original real 3x3 Yukawa matrices for initial A(-4)/D(-2) fits. D(-2) additionally **allows** `y2Conj` and an S2 holomorphic mass term; whether those are held at zero or included as nonzero initial parameters must be explicit and consistently propagated. No fitted result has yet been produced under this new complete evaluator.

## Existing components to reuse, rather than recreate

- `Lagrangian/model/LagrangianBuilder.wl`, `Lagrangian/RunFullGeneralScalarBuild.wl`: general native Matchete Lagrangian and interaction checks.
- `RGE/running/rgbeta/RunGeneralT3QuadraticRGBeta.wl`, `scripts/RunT3CompleteGeneralScalarSymbolic.py`, `scripts/RunT3FullScalarContinuous.py`: general one-loop renormalisable RG exports and numerical evolution.
- `Lagrangian/RunMatching.wl`: **already contains** `RunSequentialT3Matching`, which matches successive stages, propagates inherited one-loop corrections through later tree matching at order `hbar`, and uses a fresh kernel for subsequent thresholds. This is implemented symbolic matching logic, **not yet a proven integrated general numerical pipeline**.
- `Numerical/parallel_t3/wilson/adm.py`: existing full `psi2phi2` candidate EFT1 Wilson evolution.
- `Numerical/parallel_t3/threshold/matching.py`, `Numerical/parallel_t3/threshold/pipeline.py`: existing joint scalar threshold and SMEFT numerical infrastructure, currently with a limited bubble-based one-loop matching prescription.
- `scripts/RunT3ReducedScalarSobolComparison.py` and the local v3 `scripts/RunT3ThesisTwoModelPilot.py`: existing 18-real-Yukawa Sobol/least-squares search and numerical evaluation caching. The v3 pilot **still computes the restricted analytic triangle**, not the requested general one-loop matched result.

## Batch 1 — general running inputs and numerical interfaces

**Work:** Extend the existing A(-4)/D(-2) runnable input/state contract to consume *all applicable* UV/EFT1 general beta exports, including scalar-mass and scalar-quartic coordinates, with a same-point MF threshold handoff. Use cached symbolic exports rather than launching Wolfram for each Yukawa evaluation. Preserve model-specific operator and tensor labels; include the permitted extra Yukawa when applicable.

**Pass condition:** At one benchmark point per model, output numerical general-coupling trajectories at `muUV`, `MF`, and `MS`, without silently dropping a coupling from the declared model inventory. This does **not** yet constitute a neutrino-mass result.

## Batch 2 — F-only general one-loop matching and EFT1 Wilson RGE

**Work:** Invoke/reuse the sequential Matchete first stage with all fermions heavy and both scalars light. Export the tree and relevant one-loop EFT1 operator terms; translate Matchete's native labels and normalisations to the numerical Wilson tensor. Insert a clear finite/pole subtraction prescription and the first-stage one-loop coefficients. Feed the matched operator coefficients to the existing candidate full-tensor EFT1 RGE. Do not confuse a single analytic triangle with a complete source-derived one-loop boundary.

**Pass condition:** Each model obtains a concrete, source-derived, finite numerical EFT1 Wilson state at `MF`, runs it to `MS`, and records all included and excluded operator terms. First achieve this for A(-4), then D(-2).

## Batch 3 — joint general scalar matching, SMEFT, masses and mixing

**Work:** Propagate the general scalar masses (including D(-2)'s possible self-bilinear split), relevant scalar interactions and Wilsons to the scalar threshold. Reuse Matchete's second-stage machinery and numerical mass-eigenbasis conversion for **both scalars removed at one threshold**. Include the relevant one-loop second-threshold contributions rather than just the existing scalar bubble if claiming full one-loop. Run the SM+Weinberg equation to 100 GeV; compute neutrino observables.

**Pass condition:** One **non-fitted** end-to-end A(-4) and D(-2) benchmark per model yields numerical `C5(MS)`, `C5(100 GeV)`, Takagi masses, mixing, and a machine-readable included-terms ledger. No requirement to demonstrate an independently certified physical coefficient for the user's preliminary scope.

## Batch 4 — full-general Sobol refit and two-model report

**Work:** Plug the completed general evaluator into the existing Sobol/least-squares optimiser; begin from thesis Appendix B Yukawas, cache symbolic expressions and numerical trials, and search the 18 real original Yukawa coordinates. The local least-squares refinement may use a sensitivity-selected subset if it reaches a useful chi-square; otherwise expand to all 18. Keep any permitted additional D(-2) Yukawa input settings fixed and explicit unless a later request expands the fit.

**Pass condition:** Save the best candidate Yukawa matrices, five NuFIT residuals, approximate chi-square, mass eigenvalues, splittings, PMNS angles, and the full-general evaluator/source labels for **both** A(-4) and D(-2). A low chi-square is optimisation success, not independent experimental validation.

## Implementation effort and hard boundary

Approximately **four major batches**; **12–20 files affected** is a rough planning estimate, not an inspected file-change list. The most uncertain work is not the RG ODE or optimiser: it is the native general one-loop **operator-to-numerical matching interface**, including mass-basis treatment and the extraction of second-threshold one-loop contributions. Work should stop and report the specific missing term if the existing symbolic output is insufficient. Do **not** silently substitute the restricted analytic hard triangle or scalar-only bubble and label it full-general.

Historical frozen-scalar or one-scale results remain separate datasets. Nothing in this plan deletes historical material or modifies GitHub.
