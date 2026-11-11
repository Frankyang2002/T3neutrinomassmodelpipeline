# T3 Thesis — coupled-IVP physics completion handover (11 November 2026)

## Objective

Use **only** the existing fully coupled IVP implementation for the final Chapter 5 calculation. First resolve two physics gaps: **(1) fermion-threshold one-loop matching**, and **(2) complete, consistently normalised scalar interactions and joint scalar matching**. Do not start a new fit campaign or rewrite the thesis until those changes are tested on A(-4), then D(-2), then extended across all 16 models. Prefer a concrete working benchmark over broad new scope.

**Local project:** `C:\thesistemp\T3neutrinomassmodelpipeline`  
**Read-only remote:** `https://github.com/Frankyang2002/T3neutrinomassmodelpipeline`  
**Latest checked public main:** `35584362f537584d76100f9e44df0322202a5443` (11 Nov 2026, 11:35:46 UTC). Verify the latest source again before making changes. Never write to GitHub.

**Thesis reference:** uploaded `T3_Honours_Thesis (37).pdf`. Chapter 3 currently describes simultaneous one-loop renormalisable/Wilson RGEs and a finite two-scalar threshold; Chapter 4 uses UV=1e7 GeV, F=1e5 GeV, S=1e3 GeV, low=100 GeV; Chapter 5 compares representations, stage-freeze running, and coupling sensitivity.

## Required EFT ordering — never change

UV (SM + all three heavy fermions + both scalar multiplets) 
→ integrate out **ALL fermions at once** at `MF = 1e5 GeV`
→ EFT1 (SM + **BOTH** BSM scalar multiplets + retained dimension-five operators), with coupled running
→ integrate out **BOTH scalar multiplets at once** at `MS = 1e3 GeV`
→ SMEFT (SM + Weinberg operator) → `mu_low = 100 GeV` → `mν = -v² C5 / 2` and Takagi observables.

No scalar-only intermediate EFT. 

## Selected production numerical implementation

- `scripts/RunT3SimultaneousEFT1A4.py` — A(-4).
- `scripts/RunT3SimultaneousEFT1All16.py` — 13 other nonmixed models.
- `scripts/RunT3SimultaneousEFT1MixedBC.py` — B(-1), C(-1) with fixed original Cartesian scalar basis.
- `scripts/RunT3SimultaneousEFT1StageFreezeHeatmaps.py` — reports from existing saved cases.
- Core: `RGE/running/intermediate/GeneralEFT1C5Mixing.py`; `Numerical/running/SMWeinbergEvolution.py`.

The EFT1 integrator is one `solve_ivp(..., method='DOP853')` over state `[renormalisable RGBeta variables, Re C, Im C]`. Both beta functions use the **same instantaneous** running state:

\[
\frac{d}{d\ln\mu}\binom{r}{C} = \frac1{16\pi^2}\binom{\beta_r^{(1)}(r)}{\beta_C^{(1)}(C,r)}.
\]

Its retained dimension-five basis is the nine-leptonic-Weyl `ψ²φ²` tensor, **not** the complete EFT operator inventory. The current simultaneous scripts explicitly mark hard fermion one-loop matching absent and non-HHSS scalar tensors/normalisations provisional. **Do not replace the coupled IVP with the strict-one-loop fixed-order `MF`-frozen implementation or new experimental Batch 17 branch.** Historical conditional results are not final coupled-IVP results.

## Saved coupled-IVP files — local existence confirmed by user

1. `output/studies/general_scalar_extended/experimental_simultaneous_eft1_seven_term/simultaneous_A_m4_v2/A_m4/cases.json` (1 model)
2. `output/studies/general_scalar_extended/experimental_simultaneous_eft1_general_nonmixed/simultaneous_nonmixed13_v1/<model>/cases.json` (13 models)
3. `output/studies/general_scalar_extended/experimental_simultaneous_eft1_mixed_BC/simultaneous_mixed_BC_v1/<model>/cases.json` (2 models)

16/16 files **exist**. Contents, solver status, finite values, source provenance and matching assumptions **have not yet been validated**. Existing source/result files must not be overwritten. Confirm `1111`, `0111`, `1001`, `1110` for every model. Bit order `U E W S` = UV renormalisable, EFT1 renormalisable, EFT1 Wilson, SMEFT. `1001` freezes BOTH EFT1 components. Finite scalar bubble is retained in all cases.

## Priority A: complete relevant fermion-threshold Wilson matching

**Problem:** Coupled IVP uses primarily the source-resolved tree LLSS boundary; independent *renormalised* complete hard one-loop fermion threshold is not fed into it. Matchete source output is not automatically the finite physical Wilson boundary. Critical conversion `GammaCC`/ordered spinor → physical Weyl tensor sign and symmetry normalisation remains uncertified.

**Implement and verify:**

1. Inventory actual **full** Matchete source outputs at `MF` for A(-4) first, D(-2) second. Classify relevant tree and one-loop operators by gauge invariant and lepton-number-violating structures, including LLSS and LLHH, with source IDs; do not silently discard unclassified sectors.
2. Derive or establish **signed spinor/operator and scalar component maps** into the coupled-IVP `C[i,j,a,b]` convention. Keep uncertainty explicit; do not assume factor `1` is verified.
3. Form **finite renormalised** one-loop matching with full UV-hard minus EFT1-soft contribution and appropriate counterterms in one consistent scheme. Do **not** directly inject raw Matchete poles or raw hard amplitudes into boundary.
4. Decide whether any generated operator classes outside retained `ψ²φ²` have nonzero contribution to `C5` at the claimed order; if they do, extend the retained transport only as necessary, or state an explicit truncation. Check no double counting with RGE logarithms or scalar threshold.
5. Feed corrected Wilson boundary to the **same simultaneous** `DOP853` solver and compare tree-only vs matched results at the identical benchmark.

**Acceptance:** exact MF tensor basis and symmetries; finite and reproducible renormalised boundary; dimension `GeV^-1`; source traceability; demonstrable single inclusion of one-loop effects; finite symmetric `C5(muS)` and low-energy observables. Distinguish numerical completion from independent proof of physical correctness.

Relevant existing reference branches are *sources to reuse or audit*, not alternative production solvers: `Lagrangian/RunGeneralT3FermionFirstEFT1OneLoop.wl`, `Numerical/parallel_t3/wilson/fermion_hard_triangle_batch1.py`, `Numerical/parallel_t3/complete/two_model_hard_general_batch10.py`, `Numerical/joint_scalar_matching/fermion_boundary_bridge.py`. Inspect current paths/versions before reusing them.

## Priority B: scalar potential, tensor dictionaries and joint threshold

**Problem:** The selected simultaneous runner supplies provisional non-HHSS quartic tensor ansätze, and the native Matchete → RGBeta invariant and component conventions are not fully identified. Conditional scalar matching includes a finite joint bubble but is not demonstrated to be the complete one-loop scalar threshold.

**Implement and verify:**

1. For each selected model, identify allowed independent scalar quartics and quadratic/bilinear terms in the current general scalar registry and their actual RGBeta beta functions. Avoid treating every scalar coupling as a universal single `0.01` tensor coefficient.
2. Map the same potential into the numerical real-scalar `λ_abcd` tensor and `HHSS` projection, preserving all signs, factorial/symmetry factors, CG bases and complex phase conventions. Check using full polynomial/tensor reconstruction, not only one chosen component.
3. Update the coupled-IVP quartic context to use the **running, source-consistent** scalar tensor at each RHS evaluation, without changing solver structure.
4. For self-bilinear and mixed models, form the complete scalar mass Hessian at the endpoint, use consistent field rotations. For B(-1)/C(-1), use original Cartesian scalar basis for EFT1 running and transform appropriately at matching; do not freeze the running mass eigenbasis.
5. Evaluate the **joint** scalar-threshold finite terms and determine whether additional operator channels/topologies must be included at the selected accuracy. Double counting between `C_LLHH(muS)` and threshold finite terms is forbidden.

**Acceptance:** quartic tensor symmetry and gauge invariance; named native↔RGBeta coordinate map (with explicit unresolved cases); no omitted nonzero allowed couplings in the chosen benchmark; matching mass spectrum agrees with evolved scalar masses; both scalars integrated out in one threshold; finite symmetric low-energy `C5`.

Relevant existing source: `RGE/general/generated/general_rgbeta_quartic_registry.json`, `Numerical/joint_scalar_matching/`, `Numerical/lo_certification/physical_normalisation/`, `Numerical/parallel_t3/complete/two_model_scalar_coordinate_batch16.py`. Inspect latest files; Batch 16 is *an audit*, not physical certification.

## Execution order and time management

1. **Archive historical branches locally** using `scripts/ArchiveLegacyT3Pipeline.ps1` (provided separately), keeping current coupled scripts/outputs untouched. This is a selective archival snapshot; it is not a stand-alone restored repository. Check the generated ZIP and `ARCHIVE_MANIFEST.json` before cleanup; **do not delete anything automatically**.
2. Validate the saved coupled-IVP `cases.json` files, including input provenance and baseline `1111` before changes; preserve these outputs for before/after comparison.
3. Implement Priority A on **A(-4)**. Establish a reproducible baseline and corrections. Then test D(-2).
4. Implement Priority B for the same two models and verify low-energy `C5`, mass-scale and group/phase consistency.
5. Extend the validated maps/physics inputs to the remaining models and the mixed B/C branch; do not assume an A/D mapping is universal.
6. Recompute/refresh the 16-model coupled-IVP outputs and fitted input compatibility. Then produce Chapter 5 §5.1 representation transfer, §5.2 four-bit stage-freeze heatmaps, §5.3 T1/T2/T3 compensated sensitivity.

Do **not** claim the whole calculation is physically complete merely because its tests, model count or `chi²` pass. For an honours rough-results phase, keep a documented conditional approximation if the hard matching proves disproportionate, but never silently change approximation mid-campaign.

## User's code workflow

- Work **only locally**, never alter GitHub.
- Read the current GitHub implementation before proposing changes; do not guess local uncommitted content.
- Changes require **FULL replacement files**, not patches; when multiple files are changed, provide **one ZIP** preserving project-relative paths.
- Give intended project paths and exact PowerShell commands.
- Do not overwrite historical numerical results. Prefer new run IDs.
- Focus on the two specified physics gaps, preserve existing conventions unless a bug is verified, avoid scope expansion.
- Distinguish code verified, result verified, inferred, and merely proposed.

## First action in new conversation

Request or inspect local `cases.json` output data and confirm all 16 `1111` results. Then identify the exact Matchete fermion-source export consumed by the A(-4) simultaneous runner and audit its Wilson boundary. **Do not refit Yukawas or rerun all 16 models yet.**
