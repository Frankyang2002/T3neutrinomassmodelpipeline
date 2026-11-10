# Workstream A — Wilson physics / normalization (conditional implementation)

<!-- T3_CURRENT_SCOPE_BEGIN: 2026-11-10 -->

**Document scope:** Current pilot uses the candidate Wilson ADM; no complete general one-loop matching input is connected.

**Current implementation status (10 November 2026):** The requested A(-4)/D(-2) **full-general, complete-one-loop, fermion-first -> joint-scalar -> SMEFT Yukawa fit is not yet an executable end-to-end numerical workflow**. The repository already implements general Matchete/RGBeta components and an actual symbolic `RunSequentialT3Matching` stage engine, but its general Wilson matching outputs are not fully connected to the two-model numerical fit. The current thesis Sobol pilot is **restricted/conditional** (analytic fermion triangle and selected scalar bubble), even if it reaches a low chi-square. See [implementation status](../../../docs/IMPLEMENTATION_STATUS.md), [A(-4)/D(-2) feasibility](../../../docs/T3_A4_D2_GENERAL_PIPELINE_FEASIBILITY_2026_11_10.md) and [implementation batches](../../../docs/T3_A4_D2_FULL_GENERAL_IMPLEMENTATION_BATCHES.md). Earlier numerical results and verification claims below refer to their **dated, historical configurations** unless an individual run manifest establishes otherwise.

<!-- T3_CURRENT_SCOPE_END -->
**Implementation target:** `Numerical/parallel_t3/wilson/` in `C:\thesistemp\T3neutrinomassmodelpipeline`. No GitHub modifications. Public `main` inspected before this batch: `a944e0030b8ef24b23342d99e0d8f33cef3cffd5` (2026-11-09). Local uncommitted files are not known. The B and C interfaces in the uploaded handovers are *proposed*, not yet established implementation interfaces. The B->A input schema here is explicit and fail-closed.

## Tensor convention

EFT1 has SM Weyl fermions `Q,uc,dc,L,ec`, three families (45 complex left-Weyl components), and `H,S1,S2` with real scalar axis length `ns=4+2*dS1+2*dS2`. The scalar real order is `H[1,R],H[1,I],H[2,R],H[2,I],S1[1,R],S1[1,I],...,S2[dS2,I]`. A complex scalar component is `z=(R+i I)/sqrt(2)`. The SU(2) components of each irrep are highest-weight to lowest-weight. All tensors use zero-based Python indices and named metadata; no positional aliases to Matchete.

```
L5 = C[i,j,a,b] psi[i] psi[j] phi[a] phi[b] + h.c.   (no factorial)
V4 = lambda4[a,b,c,d] phi[a] phi[b] phi[c] phi[d]/4!
L5^W = 1/2*C5[p,q]*(L[p].H)*(L[q].H) + h.c.
C[nu_p,nu_q,H[2,R],H[2,R]] = C5[p,q]/4
```

`C` has shape `(45,45,ns,ns)`, `complex128` and is separately symmetric under `i<->j` and `a<->b`. `C5` is `complex128` symmetric `(3,3)`; it is embedded as the *entire* Weinberg gauge-singlet polynomial (charged and imaginary Higgs components too). `m_nu = -(v^2/2)*C5` is a **convention only** here, not a neutrino-mass result: scalar matching is not performed by A.

`ModelBasis.generators()` uses the A(-4) source's SU(2), U(1) and SU(3) left-Weyl, real-scalar generators, extended to each irrep. The U(1) generator is the *non-GUT-normalized* `gY` convention. `SU3` is the fundamental for Q and antifundamental for `uc,dc`. The two scalars stay until the same MS. No heavy F coordinate occurs in EFT1.

## Dimension-five closure assumption and ADM status

For the stated unbroken Z2 and SM chiral Weyl content, the **physical** non-`psi2phi2` dimension-five classes are excluded by representation counting:

- `phi^5`: Z2 parity requires an even number of Z2-odd heavy scalars, hence an odd number of Higgs doublets; A/D/E heavy representations have integer spin, and B/C heavy doublets occur in even number, so the total number of SU(2) half-integer factors is odd and cannot make a singlet.
- `phi F^2`: a single Z2-even gauge-singlet scalar is absent; H is an SU(2) doublet.
- `psi^2 F`: no pair of unbroken-SM left-Weyl fermions makes a gauge-singlet bilinear without inserting Higgs fields. Heavy F is absent in EFT1.

These are analytical representation checks, **not** a native automatic Hilbert-series enumeration. EOM-redundant off-shell operators are not independent physical running coordinates. The full `psi2phi2` tensor has all fermion-sector components available and does not assume a projected LLSS+LLHH subsystem is closed. Nonzero C5(MF), its self-running, LLSS->C5, reverse C5->LLSS, and all other same-class tensor mixings use one beta and one solver.

Three explicit candidate ADMs:

| `adm_variant` | Status |
|---|---|
| `fonseca_eq22_six_term` | **Default conditional candidate**; drops the crossed Yukawa graph; six channels agree with independent *gaugeless* Fonseca–Olgoso–Santiago Eq. (22) as tested in existing repository. Full gauge completion not physically certified. |
| `ward_covariant_candidate` | C-conjugated crossed Yukawa plus transposed fermion Casimir; passes general complex Weyl basis covariance, but crossed loop coefficient not certified. |
| `source_v1_unstarred` | Faithful arXiv:2512.16890v1 Eq. (4.85) crossed term; known to violate general complex Weyl basis covariance. Audit mode only, never a physical choice. |

For a left-Weyl change `psi_old=U psi_new`, the Wilson tensor transforms as `C_new=U U C_old`. The unstarred crossed contribution `y*y*C` does not have this transformation law; `y*y*C*` does. The fermion Casimir must be oriented as `C2(F)^T` for these coefficient indices. Tests directly show the unstarred variant fails nontrivial unitary transformations, whereas the other two pass. **Covariance does not prove the crossed one-loop counterterm coefficient**, so no ADM is physically certified. Source comparisons: GitHub `GeneralEFT1SevenTermFlow.py`, `GeneralEFT1CrossedWardCandidate.py`, `GeneralEFT1SixTermReferenceCandidate.py`, `GeneralEFT1IndependentOneLoopReference.py`, `GeneralEFT1OneLoopTopologyAudit.py`; arXiv:2512.16890v1 Eq. (4.85); Fonseca et al. JHEP07(2026)094 Eq. (22). JHEP09(2026)159 is the published Aebischer et al. version; its equation has not been independently reconciled in this batch.

The code calls the **existing** matrix-free repository beta implementations when installed, with a standalone algebraic mirror for isolated tests. Every output records which backend was actually executed. No published formula is certified by the successful ODE solution.

## Workstream B -> Workstream A file interface

**Both manifests** are JSON dictionaries containing:

```json
{
  "model": "A_-4",
  "basis_id": "T3_RGBeta_reference_Weyl45_real_interleaved_v1",
  "fermion_labels": ["the complete exact 45-label array from build_basis(model)"],
  "scalar_labels": ["the complete exact scalar-label array from build_basis(model)"],
  "full_SM_three_generation_Weyl_basis": true,
  "MF_gev": 100000.0,
  "MS_gev": 1000.0,
  "npz_sha256": "hex SHA256 of corresponding NPZ",
  "source_sha256": {"project/relative/source_file": "hex SHA256 of source file"}
}
```

**`background.npz`:** `mu_gev[N]`, `Y[N,45,45,ns]`, `lambda4[N,ns,ns,ns,ns]`, `g2[N]`, `gY[N]`, `g3[N]`, with `mu` endpoints exactly MF and MS (either order). Complex Yukawa and quartic tensors are **canonical and normalised**, not raw RGBeta parameter vectors. Interpolation is linear in `ln(mu)` for actual renormalisable tensor samples; `g2/gY/g3` are interpolated similarly. Optional arrays `generator_SU2_fermion[3,45,45]`, `generator_SU2_scalar[3,ns,ns]`, similarly `U1` and `SU3`, are checked against the reference generators. If B's field coordinates differ, use the explicit unitary/orthogonal transform helper `transforms.py` **and establish its provenance first**; label mismatch otherwise blocks.

**`boundary.npz`:** Either `C_full[45,45,ns,ns]` (set `C_full_is_canonical:true` in boundary manifest), or `C_LLSS_formal_muF[3,3,2,2,ns-4,ns-4]` plus `heavy_source_from_canonical[ns-4,ns-4]` real orthogonal. Optional `C5_hard[3,3]` can be nonzero and is *added* as the full Weinberg polynomial. Do **not** include the hard contribution twice in `C_full` and `C5_hard`. Boundary manifest must set `operator_basis:"psi2phi2_no_factorial"`. For formal Matchete LLSS specify `gamma_cc_to_weyl_signed_factor` (numeric factor). This factor and scalar map remain user-provided **assertions, not established physical normalization**. `--assume-formal-identity` explicitly assumes missing heavy map identity and missing spinor factor +1, marked conditional. It does not bypass source hashes.

A(-4) additionally needs `native_source_path` and `corrected_Aminus4_portal_evidence:true` in the B boundary manifest. Old sources from `fermion_first_all16_v1` or `eft1_stageB_saved_a_m4` are rejected. This does not guarantee native regeneration: B must supply source hashes and fresh evidence. A checks that the current `ScalarPotential.wl` retains `T3HiggsPortalHermitianPhase`, the `I` branch, and `phase` **inside** `PlusHc`. No A edits of this file are included.

Source hashes are matched against the local repo before a real-input run. The current GitHub revision is **not** assumed to equal the user's uncommitted working copy. Every `.npz` file is read with `allow_pickle=False`; all axes and symmetry conditions are checked. Unknown forms/fields/normalizations block rather than being silently set to zero.

**Outputs:** `<out>/<model encoded as A_m4>/canonical_wilson_MF.npz`, `wilson_at_MS.npz`, `manifest.json`. At MF `C_full`, `C5_EFT1_MF`, `beta_C_at_MF_16pi2`; at MS `C_full`, `C_LLSS_at_MS` (the full 6-lepton-weak-index x heavy-real projection), `C5_EFT1_MS`, `C5_EFT1_MF`, `beta_C_at_MS_16pi2` plus a *diagnostic only* linearized leading log. The full tensor `C_full` and `C5_EFT1_MS` are for C's **one joint scalar threshold**, and C must **not** add the direct leading log a second time. The manifest includes exact bases, symmetry/Weinberg projection residuals, ADM choice and individual output SHA256. No neutrino masses and no SMEFT C5 are produced by A.

## All-16 conditional numerical fixture

`synthetic_invariants.py` solves the exact SU(2) generator Ward nullspace for
`L L S1† S2` and `H H S1 S2† + h.c.` with the representation dimensions of
**each of the 16 models**. This is independent of any unverified Matchete/RGBeta
Clebsch normalization: it is purely a fixed highest-weight reference-frame
polynomial with arbitrary synthetic strengths, NOT a physical model coupling.
The source has three complex lepton flavours and is run together with an
independent nonzero Weinberg hard boundary. The corresponding portal mixes
LLSS into LLHH in every synthetic model; a Weinberg-only source feeds LLSS
through the same full RGE. The synthetic input deliberately takes all
renormalisable Yukawas and gauge couplings zero to isolate quartic mixing.
These zero couplings are *only* synthetic fixture choices, never defaults for
real Workstream B inputs.

On 2026-11-09 the local isolated test run passed **12/12**; 16/16 one-step
qualitative full-tensor trajectories completed with nonzero induced Weinberg.
The generated norms were about `3.56e-13 GeV^-1` from arbitrarily chosen
`~2e-10 GeV^-1` LLSS, `0.02` portal and nonzero hard `C5`; these numbers are
neither physical T3 neutrino masses nor normalization benchmarks.

## Matchete / RGBeta dictionary

`RunT3ParallelATensorDictionary.py` extends the existing `T3Batch1TensorMap.Polynomial` numerical polynomial engine. Unlike direct symbol-by-symbol matching it **sums reused native slots** into one coefficient column and reports separate Matchete/RGBeta matrix ranks, singular values, fit/holdout residuals and the full numeric `RGBeta coefficients = A @ Matchete coefficients` map. It requires explicit numerical polynomial components, coupling prefactors, H.c., phase and potential/Lagrangian signs, **and** an independently established common scalar component frame. The old `stage22-native-rgbeta-declarations-v1` JSON cannot be fed to it. A numerically exact *provided subspace* is not a complete physical scalar dictionary unless all ordinary and additional quartics and all mass/Yukawa maps are demonstrated.

`RunT3ParallelAExportRGBetaComponents.wl` reads `$quartics` and attempts native evaluation of the exact invariant functions at explicit SU(2) component indices for all 16, including the new and reused symbols, recording unresolved native CG expressions rather than guessing them. It deliberately marks the **RGBeta quartic coupling-prefactor extraction and common component frame blocked** until their native action normalization is confirmed; a successful component probe is not a Matchete/RGBeta dictionary. `ScalarPotential.wl` is left untouched. The A(-4) anti-Hermitian `lambdaH2Inv1` portal must retain `I` inside `PlusHc`; the independent `Inv2` remains active. The supported `transforms.py` can subsequently transform quartics, scalar and Weyl masses, Yukawas and Wilson coefficients after *proven* field matrices are supplied; no numerical mass/Yukawa identification is invented.

## Exact Windows commands (PowerShell, from repo root)

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
# First extract the delivered ZIP into this directory using Expand-Archive -Force.
git status --short
git rev-parse HEAD
.\.venv\Scripts\python.exe -m unittest discover -s tests -p 'test_t3_parallel_a*.py' -v
.\.venv\Scripts\python.exe -m scripts.RunT3ParallelA --synthetic A_-4 --steps 3 --out output\parallel_t3\wilson\a_synthetic_v1
.\.venv\Scripts\python.exe -m scripts.RunT3ParallelA --synthetic-all16 --steps 1 --out output\parallel_t3\wilson\all16_synthetic_v1
```

Real B inputs (when available) should instead use:

```powershell
.\.venv\Scripts\python.exe -m scripts.RunT3ParallelA --model A_-4 `
  --root . --adm-variant fonseca_eq22_six_term --solver DOP853 --steps 12 `
  --background-npz 'output\parallel_t3\background\A_m4\background.npz' `
  --background-json 'output\parallel_t3\background\A_m4\background_manifest.json' `
  --boundary-npz 'output\parallel_t3\background\A_m4\fermion_boundary_formal.npz' `
  --boundary-json 'output\parallel_t3\background\A_m4\fermion_boundary_manifest.json' `
  --out 'output\parallel_t3\wilson\B_connected_v1'
```

This call will **block** until B's manifest matches the documented schema. Change paths as needed; never reuse old A(-4) formal numeric sources. Native Wolfram diagnostics (not executed in this batch):

```powershell
$env:T3_EXPERIMENTAL_GENERAL_FERMIONS='1'
$env:T3_EXPERIMENTAL_GENERAL_QUARTICS='1'
$env:T3_EXPERIMENTAL_GENERAL_QUADRATICS='1'
$env:T3_EXPERIMENTAL_FLAVOR_MASSES='1'
wolframscript -file scripts\ExportGeneralT3ScalarTensorMatchete.wl `
  output\parallel_t3\wilson\matchete_native_tensors.json
wolframscript -file scripts\RunT3ParallelAExportRGBetaComponents.wl `
  output\parallel_t3\wilson\rgbeta_native_components_probe.json
```

This exporter is a new native probe; its first execution may produce actionable `Blocked_...` statuses if RGBeta does not evaluate a CG tensor numerically. Never treat the older `ExportGeneralT3ScalarTensorRGBetaDeclarations.wl` as a numerical tensor exporter. For signed/native polynomial map, separately provide the `t3-parallel-a-full-native-polynomials-v1` normalized component JSONs and call:

```powershell
.\.venv\Scripts\python.exe -m scripts.RunT3ParallelATensorDictionary `
  --model A_-4 --matchete-components 'path\to\verified_matchete_polynomials.json' `
  --rgbeta-components 'path\to\verified_rgbeta_polynomials.json' `
  --output 'output\parallel_t3\wilson\A_m4_dictionary_diagnostic.json'
```

Do not run B/C integration or scalar threshold in A.
