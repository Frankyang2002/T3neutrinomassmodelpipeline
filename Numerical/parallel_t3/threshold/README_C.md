# Workstream C — simultaneous scalar matching and SMEFT (candidate v2)

**C v2 real A/B integration and strict Wilson power-counting instructions:** see [`README_C_REAL_AB_v2.md`](README_C_REAL_AB_v2.md). This file also retains the original v1 numerical interface description. No physical complete-one-loop claim is made.

## Ownership and source baseline

All delivered files are confined to `Numerical/parallel_t3/threshold/`, `scripts/RunT3ParallelC*` and `tests/test_t3_parallel_c*`. Read-only source inspected on **9 Nov 2026**: `Frankyang2002/T3neutrinomassmodelpipeline`, `main=a944e0030b8ef24b23342d99e0d8f33cef3cffd5`. Local uncommitted code is inaccessible here; compare local code before installing. **No GitHub writes**.

The C implementation has no A/B imports and never mutates their files. It calls the **existing** `Numerical.running.SMWeinbergEvolution.SMWeinbergInitialConditions/evolve_sm_weinberg` to run the SM+Weinberg one-loop RGEs.

## Actual numerical EFT order

`UV -> all F out simultaneously at MF -> EFT1 (both S active; all required dim-5 Wilsons evolved by A) -> single simultaneous removal of every real heavy-scalar mass eigenstate at MS -> SMEFT (already evolved C5 + finite bubble + separately provided independent hard terms) -> SM one-loop C5(mu_low).`

No scalar-by-scalar decoupling. No pre-existing LLSS->C5 leading log is added. `C5_EFT1(MS)` is extracted from the **FULL** evolved tensor, not set to zero. Both scalars' general mass eigenstates are diagonalized together, with `n(n+1)/2` unordered scalar pairs and multiplicity 2 for distinct eigenstates.

### Implemented physics / convention

Source master convention:

- Real scalar components: `H=(h0+i h1,h2+i h3)/sqrt(2)`; `S_i=(R_i+i I_i)/sqrt(2)`.
- `L5=C[i,j,a,b] psi_i psi_j phi_a phi_b + h.c.`, **no factorial**, scalar and Weyl pairs symmetric.
- `L5^W=(1/2) C5[p,q] (L_p epsilon H)(L_q epsilon H)+h.c.`. Consequently `C[nu_p,nu_q,h2,h2] = C5[p,q]/4` (zero-based `h2` = neutral real Higgs component). This **differs by a factor two** from `joint_scalar_matching/core.py:weinberg_tensor()` and is intentionally **not** copied blindly.
- `V4=lambda[a,b,c,d] phi_a phi_b phi_c phi_d/4!`, **real, fully symmetric** tensor.
- Conditional heavy-bubble contribution, at one common `MS`, with all eigenstates:

  `Delta C[pqijAB] = -1/(2*16*pi^2) sum_(a,b) C[pqijab](MS) lambda[ABab](MS) B0bar(m_a^2,m_b^2;MS)`.

  `B0bar = 1 - [ma2 ln(ma2/MS^2) - mb2 ln(mb2/MS^2)]/(ma2-mb2)`, with smooth degenerate limit `-ln(m2/MS^2)`.
- Full quartic input is **retained** and checked, including multiple independent portal structures, phase-sensitive contractions, scalar self-couplings and cross-scalar quartics. For this single `LLSS + HHSS` bubble only the **full tensor's `HHSS` projection** enters. `SSSS` does not independently produce two external Higgs legs at this one-loop insertion; other operators/diagrams must be matched separately.
- `strict_one_loop` projects the *already run but explicitly one-loop-truncated* A full Wilson tensor `C_psi2phi2_ms_strict_1loop` onto C5(MS), and contracts the bubble with A-supplied **order-zero** `C_psi2phi2_tree_mf`. Both are required in C v2; simply using the fully RGE-improved C5 would not be strictly one-loop in the Wilson sector. The A/B common field-coordinate mapping must be supplied. The B mass/quartic endpoints at MS remain unexpanded running parameters.
- `rge_improved` contracts the bubble with evolved `C_psi2phi2_ms`, including higher-order terms in the fixed-order expansion.
- `m_nu(mu_low) = -(v^2/2) C5(mu_low)` in GeV. The three masses output are sorted Takagi singular values (`numpy.linalg.svd`) converted to eV; **no pole corrections**.

### Mass-basis adapters

`models.from_parameters(model, mS1Sq, mS2Sq, mu_sq=...)` supports every production model. Common canonical real basis is interleaved `(S1R,S1I,...,S2R,S2I,...)`. It returns the full real mass Hessian and explicit `R_canonical_from_native`. The master matching also accepts **any independently provided positive real-symmetric n×n mass Hessian**, so nondegenerate masses and permitted cross-sector mixings are not artificially restricted to these analytic parameterisations.

- `B_-1/C_-1`: two hypercharge `-1/2` doublets `Phi=(S1, epsilon S2*)`; `Phi=U (X,Y)`, `U† [[m1²,mu],[mu*,m2²]] U = diag(mX²,mY²)`. Convert this complex unitary map to a real **8×8 orthogonal** matrix. The `-epsilon` map is included. Both doublets' eight real eigenstates decouple **together** at MS.
- `A_-2,A_+0,D_-2,D_+0,E_-2,E_+0`: `S[k]=(Minus[k]+i Plus[k])/sqrt(2)` in fixed native Cartesian SU2 labels (no spurious spherical `W`), with `V2=m² S†S-[mu (S†)^2+h.c.]`. Hessian per component `[[m²-2Re(mu),-2Im(mu)],[-2Im(mu),m²+2Re(mu)]]`; eigenvalues `m²±2|mu|`. Imaginary mixing is retained: **complex mu is never silently made real**.
- All other points: unbroken quadratic mass Hessian with independent S1/S2 masses. Full mass matrix from B supersedes this reduced analytic constructor.

`RunT3ParallelCMassBasis.py` exports `R_canonical_from_native`, eigenbasis, native and canonical Hessians to B. It does **not** assert the A/B CG/component dictionary is certified. A/B must provide `R_wilson_from_background` for the final C match, even when that matrix is identity.

## Required A/B interface (v1 adapter contracts)

C consumes *two* sibling JSON/NPZ pairs for each model. This is an explicit **C input adapter**, not a claim that current in-progress A/B output already uses these keys. A/B can export/rename their actual tensors to this format without changing C. `np.load(..., allow_pickle=False)`, no pickles.

**A: `wilson.json`** (sibling `wilson_at_MS.npz`):

```json
{
  "format": "t3-parallel-c-wilson-input-v1",
  "tensor_npz": "wilson_at_MS.npz",
  "model": "A_-4",
  "matching_context_id": "same-state-id-on-both-sides",
  "matched_seed_id": "same-F-threshold-seed-id",
  "lagrangian_source_id": "same-source-hash-or-immutable-id",
  "mu_f_gev": 100000.0,
  "mu_s_gev": 1000.0,
  "units": "GeV",
  "scheme": "MSbar",
  "canonical_wilson_convention": "C_ijab_full_sum_no_factorial_C5_1over2_H_R_interleaved_v1",
  "wilson_evolved_to_ms": true,
  "includes_eft1_weinberg": true,
  "legacy_direct_llss_to_c5_log_added": false,
  "scalar_basis_id": "A_-4:canonical_interleaved",
  "lepton_indices": [[0,1],[2,3],[4,5]],
  "higgs_source_indices": [0,1,2,3],
  "heavy_source_indices": [4,5,6,7,8,9,10,11],
  "synthetic_fixture": false,
  "wilson_adm_certified": false,
  "wilson_dim5_closure_certified": false,
  "signed_GammaCC_to_Weyl_certified": false
}
```

Indices are **explicit zero-based** positions into A's actual full tensor; illustrative values are not assumed by C. A NPZ arrays:

| Key | Shape / meaning |
| --- | --- |
| `C_psi2phi2_ms` | `(Nw,Nw,4+n,4+n)` **full** `complex128` evolved Wilson at MS, all Weyl species retained. |
| `C5_EFT1_ms` | Optional `(3,3)` `complex128`; verified against the Weinberg **component extracted from full C**, never added separately. |
| `C_psi2phi2_tree_mf` | Same full shape, order-zero Wilson at MF, **required in strict mode**. |
| `C_psi2phi2_ms_strict_1loop` | Same full shape; the one-loop-TRUNCATED evolved full Wilson at MS; required in strict mode with manifest `strict_1loop_endpoint_at_ms: true`, `strict_1loop_origin: "..."`. |
| `C5_independent_finite_at_MS` | Optional `(3,3)` independent finite term, **not** a leading log; requires manifest `finite_hard_ms_kind="independent_finite_only"`. |

The A convention string must only be declared after its actual spinor and real-field map is understood. Using the test fixture format label does not certify the upstream map.

**B: `background.json`** (sibling `background.npz`):

```json
{
  "format": "t3-parallel-c-background-input-v1",
  "tensor_npz": "background.npz",
  "model": "A_-4",
  "matching_context_id": "same-state-id-on-both-sides",
  "matched_seed_id": "same-F-threshold-seed-id",
  "lagrangian_source_id": "same-source-hash-or-immutable-id",
  "mu_f_gev": 100000.0,
  "mu_s_gev": 1000.0,
  "units": "GeV",
  "scheme": "MSbar",
  "quartic_potential_convention": "V4_lambda_real_symmetric_over_4factorial_v1",
  "quartic_inventory_status": "full_tensor",
  "wilson_scalar_basis_id": "A_-4:canonical_interleaved",
  "background_scalar_basis_id": "A_-4:background_basis_id",
  "single_common_scalar_threshold": true,
  "scalar_mass_origin": "unbroken_electroweak",
  "synthetic_fixture": false,
  "sm_couplings_at_ms": {"gY": 0.36, "g2": 0.65, "g3": 1.1, "lambdaH": 0.26},
  "quartic_channel_names": ["all_native_invariants_already_mapped"],
  "scalar_operator_inventory_certified": false,
  "native_matchete_rgbeta_normalization_certified": false
}
```

B NPZ arrays:

| Key | Shape / meaning |
| --- | --- |
| `mass_squared_gev2_ms` | `(n,n)` real-symmetric heavy-scalar Hessian in B background basis. |
| `R_wilson_from_background` | `(n,n)` real orthogonal; **`s_A = R s_B`**, so `C_B = R^T C_A R` on the two scalar legs. Required even if identity. |
| `quartic_real_full_ms` | `(4+n,)*4` fully symmetric **real** V4 tensor at MS, in B basis, Higgs indices first. **No reduced portal accepted**. |
| `quartic_real_channels_ms` | Optional `(k,4+n,4+n,4+n,4+n)` channels summing exactly to full quartic; `quartic_channel_names` ordered. |
| `ye`, `yu`, `yd` | Full complex `(3,3)` SM Yukawas at MS; runner uses the existing SM conventions. |
| `scalar_gauge_generators_full` | Optional `(ng,4+n,4+n)` imaginary Hermitian gauge generators in B real-field basis, Higgs first. Computes mass and quartic Ward residuals; claiming `gauge_ward_certified:true` without this array is rejected. |

For the 16 models, `n=2*(dS1+dS2)` is 8 or 12 real heavy components. No field is removed from the tensors before `MS`.

**Provenance gates**: same model, `matching_context_id`, `matched_seed_id`, `lagrangian_source_id`, MF/MS, scalar-basis link, canonical conventions, unit and scheme. Different old/new A(−4) portal seeds therefore fail **if A/B provenance accurately represents the input states**. C cannot independently reproduce an unreported Matchete source hash.

## Outputs

`output/parallel_t3/threshold/<run_id>/<model_slug>/`:

- `scalar_mass_basis.npz`: input mass Hessian, `R_wilson_from_background`, eigenvalues and orthogonal mass eigenvectors.
- `C5_smeft_at_MS.npz`: selected-power-counting `C5_EFT1_MS`, also `C5_EFT1_MS_rge_improved_reference`, `C5_scalar_bubble_MS`, `C5_independent_finite_MS`, their single sum `C5_SMEFT_MS`.
- `C5_smeft_low.npz`: `C5_SMEFT_low`, `m_nu_GeV`, sorted `mass_singular_values_eV`.
- `C5_per_quartic_channel.npz` if named full-tensor invariant channels supplied.
- `manifest.json`: input/output SHA-256 digests, scale, scheme, basis, term accounting, pair-level contributions, missing terms and certification flags.
- Top-level run `summary.json`: all requested models and individual successes/blockades.

The JSON result explicitly records `NumericalCandidateNotCertified`. `--match-only` is for matching diagnostics and sets `SMEFT_running.performed=false`; it uses `C5(MS)` as an *unchanged* low-scale proxy, not as real evolved C5.

## Windows PowerShell commands

From `C:\thesistemp\T3neutrinomassmodelpipeline`, after downloading the ZIP:

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "$env:USERPROFILE\Downloads\T3_Workstream_C_RealAB_Integration_v2.zip" -DestinationPath . -Force
.\.venv\Scripts\python.exe -m pytest -q tests\test_t3_parallel_c_threshold.py tests\test_t3_parallel_c_real_integration.py
```

Synthetic A(−4) full matching **and repository SMEFT ODE**:

```powershell
.\.venv\Scripts\python.exe scripts\RunT3ParallelCThreshold.py `
  --make-fixture A_-4 --mode rge_improved --mu-low 173 --vev 246.22 `
  --fixture-dir output\parallel_t3\threshold\fixtures `
  --output output\parallel_t3\threshold\synthetic_A4 --overwrite
```

Strict one-loop counting (requires A's **order-zero** and **one-loop-truncated** full Wilson endpoints; C v2 synthetic fixture supplies both):

```powershell
.\.venv\Scripts\python.exe scripts\RunT3ParallelCThreshold.py `
  --make-fixture A_-4 --mode strict_one_loop --mu-low 173 `
  --output output\parallel_t3\threshold\synthetic_A4_strict --overwrite
```

All-16 **synthetic algebra** test, including six self-split and two mixed doublet examples:

```powershell
.\.venv\Scripts\python.exe scripts\RunT3ParallelCThreshold.py `
  --make-fixtures-all16 --mode rge_improved `
  --output output\parallel_t3\threshold\synthetic_all16 --overwrite
```

Genuine A/B **exported** endpoint data, after the respective workstreams have produced the C v1 manifests:

```powershell
.\.venv\Scripts\python.exe scripts\RunT3ParallelCThreshold.py `
  --wilson output\parallel_t3\wilson\A_m4\wilson.json `
  --background output\parallel_t3\background\A_m4\background.json `
  --mode rge_improved --mu-low 173 `
  --output output\parallel_t3\threshold\real_candidate_v1
```

All-16 A/B prepared directories `output\parallel_t3\AB_exports\<slug>\wilson.json` and `background.json`:

```powershell
.\.venv\Scripts\python.exe scripts\RunT3ParallelCThreshold.py `
  --all16-input-root output\parallel_t3\AB_exports `
  --mode rge_improved --output output\parallel_t3\threshold\all16_candidate_v1
```

Export mixed doublet basis for B:

```powershell
.\.venv\Scripts\python.exe scripts\RunT3ParallelCMassBasis.py `
  --model B_-1 --mS1sq 1900000 --mS2sq 2600000 `
  --mu-real 140000 --mu-imag 110000 `
  --output output\parallel_t3\threshold\basis_B_m1
```

Export full B-supplied nondegenerate real Hessian (`--mass-key` defaults to `mass_squared_gev2_ms`):

```powershell
.\.venv\Scripts\python.exe scripts\RunT3ParallelCMassBasis.py `
  --model A_-4 --background-npz output\parallel_t3\background\A_m4\background.npz `
  --output output\parallel_t3\threshold\basis_A_m4
```

## Explicit remaining physics uncertainties

1. **A**: signed Matchete GammaCC→canonical Weyl polynomial factor, exact full real-field index map, correct complete one-loop EFT1 ADM and other dimension-five classes. A's full tensor must already include genuine Weinberg running/mixing at MS.
2. **B**: native Matchete↔RGBeta full symmetry-allowed scalar/mass/quartic numerical dictionary. **A(−4) anti-Hermitian `lambdaH2Inv1` must remain `I` inside `PlusHc`**. Do not mix any old pre-fix A(−4) numerical state with the corrected source. For special cases B must provide an explicitly aligned native scalar field basis/rotation; C's analytic maps are not automatic certification of that alignment.
3. **C**: only the project-convention scalar bubble and separately provided independent finite hard C5 are evaluated. Other one-loop scalar matching diagrams, wavefunction/EOM/broken-phase effects, and possible contributions from dimension-five operator classes other than the selected `psi²phi²` Wilson sector are **not** computed or proved absent. Consequently **this does not assert complete one-loop MS matching**.
4. Strict counting projects the A-supplied one-loop-truncated full tensor at MS for C5, contracts MF order-zero LLSS with the common scalar bubble, and does not reinsert A's C5 leading log. B mass/quartic inputs remain renormalized at MS; their running is not expanded in this calculation and can formally introduce higher-order effects. Full fixed-order certification still requires external A/B perturbative provenance.
5. Gauge invariance: optional numerical Ward checks require A/B gauge generators. Unbroken mass/quartic tensors without them carry `not_supplied_not_certified`. EWSB-resummed masses may generate SU2-non-Weinberg remnants; matching rejects those by default. `--allow-projection-residual` outputs only a **conditional projection**, not a physical SMEFT Wilson coefficient.
6. One-loop SMEFT ODE uses the existing project's SM coupling and vev convention; low-scale matrix/singular values omit pole matching and are not physical mass measurements.

### Tests and sample numerical output

Original v1 sandbox: **26 passed, 1 skipped** (the skipped test imports the repository SMEFT runner unavailable in the isolated archive). Every model's synthetic `--match-only --mode strict_one_loop` numerical path passed. In a **temporary test mirror of the inspected repository SM beta equations**, the synthetic A(−4) `rge_improved` run integrated successfully (`nfev=38`) and yielded illustrative `m_nu` singular values approximately `(0.0141, 0.0178, 0.0202)` eV at 173 GeV. These are synthetic, NOT a T3 prediction or an execution of the user's local runner. The real repository runner is invoked by the actual delivered CLI and is tested when installed locally.

### C v2 execution status

Isolated v2 suite: **53 passed, 1 skipped** (`SMWeinbergEvolution` absent from isolated ZIP tree). All 16 C synthetic fixtures were separately executed in BOTH `strict_one_loop` and `rge_improved` `--match-only` modes through the v2 real-integration CLI (16/16 complete; 32 conditional match results). This is a contract test only, not real A/B integration. See the v2 readme for strict-full-Wilson requirement and real-input preflight.
