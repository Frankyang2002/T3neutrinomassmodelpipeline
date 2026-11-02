# T3 leading-one-loop completeness — Workstream 3

## Evidence and scope

Read-only GitHub repository `Frankyang2002/T3neutrinomassmodelpipeline`, `main` at commit `2c29e5ba2c0350b7e9be292e02ceb8bc7c52c0d3` (2026-11-02). Only paths in this Workstream 3 ownership are written. This package does not modify or rerun the existing fermion-tree, EFT1 seven-term, scalar-bubble or independent direct-UV calculations. The provided baseline ledger is source-inspection evidence, NOT an all-16 numerical evaluation. Run the audit locally to inspect the actual shared checkout. Other chats' files may differ from the read-only snapshot.

**Strict EFT sequence:**

`UV -> ALL BSM heavy fermions at ONE mu_F -> EFT1(S1,S2,H,SM fermions) -> BOTH BSM scalar multiplets at ONE mu_S -> SMEFT`.

No intermediate EFT with one heavy fermion or one heavy scalar; mass eigenstate sums do **not** constitute multiple thresholds.

## Source-backed findings

* `scripts/RunGeneralT3FermionFirstEFT1Tree.py` constructs formal all-fermion LLSS boundaries. The release bridge (`Numerical/joint_scalar_matching/fermion_boundary_bridge.py`) explicitly declares the **signed GammaCC-to-Weyl coefficient unknown**; formal per-operator extraction is not an independently normalized physical Wilson coefficient.
* `RGE/running/intermediate/GeneralEFT1C5Mixing.py` contains seven named master-RGE source contributions: Yukawa wavefunction, scalar pair, mixed Yukawa-gauge, crossed Yukawa, conjugate-coefficient Yukawa, scalar anomalous dimension and fermion anomalous dimension. Its required certificate list explicitly includes the independently unproven complete dim-5 inventory and signed normalization. Seven-term machinery is **not** seven-term verified physical input.
* `Numerical/joint_scalar_matching/core.py` implements the jointly diagonalized scalar bubble with coefficient `-C_LLSS * lambda_HHSS * B0bar /(2*16*pi²)`. Every pair and mass splitting is retained; an upstream full renormalized endpoint is not supplied for every model.
* `scripts/RunGeneralT3EFT1HardPoleAudit.py` factorizes B_-1/C_-1 raw LLHH poles and logs, but its own report says `physical_majorana_c5_computed=false`, `operator_hard_finite_at_muF_is_scheme_certified=false`, and `independent_counterterm_pole_check_performed=false`. It does not justify setting hard C5 to zero for the other 14 models.
* `Numerical/direct_one_loop/core.py` implements a fully complex UV Weyl/real-scalar zero-momentum one-loop mass sum; its `README.md` declares the automatic Matchete-to-Weyl vertex map provisional. A direct UV benchmark is an independent comparator, **not** an EFT matching counterterm that may be inserted into C5.
* **Verified UV Lagrangian scope:** `Lagrangian/interactions/T3Topology.wl` explicitly gates opposite-complex-scalar **extra Yukawas** behind `T3_EXPERIMENTAL_GENERAL_FERMIONS=1`, for those models with one hypercharge-zero heavy scalar. From `Lagrangian/model/T3ModelCatalog.wl`, `Y(S1)=alpha/2`, `Y(S2)=(alpha+2)/2`, so this additional coupling exists for `A_-2`, `A_+0`, `D_-2`, `D_+0`, `E_-2`, `E_+0` only. It can change the tree LLSS boundary and neutrino mass at the *same* one-loop order; its exclusion is NOT higher-loop suppression. `Lagrangian/model/LagrangianBuilder.wl` also flags that general heavy-fermion flavour masses are not certified.
* `Lagrangian/interactions/GeneralScalarOperators.wl` and `Lagrangian/interactions/ScalarPotential.wl` add scalar operators with opt-in flags; scalar potential completeness of a particular pipeline run cannot be inferred from their presence. No automated all-16 enumeration of every symmetry-allowed dim-5 EFT1 monomial is verified by the inspected modules. Thus the proposition "tree LLSS alone closes the EFT1 operator basis" remains unresolved.

## Leading-one-loop power counting

Assign `hbar` to every loop and treat quartic, Yukawa and mass couplings at `hbar^0`. There is no allowed ordinary tree neutrino mass in the exact inert-Z2 models. At first nonzero order,

```
C_LLSS(muF)    ~ Y Y / M_F                              O(hbar^0)
C5_hard(muF)   ~ lambda_HHSS * Y Y / (16pi² M_F)        O(hbar^1)
Delta C5_EFT1  ~ beta_C5[C_LLSS] * ln(muS/muF)/(16pi²) O(hbar^1)
Delta C5_muS   ~ - C_LLSS * lambda * B0bar/(32pi²)     O(hbar^1)
```

There may be additional independent hard amplitudes and operators. Additional opt-in conjugate-neutral-scalar Yukawas, when hypercharge allows them, belong at **leading one loop**, not in a later correction. One-loop corrections to the tree LLSS Wilson itself, multiplied by a further scalar bubble, are **two loops** in C5. One-loop running of couplings inside an O(hbar) source is also two loops at fixed order. Crucially, terms suppressed by `mS²/MF²` in a **single UV loop** are still **one loop**, although they enter the heavy-mass EFT as dimension 7 and above. The code retains and isolates the complete one-loop scalar-portal triangle remainder; it does not insert it spuriously into dimension-five SMEFT.

There is exactly one SM Weinberg operator (up to flavour indices) at dimension five; EFT1 can have multiple independent gauge contractions of `LL S_i S_j` because the BSM scalars remain dynamical. A field content and SU(2)xU(1)xZ2 operator-enumeration proof is still needed before declaring the LLSS basis closed. `LLHS` with exactly one inert-odd scalar vanishes by exact Z2; `LLHH` is Z2 even, **allowed** at one loop and must not be classified as symmetry-zero. `LLSS` operators generated at one loop normally enter neutrino mass only with an additional loop, assuming no odd scalar VEV.

## New conditional hard one-loop scalar-portal term

The new reusable kernel `hard_matching.py` requires a **separate, canonical real-component tree LLSS tensor for every heavy fermion mass eigenstate** `C[k,flavour1,flavour2,SU2_1,SU2_2,scalar_a,scalar_b]`. It is **not** valid to assign the *total* LLSS tensor to every heavy fermion; nor can one extract the individual fermion masses from a generic total C. All Mk are included together at **one** muF.

Let `R^T M_S² R = diag(m_a²)` and rotate `C_k` and `lambda` by this **one** orthogonal matrix for both scalar multiplets. The projective amplitude is defined using the existing real-Higgs Weinberg tensor `T_ijAB` and `LOOP=16pi²`. The matched scalar-portal hard LLHH boundary is

\[
 {\cal A}_{\mathrm{hard}}(\mu_F) = \frac{1}{2\mathrm{LOOP}}\sum_{k,a,b} C_{k,ab}\lambda_{ab}
\left[1+\ln\frac{\mu_F^2}{M_k^2}\right].
\]

The corresponding **one** scalar-source EFT1 log and **one** simultaneous finite scalar matching are

\[
 {\cal A}_{\mathrm{RG,scalar}}=\frac{\ln(\mu_S/\mu_F)}{\mathrm{LOOP}}\sum_{k,a,b}C_{k,ab}\lambda_{ab},\quad
 {\cal A}_{\mathrm{scalar}}(\mu_S)=-\frac{1}{2\mathrm{LOOP}}\sum_{k,a,b}C_{k,ab}\lambda_{ab}
\bar B_0(0;m_a^2,m_b^2;\mu_S).
\]

All amplitudes above retain their weak/Higgs/flavour tensor indices; they are projected onto C5 **only after summing**. The `1` in the hard expression is not an arbitrary fit: it comes from the difference between the zero-momentum **three-propagator full UV triangle** and the two-scalar EFT1 bubble, at leading order in `m_S²/M_F²` under the **declared** local real-scalar Wilson normalization. It is contingent on the input matching these conventions; it is not a complete model-independent renormalization proof for all other UV hard diagrams.

At any fixed `muF,muS`, the code additionally evaluates

\[
 {\cal A}_{\rm full\;triangle}=\frac{1}{2\mathrm{LOOP}}\sum_{k,a,b} C_{k,ab}\lambda_{ab}\,
 K(M_k^2,m_a^2,m_b^2),\quad
 K=\frac{F(m_a^2/M_k^2)-F(m_b^2/M_k^2)}{(m_a^2-m_b^2)/M_k^2},\quad
 F(r)=\frac{r\ln r}{1-r}.
\]

Thus `A_full - (A_hard + A_RG,scalar + A_scalar)` is the **one-loop higher-power remainder** at zero momentum for this topology. It is nonzero at finite `mS/MF`, is independent of the chosen common matching scales, and is not double-counted as a separate physical diagram. For coincident scalar masses `a=b=m²`, with `r=m²/M²`, the exact one-loop expansion is `K=(1-r+ln r)/(1-r)² = 1+ln r + r(1+2 ln r)+O(r² ln r)`; this explicitly shows the one-loop, higher-dimension power remainder. The exact coincident `M²=a=b` kernel is `-1/2`. Stable coincident scalar-mass limits are supported. If unbroken gauge invariance fails in an EWSB-split input, `projection_passed=false` and the result **cannot** be treated as a physical Weinberg coefficient. Individual heavy-fermion/eigenstate pair tensors, including suppressed terms, and independent labelled quartic channel projections are exported; individual pair projections are diagnostic and need not be gauge invariant by themselves.

### Scheme, overlaps and cancellations

* The definition of `B0bar` is MSbar; `1+ln(muF²/Mk²)` is **in the compatible MSbar definition and local operator convention**. Changing the subtraction convention shifts finite hard and scalar pieces while leaving the matched observable invariant only after a full counterterm calculation.
* At fixed source `C` and quartic, `d A_scalar / d ln muS = -sum(C lambda)/LOOP`; the scalar-source `beta(C5)` cancels it. The hard scale derivative cancels the change of the explicit scalar RG-log partition. This is an algebraic partial-sector check, not a full physical RGE certificate.
* Do **not** add the new `eft1_scalar_log` to an existing **seven-term EFT1 transport** that already includes the scalar source. In that case supply its complete transported C5 instead, and use only the new hard term once as the missing muF boundary, with the already existing joint scalar threshold once.
* Do **not** append the full UV triangle to hard + RG + scalar matching. The exact triangle is a comparator (and a higher-dimension one-loop remainder), not another independent diagram.
* The source-level `C_k` could depend on the phase/sign of the Matchete `GammaCC -> Weyl` dictionary, including vectorlike/Majorana conventions. The current source does not independently settle that sign/normalization. Therefore **the numerical kernel is implemented, but all-16 physical hard C5 remains unresolved**.
* Universal one-quartic topologies are not a proof that other hard diagrams, extra EFT1 dim-5 structures, finite wavefunction/counterterm pieces, and non-scalar mixing are absent. Never label the partial result physics-complete.

## Contract for optional per-model numerical inputs

Supply `<hard-input-dir>/A_m4.json`, `A_m4.npz`, etc. For instance, `A_m4.json` must have:

```json
{
  "model": "A_-4",
  "tensor_convention": "T3_REAL_LLSS_HHSS_WEINBERG_v1",
  "scheme": "MSbar",
  "mu_f_gev": 100000.0,
  "mu_s_gev": 1000.0,
  "tensor_npz": "A_m4.npz",
  "source_normalization_certified": false,
  "quartic_channel_names": ["physical_source_portal_1", "physical_source_portal_2"],
  "provenance": {"source": "PUT_VERIFIED_LOCAL_SOURCE_FILE_AND_CERTIFICATE_HERE"}
}
```

The numerical values above demonstrate **the schema only**. They are not fitted benchmark inputs. `A_m4.npz` must contain:

| NPZ entry | Shape | Meaning |
|---|---|---|
| `C_LLSS_by_fermion` | `(K,3,3,2,2,n,n)` complex | Separate tree Wilson tensor for **each mass eigenstate**; all heavy species |
| `C_LLSS_total_muF` | `(3,3,2,2,n,n)` complex | Must equal sum over the K sources (explicitly verified) |
| `fermion_masses_gev` | `(K,)` positive | Mass eigenvalues of all heavy BSM fermions |
| `lambda_HHSS` | `(4,4,n,n)` real | Full mixed Higgs/scalar quartic at the relevant fixed order |
| `scalar_mass_squared_gev2` | `(n,n)` real symmetric | Both multiplets at once, full mixing and self-bilinears |
| `lambda_HHSS_channels` | `(J,4,4,n,n)` real, optional | Independent labelled source contributions whose sum equals the full quartic |

Use `n=8` for A/B/C/D or `n=12` for E, with scalar ordering confirmed from the actual source. The code refuses missing per-species sums or mismatched quartet channels. No unknown normalizations, masses or extra scalar couplings are assigned default zeros or guessed values. Additional full seven-term source inputs belong to their separate workstream; the new calculator does not replace that machinery.

## PowerShell installation and execution

From the repository root:

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "$HOME\Downloads\T3_LO_Completeness_Workstream3.zip" -DestinationPath . -Force
python -m unittest discover -s tests -p 'test_lo_completeness_*.py' -v
python -m scripts.RunLOCompletenessAudit --all --check-only `
  --output-dir 'output\general_observable\lo_certification_completeness\check_20261102'
python -m scripts.RunLOCompletenessAudit --all `
  --output-dir 'output\general_observable\lo_certification_completeness\local_audit_01'
Get-Content 'output\general_observable\lo_certification_completeness\local_audit_01\summary.json' -Raw
```

If and only if **source-verified per-heavy-fermion** input contracts exist, request the conditional all-16 matrix calculation (missing models are reported as `MISSING_INPUT`, not assigned zero):

```powershell
python -m scripts.RunLOCompletenessAudit --all `
  --hard-input-dir 'input\lo_hard_conditional' `
  --output-dir 'output\general_observable\lo_certification_completeness\conditional_hard_01'
```

The runner requires **a fresh output directory** and refuses overwrites or output paths outside Workstream 3. Its ledger uses exactly `INCLUDED`, `ZERO_BY_SYMMETRY`, `HIGHER_ORDER`, `MISSING`, `UNRESOLVED`. Its optional per-model numeric status is a **separate** readiness field, not a ledger status. The package includes a baseline **read-only-source snapshot** under `output/general_observable/lo_certification_completeness/source_snapshot_20261102/`; run the local audit to record checked-out source-presence and git blob comparisons.

## Outstanding blocks to full one-loop certification

1. Enable and verify every symmetry-allowed extra neutral-scalar conjugate Yukawa (six models), verify general heavy-fermion flavour masses, and certify the full opt-in UV scalar quartic/mass basis.
2. An independently verified signed and absolutely normalized `GammaCC -> Weyl` map plus mass-eigenstate `C_LLSS_by_fermion` for every model, including vectorlike/mixed chiralities.
3. Exhaustive symmetry/SU(2)/hypercharge operator enumeration through EFT1 dimension five, checking possible lepton bilinears and redundant/EOM-related operators, with a cross-check of every generated UV term.
4. A complete counterterm/subtraction and field-renormalization calculation for all one-loop hard muF LLHH diagrams. The B/C raw poles/logs are insufficient.
5. Full **renormalized**, independently normalized seven-term EFT1 C5 logs and a simultaneous muS endpoint for both scalars with the same tensor convention; independent cancellation against scalar-threshold logarithms.
6. Comparison of the physical all-16 matched mass matrices to the exact UV direct loop (with the source-derived complex vertex dictionary), including the suppressed higher-dimension one-loop remainder when the EFT hierarchy is not extreme.

Until these blocks are independently resolved, the deliverable is **implemented and algebraically tested but not a physics-complete leading-one-loop T3 pipeline**.

## Source-backed formal per-heavy-flavour input builder (Workstream 3 follow-up)

Additional owned files:

- `Numerical/lo_certification/completeness/source_bridge.py`
- `scripts/BuildLOHardInputs.py`
- `tests/test_lo_completeness_source_bridge.py`

**Verified public-source basis for this stage:** the November 2, 2026 GitHub
`main` tree `2c29e5ba2c0350b7e9be292e02ceb8bc7c52c0d3` has all 16
`source_021dcc1e8fe8ba0f/<tag>/boundary.json` and all 16
`all16_v1/<tag>/running.json`. All 16 published formal source tensors have
lepton-flavour-diagonal (`p=q`) support, and the corresponding saved UV states
have diagonal `y1`, `y2` and `MF` at the common muF. Those facts support a
**specific benchmark-only formal decomposition** into `k=1,2,3` flavour
blocks, not a general mass-eigenstate reconstruction. The algorithm verifies
rather than assumes these conditions on each local file and checks that all PL
source contributions sum to the recorded tree boundary. No spurious term is
introduced from a PR Hermitian conjugate or scalar exchange-orbit count.

### Ready-now commands (PowerShell from repository root)

```powershell
python -m unittest discover -s tests -p 'test_lo_completeness_*.py' -v
python -m scripts.BuildLOHardInputs `
  --check-only `
  --output-dir 'output\general_observable\lo_certification_completeness\source_formal_01'
python -m scripts.BuildLOHardInputs `
  --output-dir 'output\general_observable\lo_certification_completeness\source_formal_01'
```

Creates (under a new output folder) all 16:

```text
formal_fermion_by_flavour/A_m4.json
formal_fermion_by_flavour/A_m4.npz
...
formal_fermion_by_flavour/E_p2.json
formal_fermion_by_flavour/E_p2.npz
summary.json
```

NPZ arrays are specifically named
`C_LLSS_by_fermion_formal_per_NGamma`,
`C_LLSS_total_formal_per_NGamma`, and `fermion_masses_gev`. They are
**not** files intended for `input/lo_hard_conditional`, which is absent,
because the scalar tensor and signed Weyl normalization are still missing.
All statuses remain non-physical and are not accepted as independently
normalized full-one-loop C5 matching.

### Optional conditional kernel input generation

The same runner can write complete NPZ and metadata pairs *inside its owned
output folder* if an independently prepared scalar source is available:

```powershell
python -m scripts.BuildLOHardInputs `
  --scalar-input-dir 'PATH_TO_EXPLICIT_ALL16_SCALAR_BUNDLES' `
  --output-dir 'output\general_observable\lo_certification_completeness\conditional_inputs_01'
python -m scripts.RunLOCompletenessAudit `
  --all `
  --hard-input-dir 'output\general_observable\lo_certification_completeness\conditional_inputs_01\conditional_inputs' `
  --output-dir 'output\general_observable\lo_certification_completeness\conditional_hard_02'
```

Each scalar bundle must be `<tag>.json` plus the named adjacent `.npz`
for the *same source basis*. JSON keys: `model`, `tensor_npz`,
`tensor_convention="T3_REAL_LLSS_HHSS_WEINBERG_v1"`, `scheme="MSbar"`,
`mu_s_gev`, `source_scalar_real_labels` (exact labels from formal
output), `quartic_channel_names` (non-empty and unique),
`complete_scalar_HHSS_channels_supplied=true`, and a non-empty `provenance`
object. NPZ arrays: `lambda_HHSS` with shape `(4,4,n,n)`,
`lambda_HHSS_channels` with shape `(k,4,4,n,n)` (every independent
channel), and `scalar_mass_squared_gev2` with shape `(n,n)`;
`n=8` (A-D), `n=12` (E). The runner verifies sum over independent
quartic channels and a strict full Weinberg projection. It **does not**
rotate from an unspecified mass basis or guess scalar masses/portals.

The optional output numerically represents **C5 per unknown signed
GammaCC-to-Weyl factor** in the specified diagonal benchmark, not physical
C5, and the direct one-loop cross-check remains independently uncertified.
Unprovided scalar inputs are *not* treated as zero. Missing/general Yukawa
branches, hard counterterms, complete dimension-five operator closure and
non-diagonal heavy-flavour decomposition remain unresolved.

There are still exactly **one common muF and one common muS**; heavy-mass
eigenvalues label summed propagator contributions but are never used as
sequential EFT thresholds. Both heavy scalar multiplets remain in EFT1 until
integrated together at muS. No original upstream files, shared output or
GitHub data are overwritten.

## Workstream 3 scalar source bridge — full all-16 exploratory batch (2026-11-03)

**Files added:**

- `Numerical/lo_certification/completeness/scalar_bridge.py` — checks the *existing*
  `Numerical/joint_scalar_matching` native-source scalar backgrounds; pulls **both**
  heavy-scalar mass and every named HHSS channel into the fermion LLSS source basis
  using the same real orthogonal map; retains independent suppressed channels.
- `scripts/BuildLOScalarHardBatch.py` — all-16 production of source-aligned scalar
  inputs and, **only where the full Weinberg projection passes**, numerical hard
  scalar-portal matching diagnostics. Reuses existing v4, v5 and v6 code, **read only**.
- `tests/test_lo_completeness_scalar_bridge.py` — tests including a true synthetic
  bridge → per-fermion hard-matching numerical calculation.

Read-only upstream paths:

```
output/general_observable/lo_fermion_boundary/releases/source_021dcc1e8fe8ba0f/
output/general_observable/all16_v1/<model>/running.json
output/general_observable/lo_joint_scalar_matching/source_quartic_v4/native_quartic_channels/
output/full16_symbolic/T3_<family>_alpha_<alpha>/data/eft1_after_F_scalar_quartic_seed.json
Numerical/joint_scalar_matching/{source_quartic_v4,numerical_benchmark_v5,special_mass_v6,large_batch}.py
```

**Explicitly exploratory:** the pre-existing v5/v6 scalar builders use a candidate
native Matchete → RGBeta quartic parameter-name/normalization map, and the eight
mixed/split scalar models identify muF source scalar fields with muS prototype
fields. Their correctness remains unverified. The builder requires explicit
opt-in flags and never labels a conditional result physical. Full general UV
quartic closure, extra Yukawas, the signed GammaCC/Weyl factor, independent
hard counterterm matching, and full seven-term transport remain outstanding.

**Source audit, without writing:**

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
python -m unittest discover -s tests -p 'test_lo_completeness_*.py' -v
python -m scripts.BuildLOScalarHardBatch `
  --check-only `
  --output-dir 'output\general_observable\lo_certification_completeness\scalar_hard_01'
```

**Source-backed exploratory scalar + hard one-loop portal run for ALL 16:**

```powershell
python -m scripts.BuildLOScalarHardBatch `
  --allow-exploratory-map `
  --allow-exploratory-special-bases `
  --output-dir 'output\general_observable\lo_certification_completeness\scalar_hard_01'
```

Results are written ONLY below the new owned `scalar_hard_01` directory:

```
scalar_hard_01/
  summary.json                               # all-16 statuses and exact blockers
  scalar_backgrounds/A_m4.json, A_m4.npz    # existing v4/v5/v6 builder output
  conditional_inputs/A_m4.json, A_m4.npz    # only if projection passes
  conditional_hard/A_m4.json                 # only if numerical candidate passes
  projection_rejected/A_m2.json              # diagnostic ONLY if Ward projection fails
  ...
```

An existing result may instead be reused **read-only** without any automatic
scalar source regeneration:

```powershell
python -m scripts.BuildLOScalarHardBatch `
  --existing-background-dir 'output\general_observable\lo_joint_scalar_matching\exploratory_all16_v6\numeric_backgrounds' `
  --output-dir 'output\general_observable\lo_certification_completeness\scalar_hard_from_v6_01'
```

A fresh output path is mandatory for every run. No GitHub writes, no scalar
threshold splitting, and no individual fermion threshold running occur.

The upstream v6 snapshot on GitHub (`main`, inspected 2026-11-03) contains
12 conditional scalar-C5 models and **four full-Weinberg Ward failures**:
`A_-2`, `D_+0`, `E_-2`, `E_+0`. They are not assigned a C5 or silently
accepted. The new hard-portal projection may reject additional models because
its per-heavy-flavour kernels are not identical to the old scalar-only bubble.
Record the actual generated `summary.json` rather than presuming a count.

**Interpretation of hard output:** `hard_muF` is a leading-power *conditional*
MSbar hard-portal candidate; `eft1_scalar_log` is the corresponding **single**
logarithmic partition; `scalar_threshold_muS` uses **one** joint scalar bubble;
`dim5_sum` is the sum of those three pieces; `full_triangle` is an independent
fixed-order one-loop mass-ratio kernel; `higher_dimension_remainder` is one-loop
but of higher order in the heavy-fermion expansion. These are not six
independent physical mass contributions and must not be added on top of one
another. The other six EFT1 RGE terms and other hard topologies are **not**
claimed complete. None of these outputs is a certified physical mass matrix.

## Workstream 3: conditional hard-portal numerical consistency + external evidence gate

New full files (exclusive owned paths):

- `Numerical/lo_certification/completeness/consistency_audit.py`
- `scripts/AuditLOConditionalConsistency.py`
- `tests/test_lo_completeness_consistency.py`

The reader consumes the **already generated** `scalar_hard_01` and does not modify it. Every model's saved `C5` components are independently recomputed with the existing hard kernel and source input NPZ, after checking that the adjacent scalar background files have the SHA-256 recorded in the conditional input. The check rejects a mismatched model, missing channels, changed saved matrices, different scalar/fermion thresholds or an unexpected physical-certification claim. It is benchmark-specific and remains per-unknown-signed `N_GammaCC`.

The diagnostic verifies the declared **partial scalar-portal topology**:

1. `hard_muF + eft1_scalar_log + scalar_threshold_muS = dim5_sum`, as reconstructed from unchanged source input.
2. With `C_LLSS,k`, quartics, and masses held constant, independently varying `muF` and `muS` leaves this combined partition unchanged to relative `1e-9` (algebraic matching-scale consistency, **not** the full seven-term running of actual coupling trajectories).
3. For the same analytic zero-momentum full triangle, scaling the scalar *mass squares* by `1/4` and `1/16` produces a higher-dimensional remainder with approximate exponent near one in `mS²/MF²`. This is a **power-counting stress test**, not an independent UV amplitude computation.
4. Weinberg residuals are checked separately for the retained hard, logarithm, scalar finite, dim5 sum and full triangle. No hard-only Ward result is allowed to erase an earlier separate scalar-threshold Ward failure.
5. Read-only Workstream 2 evidence `output/general_observable/lo_certification_couplings/ws2_portal_normalization_v1/summary.json` is incorporated. This audit records missing native-source portal directions and as-yet unnormalized additional HHSS candidates, even if all numerical source channels are present in the old v4 cache. Source-channel completeness **is not general UV portal completeness**.
6. Read-only previous v6 scalar matching summary is incorporated. The older projection-failing A_-2, D_+0, E_-2 and E_+0 models (as recorded in v6) remain a separate consistency blocker even when new hard-sector projection passes. The optional Workstream C GammaCC calibration report, if present, is indexed as a **conditional** Clifford/Grassmann result, never promoted into a signed physical Wilson certificate.

All output fields explicitly retain `physical_C5_certified=false` and `physical_full_one_loop_certified=false`. Do **not** add this partial scalar-source log to a seven-term EFT1 RGE solution which already includes the scalar source. Neither the synthetic UV triangle nor its own truncated expansion is an independent one-loop full UV matching calculation.

From the local project root in PowerShell:

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "$HOME\Downloads\T3_LO_Completeness_Workstream3_ConsistencyGate.zip" -DestinationPath . -Force
python -m unittest discover -s tests -p 'test_lo_completeness_*.py' -v
python -m scripts.AuditLOConditionalConsistency --check-only
python -m scripts.AuditLOConditionalConsistency --output-dir 'output\general_observable\lo_certification_completeness\conditional_consistency_01'
Get-Content 'output\general_observable\lo_certification_completeness\conditional_consistency_01\summary.json' -Raw
Import-Csv 'output\general_observable\lo_certification_completeness\conditional_consistency_01\numerical_table.csv' | Format-Table -AutoSize
```

Outputs are `summary.json`, `numerical_table.csv`, and `blockers.json` within the requested **new** Workstream 3 output directory. The output path is refused if already present; choose `conditional_consistency_02` on a rerun. If an upstream optional evidence report is absent, its status remains `NOT_PROVIDED`, not treated as a pass. An incomplete numerical input or failed consistency check returns nonzero and a fail-closed summary.

**Physical interpretation:** The latest observed all-16 portal-span audit has nontrivial extra symmetry-allowed channels; the previous 16/16 conditional hard portal successes therefore do *not* close the general one-loop T3 problem. The independent signed GammaCC-to-Weyl calibration and full counterterm/operator inventory still belong to their original workstreams. Workstream 3 does not change those workers or their data.
