# Research General T3 Tensors And Quartics

> Historical research record. This document combines **original source
> Markdown** under original-path headings. Local Markdown links to moved
> files are updated; exact original bytes are in the external archive.
> This is not a new certification
> of physics, a rewritten conclusion, or a new computational workflow.
> Original bytes and SHA-256 values are stored in the external archive.

## Contents

- [General T3: Source sign and conditional quartic action identities](#source-general-t3-action-identity-md) — `docs/GENERAL_T3_ACTION_IDENTITY.md`
- [General T3: four shared B/C RGBeta quartics — local correction](#source-general-t3-bc-local-correction-md) — `docs/GENERAL_T3_BC_LOCAL_CORRECTION.md`
- [General T3: independent SU(2) audit of four remaining reused B/C quartics](#source-general-t3-bc-reused-covariance-md) — `docs/GENERAL_T3_BC_REUSED_COVARIANCE.md`
- [General T3 scalar EFT1 — first narrow generalisation correction](#source-general-t3-eft1-quartic-fix-notes-md) — `docs/GENERAL_T3_EFT1_QUARTIC_FIX_NOTES.md`
- [Stage 20 RGBeta v1.2.0 compatibility and candidate-gauge audit](#source-general-t3-general-quartic-rgbeta-compat-fix-md) — `docs/GENERAL_T3_GENERAL_QUARTIC_RGBETA_COMPAT_FIX.md`
- [General T3: RGBeta quartic projection fix](#source-general-t3-general-quartic-rgbeta-projected-betas-md) — `docs/GENERAL_T3_GENERAL_QUARTIC_RGBETA_PROJECTED_BETAS.md`
- [General T3: Native RGBeta complex-quartic assembly probe (read-only)](#source-general-t3-native-quartic-action-md) — `docs/GENERAL_T3_NATIVE_QUARTIC_ACTION.md`
- [General T3: native quartic prefactors and conditional signed coupling maps](#source-general-t3-quartic-prefactors-md) — `docs/GENERAL_T3_QUARTIC_PREFACTORS.md`
- [General T3: Native RGBeta SU(2) `S2` and triplet conjugation structure](#source-general-t3-rgbeta-s2-structure-md) — `docs/GENERAL_T3_RGBETA_S2_STRUCTURE.md`
- [General T3: SU(2) triplet component-basis witness](#source-general-t3-triplet-intertwiner-md) — `docs/GENERAL_T3_TRIPLET_INTERTWINER.md`

---

<a id="source-general-t3-action-identity-md"></a>

## Original document: `docs/GENERAL_T3_ACTION_IDENTITY.md`

**Original SHA-256:** `9a7f97385f8e8fa3dc8d8f9e889745a64e4d54e3b7f056a499af6469d1b95113`

**Status:** Historical record; check the current source and main
documentation before executing commands reproduced below.

<details>
<summary>Expand original: General T3: Source sign and conditional quartic action identities</summary>

# General T3: Source sign and conditional quartic action identities

## Scope and result

This stage is **local/read-only** and does **not** change RGBeta production definitions or promote physical-coupling certification. It uses the user's **fresh post-Stage24** Stage22/23 evidence and successful **live Stage25** native assembly report, then checks the quartic *polynomials*, including repeated identical scalar fields and Hermitian-conjugate contributions, in the **provisional common component labels** used by Stage22.

The result is stronger than matching the ordered invariant tensors alone. It is still **not** an independently established global Matchete/RGBeta Cartesian real-field intertwiner; the six `CalibratedTripletCGCandidate` cases require that separate verification.

## Relative scalar-action sign: published and project source

RGBeta's formal convention is explicit in A. E. Thomsen, *Introducing RGBeta*, EPJC **81** (2021) 408 (DOI [10.1140/epjc/s10052-021-09142-4](https://doi.org/10.1140/epjc/s10052-021-09142-4)):

- Eq. (2.3): the model-builder Lagrangian contains `- Σ lambda (phi phi phi phi)`.
- Eq. (2.9): the complex field is `phi = Delta.Phi/sqrt(2)`.
- Eq. (2.13): the quartic piece in the generic real-field Lagrangian has `- 1/4 (lambda Delta^4 Phi^4)`.
- Eq. (2.15): the real symmetric tensor uses `1/4` times its 24-permutation sum.

RGBeta's `AddQuartic` native implementation stores, for four complex scalar legs,

\[
  \Lambda\text{-normalisation} = 4!\,(1/\sqrt2)^4 = 6.
\]

Combining the generic `1/4!` quartic action, the native factor and the complex-field projection gives

\[
  \frac{6}{4!}(\sqrt{2})^4 = \frac14\times4 = 1.
\]

Thus the **RGBeta** contribution in model-builder complex fields is **`-lambda * invariant * fields`**, not `-6*lambda` and not `-lambda/4`. The conjugate interaction is present as a separate term for `SelfConjugate -> False`.

In this project's `Lagrangian/model/LagrangianBuilder.wl`, the BSM interactions enter the action through `LInt = Total[Lookup[valid, "Expression", {}]]`. `scripts/ExportGeneralT3ScalarTensorMatchete.wl` records these new quartics as `+lambda[] times contracted fields in PlusHc[contribution]` with prefactor 1. Stage26 verifies these exact source strings locally, and **fails closed** if they change. This is a **source-level relative-sign argument**, not a new live Matchete-vs-RGBeta action evaluation.

Consequently, under the **same physical field-component basis** and identical holomorphic orientation, a Stage22 tensor block `R=M A` requires

\[
  r = -A^{-1}m,
\]

where `m` and `r` are Matchete and RGBeta coupling column vectors. Do not insert this map in production until the global component basis is established.

## Stage26 implementation

`python -m scripts.QuarticActionIdentity`:

1. Recomputes Stage23 and Stage25 diagnostics from their native inputs, compares them with stored JSON, and fails on stale evidence. In particular it needs the **live Stage25** RGBeta assembly file and its 50-slot/46-block result.
2. Checks the local Matchete builder and Stage22 exporter source conventions.
3. Constructs an **exact commuting polynomial** from every Matchete component tensor and an **independently reconstructed** RGBeta invariant tensor in the Stage22 provisional component basis. All repeated-field monomials are added, retaining exact coefficients (`sympy`, no floats).
4. Explicitly adds Hermitian conjugates using the formal species swaps `H↔Hdag`, `S1↔S1dag`, `S2↔S2dag` and conjugates coefficients.
5. Substitutes the **conditional** signed Stage23 coupling matrices and checks the *entire* polynomial identity `L_M + V_R = 0` exactly, plus negative controls for a wrong sign or a spurious factor `1/4`.
6. Requires exactly **50 slots / 46 blocks** across 16 models and refuses non-post-Stage24 reindexing classes.

**Important limitation:** Formal conjugation in the provisional component labels is not proof that the underlying native triplet components are physically identical or correctly normalised. The source sign is well documented; the physical coupling map is **still UNCERTIFIED**. No JSON flag in this stage claims otherwise.

## Files

- `scripts/QuarticActionIdentity.py` — complete new diagnostic
- `tests/python/test_general_t3_quartic_action_identity.py` — regression tests and optional live integration
- `docs/GENERAL_T3_ACTION_IDENTITY.md` — this document

The output report goes to:

`output/general_t3_stage22/stage26_conditional_quartic_action.json`

## PowerShell commands (project root)

```powershell
Expand-Archive -LiteralPath .\T3_General_Scalar_Stage26_ActionIdentity.zip `
    -DestinationPath . -Force

python -m pytest -q tests\python\test_general_t3_quartic_action_identity.py

python -m scripts.QuarticActionIdentity
```

Requires successful output from the previous commands:

- `python -m scripts.RGBetaAlignedGeometry`
- `python -m scripts.QuarticPrefactorAudit`
- `python -m scripts.NativeQuarticActionAudit`

**Expected output after the user's successful Stage25 live run:** 16 models, **50 quartic slots**, **46 complete conditional action blocks**, with source-supported sign and **physical coupling mapping UNCERTIFIED**.

## Next gate

Independent **physical component-basis** verification of the six calibrated triplet structures. In particular the calibrated `sqrt(2)/3^(1/4)` tensor coefficient must not be mistaken for a proven, kinetic-term-preserving unitary field redefinition. Only after the global field identification passes can the signed coupling matrices be used in the numerical RGBeta pipeline.

</details>

---

<a id="source-general-t3-bc-local-correction-md"></a>

## Original document: `docs/GENERAL_T3_BC_LOCAL_CORRECTION.md`

**Original SHA-256:** `ba2d8c07f7cc3367bb09be512dad2e9300898d63c9111144945db1deae89fd8e`

**Status:** Historical record; check the current source and main
documentation before executing commands reproduced below.

<details>
<summary>Expand original: General T3: four shared B/C RGBeta quartics — local correction</summary>

# General T3: four shared B/C RGBeta quartics — local correction

**Scope:** repair only four verified, native SU(2)-noncovariant `AddQuartic` invariant definitions in `RGE/running/rgbeta/T3RGBetaModel.wl`. This file is the full replacement of current GitHub `main` blob SHA `1711bd0c8e639347e11494013bd5b9863df0fa11`, plus precisely four expression changes. No GitHub edit.

## Inputs and evidence

The user's Stage24 local run reported **13 Python tests passed**, with **8/8 native B/C covariance failures** and **8/8 proposed SU(2)-covariant exact Matchete matches**. An independent live RGBeta export confirmed 2 models, 8 expected native failures, 8 proposed passes, and exact Wolfram/Python component agreement. This tests the standard `U`/`Conjugate[U]` field representations, not signed physical coupling normalisation.

## Exact replacements

Each is written for its existing registered field order; each `SelfConjugate -> False` stays unchanged. Expressions below use `D(a,b)=del[SU2L@fund,#a,#b]`, `E(a,b)=eps[SU2L@fund,#a,#b]`.

- `lambdaS1barS2S2bar2`, `{Bar[S1], S2, Bar[S2], Bar[S2]}`: `(E(1,4) D(2,3) + E(1,3) D(2,4))/2`.
- `lambdaS1bar2S2bar2`, `{Bar[S1], Bar[S1], Bar[S2], Bar[S2]}`: `(E(1,4) E(2,3) + E(1,3) E(2,4))/2`.
- `lambdaS1S1bar2S2bar`, `{S1, Bar[S1], Bar[S1], Bar[S2]}`: `(D(1,3) E(2,4) + D(1,2) E(3,4))/2`.
- `lambdaHHdagS1barS1bar`, `{Bar[H], Bar[H], Bar[S1], Bar[S1]}`: `(E(1,4) E(2,3) + E(1,3) E(2,4))/2`.

No other `AddQuartic`, coupling name, or potential prefactor is altered. The earlier `lambdaHHdagS1barS2barCross` correction is preserved.

## Diagnostic compatibility

`scripts/RGBetaAlignedGeometry.py` now recognises both the historical `J`-reindexed tensor candidates and the *fully* repaired Stage24 geometry profile. For the repaired profile, each of these eight slots has an **exact direct tensor equality** `R=M`, without a per-block `J`.

`scripts/QuarticPrefactorAudit.py` accepts the old and repaired **whole-model** classification profiles. It still emits a provisional signed matrix and retains `signed_physical_coupling_map_certified = false`. The shared potential-sign assumption, live real-scalar tensor normalisation and global field convention still require independent verification.

The old `BCQuarticRemainderCovariance.py`, `ExportGeneralT3BCQuarticNativeCovariance.wl`, and `VerifyGeneralT3BCQuarticNativeInvariants.py` are historical **pre-fix** diagnostics; do not rerun them against repaired native exports. Use the new post-fix exporter and verifier below.

## Windows PowerShell (run from project root)

Make a backup and install all complete replacement/new files:

```powershell
Copy-Item .\RGE\running\rgbeta\T3RGBetaModel.wl .\RGE\running\rgbeta\T3RGBetaModel.wl.pre-stage24.bak
Expand-Archive -LiteralPath .\T3_General_Scalar_Stage24_BC_LocalFix.zip -DestinationPath . -Force
```

Check the Python regression suite:

```powershell
python -m pytest -q tests\python\test_general_t3_bc_postfix.py tests\python\test_general_t3_prefactors.py tests\python\test_general_t3_tensor_aligned_geometry.py
```

Recheck the **live registered** quartics and compare all 8 SU(2) tensors with Matchete:

```powershell
wolframscript -file .\scripts\ExportGeneralT3BCQuarticCovarianceCorrection.wl `
    output/general_t3_stage22/stage24_bc_postfix_wolfram_native_covariance.json
python -m scripts.VerifyGeneralT3BCQuarticCovarianceCorrection
```

Refresh the complete Stage22/23 native diagnostic evidence, in this order:

```powershell
powershell -ExecutionPolicy Bypass -File .\scripts\RunGeneralT3ScalarTensorCertification.ps1 -ExportRGBetaDeclarations -NativeDiagnosticsOnly
python -m scripts.RGBetaAlignedGeometry
python -m scripts.QuarticPrefactorAudit
```

Expected Stage22 **conditional** geometry classes after the full correction:

```text
DirectGeometryCandidate:                 32
CorrectedBCReusedInvariantCandidate:      8
CalibratedTripletCGCandidate:             6
CorrectedReusedInvariantCandidate:        4
slots:                                   50
blocks:                                  46
```

The Stage23 signed coupling-map count remains 46. Neither the live SU(2) verification nor these exact component matches certify the physical Matchete–RGBeta coupling map; all certification flags remain false.

## Synthetic checks performed before packaging

- Source block fingerprint: only the four previously `delS2 * delS2` contractions replaced; previous 110101 repair preserved.
- Exact SU(2) test: proposed tensors pass all three generators.
- Synthetic Stage22 post-fix registry fixture: 32 + 8 + 6 + 4 = 50 slots, 46 blocks.
- Synthetic Stage23 prefactor fixture: 50/50 native factor checks, 46 conditional signed matrices.
- Python post-fix verifier regression: accepts exact in-memory component tensors and rejects corrupted SU(2) residuals, incorrect components, wrong Matchete tensors and old-format exports.

**No live post-fix RGBeta run was performed in the packaging environment.**

</details>

---

<a id="source-general-t3-bc-reused-covariance-md"></a>

## Original document: `docs/GENERAL_T3_BC_REUSED_COVARIANCE.md`

**Original SHA-256:** `fb97e05700b6c19e8bbcd115f7f4a0dfdb55f06987aad99236068a6c32e09934`

**Status:** Historical record; check the current source and main
documentation before executing commands reproduced below.

<details>
<summary>Expand original: General T3: independent SU(2) audit of four remaining reused B/C quartics</summary>

# General T3: independent SU(2) audit of four remaining reused B/C quartics

## Why this stage is needed

Stage 23 returned 50/50 native RGBeta prefactors, 46/46 full-rank **conditional** coupling matrices, and `physical_coupling_mapping: UNCERTIFIED`. During the next global field-basis check, no single conjugation-consistent `I` or epsilon-`J` transformation on the three complex doublets `(H, S1, S2)` aligned all six B/C quartic blocks with Matchete simultaneously. The maximum was five of six blocks. The reported **eight pseudoreal-reindexed slots** are only per-block tensor identities, not a verified global field basis.

The actual **registered** B/C RGBeta expressions (from the Stage 22 native Wolfram exports) reveal that four *other* reused coupling definitions, each used in **both** `B_-1` and `C_-1`, fail the explicit SU(2) generator test under the same fundamental/conjugate convention that verified the corrected `110101` invariant.

## Exact candidate contractions

For each `AddQuartic`, component indices `1,2,3,4` refer to the ordered fields in the declaration; `eps` is the antisymmetric SU(2) fundamental invariant with `eps[1,2]=+1`, and `del` is a Kronecker delta. The replacements below are **proposed**: they are not yet installed in `RGE/running/rgbeta/T3RGBetaModel.wl`.

| Registered coupling | RGBeta fields | Candidate `GroupInvariant` function (schematic component form) |
|---|---|---|
| `lambdaS1barS2S2bar2` | `(S1†, S2, S2†, S2†)` | `(ε₁₄δ₂₃ + ε₁₃δ₂₄)/2` |
| `lambdaS1bar2S2bar2` | `(S1†, S1†, S2†, S2†)` | `(ε₁₄ε₂₃ + ε₁₃ε₂₄)/2` |
| `lambdaS1S1bar2S2bar` | `(S1, S1†, S1†, S2†)` | `(δ₁₃ε₂₄ + δ₁₂ε₃₄)/2` |
| `lambdaHHdagS1barS1bar` | `(H†, H†, S1†, S1†)` | `(ε₁₄ε₂₃ + ε₁₃ε₂₄)/2` |

The full native Mathematica expressions are in `scripts/ExportGeneralT3BCQuarticNativeCovariance.wl`; the explicit index-by-index Python expressions are in `scripts/BCQuarticRemainderCovariance.py`.

For the two B/C models, the **Python exact calculation** establishes:

- Eight original source-extracted SU(2) tensors have nonzero Lie-generator residuals; all eight proposed tensors have zero residual for `T1`, `T2`, and `T3`.
- Each candidate tensor equals the corresponding native Matchete 1D tensor **component by component**, without an additional global `J` reindexing: `R=M`, `A=[1]`.
- Each candidate agrees with the *previous conditional* per-block `J` witness. This does not verify that `J` is a globally valid field redefinition.
- A restricted scan of all `2³ = 8` common transformations that apply `I` or `J` to both members of each doublet/conjugate pair finds **no** global alignment of all six B/C blocks with the unmodified four reused invariants. This scan is not a proof about arbitrary unitary basis transformations.
- These are geometric and covariance results only. **No physical coupling map, scalar-potential sign, or beta function is certified.**

## Evidence integrity

The Python script calls the existing `QuarticPrefactorAudit.run()` first; this checks 16 models, 50 native slots, 46 blocks, registered RGBeta prefactors, Matchete tensor ranks, and the aligned geometry SHA-256 input hashes. It then checks the exact names, fields, and source expressions of the eight remaining B/C reused slots. It **fails closed** if the native expressions are stale or have already been changed in an unreviewed way.

The independent Wolfram script builds the real RGBeta `B_-1` and `C_-1` models, loads the generated extension and reads their **registered** `$quartics` invariant functions using `RGBetaNativeAccess.wl`, then evaluates both the native and candidate contractions under the three generators. The component heads `del`/`eps` are converted before substituting numerical indices to avoid RGBeta's `del[rep,a,a] = Dim[rep]` dummy-index trace convention being misapplied to component labels.

The live exporter has **not yet been run in this environment**. After it runs, `VerifyGeneralT3BCQuarticNativeInvariants.py` checks all sixteen original/proposed component arrays (two models × four couplings × two versions) and all covariance statuses against the independent Python evaluation.

## Files in this ZIP (complete new files)

- `scripts/BCQuarticRemainderCovariance.py`
- `scripts/ExportGeneralT3BCQuarticNativeCovariance.wl`
- `scripts/VerifyGeneralT3BCQuarticNativeInvariants.py`
- `tests/python/test_general_t3_bc_covariance.py`
- `docs/GENERAL_T3_BC_REUSED_COVARIANCE.md`

**No production file, Stage 22 script, Stage 23 script, or GitHub repository is changed by this package.**

## PowerShell commands (repository root)

```powershell
Expand-Archive -LiteralPath .\T3_General_Scalar_Stage24_BC_Covariance.zip -DestinationPath . -Force

python -m pytest -q tests\python\test_general_t3_bc_covariance.py
python -m scripts.BCQuarticRemainderCovariance

wolframscript -file .\scripts\ExportGeneralT3BCQuarticNativeCovariance.wl `
    output/general_t3_stage22/stage24_bc_wolfram_native_covariance.json

python -m scripts.VerifyGeneralT3BCQuarticNativeInvariants
```

Expected Python result: 8/8 native failures, 8/8 covariant proposed tensors, no restricted global `I/J` field basis, physical mapping **UNCERTIFIED**. The live verifier should report `native_failed_expected: 8`, `proposed_su2_pass: 8`, `exact_python_wolfram_components_agree: True`. If a check fails, **do not** apply the candidates as production changes. Please provide the live output JSON for source review.

## Next code change, after successful live verification

Modify only the four shared `GroupInvariant` functions in the **complete local** `RGE/running/rgbeta/T3RGBetaModel.wl` replacement. Then regenerate native RGBeta declarations/invariants and update Stage 22–23 geometry and prefactor diagnostics to classify these eight slots as *direct corrected-invariant candidates*, not pseudoreal reindexings. Rerun the local tests and keep physical normalisation uncertified pending a genuinely native common real-scalar action comparison.

</details>

---

<a id="source-general-t3-eft1-quartic-fix-notes-md"></a>

## Original document: `docs/GENERAL_T3_EFT1_QUARTIC_FIX_NOTES.md`

**Original SHA-256:** `bbff55c61ae582c47d2f0f1e3f1c056980cbad49145f9341d4218b3107de151f`

**Status:** Historical record; check the current source and main
documentation before executing commands reproduced below.

<details>
<summary>Expand original: General T3 scalar EFT1 — first narrow generalisation correction</summary>

# General T3 scalar EFT1 — first narrow generalisation correction

## Verified source boundary

Source inspected: GitHub `main` at the time of packaging (read-only).
Files inspected: `RGE/general/generated/general_rgbeta_quartic_registry.json`,
`Lagrangian/interactions/{GeneralScalarOperators,GeneralScalarVerifiedQuartics,ScalarPotential}.wl`,
`RGE/running/rgbeta/{T3RGBetaModel,RunGeneralT3RGBeta}.wl`,
`Numerical/experimental/GeneralT3ScalarRGState.py`,
`scripts/{FermionThresholdHandoff,IntermediateEFTRunning,ExportGeneralT3EFT1OneLoop}.py/.wl`,
and documentation for Stages 22–29 and 34–36.

The public source contains all 16 generated RGBeta extension files, the
38 new coupling declarations, and all 12 B/C reuse records.
The ordinary Stage34 exporter **explicitly omits** the new Stage29 quartics.

## First correction

Install **only** `scripts/RunGeneralT3EFT1Quartics.py`. This separate
experimental entrypoint uses the current `Stage33` EFT1 boundary and its
38 new quartic sidecars, converts them to the existing
`GeneralT3ScalarRGState`, regenerates **EFT1** full-projected RGBeta beta
functions using the present local `T3RGBetaModel.wl`, and numerically evolves
all baseline, reused, and new quartics from 1e5 to 1e3 GeV.

Pre-existing `Stage34` reports and `Stage35–36` neutrino-mass processing remain
unchanged. This is **renormalisable EFT1 only** (not full EFT including
higher-dimensional Wilson-coefficient running).

## Essential limitations

- Native RGBeta/Matchete physical component/phase/potential conventions are
  **not certified**. The initial new couplings are used in the explicitly
  convention-labelled RGBeta basis already present in Stage33.
- Eight allowed scalar quadratic structures remain absent from scalar mass
  matrices. Six are neutral singlet/triplet holomorphic bilinears; two B/C
  structures are inter-doublet mixing terms. In a mass-eigenstate basis some
  off-diagonal mass parameters are removable via field rotation, but this
  changes Yukawa/quartic tensors. This does **not** justify simply dropping
  interactions from a fully general fixed-basis Lagrangian.
- Fermion-threshold finite matching and the resulting higher-dimensional
  Wilson coefficients are **not** recalculated.
- The new quartic flow is **not** inserted into scalar threshold matching,
  `C5`, or neutrino mass observables.
- Do not feed the output into Stage35 or Stage36 and describe it as a fully
  general T3 result.
- Older `output/general_t3_stage20/*_eft1_rgbeta.json` must **not** be reused
  blindly; the runner regenerates expressions after the later B/C covariance fixes.

## Installation, PowerShell (project root)

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath C:\path\to\T3_Generalisation_EFT1_Fix.zip -DestinationPath . -Force
python -m py_compile .\scripts\RunGeneralT3EFT1Quartics.py
python -m scripts.RunGeneralT3EFT1Quartics --model B_+1
```

For all 16 once the single-model check succeeds:

```powershell
python -m scripts.RunGeneralT3EFT1Quartics --all
```

`wolframscript` must be on `PATH` with RGBeta installed. The Stage33 files
and hashed Stage30/31 source artifacts must exist at their current standard
location. A missing dependency or mismatched RGE-key set causes a failure,
not substitution of zero beta functions.

## Output project locations

- `output/general_t3_eft1_quartics/<MODEL>_eft1_live_full_rgbeta.json`
- `output/general_t3_eft1_quartics/<MODEL>_eft1_full_quartic_running.json`

The `<MODEL>` in filenames uses `B_p1`, `A_p0`, `E_m2`, etc.

`docs/general_t3_operator_gap_inventory.csv` enumerates 38 new, 12 reused,
and eight blocked quadratic model-specific slots. "Additional quartic slots"
are model-specific: same named symbolic contraction in two different models
is counted twice. The 12 reuse records are six in B_-1 and six in C_-1,
not 12 distinct globally named RGBeta couplings.

</details>

---

<a id="source-general-t3-general-quartic-rgbeta-compat-fix-md"></a>

## Original document: `docs/GENERAL_T3_GENERAL_QUARTIC_RGBETA_COMPAT_FIX.md`

**Original SHA-256:** `20026d126de7143a42f2d699594316cc2fad5da4bd64add92726d075af0a8b89`

**Status:** Historical record; check the current source and main
documentation before executing commands reproduced below.

<details>
<summary>Expand original: Stage 20 RGBeta v1.2.0 compatibility and candidate-gauge audit</summary>

# Stage 20 RGBeta v1.2.0 compatibility and candidate-gauge audit

## Failure identified

The local Wolfram log reported:

    General::invalidopt: Option CheckInvariance for function AddQuartic received invalid value True.
    General scalar quartic registration failed

The 16 generated `AddQuartic` expressions previously set `CheckInvariance ->
(Environment["T3_GENERAL_RGBETA_VERIFY_GAUGE"] === "1")`. When the user selected
`-CheckGauge`, that option evaluated to `True` and RGBeta v1.2.0 rejected the
registration. The Wolfram process therefore failed **before** computing betas.

## Repair

- The generator no longer passes the unsupported `CheckInvariance` option.
  The option uses RGBeta's default; it is not silently re-enabled.
- The `-CheckGauge` flag now performs an independent, explicit **structural**
  audit of all 38 *new* candidate quartic operators: U(1)Y neutral charge,
  Z2 parity, conjugated field orientation and SU(2) candidate multiplicity.
  The audit **is not RGBeta's native Ward-identity gauge check**, and the
  output records `native_rgbeta_checkinvariance_used=false`.
- The generated extra B/C(-1) `H Hdag S1dag S2dag` contraction was corrected
  from `eps(H,S2dag)*eps(Hdag,S1dag)` to
  `del(S2dag,H)*eps(Hdag,S1dag)` in the standard dual-doublet convention.
  A regression compares the two tensors under finite SU(2) rotations.
- The existing six B/C(-1) RGBeta couplings are **not modified**. Their
  independent exact tensor-normalisation checks remain outstanding.
- The runner retains the prior CLI fixes (`m1` for negative alpha; comma-separated
  model selection) and captures complete Wolfram logs.

## Local tests

    python -m pytest -q tests/python/test_general_t3_general_quartic_r_g_beta.py tests/python/test_general_t3_general_quartic_quadratics.py

On the isolated packaging filesystem: 35 passed, 2 skipped (production-only
numerical-state tests). WolframScript is not installed in this environment;
no RGBeta beta-export success is claimed.

## Local Windows command

    powershell -ExecutionPolicy Bypass -File .\scripts\RunGeneralT3GeneralQuarticRGBeta.ps1 -ExportBetas -CheckGauge -Models B_-1,C_-1

The JSON output metadata is candidate-only. Do **not** use these couplings as
certified Matchete parameters or as complete physical neutrino-mass predictions.

</details>

---

<a id="source-general-t3-general-quartic-rgbeta-projected-betas-md"></a>

## Original document: `docs/GENERAL_T3_GENERAL_QUARTIC_RGBETA_PROJECTED_BETAS.md`

**Original SHA-256:** `6f2e68a75498ccc4c65884841cf52b5a1e7b613290db089445df2625922a7436`

**Status:** Historical record; check the current source and main
documentation before executing commands reproduced below.

<details>
<summary>Expand original: General T3: RGBeta quartic projection fix</summary>

# General T3: RGBeta quartic projection fix

## Observed local failure

B_-1 UV: `BetaTerm failed for lambdaT3GeneralQ110101I1` (exit 8).
This coupling shares its field list `{H, Bar[H], Bar[S1], Bar[S2]}`
with the existing `lambdaHHdagS1barS2barCross`, but uses a second candidate
SU(2) contraction.

RGBeta's `Package/Betas.m` implements `CheckPotentialCouplingMixing` and warns
that calling `BetaTerm[oneCoupling, 1]` when two couplings share field content
can return an unpredictable mixture. Its supported alternative is
`BetaTerm[Quartic, 1]`, which uses `ProjectionToUnmixedBetas` to invert the
quartic projection matrix.

## Code change

`RGE/running/rgbeta/RunGeneralT3RGBeta.wl` now calls `BetaTerm[Quartic,1]`
**once** and assigns its *projected* entries to **all** existing quartic
keys as well as new generated couplings. It checks for missing or unevaluated
betas, refuses duplicate names, and records
`GeneralScalarQuarticBetaBasis=RGBetaFullQuarticProjection` in JSON metadata.
Gauge, Yukawa and scalar mass beta entries retain their previous paths.

This avoids a logically inconsistent calculation that would use an unmixed new
quartic alongside mixed old quartic betas. The RGBeta projection can still
fail if the generated contractions are redundant in RGBeta's actual tensor
basis. That must be diagnosed on the local Wolfram installation; the Python
tests alone do NOT validate physical independence, exact normalisations, or
native gauge invariance.

## Test and run

```powershell
python -m pytest -q tests\python\test_general_t3_general_quartic_r_g_beta.py tests\python\test_general_t3_general_quartic_quadratics.py
powershell -ExecutionPolicy Bypass -File .\scripts\RunGeneralT3GeneralQuarticRGBeta.ps1 -ExportBetas -CheckGauge -Models B_-1,C_-1
```

If RGBeta reports `ProjectionToUnmixedBetas::singular`, treat it as a genuine
possible redundant/invalid tensor-basis issue; don't mask it by reverting to
`BetaTerm[symbol,1]`.

</details>

---

<a id="source-general-t3-native-quartic-action-md"></a>

## Original document: `docs/GENERAL_T3_NATIVE_QUARTIC_ACTION.md`

**Original SHA-256:** `717e571b3831adc61820dc11269073964e33056c7bafbdd5bf7d7c5efa97744a`

**Status:** Historical record; check the current source and main
documentation before executing commands reproduced below.

<details>
<summary>Expand original: General T3: Native RGBeta complex-quartic assembly probe (read-only)</summary>

# General T3: Native RGBeta complex-quartic assembly probe (read-only)

## Scope and result status

Stage 24 repaired the eight remaining reused B/C invariant slots and its live
native post-fix audit reported eight exact matches. Stage 23 gave 50/50 native
quartic prefactors and 46 *conditional* signed maps. **None of these is a
certification of the physical Matchete–RGBeta potential mapping.**

Stage 25A checks an independent part of the missing normalisation step: the
**actual RGBeta runtime quartic tensor assembly**, not merely the stored
`Coupling` and `UniqueArrangements` strings. It also checks RGBeta's complex
scalar species index map and scalar contraction metric. This is a necessary
cross-check, but the exporter intentionally does **not** claim to have obtained
a Cartesian real-scalar component tensor or established a potential sign.

## Source conventions checked

Checked against the upstream `aethomsen/RGBeta` `Package/FieldsAndCouplings.m`
source (current repository at preparation):

- `UpdateFieldIndexMap` (around lines 130–155) represents each complex scalar
  `phi` and `Bar[phi]` as separate index entries, with scalar contraction
  metric 2 for the conjugate pair.
- `AddQuartic` (around lines 751–787) constructs
  `normalization=24 (1/sqrt(2))^4=6` for four complex slots and calls
  `CouplingPermutations`.
- `Package/GroupsAndIndices.m` (lines 275–277) contains
  `Matrix[m_][] = m`.  Thus the stored function
  `6*Matrix[lambda][##1] &` evaluates at zero flavour indices to the
  symbolic coupling **`6*lambda`**, not the unevaluated
  `6*Matrix[lambda][]`.  Stage 25 checks this exact expression in Wolfram
  and independently checks its serialisation in Python.  This change fixes
  only the runtime-format mismatch; it does not modify physical normalisation.
- `AvgPermutations` / `CouplingTensorMap` (around lines 925–943)
  populate the tensor using each unique field-index permutation with
  weight `1/N`.
- `Lam` (around lines 1035–1046) builds its quartic tensor by using
  `CouplingTensorMap` and explicitly including the conjugate term when
  `SelfConjugate=False`.

The live exporter uses the *runtime* `CouplingTensorMap`, not a reimplementation
substituted in its place. For each of the 50 Stage 22 native declarations, it:

1. Rebuilds the relevant RGBeta model and loads the actual generated extra
   quartic declarations (16 models).
2. Confirms that the registered `Fields`, `Invariant`, `Coupling`, and
   `UniqueArrangements` match the refreshed Stage 22 declarations.
3. Evaluates `CouplingTensorMap` and independently reconstructs its ordered
   rules from the native permutation list; compares each LHS and full symbolic
   RHS (including field-index permutation and SU(2) tensor factors).
4. Checks that all four scalar slots have a conjugate partner in the native
   field-index map and that each pair has metric component 2.
5. Checks the exact native complex-quartic coupling factor.

The Python verifier additionally recalculates Stage 23 (including the Stage 22
source-file SHA256 consistency checks), compares it to the saved Stage 23 report,
and checks all 50 live permutation counts across the 46 conditional blocks.

## What this does NOT show

`$fieldIndexMap` uses complex field species and conjugate species as independent
bookkeeping entries; it does **not** directly label the Cartesian real and
imaginary components. In particular, the factor `6` in `Lam` alone does not
prove the sign or physical coefficient in the Matchete Lagrangian.

The following are still uncertified:

1. A single globally consistent Matchete–RGBeta SU(2) component intertwiner,
   particularly the independent triplet basis and relative phases.
2. The physical action sign connecting a Matchete `+lambda * operator` in
   `LUV` to RGBeta's real/quartic action convention. Stage 23's
   `lambda_R=-A^{-1} lambda_M` remains an **assumption-dependent candidate**.
3. Direct comparison of the **entire native quartic action**, with the
   Hermitian-conjugate contribution and identical-field combinatorics, after
   converting complex scalar species to an independently verified common
   Cartesian basis. Stage 25A tests RGBeta's assembly engine but not the
   `Lam` contraction against a separately exported physical Matchete action.
4. The eight previously blocked scalar quadratic entries.

No production RGE, vertex definition, numerical benchmark, or certification
flag is modified by this package.

## Installation (PowerShell from repository root)

Place `T3_General_Scalar_Stage25_NativeAssembly.zip` in
`C:\thesistemp\T3neutrinomassmodelpipeline` and run:

```powershell
Expand-Archive -LiteralPath .\T3_General_Scalar_Stage25_NativeAssembly.zip `
    -DestinationPath . -Force

python -m pytest -q tests\python\test_general_t3_native_quartic_action.py

wolframscript -file .\scripts\ExportGeneralT3NativeQuarticAssembly.wl `
    output/general_t3_stage22/stage25_live_native_quartic_assembly.json

python -m scripts.NativeQuarticActionAudit
```

Prerequisites: regenerated Stage 22 native declarations and aligned geometry,
Stage 23 prefactor report after Stage 24 corrections; `RGBeta` available in
`wolframscript`. If these have not been run, regenerate before invoking Stage 25:

```powershell
powershell -ExecutionPolicy Bypass -File .\scripts\RunGeneralT3ScalarTensorCertification.ps1 `
    -ExportRGBetaDeclarations -NativeDiagnosticsOnly
python -m scripts.RGBetaAlignedGeometry
python -m scripts.QuarticPrefactorAudit
```

Successful Stage 25 live run writes:

- `output/general_t3_stage22/stage25_live_native_quartic_assembly.json`
- `output/general_t3_stage22/stage25_native_quartic_assembly_diagnostics.json`

Expected verified coverage: **16 models, 50 slots, 46 blocks**. The total
number of native permutation rules is measured at runtime; it is intentionally
not predicted here.

If the Wolfram run fails, retain the original error output and share it;
do not interpret the lack of a JSON report as a successful certification.
The Wolfram exporter was not executable in the preparation environment.

## Next numerical decision

If the live assembly passes, the remaining blocking substep is to compare a
minimal set of independently oriented, Hermitian-conjugated quartic *actions*
with the same component basis and fix the relative sign. Apply the mapping to
the RGE only if that explicit comparison is successful. For rough benchmarks,
limit this check to the couplings actually used at the target point rather
than expanding to unnecessary phenomenology.

</details>

---

<a id="source-general-t3-quartic-prefactors-md"></a>

## Original document: `docs/GENERAL_T3_QUARTIC_PREFACTORS.md`

**Original SHA-256:** `f0f0e6f2a869a0de4add567f068cdb84551de30b627a542fdda379c1ec650a17`

**Status:** Historical record; check the current source and main
documentation before executing commands reproduced below.

<details>
<summary>Expand original: General T3: native quartic prefactors and conditional signed coupling maps</summary>

# General T3: native quartic prefactors and conditional signed coupling maps

## Scope

This is a **read-only diagnostic** on the Stage 22 post-fix native evidence. It does **not** modify `T3RGBetaModel.wl`, beta functions, numerical inputs, or any Stage 22 certification flags. The result is a *conditional signed map*, not a certified physical coupling transformation. The calculation fails closed if the `B_-1`/`C_-1` corrected `110101` invariant is absent or if the Stage 22 exports are stale/incomplete.

## Source-verifiable RGBeta factors

The RGBeta `Package/FieldsAndCouplings.m` implementation (upstream `aethomsen/RGBeta`) uses, in `AddQuartic`:

```wolfram
normalization = 24;
Do[
  If[!$scalars[phi, SelfConjugate] || Head @ phi === Bar,
    normalization *= 1/Sqrt[2]
  ],
  {phi, {phi1, phi2, phi3, phi4}}
];
lam = Evaluate[normalization Matrix[coupling] @ ##] &;
```

Consequently, for four complex scalar slots, the registered quartic tensor coefficient is

\[
n_{\mathrm{RGBeta}} = 4!\left(\frac{1}{\sqrt2}\right)^4 = 6.
\]

The native `rgbeta_native_declarations.json` and `rgbeta_native_invariants.json` files are checked for the **actual registered** `6*Matrix[lambda][##1]&` on **all 50 slots**, including the reuse cases. `SelfConjugate` must be `False` and the two independent exports must agree on their invariant functions and permutation lists.

In the same upstream source, `CouplingTensorMap` calls `AvgPermutations`, assigning each of the `N` `UniqueArrangements` terms weight `1/N`. With the standard real-field quartic contraction `\Lambda_{abcd}\phi_a\phi_b\phi_c\phi_d/4!`, their total coefficient per unsymmetrised complex-field orientation is

\[
\frac1{4!}\sum_{p=1}^N\frac6N=\frac14
=\left(\frac1{\sqrt2}\right)^4.
\]

The matching `1/4` is a consistency check of the **complex-to-real component normalisation**. In particular, **the visible `6*Matrix[...]` is not an extra factor of six in the physical scalar potential**. This check does not independently identify RGBeta's live real-scalar basis with Matchete's basis and does not establish the physical relative sign.

The Stage 22 Matchete tensor exporter marks all 50 quartic contributions with `potential_prefactor="1"`, and the Matchete `Lagrangian/interactions/ScalarPotential.wl` constructs the scalar interactions as Lagrangian expressions, with `PlusHc` for non-Hermitian terms. This annotation is source evidence, not a proof that all component conventions have been globally identified.

## Exact conditional coupling map

For each field-content block, the Stage 22 aligned diagnostic has column-basis equality

\[
R=MA,\quad R_j=\sum_i M_iA_{ij}.
\]

Let `m` and `r` denote column vectors of Matchete and RGBeta couplings. Under the **additional hypotheses** that the two programs use exactly the same physical fields and that the Matchete quartic term appears with `+m` in the Lagrangian whereas RGBeta uses `-V(r)`, equality requires

\[
M m=-R r\quad\Longrightarrow\quad r=-A^{-1}m.
\]

The script calculates and checks `A(-A^{-1})=-I` **exactly** in each of 46 full-rank blocks, retaining exact radicals. The sign `-` is explicitly a *conditional physical-convention assumption*; it is not inferred from the number `6`.

Examples (columns/rows in the order contained in each block's JSON record):

- Direct 1D block, `A=[1]`: `r=-m` **if the hypotheses hold**.
- `B_-1`/`C_-1`, corrected `110101` block:

  \[
  A=\begin{pmatrix}\tfrac12&1\\-\tfrac{\sqrt3}2&0\end{pmatrix},\qquad
  -A^{-1}=\begin{pmatrix}0&\tfrac{2\sqrt3}{3}\\-1&-\tfrac{\sqrt3}{3}\end{pmatrix}.
  \]

- `E_+0`/`E_-2` 2D triplet-only blocks:

  \[
  A=\begin{pmatrix}\tfrac12&1\\\tfrac{\sqrt{33}}6&0\end{pmatrix},\qquad
  -A^{-1}=\begin{pmatrix}0&-\tfrac{2\sqrt{33}}{11}\\-1&\tfrac{\sqrt{33}}{11}\end{pmatrix}.
  \]

All these matrices remain **candidates**; they must **not** be silently inserted into the production numerical or RGE pipeline.

## Input integrity and strict failure conditions

The tool reads the user's **fresh post-fix**:

- `output/general_t3_stage22/matchete_tensors.json`
- `output/general_t3_stage22/rgbeta_native_declarations.json`
- `output/general_t3_stage22/rgbeta_native_invariants.json`
- `output/general_t3_stage22/stage22_aligned_geometry_diagnostics.json`

It checks the aligned report's SHA-256 hashes against all three native inputs, total 16 models/50 slots/46 blocks, the four post-fix geometry classes (32 direct, 8 pseudoreal, 6 triplet-CG, 4 corrected B/C), exact invariant function agreement, registered normalisation factor, valid unique native field permutations, geometric component residual zero, invertibility, and full slot coverage.

Do **not** use a pre-fix `stage22_aligned_geometry_diagnostics.json`. The Python tests can synthesize a corrected fixture from older saved evidence solely for regression testing; **synthetic fixtures are never treated as live exports**.

## PowerShell (repository root)

Install `T3_General_Scalar_Stage23_Prefactors.zip` locally:

```powershell
Expand-Archive -LiteralPath .\T3_General_Scalar_Stage23_Prefactors.zip -DestinationPath . -Force

python -m pytest -q tests\python\test_general_t3_prefactors.py

python -m scripts.QuarticPrefactorAudit
```

The last command requires the actual corrected native exports already generated during Stage 22, including `stage22_aligned_geometry_diagnostics.json`. The script writes:

`output/general_t3_stage22/stage23_quartic_prefactor_candidates.json`

Expected message: `50/50` native prefactors, `46/46` exact conditional map blocks, **UNCERTIFIED** physical mapping. Existing Stage 22 certification remains false.

## Remaining certification gate

A physical certificate requires an independently verified common **real-scalar field-component basis**, including the `B/C` pseudoreal reidentification and triplet CG intertwiners, and an explicit live-RGBeta/Matchete scalar-action comparison fixing the **relative sign** and `PlusHc` counting. Until then:

- `source_prefactors_verified=true`: the **specified source factor/averaging identities** and exported native factors passed;
- `global_su2_component_basis_verified=false`;
- `live_real_scalar_quartic_action_verified=false`;
- `signed_physical_coupling_map_certified=false`.

</details>

---

<a id="source-general-t3-rgbeta-s2-structure-md"></a>

## Original document: `docs/GENERAL_T3_RGBETA_S2_STRUCTURE.md`

**Original SHA-256:** `1daa2947918bc07429cf5cfc74ec898dc4ede7fea093b38d1550ab140ac70802`

**Status:** Historical record; check the current source and main
documentation before executing commands reproduced below.

<details>
<summary>Expand original: General T3: Native RGBeta SU(2) `S2` and triplet conjugation structure</summary>

# General T3: Native RGBeta SU(2) `S2` and triplet conjugation structure

**Scope:** Read-only, fail-closed live RGBeta group-structure probe, plus exact Matchete triplet-metric audit. All output remains **UNCERTIFIED** for physical inter-package component-basis identification and quartic coupling normalisation. No production edits and no GitHub edits.

## Why a numerical RGBeta triplet matrix is not directly exportable

Inspected RGBeta v1.2.0 `Package/GroupsAndIndices.m` defines `tGen[rep,A,a,b]` as an **abstract representation generator** and `delS2[group,a,i,j]` as an invariant between a symmetric-square index and two fundamental indices. `RefineGroupStructures` expresses the `S2` generator in terms of fundamental **symbolic** generators. The implementation does *not* select numeric Pauli generator matrices or three explicit triplet component coordinates. `Package/FieldsAndCouplings.m` associates complex fields with a scalar/conjugate-scalar index pair, not with explicit Cartesian real/imaginary components.

A numeric triplet generator constructed from Pauli matrices and a conventional symmetric-product basis is therefore a **chosen reference representation**, not an independent native RGBeta generator-component measurement. The existence of a unitary equivalence alone does not select a unique field-phase or barred-field convention. Do not label that equivalence a physical coupling certificate.

## Independent native checks made by the Wolfram script

For `A_+0`, `A_-2`, `D_+0`, `D_-2`, `E_+0`, `E_-2`, after rebuilding each registered model, the exporter requires:

- Triplets use the actual registered `SU2L[S2]` gauge representation.
- Native dimension `Dim[S2] = 3`, quadratic Casimir `Casimir2[S2] = 2`, and trace normalisation `TraceNormalization[S2] = 2`.
- The group-level conjugate rule `tGen[Bar[rep],A,i,j] = -tGen[rep,A,j,i]` **as evaluated in the installed RGBeta runtime**.
- Native symmetric-doublet orthogonality `delS2[SU2L,a,i,j] delS2[SU2L,b,i,j] = del[SU2L[S2],a,b]`.
- The native complex-triplet/conjugate-triplet scalar contraction metric has both corresponding entries equal to 2, consistent with Stage 25.
- The native generator remains **symbolic** rather than a numerical 3×3 component matrix.

Any mismatch aborts and no JSON evidence is produced. No synthetic result may substitute for the user's actual Wolfram export.

## Exact Matchete triplet invariant form

From Stage 27's independently extracted normalised CG matrices \(C_a\) and \(\epsilon_{12}=1\), form

\[
\Omega_{ab}=\sum_{i,j,k,l} C_a^{ij}C_b^{kl}\,\epsilon_{ik}\epsilon_{jl}.
\]

For the project native Matchete CG data, **\(\Omega=-I_3\)**. The Stage 28 Python audit recomputes and checks the exact matrix and its invariance:

\[
(T^a_{\mathrm{Mat}})^T\Omega+\Omega T^a_{\mathrm{Mat}}=0,
\qquad \sum_a(T^a_{\mathrm{Mat}})^2=2I_3.
\]

This is a consistency check on the Matchete CG and a *source-compatible* abstract `S2` group structure. It does not determine the Matchete/RGBeta triplet **component coordinates**, nor does it prove the field conjugation map. The previously reported \(UU^T=\mathrm{diag}(1,1,-1)\) is not by itself the physical `Bar` map.

## Local commands

From the repository root, after extracting the ZIP:

```powershell
python -m pytest -q tests\python\test_general_t3_triplet_native_structure.py

wolframscript -file .\scripts\ExportGeneralT3NativeTripletStructure.wl `
    output/general_t3_stage22/stage28_live_rgbeta_triplet_structure.json
if ($LASTEXITCODE -ne 0) { throw "Stage28 native probe failed" }

python -m scripts.RGBetaTripletNativeStructureAudit
python -m pytest -q tests\python\test_general_t3_triplet_native_structure.py
```

Native input `output/general_t3_stage22/stage28_live_rgbeta_triplet_structure.json`; Python diagnostic output `output/general_t3_stage22/stage28_triplet_native_structure_audit.json`. The Python audit re-runs Stage 27 and Stage 26 checks before trusting the saved witness.

## Practical next step

There are two distinct goals: (a) show the RGBeta native `S2` algebra matches spin 1; this stage can verify it, and (b) fix a concrete component coordinate and physical field conjugation convention across Matchete and RGBeta; the abstract RGBeta group interface **cannot uniquely certify (b)**. To reach rough numerical RGE benchmarks efficiently, the next implementation should **explicitly declare** one shared component convention (rather than claim it was extracted from RGBeta), calculate corresponding couplings with tracked convention choice, and perform a basic basis-independent benchmark cross-check. A physically certified map would still require an independent cross-package action/field check.

</details>

---

<a id="source-general-t3-triplet-intertwiner-md"></a>

## Original document: `docs/GENERAL_T3_TRIPLET_INTERTWINER.md`

**Original SHA-256:** `70352009e4a3394f19c2b78420495a7b809a034c240ee845fa5f03160ebf0995`

**Status:** Historical record; check the current source and main
documentation before executing commands reproduced below.

<details>
<summary>Expand original: General T3: SU(2) triplet component-basis witness</summary>

# General T3: SU(2) triplet component-basis witness

**Scope:** New, read-only source-informed diagnostic for the six `CalibratedTripletCGCandidate` scalar quartic slots. No production RGBeta, Matchete, or running parameters are altered. Not a physical coupling certificate.

## Provenance and prerequisites

- Stage 22 native Matchete `matchete_tensors.json` (50 slots, 16 models) and native RGBeta `rgbeta_native_invariants.json` (50 slots).
- Stage 25 live `CouplingTensorMap` audit, then Stage 26 signed-action/HC report. The command-line entry point calls the Stage 26 audit afresh and compares its **entire** result with the saved Stage 26 JSON; missing or stale inputs cause an error.
- Inspected source: `Lagrangian/model/SU2Invariants.wl`, `Lagrangian/interactions/ScalarPotential.wl`, `RGE/running/rgbeta/T3RGBetaModel.wl`, and `scripts/RGBetaAlignedGeometry.py` from the project, plus RGBeta v1.2.0 `Package/FieldsAndCouplings.m`.

## Independent A/D calibration

The Matchete native quartic tensors for `A_+0`, `A_-2`, `D_+0`, and `D_-2` contain two barred Higgs doublets and exactly one SU(2) triplet (the other BSM scalar being a singlet). Extract the **native Matchete** 2×2 matrices \(Q_a(i,j)\) and fix the positive norm \(N\):

\[
N^2=\sum_{ij}|Q_a(i,j)|^2=\frac{2}{\sqrt3},\qquad C_a=\frac{Q_a}{N}.
\]

All four models must produce *the same* \(C_a\), without independent rephasings. On the recorded native data:

\[
C_1=\frac{1}{\sqrt2}\begin{pmatrix}1&0\\0&-1\end{pmatrix},\quad
C_2=\frac{1}{\sqrt2}\begin{pmatrix}0&-1\\-1&0\end{pmatrix},\quad
C_3=\frac{-i}{\sqrt2}\begin{pmatrix}1&0\\0&1\end{pmatrix}.
\]

These form an orthonormal basis of \(\mathrm{Sym}^2\mathbb C^2\). They are **read from Matchete component tensors**, not accepted from Stage 22's previously calibrated `CG` constant.

## Independently predict the two E-series tensors

Insert that single \(C\) in the contraction defining `T3RGBetaS2TripleInvariant`:

\[
F_{abc}=\sum_{ijklmn}C_a^{ij}C_b^{kl}C_c^{mn}
\epsilon_{ik}\epsilon_{jm}\epsilon_{ln},\quad
R_{pqbc}=\sum_a C_a^{pq}F_{abc}.
\]

Compare against every *native* Matchete tensor component for `E_+0` and `E_-2`. The verified tensor-only scales \(R=kM\) are:

| Four simple A/D slots | Two E slots |
|---|---|
| \(k=3^{1/4}/\sqrt2\) | \(k=1/\sqrt2\) |

Different \(k\)'s are **tensor normalisations, not running-coupling transformations**. The code rejects a changed native contraction, sign, phase, component, or source scaling.

## Exact SU(2) generator intertwining

Choose an *explicit reference*, **not a claimed native RGBeta basis**:

\[
B=(|11\rangle,\, (|12\rangle+|21\rangle)/\sqrt2,\,|22\rangle),\qquad
U_{ab}=\langle B_b|C_a\rangle.
\]

With \(t^a=\sigma^a/2\) acting on the fundamental, the two barred doublets use \(\bar t^a=-t^{a*}\). In the symmetric-product reference basis:

\[
(T_B^a)_{bc}=\langle B_b|\bar t^a\otimes1+1\otimes\bar t^a|B_c\rangle,
\qquad T_C^a=U T_B^a U^\dagger.
\]

The audit checks all three exact Hermitian generators, the commutators \([T^1,T^2]=iT^3\) and cycles, the spin-1 Casimir \(\sum_a(T^a)^2=2I\), \(U^\dagger U=I\), and \(T_C^aU=UT_B^a\).

The resulting unitary candidate and conjugation diagnostic are

\[
U=\begin{pmatrix}
1/\sqrt2&0&-1/\sqrt2\\
0&-1&0\\
-i/\sqrt2&0&-i/\sqrt2
\end{pmatrix},\qquad UU^T=\operatorname{diag}(1,1,-1).
\]

**Critical distinction:** \(UU^T\) is a basis-conjugation diagnostic. It does **not** establish the action of Matchete `Bar` or RGBeta `Bar` on physical spin-1 components. Nor have we independently extracted Matchete and RGBeta native triplet generator component matrices; the two generator sets above come from an explicit mathematical reference and the Matchete-inferred \(C\). Therefore the generator intertwining alone cannot certify a physical inter-package basis change.

## Files and outputs

New files only:
- `scripts/TripletIntertwiner.py`
- `tests/python/test_general_t3_triplet_intertwiner.py`
- `docs/GENERAL_T3_TRIPLET_INTERTWINER.md`

From the repository root in PowerShell (after unpacking):

```powershell
python -m pytest -q tests\python\test_general_t3_triplet_intertwiner.py
python -m scripts.TripletIntertwiner
```

Report: `output/general_t3_stage22/stage27_triplet_generator_witness.json`.

## Certification status / next numerical step

This tests whether one **Matchete-derived triplet CG basis** works for all six cases and whether its spin-1 generator representation is consistent. It does **not** measure native RGBeta `S2` component matrices, Matchete native generator components, or the cross-package field conjugation rule. Consequently:

- `global_matchete_rgbeta_cartesian_basis_certified = false`;
- `physical_coupling_mapping_certified = false`.

Before promoting the 46 Stage 23/26 conditional signed maps to physical RG input, independently extract and identify the native component representation and barred-field conjugation on **both sides**. A candidate \(U\) without these checks is not a certificate. No further arbitrary per-block reindexing is justified.

</details>

---
