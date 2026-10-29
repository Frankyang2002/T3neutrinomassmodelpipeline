# Research Matchete And Mass Eigen

> Historical research record. This document combines **original source
> Markdown** under original-path headings. Local Markdown links to moved
> files are updated; exact original bytes are in the external archive.
> This is not a new certification
> of physics, a rewritten conclusion, or a new computational workflow.
> Original bytes and SHA-256 values are stored in the external archive.

## Contents

- [Matchete mass-eigenstate source coupling-presence fix](#source-t3-mass-eigen-coupling-detection-fix-md) — `docs/T3_MASS_EIGEN_COUPLING_DETECTION_FIX.md`
- [B/C(-1) — native Matchete mass eigenfields and matching preparation](#source-t3-mass-eigen-matching-preparation-md) — `docs/T3_MASS_EIGEN_MATCHING_PREPARATION.md`
- [Matchete SM field-label collision fix](#source-t3-mass-eigen-matchete-sm-field-symbol-fix-md) — `docs/T3_MASS_EIGEN_MATchete_SM_FIELD_SYMBOL_FIX.md`
- [T3 mass-eigenstate provenance path correction](#source-t3-mass-eigen-provenance-fix-md) — `docs/T3_MASS_EIGEN_PROVENANCE_FIX.md`
- [Mixed-doublet Matchete live-stage diagnostics](#source-t3-mass-eigen-wolfram-stage-diagnostics-md) — `docs/T3_MASS_EIGEN_WOLFRAM_STAGE_DIAGNOSTICS.md`
- [T3: Matchete scalar operator audit with installed-version baseline calibration](#source-t3-matchete-detailed-mass-basis-audit-md) — `docs/T3_MATCHETE_DETAILED_MASS_BASIS_AUDIT.md`
- [Matchete mass-eigenstate: typed-field transport and baseline isolation](#source-t3-matchete-eigenfield-typed-index-and-baseline-isolation-md) — `docs/T3_MATCHETE_EIGENFIELD_TYPED_INDEX_AND_BASELINE_ISOLATION.md`
- [T3 mixed scalar native Matchete: internal objects compatibility fix](#source-t3-matchete-internal-objects-fix-md) — `docs/T3_MATCHETE_INTERNAL_OBJECTS_FIX.md`
- [Native mass-eigenstate Matchete: both S2 Higgs-portal CG tensors](#source-t3-matchete-native-h2-portal-cg-fix-md) — `docs/T3_MATCHETE_NATIVE_H2_PORTAL_CG_FIX.md`
- [Native SU(2) mass-basis tensors for S2-dependent baseline scalar quartics](#source-t3-matchete-native-s2-scalar-quartics-md) — `docs/T3_MATCHETE_NATIVE_S2_SCALAR_QUARTICS.md`
- [Native HHXX topology CG reconstruction — experimental B_-1/C_-1](#source-t3-matchete-native-topology-cg-fix-md) — `docs/T3_MATCHETE_NATIVE_TOPOLOGY_CG_FIX.md`
- [B_-1 / C_-1 native Matchete y1 and additional-quartic recoupling](#source-t3-matchete-native-yukawa1-and-extra-quartics-md) — `docs/T3_MATCHETE_NATIVE_YUKAWA1_AND_EXTRA_QUARTICS.md`
- [Native Yukawa2 SU(2) recoupling (experimental)](#source-t3-matchete-native-yukawa2-cg-fix-md) — `docs/T3_MATCHETE_NATIVE_YUKAWA2_CG_FIX.md`
- [T3 B_-1 / C_-1: exact symbolic mixing in native Matchete check](#source-t3-matchete-symbolic-mixing-invariance-diagnostic-md) — `docs/T3_MATCHETE_SYMBOLIC_MIXING_INVARIANCE_DIAGNOSTIC.md`
- [Typed pseudoreal-bridge indices: B_-1 / C_-1](#source-t3-matchete-typed-bridge-indices-fix-md) — `docs/T3_MATCHETE_TYPED_BRIDGE_INDICES_FIX.md`
- [T3 Matchete mixed-heavy-mass audit fix](#source-t3-matchete-mixed-mass-basis-md) — `docs/T3_MATCHEte_MIXED_MASS_BASIS.md`
- [T3 mixed scalar mass-basis batch — experimental B(-1)/C(-1)](#source-t3-mixed-scalar-mass-basis-batch-md) — `docs/T3_MIXED_SCALAR_MASS_BASIS_BATCH.md`

---

<a id="source-t3-mass-eigen-coupling-detection-fix-md"></a>

## Original document: `docs/T3_MASS_EIGEN_COUPLING_DETECTION_FIX.md`

**Original SHA-256:** `5475113ae5f6951e30b736825e046ecf3d0ccb711d82f39f54644dcfb44f966e`

**Status:** Historical record; check the current source and main
documentation before executing commands reproduced below.

<details>
<summary>Expand original: Matchete mass-eigenstate source coupling-presence fix</summary>

# Matchete mass-eigenstate source coupling-presence fix

## Evidence from the local Windows run

The `B_-1` and `C_-1` live scripts reached `original-interaction-assembly-started`
after Matchete generated all seven extra quartic candidates. Both failed with
`AdditionalQuarticCouplingsAbsentInOriginalSource`, *before* eigenfield
registration or any heavy-mass-basis validation.

The previous dynamic check

```wl
With[{coupling = symbol}, !FreeQ[expr, HoldPattern[coupling[___]]]]
```

was evaluated independently in a Wolfram kernel and falsely returned `False`
for a present applied coupling head (`aa[]`). Thus the reported failure was a
false-negative detector, not evidence of a missing physical interaction.

## Implementation

Only `Lagrangian/RunGeneralT3MassEigenMatchete.wl` has physics-adjacent changes.
The helper now obtains all **applied expression heads** via `Cases[..., HoldPattern[head_Symbol[___]] :> HoldComplete[head], ...]` and checks whether the **exact symbolic coupling head** occurs. It distinguishes different symbol contexts, unlike a comparison of bare printed names. An early in-kernel positive/negative control checks this mechanism on each run. When a coupling really is absent, the script still aborts, now listing the missing symbols.

No numerical boundary, RG equation, scalar field transformation, Yukawa, quartic
expression, CG registration or Matchete pass criterion was altered. The existing
SHA-256 input provenance checks remain unchanged.

## Exact project paths

- `Lagrangian/RunGeneralT3MassEigenMatchete.wl` — full replacement
- `tests/python/test_general_t3_mass_eigen_matchete.py` — full replacement
- `docs/T3_MASS_EIGEN_COUPLING_DETECTION_FIX.md` — new note

## Install and run (PowerShell)

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath 'C:\path\to\T3_MassEigen_Coupling_Head_Fix.zip' -DestinationPath . -Force
python -m pytest -q tests\python\test_general_t3_mass_eigen_matchete.py tests\python\test_general_t3_mixed_scalar_mass_basis.py
python -m scripts.RunGeneralT3MassEigenMatchete --model B_-1
python -m scripts.RunGeneralT3MassEigenMatchete --model C_-1
python -m scripts.RunGeneralT3MassEigenMatchete --all
```

Offline tests check source contracts and Python computations; a standalone
Wolfram kernel reproduced both the old false-negative and the replacement's
correct match. Full live Matchete execution on the user's installation remains
necessary. No physical threshold matching or Matchete/RGBeta normalization is
certified by this change.

</details>

---

<a id="source-t3-mass-eigen-matching-preparation-md"></a>

## Original document: `docs/T3_MASS_EIGEN_MATCHING_PREPARATION.md`

**Original SHA-256:** `a183f11901ca21bd76c74c4dc22ca616e7cf605d9a7ae1f285ea6b313b2ff528`

**Status:** Historical record; check the current source and main
documentation before executing commands reproduced below.

<details>
<summary>Expand original: B/C(-1) — native Matchete mass eigenfields and matching preparation</summary>

# B/C(-1) — native Matchete mass eigenfields and matching preparation

This experimental batch follows the successful `RunGeneralT3MixedScalarMassBasis --all`
outputs. It does **not** modify production UV/EFT states, RGBeta, or threshold matching.

## Physics and representation

For `B_-1` and `C_-1` the original scalar charges are
`Y(S1)=-1/2`, `Y(S2)=+1/2`. The numerical convention is

\[
\epsilon=\begin{pmatrix}0&1\\-1&0\end{pmatrix},\quad
\Phi=\begin{pmatrix}S_1\\\epsilon S_2^*\end{pmatrix}
=U\begin{pmatrix}X\\Y\end{pmatrix},\qquad
U^\dagger A U=\mathrm{diag}(m_X^2,m_Y^2).
\]

`X` and `Y` are **distinct heavy complex SU(2) doublets**, both with `Y=-1/2`.
The Matchete heavy free Lagrangian contains only their diagonal masses `MSX`,
`MSY`. The original off-diagonal quadratic term is represented in the numerical
matrix `A`; it is **not** reinserted into the diagonal Matchete Lagrangian.

The Wolfram script loads the existing ordinary model and all seven extra Matchete
quartic candidates, constructs two native pseudoreal SU(2) bridge tensors, and
rewrites all original **eleven baseline interaction candidates** and **seven extra
quartic candidates**. The resulting Matchete Lagrangian contains the SM, heavy
T3 fermion, heavy X and Y, both Yukawas, and all quartic structures. This is a
**UV interaction-basis verification**: the provided numerical rotation comes
from the scalar-only EFT1 endpoint; the Yukawa vertices are not claimed to be
active in that EFT1.

The script fails if original scalar fields remain, any source vertex is absent,
`CheckLagrangian` fails, or the detailed `HeavyMassBasis` check is not `True`.

## Key limitations

- Both native Matchete bridge Clebsches are individually checked to be
  proportional to the antisymmetric epsilon tensor, aligned to the fixed
  numerical epsilon component convention and verified unitary. This is a
  component-order check *inside Matchete*, not an independent physical
  Matchete-to-RGBeta source-tensor map or quartic normalization certificate.
- The numerical unitary matrix is used as an interaction-coordinate choice.
  Matchete's free scalar masses `MSX`, `MSY` are independent symbolic couplings;
  they are not inferred to reproduce a perturbative EFT1 mass RG trajectory.
- No `Match[...]`, scalar threshold Wilson matching, or physical neutrino-mass
  prediction is performed. The matching flags remain false.
- The scalar-only EFT1 sample at 1000 GeV is below **both** benchmark scalar
  thresholds (~1030 and ~1185 GeV). Proper sequential matching must first run
  parameters to the *heavier* matching scale, integrate out Y, run the remaining
  EFT, and match X at its own threshold. Using the 1000 GeV endpoint directly
  as either threshold boundary is not justified.
- This installation contains no Wolfram kernel. Native `CheckLagrangian`
  success can only be established by running on your Windows installation.
- If the full numerical source has changed since the saved mass-basis JSON,
  preparation fails rather than quietly reusing stale coefficients.

## Install (PowerShell, repository root)

```powershell
Expand-Archive -LiteralPath "C:\path\to\T3_Native_Mass_Eigen_Matchete_Batch.zip" -DestinationPath . -Force
python -m pytest -q tests\python\test_general_t3_mass_eigen_matchete.py
python -m scripts.RunGeneralT3MassEigenMatchete --all --prepare-only
python -m scripts.RunGeneralT3MassEigenMatchete --model B_-1
python -m scripts.RunGeneralT3MassEigenMatchete --model C_-1
python -m scripts.RunGeneralT3MassEigenMatchete --all
```

Output directory: `output/general_t3_mass_eigen_matchete/`.

Each model gets a `*_parameters.json`, a Wolfram `*_matchete.json`, and a full
`*_matchete.log`. Summary JSON includes `physical_matching_performed: false`.
A non-zero exit from Wolfram or any failed heavy-basis diagnostic is an actual
blocker; use the log, rather than assuming `HeavyMassBasis` passed.

## Source status

Checked against the public read-only GitHub repository on 2026-10-27:
`Lagrangian/model/{T3Fields,LagrangianBuilder,SU2Invariants}.wl`,
`Lagrangian/interactions/{T3Topology,GeneralScalarOperators}.wl`, and
`RGE/general/generated/general_rgbeta_quartic_registry.json`.
The two prior mixed-basis Python modules and saved JSON schemas were inspected
from the delivered mass-basis ZIP. GitHub was not modified.

</details>

---

<a id="source-t3-mass-eigen-matchete-sm-field-symbol-fix-md"></a>

## Original document: `docs/T3_MASS_EIGEN_MATchete_SM_FIELD_SYMBOL_FIX.md`

**Original SHA-256:** `ef2c0817d3e41c2008161325970744c11f0ce40a28b0c7640e329f47bb15bb4a`

**Status:** Historical record; check the current source and main
documentation before executing commands reproduced below.

<details>
<summary>Expand original: Matchete SM field-label collision fix</summary>

# Matchete SM field-label collision fix

## Verified live failure

The B_-1 run reached `baseline-build-started`, then `LoadModel["SM"]`
raised `DefineField::FieldLabel` with the **numerical 2x2 U matrix** as the
field label. The original Wolfram runner assigned `u = ...` before loading
Matchete's Standard Model, overwriting Matchete's standard up-quark label `u`.
This is a code namespace collision, not a scalar mass-basis physics failure.

## Correction

In `Lagrangian/RunGeneralT3MassEigenMatchete.wl`, all 15 uses of the global
`u` matrix variable are renamed to `$T3MEUnitaryMatrix`. The JSON input key
`unitary_U`, formulas for unitarity/diagonalisation and the full Yukawa/quartic
field substitution are otherwise unchanged. Existing numerical outputs remain
valid. The diagnostic checkpoint runner and Python provenance guard are kept.

## Installation and exact commands (PowerShell)

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath 'C:\path\to\T3_MassEigen_Matchete_SM_u_Symbol_Fix.zip' -DestinationPath . -Force
python -m pytest -q tests\python\test_general_t3_mass_eigen_matchete.py tests\python\test_general_t3_mixed_scalar_mass_basis.py
python -m scripts.RunGeneralT3MassEigenMatchete --model B_-1
python -m scripts.RunGeneralT3MassEigenMatchete --model C_-1
python -m scripts.RunGeneralT3MassEigenMatchete --all
```

If another live error appears, inspect `output/general_t3_mass_eigen_matchete/B_m1_final_matchete.progress.json`
and the matching `.log`; do not interpret the offline tests as a live
Matchete validation. No Wilson matching is performed by this fix.

</details>

---

<a id="source-t3-mass-eigen-provenance-fix-md"></a>

## Original document: `docs/T3_MASS_EIGEN_PROVENANCE_FIX.md`

**Original SHA-256:** `d7dffea3386d2114a0c8231264d5bf0bdfde2e9adb9df883fedcd5ddec74cb40`

**Status:** Historical record; check the current source and main
documentation before executing commands reproduced below.

<details>
<summary>Expand original: T3 mass-eigenstate provenance path correction</summary>

# T3 mass-eigenstate provenance path correction

## Cause

`scripts/RunGeneralT3MixedScalarMassBasis.py` originally serialized provenance
keys via `str(p.relative_to(ROOT))`. On Windows this produced backslash-separated
keys, whereas `scripts/RunGeneralT3MassEigenMatchete.py` looked up fixed
forward-slash keys and falsely reported a changed registry.

## Fix

- Rotation JSON producer now writes portable POSIX-style relative keys.
- Matchete runner normalizes existing Windows and POSIX path keys before
  comparing each SHA-256 against the current file bytes.
- Missing source records, altered digests, traversal paths, malformed values,
  and duplicate aliases still fail closed. No bypass option is introduced.
- Existing rotation results on disk remain readable without regeneration.

## Installation

Unzip from repository root with `Expand-Archive ... -DestinationPath . -Force`.
All files in this archive are complete source files, not patches.

```powershell
python -m pytest -q tests\python\test_general_t3_mass_eigen_matchete.py tests\python\test_general_t3_mixed_scalar_mass_basis.py
python -m scripts.RunGeneralT3MassEigenMatchete --all --prepare-only
python -m scripts.RunGeneralT3MassEigenMatchete --model B_-1
python -m scripts.RunGeneralT3MassEigenMatchete --model C_-1
python -m scripts.RunGeneralT3MassEigenMatchete --all
```

This only fixes provenance; live Matchete heavy-basis validation and Wilson
threshold matching remain unverified until actually executed.

</details>

---

<a id="source-t3-mass-eigen-wolfram-stage-diagnostics-md"></a>

## Original document: `docs/T3_MASS_EIGEN_WOLFRAM_STAGE_DIAGNOSTICS.md`

**Original SHA-256:** `7182b6c159c30a6235933a1509a661799288f9b5d78d88b074715fe6a9328a92`

**Status:** Historical record; check the current source and main
documentation before executing commands reproduced below.

<details>
<summary>Expand original: Mixed-doublet Matchete live-stage diagnostics</summary>

# Mixed-doublet Matchete live-stage diagnostics

## Reason for replacement

The live `B_-1` eigenbasis run reached Wolfram but returned only
`UncaughtMatcheteError`, exit 5. In the previous Wolfram file that status
is written when `UsingFrontEnd[Catch[T3MEPerform[], ...]]` returns `$Failed`,
`$Aborted`, or a non-association. The result did not preserve the failing
stage. Thus **there is no evidence yet whether `HeavyMassBasis` is valid.**
The provenance mismatch and the numerical mass rotation are already solved.

## Changes (experimental only)

- The existing `Lagrangian/RunGeneralT3MassEigenMatchete.wl` is replaced in
  full. It writes a separate stage-progress JSON before each high-risk step
  and emits `T3ME_STAGE: ...` in stdout. `UncaughtMatcheteError` now includes
  the last stage, a short stage history, the kernel result and message list.
- The existing `scripts/RunGeneralT3MassEigenMatchete.py` is replaced in full.
  It deletes stale progress files and reports the checkpoint and most recent
  output when live Wolfram exits unsuccessfully. SHA-256 provenance checks
  are unchanged and Windows/POSIX key normalization is preserved.
- `tests/python/test_general_t3_mass_eigen_matchete.py` is replaced in full.
  The tests enforce the checkpoint contract and mock an uncaught Matchete
  failure with/without progress data.
- This note is new.

The Lagrangian, scalar rotation, source inventory, CG conventions, and
`CheckLagrangian` pass/fail criteria are **unchanged**. This package does not
fix a yet-unidentified native Matchete issue; it exposes where it fails.

## Install (PowerShell, project root)

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "C:\path\to\T3_MassEigen_Wolfram_Stage_Diagnostics.zip" -DestinationPath . -Force
python -m pytest -q tests\python\test_general_t3_mass_eigen_matchete.py tests\python\test_general_t3_mixed_scalar_mass_basis.py
python -m scripts.RunGeneralT3MassEigenMatchete --model B_-1
```

`--all` is **not necessary** while `B_-1` fails. After a failure, inspect:

```powershell
Get-Content .\output\general_t3_mass_eigen_matchete\B_m1_final_matchete.progress.json
Get-Content .\output\general_t3_mass_eigen_matchete\B_m1_final_matchete.json
Get-Content .\output\general_t3_mass_eigen_matchete\B_m1_final_matchete.log -Tail 120
```

The `last_stage` and nearby Wolfram messages identify the specific operation
that needs correction. `B_-1` and `C_-1` still require **successful** ordinary
`CheckLagrangian` and `HeavyMassBasis` before the mass-eigenstate Lagrangian
can be considered validated.

## Verification

67/67 offline Python tests passed in a workspace containing the previous
mass-rotation source and tests. No Wolfram/Matchete kernel is available in
this packaging environment: the native `CheckLagrangian` remains unverified.

</details>

---

<a id="source-t3-matchete-detailed-mass-basis-audit-md"></a>

## Original document: `docs/T3_MATCHETE_DETAILED_MASS_BASIS_AUDIT.md`

**Original SHA-256:** `a39c7e4428f474ddb6516205658f8f02ef31ae7c1a08d5921fa55399de379da7`

**Status:** Historical record; check the current source and main
documentation before executing commands reproduced below.

<details>
<summary>Expand original: T3: Matchete scalar operator audit with installed-version baseline calibration</summary>

# T3: Matchete scalar operator audit with installed-version baseline calibration

## Source of the latest failure

The 43/43 Python tests from the previous ZIP passed, but the user's **live**
Matchete `CheckLagrangian[..., DetailedOutput -> True]` returned:

- ordinary baseline `CheckLagrangian[base]` **True**, established earlier by
  `BuildT3Lagrangian`;
- explicit `Hermiticity`, `GaugeInvariant`, `CanonicallyNormalized`,
  `HeavyMassBasis`, and `AllObjectsDefined` checks **True**;
- `UncontractedIndices`, `HeavyTadpoles` and spin-chain failures **empty**;
- **`Charges`** rather than the `ChargeNeutral` key assumed by the previous parser;
- `GaugeAnomalies` and `ZeroModeMixing` carrying native values other than
  literal `True` or `{}`.

The previous parser mislabelled these last three as `FailedOrUnsupported`,
causing `BaselineDetailedCheckFailed` even in A_+0. This was a parser problem;
no new physics discrepancy was established.

## Corrected acceptance criteria

The replacement `Lagrangian/RunGeneralQuadraticMatcheteAudit.wl` performs two
actual Matchete calls for each candidate Lagrangian:

1. `CheckLagrangian[lag]` (the ordinary authoritative full check), and
2. `CheckLagrangian[lag, DetailedOutput -> True]` (per-check diagnostic data).

For the original **unextended baseline** it requires the ordinary check to be
`True` and explicit structural diagnostics to pass. It saves exact `InputForm`
values of native diagnostic entries whose pass/fail semantics have not been
established from the installed Matchete version.

For each **individual additional quartic**, **individual quadratic** and the
**combined extra scalar Lagrangian** it demands:

- expected ordinary check: `True` unless the extra B/C(-1) scalar mixing is
  present, in which case the ordinary check must be `False`;
- explicit passing values for known structural tests: Hermiticity, no open
  spin chains or uncontracted indices, canonical normalization, no heavy
  tadpoles, all objects defined, gauge invariance;
- `HeavyMassBasis == True` for ordinary cases and `False` for the mixed
  B/C(-1) quadratic cases;
- the **exact same** native opaque `InputForm` diagnostic values as the
  independently passing baseline, including `Charges`, `GaugeAnomalies`
  and `ZeroModeMixing` when exposed by that installation.

A new/unrecognised native diagnostic key, or a different opaque value,
causes a **failure**, not an automatic pass. Individual and combined reports
include both the normal Matchete result and the full recorded native values.

This is deliberately a **baseline-anchored consistency check**, not a proof
that the semantic meaning of every version-specific Matchete return object
has been decoded. In particular the B/C(-1) mixed operator still cannot be
matched in the original non-diagonal heavy mass basis.

## Files in the ZIP

- `Lagrangian/RunGeneralQuadraticMatcheteAudit.wl` — full replacement
- `scripts/RunGeneralT3QuadraticMatcheteAudit.py` — full replacement
- `tests/python/test_general_t3_matchete_mass_basis.py` — full replacement
- `docs/T3_MATCHETE_DETAILED_MASS_BASIS_AUDIT.md` — this full replacement

No production physics files or GitHub code are modified.

## Commands (PowerShell)

From `C:\thesistemp\T3neutrinomassmodelpipeline`:

```powershell
Expand-Archive -LiteralPath "C:\path\to\T3_Matchete_Baseline_Anchored_Audit_Fix.zip" -DestinationPath . -Force
python -m pytest -q tests\python\test_general_t3_matchete_mass_basis.py
python -m scripts.RunGeneralT3QuadraticMatcheteAudit --model A_+0
python -m scripts.RunGeneralT3QuadraticMatcheteAudit --model B_-1
python -m scripts.RunGeneralT3QuadraticMatcheteAudit --model C_-1
python -m scripts.RunGeneralT3QuadraticMatcheteAudit --all
```

Generated diagnostics are written to
`output/general_t3_quadratic_eft1/matchete/`.

## Certification boundary

A successful experimental audit should report 14 ordinary heavy-basis checks
and two (`B_-1`, `C_-1`) baseline-anchored operator checks with
`HeavyMassBasis` specifically blocked. No Matchete-RGBeta physical coefficient
normalization, mixed scalar mass eigenbasis, or threshold matching is
certified. Existing numerical EFT1 evolutions are not modified.

If the native `Charges` / `GaugeAnomalies` / `ZeroModeMixing` entries legitimately
vary between a baseline and an extended Lagrangian, the audit rejects the
comparison and prints both versions as `InputForm`. That difference must be
analysed rather than forcibly whitelisted.

</details>

---

<a id="source-t3-matchete-eigenfield-typed-index-and-baseline-isolation-md"></a>

## Original document: `docs/T3_MATCHETE_EIGENFIELD_TYPED_INDEX_AND_BASELINE_ISOLATION.md`

**Original SHA-256:** `5a034137bee5a9bdb598fe9a8b8f10d0c9a1bb7c2095436ef01520cf132ba15f`

**Status:** Historical record; check the current source and main
documentation before executing commands reproduced below.

<details>
<summary>Expand original: Matchete mass-eigenstate: typed-field transport and baseline isolation</summary>

# Matchete mass-eigenstate: typed-field transport and baseline isolation

**Scope:** `B_-1`, `C_-1`, experimental `RunGeneralT3MassEigenMatchete` only.

## Observed Windows live failure

After the symbolic mixing-spurion and bridge-index batches, the live run passes original model construction, all seven additional quartic source registrations, SU(2) bridge setup, heavy eigenfield registration, and eigenfield Lagrangian assembly. The *combined* `CheckLagrangian` aborts with `DefineCG::notInvariant`; the previous isolator only identifies `Baseline11Combined`, which is not an individual operator. Neither `HeavyMassBasis` nor matching is certified.

## Source-level correction

The original source's already-evaluated scalar object is `Field[NewScalar1, Scalar, {Index[i,T3BSMd2]}, {}]`. The previous rule sent the **already typed** `Index[i,T3BSMd2]` back through the public accessor `NewScalarX[...]`, risking an `Index[Index[i,T3BSMd2],T3BSMd2]` nested index. The new rule builds the native `Field[NewScalarX,Scalar,{index},{}]` using the pre-existing typed index. It does not alter the SU(2) epsilon bridge or the numerical mass rotation. The S2-side new dummy label remains passed to the public field accessor exactly once, and the bridge receives its independently typed `Index[j,T3BSMd2]`.

The Wolfram script checks the shapes of new scalar `Field` indices immediately after the interaction substitution and again after final `RelabelIndices` on the complete Lagrangian. The preflight fails on any non-`T3BSMd2` or doubly wrapped field index, and still rejects untyped bridge CG slots.

## Improved live failure isolation

If `CheckLagrangian` still fails, the runner now checks `LSM + FreeLag[F,X,Y]` by itself, then partitions the **already validated** original `built["LInt"]` by coupling symbol and CG label, and checks these disjoint baseline partitions before any of the seven additional quartics. This avoids rebuilding/duplicating `DefineCG` registrations. The partitioner verifies its parts sum back to the exact expanded original interaction Lagrangian. It prints each tested signature, and returns the first failing one in `isolated_source_failure`.

If the baseline partition check is inconclusive because a partition is not separately Hermitian, its failure must be interpreted with the message/log, **not** as a proof that the physical interaction is gauge non-invariant. `DefineCG::notInvariant` specifically indicates a tensor issue; `CheckLagrangian == False` alone may have other causes.

## Installation paths

Full replacements:

- `Lagrangian/RunGeneralT3MassEigenMatchete.wl`
- `tests/python/test_general_t3_mass_eigen_matchete.py`

New documentation:

- `docs/T3_MATCHETE_EIGENFIELD_TYPED_INDEX_AND_BASELINE_ISOLATION.md`

No change to the Python Matchete runner, numerical mass rotation, RGEs, production Lagrangian, or source hash logic. No GitHub changes.

## Commands

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath 'C:\path\to\T3_Matchete_Eigenfield_Index_and_Isolation_Fix.zip' -DestinationPath . -Force
python -m pytest -q tests\python\test_general_t3_mass_eigen_matchete.py tests\python\test_general_t3_mixed_scalar_mass_basis.py
python -m scripts.RunGeneralT3MassEigenMatchete --model B_-1
python -m scripts.RunGeneralT3MassEigenMatchete --model C_-1
```

Expected diagnostic destinations:

- `output/general_t3_mass_eigen_matchete/B_m1_final_matchete.log`
- `output/general_t3_mass_eigen_matchete/B_m1_final_matchete.json`
- `output/general_t3_mass_eigen_matchete/C_m1_final_matchete.log`
- `output/general_t3_mass_eigen_matchete/C_m1_final_matchete.json`

Offline verification: Python test suite and standalone representative Wolfram-language index-shape and expression-partition checks. **Live Matchete has not been run on the user's environment.** A real `HeavyMassBasis` pass is not claimed. External numeric U phases, Matchete–RGBeta normalization, and two-threshold matching are still uncertified.

</details>

---

<a id="source-t3-matchete-internal-objects-fix-md"></a>

## Original document: `docs/T3_MATCHETE_INTERNAL_OBJECTS_FIX.md`

**Original SHA-256:** `e7b557e0e63bd780050aa12ec3a3b771feb822ddb65295ab2f21081388aa0fe1`

**Status:** Historical record; check the current source and main
documentation before executing commands reproduced below.

<details>
<summary>Expand original: T3 mixed scalar native Matchete: internal objects compatibility fix</summary>

# T3 mixed scalar native Matchete: internal objects compatibility fix

## Trigger

Both B_-1 and C_-1 exited at `original-interaction-assembly-started` with
`AdditionalQuarticCouplingsAbsentInOriginalSource`, listing the same seven
`lambdaT3ExtraQ...` labels. Live logs already showed all seven independent
Matchete quartic candidates generated.

## Verified cause

Matchete's user-facing shortcuts are converted to internal expressions.
Its documented examples show `phi[]` represented by `Field[phi,Scalar,{},{}]`
and `y[]` represented by `Coupling[y,{},0]`. The old check looked for
`lambdaT3ExtraQ...[___]` as an applied head and therefore failed when the
label appeared as an argument to `Coupling`. The old field transformation also
matched only unevaluated `NewScalar1[i]`/`NewScalar2[i]` shortcuts.

## Changes

Full replacement `Lagrangian/RunGeneralT3MassEigenMatchete.wl`:

- Coupling-presence check searches for the exact symbol *at any depth* in the
  already evaluated Matchete expression. The kernel self-test covers both
  shorthand application and `Coupling[label,{},0]`, plus absent-label controls.
- Old-field detection searches for the exact `NewScalar1` and `NewScalar2`
  labels, including their occurrence inside internal `Field[...]` objects.
- Four scalar rotation rules now match `Field[NewScalar1,...]`,
  `Bar[Field[NewScalar1,...]]` and the analogous `NewScalar2` expressions.
- Representative Matchete field/coupling signatures are printed before the
  SU(2) bridge. If the installed version has a different conjugation object
  or index structure, the runner fails closed and reports these structures.
- Existing UV/EFT1 numerical data, mixing matrix, SU(2) bridge checks,
  physical-matching flags and source-provenance checks are not changed.

## Scope and limitations

Offline tests and a kernel-only symbolic mock check pass. No local Matchete
installation was available for a live run. This change does NOT certify the
new heavy-field mass basis, the Matchete–RGBeta normalization map or physical
Wilson-coefficient matching.

If `OriginalScalarFieldsRemainAfterRotation` occurs, inspect the printed
`Raw original Matchete signatures`; they will identify whether the installed
Matchete encodes conjugate fields differently from the expected representation.
Do not disable source-inventory or heavy-mass-basis checks.

## Installation and execution (PowerShell)

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath 'C:\path\to\T3_Matchete_Internal_Objects_Fix.zip' -DestinationPath . -Force
python -m pytest -q tests\python\test_general_t3_mass_eigen_matchete.py tests\python\test_general_t3_mixed_scalar_mass_basis.py
python -m scripts.RunGeneralT3MassEigenMatchete --model B_-1
python -m scripts.RunGeneralT3MassEigenMatchete --model C_-1
```

Generated run logs and progress files remain under
`output/general_t3_mass_eigen_matchete/`.

</details>

---

<a id="source-t3-matchete-native-h2-portal-cg-fix-md"></a>

## Original document: `docs/T3_MATCHETE_NATIVE_H2_PORTAL_CG_FIX.md`

**Original SHA-256:** `8f0b7361427a4e8a4397867db458a8902394c36383af1f817af2812224c9a20f`

**Status:** Historical record; check the current source and main
documentation before executing commands reproduced below.

<details>
<summary>Expand original: Native mass-eigenstate Matchete: both S2 Higgs-portal CG tensors</summary>

# Native mass-eigenstate Matchete: both S2 Higgs-portal CG tensors

## Live failure (27 October 2026)

Both `B_-1` and `C_-1` pass the isolated diagonal free Lagrangian, two S1 Higgs portals, the reconstructed `lambdaT3` topology and the S1 self-interaction. The first remaining failed partition is:

```
DefineCG::notInvariant: The tensor is not invariant under group transformations.
BaselinePart6 lambdaH2Inv1|CG:T3HiggsPortal2CGInv1 -> $Aborted
```

The previous direct reconstruction of `lambdaT3` **passed the live isolated Matchete test** and is retained here.

## Change

The S2 portal includes two S2 legs. Naively transforming each field leaves the old portal CG multiplied by two epsilon bridges in the resulting Matchete expression. Both independent portal structures are now reconstructed as native four-leg CGs via the exact tensor relation

```
(Hbar, H, Bar[S2], S2) -> (Hbar, H, X_or_Y, Bar[X_or_Y])
T_new[i,j,k,l] = sum_{a,b} T_old[i,j,a,b] eps[a,k] eps[b,l]
```

The old/new basis and orientation are obtained from the existing `PhysicalSU2InvariantBasis` implementation, not from guessed Clebsch coefficients. Both original and new invariant spaces have dimension 2. The runner rejects any missing rank or span. It registers `T3MENativeH2CGInv1/Inv2`, retains the real source couplings `lambdaH2Inv1/Inv2`, and uses the existing four symbolic U entries and their conjugates. It checks exact source-term reconstruction before excluding the old portal terms and adding the projected native terms. The isolated-operator diagnostic uses the same direct native portal expressions.

The success JSON now includes `native_H2_portal_epsilon_tensors` and the Python validator rejects missing or inconsistent invariant-count, rank, and exact-projection reports.

## Verification / nonclaims

- 92/92 offline Python tests passed.
- Actual Matchete `InvariantTensors` checked in a remote Wolfram kernel: old and new invariant space dimensions both 2, projected rank 2, projected tensors lie in the native two-dimensional span.
- Full native `CheckLagrangian` is **NOT VERIFIED**: the remote kernel cannot run the front-end-dependent Matchete validation. This is a proposed fix for the user's observed live failure, not a certified live pass.
- The transformation can reveal further CG problems with other S2-dependent operators. Fail-closed diagnostics remain active.
- Numerical mass eigenvalues and U, EFT1, threshold matching and Wilson coefficients are unchanged. No physical coupling normalization is claimed.

## Full replacement project paths

- `Lagrangian/RunGeneralT3MassEigenMatchete.wl`
- `scripts/RunGeneralT3MassEigenMatchete.py`
- `tests/python/test_general_t3_mass_eigen_matchete.py`
- `docs/T3_MATCHETE_NATIVE_H2_PORTAL_CG_FIX.md` (new)

## PowerShell

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "C:\path\to\T3_Matchete_Native_H2_Portal_CG_Fix.zip" -DestinationPath . -Force
python -m pytest -q tests\python\test_general_t3_mass_eigen_matchete.py tests\python\test_general_t3_mixed_scalar_mass_basis.py
python -m scripts.RunGeneralT3MassEigenMatchete --model B_-1
python -m scripts.RunGeneralT3MassEigenMatchete --model C_-1
```

</details>

---

<a id="source-t3-matchete-native-s2-scalar-quartics-md"></a>

## Original document: `docs/T3_MATCHETE_NATIVE_S2_SCALAR_QUARTICS.md`

**Original SHA-256:** `a4c049195a9387aae9d456d4ba86e56d776c6f89ef4e97bed8cb132376c49fee`

**Status:** Historical record; check the current source and main
documentation before executing commands reproduced below.

<details>
<summary>Expand original: Native SU(2) mass-basis tensors for S2-dependent baseline scalar quartics</summary>

# Native SU(2) mass-basis tensors for S2-dependent baseline scalar quartics

## Scope

Experimental `B_-1` and `C_-1` native Matchete mass-eigenbasis validation.
This package replaces only the following full project files:

- `Lagrangian/RunGeneralT3MassEigenMatchete.wl`
- `scripts/RunGeneralT3MassEigenMatchete.py`
- `tests/python/test_general_t3_mass_eigen_matchete.py`

This document is new at `docs/T3_MATCHETE_NATIVE_S2_SCALAR_QUARTICS.md`.
No GitHub files are changed. The existing numerical mass-basis JSON and EFT1 output are unchanged.

## Live failure that motivated the change

Both models passed the native reconstructed `lambdaT3` and `lambdaH2Inv1/Inv2` tests but failed `DefineCG::notInvariant` for `BaselinePart8 lambda12Inv1|CG:T3CrossScalarCGInv1`. The old four-leg CG multiplied by independent epsilon bridge CGs is not always accepted by the final Matchete invariant check.

## Implementation

The original source Lagrangian is strictly partitioned into the already-native topology and H2 portals; all remaining `lambda12Inv1`, `lambda12Inv2`, and `lambdaS2Inv1` terms are then separated with an exact algebraic reconstruction check. No source coupling is deleted without native reconstruction.

For `lambda12Inv1/Inv2`, the old tensor with legs `(S1bar, S1, S2bar, S2)` is contracted with epsilon in its last two slots to obtain the direct mass-basis invariant with legs `(Xbar, X, X, Xbar)` (where `X` denotes either eigenfield `NewScalarX` or `NewScalarY`). The two projected tensors and the new invariant basis must both have rank 2, and their union must also have rank 2.

For `lambdaS2Inv1`, all four old scalar legs are epsilon-transported from `(S2bar,S2,S2bar,S2)` to `(X,Xbar,X,Xbar)`. Both source and target are projected onto their correct symmetric-boson subspaces; the single projected tensor must lie in the single native invariant span.

The exact symbolic coefficients follow `S1 = U11 X + U12 Y`, `S2bar = - U21 epsilon X - U22 epsilon Y`, and their conjugates. For all three quartics the two/four minus signs cancel. No Matchete/RGBeta normalization map is inferred.

The actual new CG tensors are registered using `DefineCG` and the new basis's `CGReps`. The fail-closed result JSON now includes `native_cross_scalar_epsilon_tensors` and `native_S2_self_epsilon_tensor`. Python validates their counts and exact-rank flags alongside the existing topology and H2 certificates. Isolation diagnostics use the same reconstructed expressions as the final Lagrangian.

## Verified and not verified

- The Python regression suite passed (99 tests in this package's combined local suite).
- Wolfram-Language delimiters and Python syntax compiled.
- A temporary native Matchete kernel returned projected/native span ranks `(2,2)` and `(1,1)` and exact Bose symmetry for the self tensor.
- **Not verified:** full `CheckLagrangian` and `HeavyMassBasis` in the user's Matchete instance; no numerical threshold or Wilson matching has been performed.

## Run

From the repository root (PowerShell):

```powershell
Expand-Archive -LiteralPath "C:\path\to\T3_Matchete_Native_S2_Scalar_Quartics_Fix.zip" -DestinationPath . -Force
python -m pytest -q tests\python\test_general_t3_mass_eigen_matchete.py tests\python\test_general_t3_mixed_scalar_mass_basis.py
python -m scripts.RunGeneralT3MassEigenMatchete --model B_-1
python -m scripts.RunGeneralT3MassEigenMatchete --model C_-1
```

For a failure, the runner prints the first failing isolated operator. The corresponding full log is in `output/general_t3_mass_eigen_matchete`.

</details>

---

<a id="source-t3-matchete-native-topology-cg-fix-md"></a>

## Original document: `docs/T3_MATCHETE_NATIVE_TOPOLOGY_CG_FIX.md`

**Original SHA-256:** `eb55e5215ba508f2bbcb881f88fa84d3b7cf2c1b572cc1069846dd7d503086d6`

**Status:** Historical record; check the current source and main
documentation before executing commands reproduced below.

<details>
<summary>Expand original: Native HHXX topology CG reconstruction — experimental B_-1/C_-1</summary>

# Native HHXX topology CG reconstruction — experimental B_-1/C_-1

## Evidence from local Windows logs

On both mixed models, the original 11-interaction baseline is valid, and the diagonal `FreeLag` passes. The isolated transformed Higgs portal groups also pass. The first rejected transformed term is `lambdaT3|CG:T3MixCG` (`DefineCG::notInvariant`). There is no evidence in those logs that the additional seven quartic operators have failed independently.

## Physics / source change

Original topology vertex: `PlusHc[lambdaT3[] H H S1 Bar[S2] CG[T3MixCG,...]]`.

For `B_-1/C_-1`, define `Phi=(S1,epsilon Conjugate[S2])^T = U (X,Y)^T`, with `epsilon={{0,1},{-1,0}}`. The original `Bar[S2]` is `-epsilon (U21 X+U22 Y)` in this component convention.

We now *remove only the original lambdaT3-expanded terms* from the interaction expression before applying the generic field substitution. We re-express them in a **new Matchete-native four-scalar invariant**:

```
O_HH := H[i] H[j] CG[T3MENativeHHMassCG,{i,j,k,m}]
V_topology in new fields proportional to
    -lambdaT3 O_HH [ U11 U21 X[k] X[m]
                    +(U11 U22 + U12 U21) X[k] Y[m]
                    +U12 U22 Y[k] Y[m] ] + h.c.
```

The code computes the exact scalar tensor factor `alphaCG` relating the original native `H H S1 Bar[S2]` tensor contracted with `epsilon` to the new Matchete tensor. It explicitly requires full array equality and last-pair bosonic symmetry. No normalization is guessed. The four U entries remain symbolic Matchete couplings with the earlier separate numerical mapping.

The remaining 10 baseline quartics, 7 additional quartics and 2 Yukawa vertices remain in the pre-existing rotation pathway. Nothing is changed in RGBeta, EFT1 running or the numerical mass eigenvalues.

## Standalone Matchete check performed in this environment

Temporarily installed Matchete in a stateless Wolfram kernel and evaluated its native `InvariantTensors` routines for the relevant SU(2) reps. They returned one original tensor and one new tensor. Contracting the old tensor's fourth index with `epsilon` gave exactly the new tensor (factor `1`), and the new tensor was symmetric in the two scalar slots:

```
{oldInvariantCount, newInvariantCount, factor, tensorIdentity, bosonicSymmetry}
= {1,1,1,True,True}
```

The same environment **cannot run the complete `CheckLagrangian`**: its Matchete model-check procedure returned `$Failed` with `FrontEndObject::notavail`. It also has no access to the local Windows project and may use a different Matchete package version. This is a native invariant-tensor check, not an end-to-end certification.

## Validation and limitations

- The Python runner now requires the Mathematica export `native_topology_epsilon_tensor` to state that the exact tensor identity and scalar-pair symmetry were verified. Failing or missing evidence is rejected.
- The existing strict inventories, checks for eliminated old scalar fields, and actual `CheckLagrangian`/`HeavyMassBasis` requirements are unchanged.
- If a different operator now fails, the failure-isolation path uses the same newly built native topology expression and reports the next specific coupling/CG group.
- Matchete-to-RGBeta normalisation, two-scale threshold matching and physical Wilson coefficients remain **uncertified/not performed**.

## Install in project root (PowerShell)

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath 'C:\path\to\T3_Matchete_Native_Topology_CG_Fix.zip' -DestinationPath . -Force
python -m pytest -q tests\python\test_general_t3_mass_eigen_matchete.py tests\python\test_general_t3_mixed_scalar_mass_basis.py
python -m scripts.RunGeneralT3MassEigenMatchete --model B_-1
python -m scripts.RunGeneralT3MassEigenMatchete --model C_-1
```

Results/logs: `output/general_t3_mass_eigen_matchete/`. This is a diagnostic mass-basis validation, **not matching**.

</details>

---

<a id="source-t3-matchete-native-yukawa1-and-extra-quartics-md"></a>

## Original document: `docs/T3_MATCHETE_NATIVE_YUKAWA1_AND_EXTRA_QUARTICS.md`

**Original SHA-256:** `65e1a12012f841cd87b00705ac50d43c55c4a2fca1d01d421fb1ddf631f5c322`

**Status:** Historical record; check the current source and main
documentation before executing commands reproduced below.

<details>
<summary>Expand original: B_-1 / C_-1 native Matchete y1 and additional-quartic recoupling</summary>

# B_-1 / C_-1 native Matchete y1 and additional-quartic recoupling

Experimental mass-eigenstate Matchete source replacement, 2026-10-27.

## Verified live input motivating this change

The user's preceding run passed the diagonal free Lagrangian and the first eleven disjoint baseline interaction groups. `BaselinePart12 y1|CG:T3Y1CG` returned `$Failed` for both B_-1 and C_-1. The complete mass-eigenstate `CheckLagrangian` therefore remains unverified.

## Code changes

- **Rebuild `y1` as a native `PlusHc` Yukawa vertex.** The heavy scalar is S1 (not S2) before transport, so the scalar doublet index needs no epsilon: S1 -> U11 X + U12 Y. The vertex retains the original `NCM[Bar[l], CConj[F]]` structure and carries a single direct `T3MENativeY1CG` invariant. One independent SU(2) invariant for dF=1 and dF=3; exact original/new component equality is checked.
- **Rebuild all seven additional scalar quartic couplings with native four-leg CG tensors.** Original species are ordered H, Hdag, S1, S1dag, S2, S2dag. For each generated operator, use the original, boson-symmetrized physical tensor from `T3ExtendedBosonicTensorBasis`, contract its S2/S2dag tensor legs with canonical epsilon, and verify its nonzero rank lies in `PhysicalSU2InvariantBasis` for the new mass-field conjugation orientations. Register one direct CG per source coupling; expand X/Y mixtures only on the external field factors.
- **Preserve source inventory.** Losslessly separate both y1 and all seven extra-quartic source terms, including Hermitian conjugates. Reassemble the eigenstate Lagrangian from the existing direct-CG baseline terms, native y1/y2 vertices and native quartics, with no old scalar labels and with the same source coefficients.
- **Keep failure isolation.** Isolated baseline y1 and each of the seven extra quartics now use the exact same new expression as the full Lagrangian. Do not interpret an isolated check as a substitute for full heavy-mass-basis validation.
- **Fail-closed output schema.** Two new certificates, `native_yukawa1_tensor` and `native_additional_quartic_tensors`, are checked in Python, including all seven exact expected quartic keys, coupling labels, invariant ranks and epsilon-leg counts. Physical Matchete-to-RGBeta coupling normalization remains uncertified.

## Numerical basis and physical interpretation

For Phi=(S1,epsilon conjugate(S2))^T = U (X,Y)^T, both X and Y are heavy, complex SU(2) doublets with Y=-1/2. Source-leg transformations are:

```
S1       = U11 X + U12 Y
S1dag    = conj(U11) Xdag + conj(U12) Ydag
S2       = -eps [conj(U21) Xdag + conj(U22) Ydag]
S2dag    = -eps [U21 X + U22 Y]
```

The SU(2) rank/epsilon checks are mathematical group-tensor tests only; these do not certify the physical phase/normalization between RGBeta and Matchete.

## Verification performed

- Local Python regression: **111 passed** (includes the unchanged mixed-scalar numerical tests).
- Python source compile: passed.
- WL source lexical bracket/comment/string balance: passed.
- Temporary installed Matchete `InvariantTensors` check: B_-1 and C_-1 y1 each have one invariant, and for all six independent *structures* (seven couplings) the projected four-leg tensors lie in the mass-basis invariant span with rank 2, exactly. The Matchete kernel has no working front end for `CheckLagrangian`; a complete live test has **not** passed here.

## Install from project root (PowerShell)

```powershell
Expand-Archive -LiteralPath "C:\path\to\T3_Matchete_Native_Y1_And_Extra_Quartics.zip" -DestinationPath . -Force
python -m pytest -q tests\python\test_general_t3_mass_eigen_matchete.py tests\python\test_general_t3_mixed_scalar_mass_basis.py
python -m scripts.RunGeneralT3MassEigenMatchete --model B_-1
python -m scripts.RunGeneralT3MassEigenMatchete --model C_-1
```

The next run must report both ordinary `CheckLagrangian=True` and detailed `HeavyMassBasis=True` to pass. No threshold Wilson coefficients are computed and neither numerical EFT1 data nor numerical mass rotation needs regeneration.

</details>

---

<a id="source-t3-matchete-native-yukawa2-cg-fix-md"></a>

## Original document: `docs/T3_MATCHETE_NATIVE_YUKAWA2_CG_FIX.md`

**Original SHA-256:** `d476eea0ea4ed8ca2e862bd32b03764303614ecf096b89e0bb2df8cb636b9041`

**Status:** Historical record; check the current source and main
documentation before executing commands reproduced below.

<details>
<summary>Expand original: Native Yukawa2 SU(2) recoupling (experimental)</summary>

# Native Yukawa2 SU(2) recoupling (experimental)

B_-1 and C_-1 pass the free Lagrangian and first ten separated baseline
contractions, but y2|CG:T3Y2CG aborts `CheckLagrangian` with
`DefineCG::notInvariant`. The source y2 interaction contains `Bar[S2]`;
the mass eigenstate X/Y is an unconjugated doublet after the epsilon
intertwiner. Leaving original `T3Y2CG` times an epsilon bridge in the
transformed expression does not produce a native Matchete CG for those
external indices.

The new `T3MEBuildNativeYukawa2` computes, using the source's current Matchete
SU(2) basis, the old tensor and new tensor, contracts epsilon on the old
scalar slot, and checks every component of `T_old eps = alphaCG T_native`.
The sign alphaCG is **calculated**; in a temporary native Matchete tensor
kernel it was -1 for the B singlet-fermion variant and +1 for the C
triplet-fermion variant. The formula contains its own explicit - sign
from `Bar[S2] = -epsilon (epsilon Bar[S2])` as in the mass-rotation
convention. The new tensor replaces all source y2 terms, including their
Hermitian-conjugate contributions, after lossless algebraic source separation.
The y1 source tensor needs no epsilon transport and remains in the generic
rotation; it is checked by the existing per-operator validator.

## Validation

- Native Matchete `InvariantTensors`: old and new spaces one-dimensional
  for B (dF=1) and C (dF=3); entrywise tensor identity verified.
- Offline Python and source structure tests check the output certificate,
  source completeness, dimensions, and strict failure handling.
- **Full Matchete CheckLagrangian/HeavyMassBasis not yet live-verified**.
- Does not perform physical matching or certify RGBeta coupling normalization.

## Project paths

- `Lagrangian/RunGeneralT3MassEigenMatchete.wl` (replacement)
- `scripts/RunGeneralT3MassEigenMatchete.py` (replacement)
- `tests/python/test_general_t3_mass_eigen_matchete.py` (replacement)
- `docs/T3_MATCHETE_NATIVE_YUKAWA2_CG_FIX.md` (new)

</details>

---

<a id="source-t3-matchete-symbolic-mixing-invariance-diagnostic-md"></a>

## Original document: `docs/T3_MATCHETE_SYMBOLIC_MIXING_INVARIANCE_DIAGNOSTIC.md`

**Original SHA-256:** `42fab890839e7f379820105be8be1b6b23a079227b7213c18e40cf9ae6c86c21`

**Status:** Historical record; check the current source and main
documentation before executing commands reproduced below.

<details>
<summary>Expand original: T3 B_-1 / C_-1: exact symbolic mixing in native Matchete check</summary>

# T3 B_-1 / C_-1: exact symbolic mixing in native Matchete check

## Trigger

After fixing the typed bridge indices, both models reached `CheckLagrangian`
but aborted with `DefineCG::notInvariant: The tensor is not invariant under group
transformations.` The failure was observed for the assembled transformed
Lagrangian; it has **not** established which vertex, CG orientation, or numerical
operation is at fault.

## Change and rationale

The prior Wolfram runner directly inserted **machine-precision** values of the
2x2 numerical scalar mixing matrix `U` into Matchete expressions. Matchete
constructs and checks group-invariant Clebsch-Gordan tensors exactly. Floating
coefficients might interfere with that exact CG reduction; the existing
messages alone do **not** prove that this is the cause.

This experimental runner now defines four neutral complex Matchete couplings

    T3MEU11[], T3MEU12[], T3MEU21[], T3MEU22[]

and transports interactions algebraically with them and their `Bar[...]`
conjugates. The numerical `U` is *still* checked for unitarity and scalar mass
squared diagonalisation, using the unchanged input JSON. A new result field
`rotation_spurion_values_from_numerical_U` contains the exact double-precision
boundary values of all four symbolic couplings. Python checks equality against
the input U on any purported successful Matchete run.

This avoids approximate coefficients in the gauge-CG calculation. **The four
couplings are formally independent inside Matchete**: a successful structural
`CheckLagrangian` alone would not prove they obey unitarity identities or
certify the Matchete/RGBeta component basis map. The separately verified
numerical map supplies those constraints at the benchmark point only.

## Failure isolation

If the full `CheckLagrangian` still fails, the Wolfram runner first tests
the combined **11-vertex baseline**, then tests each of the seven additional
quartic operators separately. It uses already registered source CGs; it does
not regenerate them. It prints the first failing group or quartic and records
a diagnostic in
`isolated_source_failure`. This is an **investigation tool**, never a way to
turn a failed full check into a pass. If even the first individual term fails,
the bridge definition and conjugate-index orientation need inspection against
the installed Matchete group convention, not further numeric fixes.

No Wilson matching is done, no EFT1 data regenerated, no repository GitHub
content changed, and no cross-package physical tensor normalization certified.

## Files

Full replacements:

- `Lagrangian/RunGeneralT3MassEigenMatchete.wl`
- `scripts/RunGeneralT3MassEigenMatchete.py`
- `tests/python/test_general_t3_mass_eigen_matchete.py`

New document:

- `docs/T3_MATCHETE_SYMBOLIC_MIXING_INVARIANCE_DIAGNOSTIC.md`

## PowerShell

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath 'C:\path\to\T3_Matchete_Symbolic_Mixing_CG_Check.zip' -DestinationPath . -Force

python -m pytest -q tests\python\test_general_t3_mass_eigen_matchete.py tests\python\test_general_t3_mixed_scalar_mass_basis.py

python -m scripts.RunGeneralT3MassEigenMatchete --model B_-1
python -m scripts.RunGeneralT3MassEigenMatchete --model C_-1
python -m scripts.RunGeneralT3MassEigenMatchete --all
```

If there is still a failure:

```powershell
Get-Content .\output\general_t3_mass_eigen_matchete\B_m1_final_matchete.log -Tail 100
Get-Content .\output\general_t3_mass_eigen_matchete\B_m1_final_matchete.json
```

Expected offline test status for this package: 82 tests passing in the
packaging environment. **Live installed Matchete remains untested by the author
of this package.**

</details>

---

<a id="source-t3-matchete-typed-bridge-indices-fix-md"></a>

## Original document: `docs/T3_MATCHETE_TYPED_BRIDGE_INDICES_FIX.md`

**Original SHA-256:** `c38366a497323bbb73529eef9a3275ec9a25dd51ce5c5dd0bf0a4a1a5fd4c6c0`

**Status:** Historical record; check the current source and main
documentation before executing commands reproduced below.

<details>
<summary>Expand original: Typed pseudoreal-bridge indices: B_-1 / C_-1</summary>

# Typed pseudoreal-bridge indices: B_-1 / C_-1

## Observed live failure

Both transformed B_-1 and C_-1 Matchete Lagrangians reached the final check,
but the installed Matchete aborted with:

```
DefineCG::invldIndType: j$32088 is not defined as a representation.
```

The original scalar fields have already evaluated into Matchete's internal
form `Field[NewScalar2,Scalar,{Index[d$$3,T3BSMd2]},{}]`.
The bridge interaction inserted a fresh bare `j$...` in
`CG[T3MEBridgeBB|T3MEBridgeUU,{i,j}]`. Matchete cannot infer a
representation for that index when checking the result.

## Change

The two CG slots now contain indices of the same explicit type:

```
CG[T3MEBridgeBB, {i, Index[j,T3BSMd2]}] Bar[NewScalarX[j]]
CG[T3MEBridgeUU, {i, Index[j,T3BSMd2]}] NewScalarX[j]
```

The `j` in `NewScalarX[j]` is deliberately a bare *short-form field index*:
Matchete converts that field shortcut to an internal `Field[..., {Index[j,T3BSMd2]}, ...]`.
The bridge and new field thus contract the same typed dummy index.
Identical changes are made for `NewScalarY`. Original first index `i` is
already a typed `Index[...]` from the source `Field` object.

An additional preflight rejects bridge contractions with untyped or wrongly
typed CG index slots before the ordinary and detailed `CheckLagrangian` calls.
It does not certify numerical component-normalization alignment.

## Install / run

Extract the ZIP into the repository root, replacing the existing Wolfram runner
and test file; no production source or numerical outputs are changed.

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath 'C:\path\to\T3_Matchete_Typed_Bridge_Indices_Fix.zip' -DestinationPath . -Force
python -m pytest -q tests\python\test_general_t3_mass_eigen_matchete.py tests\python\test_general_t3_mixed_scalar_mass_basis.py
python -m scripts.RunGeneralT3MassEigenMatchete --model B_-1
python -m scripts.RunGeneralT3MassEigenMatchete --model C_-1
```

If a new live check fails, the progress JSON and execution log in
`output/general_t3_mass_eigen_matchete/` identify the stage and Matchete messages.

## Validation scope

The source and offline tests verify the typed-index code path. The Matchete
kernel itself was **not** available in the packaging environment. A native
`HeavyMassBasis` pass remains unverified until both live runs succeed.
No scalar-threshold Wilson matching has been performed.

</details>

---

<a id="source-t3-matchete-mixed-mass-basis-md"></a>

## Original document: `docs/T3_MATCHEte_MIXED_MASS_BASIS.md`

**Original SHA-256:** `45488df76ad13d7913bf9586413589a16e43c06f5c64b1a907d2f73253ce014b`

**Status:** Historical record; check the current source and main
documentation before executing commands reproduced below.

<details>
<summary>Expand original: T3 Matchete mixed-heavy-mass audit fix</summary>

# T3 Matchete mixed-heavy-mass audit fix

## Observed failure (B_-1)

The original experiment generated an independent gauge-invariant bilinear
`S1dag S2dag + h.c.` and added it to a *heavy-field* Lagrangian passed to
`CheckLagrangian`. Matchete reported:

`HeavyMassBasisQ::mixedfields` / `CheckLagrangian::HeavyMassBasis`.

This is a **mass-basis obstruction**. It does not by itself establish a
violation of U(1)Y, SU(2)L, or Z2. The same obstruction applies to C_-1.

The model hypercharges are

    Y(S1) = alpha/2 = -1/2,
    Y(S2) = (alpha+2)/2 = +1/2.

Their conjugate bilinear is neutral and their SU(2) doublets admit one
antisymmetric-singlet contraction. A nonzero coefficient mixes scalar
states. To use Matchete for heavy-scalar matching, one must diagonalize the
quadratic Lagrangian and consistently rotate all Yukawa and quartic vertices.

## What this batch changes

**Complete replacement** of:

- `Lagrangian/RunGeneralQuadraticMatcheteAudit.wl`
- `scripts/RunGeneralT3QuadraticMatcheteAudit.py`

**New**:

- `tests/python/test_general_t3_matchete_mass_basis.py`
- `docs/T3_MATCHEte_MIXED_MASS_BASIS.md`

For all 16 models the audit builds and checks the original heavy-field
Lagrangian. It generates all additional candidate quartics and quadratic
bilinears in the existing Matchete generator and checks the expanded
Lagrangian.

For the 14 models other than B_-1/C_-1, it validates in the original
heavy-field basis. For B_-1/C_-1, it uses the *existing* repository helper
`ReclassifyT3Fields[model, {"F"}, {"F", "S1", "S2"}]`, which classifies S1/S2
as light but keeps F heavy. It rebuilds the free Lagrangian and checks the
entire coupled interaction in the **pre-diagonalization scalar basis**.
This deliberately bypasses only the requirement that external heavy fields
enter Matchete already in their matching mass basis. It does *not* perform
matching or provide a physical mass-basis certificate.

The Python wrapper checks all JSON flags and refuses to label the two
pre-diagonalization checks as heavy mass-basis passes. There must be 14
heavy-basis plus 2 pre-diagonalization audits and 50 quartic/8 quadratic
model-specific slots. An unsuccessful light-basis Matchete run remains an
explicit failure, not a reinterpreted success.

## Commands (PowerShell, project root)

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath 'C:\path\to\T3_Matchete_Mixed_Mass_Basis_Fix.zip' -DestinationPath . -Force
python -m pytest -q tests\python\test_general_t3_matchete_mass_basis.py
python -m scripts.RunGeneralT3QuadraticMatcheteAudit --model B_-1
python -m scripts.RunGeneralT3QuadraticMatcheteAudit --model C_-1
python -m scripts.RunGeneralT3QuadraticMatcheteAudit --all
```

Results: `output/general_t3_quadratic_eft1/matchete/`; combined audit
summary: `matchete_16_model_summary.json`.

## Certainty/status

- **Verified from source**: two doublet mixed-field mass structures; production
  `T3Fields.wl` supports heavy/light reclassification; the original audit
  unconditionally added a mixed mass to a heavy-field Lagrangian.
- **Verified offline**: this replacement's Python syntax and report basis
  flags, model totals, CLI arguments, ZIP structure and tests.
- **Must be tested on your machine**: the live Matchete `UsingFrontEnd` behavior
  after reclassifying S1 and S2 as light. This environment has no Wolfram
  front end and cannot assert that Matchete accepts the reclassified model.
- **Not established**: physical Matchete/RGBeta coupling mapping,
  diagonalized heavy-scalar states, mass-matrix RG conventions,
  generalized scalar matching, or neutrino-mass predictions.

The two models' pre-diagonalization operator test **is not** evidence that
Matchete supports mixed heavy states during matching. The next physics task
is scalar mass-basis transformations (including Yukawa/quartic rotation) and
mass-eigenstate matching.

</details>

---

<a id="source-t3-mixed-scalar-mass-basis-batch-md"></a>

## Original document: `docs/T3_MIXED_SCALAR_MASS_BASIS_BATCH.md`

**Original SHA-256:** `31a4dd3a174a9ac3255e4487198419be2ca66d25a5e744ec380862aa6759ae01`

**Status:** Historical record; check the current source and main
documentation before executing commands reproduced below.

<details>
<summary>Expand original: T3 mixed scalar mass-basis batch — experimental B(-1)/C(-1)</summary>

# T3 mixed scalar mass-basis batch — experimental B(-1)/C(-1)

## Source and scope

Inspected repository `main` definitions:

- `Lagrangian/model/T3ModelCatalog.wl`: for alpha = -1, the two doublets have
  `Y(S1)=-1/2` and `Y(S2)=+1/2`.
- `RGE/running/rgbeta/T3RGBetaModel.wl`: 10 baseline scalar quartic slots
  (including SM lambdaH), 6 pre-existing structural quartics for each model,
  and 2 Yukawa scalar legs (`Bar @ S1` in `y1`, `S2` in `y2`).
- `RGE/general/generated/{B,C}_-1_extra_quartics.wl`: a seventh, independent
  RGBeta candidate `lambdaT3GeneralQ110101I1`, alongside the pre-existing
  `lambdaHHdagS1barS2barCross` contraction of the same field content.
- `Numerical/experimental/GeneralT3ScalarQuadraticMasses.py`: existing
  eight-real-mode Hessian in a specified SU(2) spherical convention.
- `scripts/RunGeneralT3QuadraticEFT1.py`: recorded Stage EFT1 data at
  `output/general_t3_quadratic_eft1/{B,C}m1_eft1_quadratic_running.json`.

No production files are overwritten. No GitHub changes. The two new Python
files operate on saved numerical data, so they do not need Wolfram for execution.

## Exact mixing convention

Set `epsilon=[[0,1],[-1,0]]`, `tilde(S2)=epsilon S2*`. For each weak-isospin
component define `Phi=(S1,tilde(S2))` and 

\[
  V_2=\Phi^\dagger
  \begin{pmatrix}m_1^2&\mu^2\\(\mu^2)^*&m_2^2\end{pmatrix}\Phi,
  \qquad \Phi=U\begin{pmatrix}X\\Y\end{pmatrix}.
\]

All `X,Y` are SU(2) **doublets of hypercharge -1/2**. By construction
`U^dagger A U=diag(m_-^2,m_+^2)` with

\[
 m_\pm^2=\frac{m_1^2+m_2^2}{2}\pm
 \sqrt{\left(\frac{m_1^2-m_2^2}{2}\right)^2+|\mu^2|^2}.
\]

The inverse species substitutions are

\[
 S_1=U_{11}X+U_{12}Y, \quad
 S_2=-\epsilon(U_{21}^*X^*+U_{22}^*Y^*).
\]

Here `U` indices are 1-based in the equations and 0-based in the JSON.
`mu2` is the source coefficient of `S1dag epsilon S2dag+h.c.`; there is NO
extra factor 1/2 for this mixed term. If both squared masses are identical and
`mu2=0`, the eigenbasis is nonunique and reported as degenerate.

## What is transported, and what is not

For **each of the 17 source quartic symbols** in B(-1) and C(-1) we record
one separate multilinear tensor transformation. The basis holds 10 ordinary
couplings, 6 pre-existing structural analogues and the 1 additional RGBeta
candidate. Both original Yukawa scalar legs are transformed too.

For a source interaction `lambda T[i1,...,in] f1[i1]...fn[in]`, every rotated
term stores an **ordered list of new fields and its full 2x2 SU(2) index map**
`M1,...,Mn`. The interaction is

\[
 \lambda\,T_{i_1\ldots i_n}
 M^{(1)}_{i_1j_1}\cdots M^{(n)}_{i_nj_n}
 f'_1{}_{j_1}\cdots f'_n{}_{j_n}.
\]

The source `T` is **preserved rather than replaced by an assumed common
CG**. This is why two different contractions with identical field lists
stay independent. The transformation is valid for **any** supplied original
SU(2) component tensor; the supplied 2x2 matrices can be used when a physical
Matchete/RGBeta common-basis tensor has been certified. The report separately
stores the corresponding source numerical RGBeta couplings. It does not
pretend that these are already physical mass-basis independent coupling
constants.

Only the scalar factor in the Yukawa vertices is transported. The fermion
and lepton SU(2) tensor factors and flavor matrices are unchanged symbols;
they cannot be numerically evaluated from scalar-only EFT1 endpoints because
F has already been integrated out there.

The **whole quadratic mass tensor**, including diagonal terms and the
holomorphic mixing term plus its hermitian conjugate, is independently
transformed. It is checked against the diagonal `U^dagger A U` and against
the existing eight-real-mode Hessian spectrum (each complex-doublet mass
squared repeated four times as a real eigenvalue).

## Install (all files in one ZIP)

Expand the ZIP into `C:\thesistemp\T3neutrinomassmodelpipeline`. The files go to:

- `Numerical/experimental/GeneralT3MixedScalarMassBasis.py` (new)
- `scripts/RunGeneralT3MixedScalarMassBasis.py` (new)
- `tests/python/test_general_t3_mixed_scalar_mass_basis.py` (new)
- `docs/T3_MIXED_SCALAR_MASS_BASIS_BATCH.md` (new)

### Commands (PowerShell, from the repository root)

```powershell
python -m pytest -q tests\python\test_general_t3_mixed_scalar_mass_basis.py
python -m scripts.RunGeneralT3MixedScalarMassBasis --model B_-1
python -m scripts.RunGeneralT3MixedScalarMassBasis --model C_-1
python -m scripts.RunGeneralT3MixedScalarMassBasis --all

# Optional: check the boundary immediately below the fermion threshold
python -m scripts.RunGeneralT3MixedScalarMassBasis --all --which initial
```

If the existing EFT1 inputs are absent, first run the already-installed
`python -m scripts.RunGeneralT3QuadraticEFT1 --all --mu2-fraction 0.01`.

Outputs (one for each model and endpoint):

- `output/general_t3_mass_basis/B_m1_final_mass_basis.json`
- `output/general_t3_mass_basis/C_m1_final_mass_basis.json`

Each report includes the mass eigenvalues, unitary rotation, source couplings,
**all** vertex-index transport terms, the exact quadratic mass tensor check,
SU(2) covariance residual, tensor-identity residuals, provenance hashes and
explicit uncertified matching flags. The printed eigenmasses are *running*
scalar mass parameters, not physical pole masses.

## Verified versus blocked

- **Verified offline:** analytic eigenvalue identity; unitary and inverse map;
  positive/negative mixing phases; SU(2) gauge covariance; mass potential
  equality for arbitrary field configurations; consistency with the existing
  B/C quadratic Hessian; formal quartic tensor transport for 17/17 source
  vertices for each model; Yukawa scalar-leg transport; source field-inventory
  consistency; fail-closed corrupt input cases.
- **Source verified:** field charge assignments, 17 RGBeta scalar couplings,
  2 original Yukawa scalar legs, 7 additional quartic structures with one
  two-dimensional contraction sector, Stage EFT1 report contract.
- **Not certified:** Matchete ↔ RGBeta signed component and coupling prefactors,
  physical field component CGs, native Matchete mass-eigenfield definitions,
  heavy-mass-basis matching, transformed Wilson coefficient matching, and
  low-energy neutrino masses. The code does **not** claim Matchete
  `HeavyMassBasis` succeeds on a redefined Lagrangian.

## Next physical task

Build certified Matchete component tensors and a physical signed coupling map;
install the two mass-eigen heavy doublets into the Matchete Lagrangian with all
transformed vertices and verify `CheckLagrangian`/`HeavyMassBasis` there. Only
then calculate scalar threshold Wilson coefficients with the two nondegenerate
mass eigenstates. **Do not** substitute the present tensor-routing report for
an actual Matchete mass-basis model or loop-matching computation.

</details>

---
