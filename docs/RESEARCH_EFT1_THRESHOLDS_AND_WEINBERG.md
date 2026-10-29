# Research Eft1 Thresholds And Weinberg

> Historical research record. This document combines **original source
> Markdown** under original-path headings. Local Markdown links to moved
> files are updated; exact original bytes are in the external archive.
> This is not a new certification
> of physics, a rewritten conclusion, or a new computational workflow.
> Original bytes and SHA-256 values are stored in the external archive.

## Contents

- [B_-1 / C_-1 scalar-threshold EOMEFT continuation fix](#source-t3-general-eft1-canonical-eom-handoff-fix-md) — `docs/T3_GENERAL_EFT1_CANONICAL_EOM_HANDOFF_FIX.md`
- [Scalar-threshold EOMEFT Hermiticity diagnostic and conservative rewrite](#source-t3-general-eft1-hermiticity-second-threshold-md) — `docs/T3_GENERAL_EFT1_HERMITICITY_SECOND_THRESHOLD.md`
- [T3 EFT1 Hermiticity: exact-sector isolation](#source-t3-general-eft1-hermiticity-sector-isolation-md) — `docs/T3_GENERAL_EFT1_HERMITICITY_SECTOR_ISOLATION.md`
- [Experimental leading-C5-only scalar-threshold mode (`B_-1`, `C_-1`)](#source-t3-general-eft1-leading-c5-scalar-threshold-md) — `docs/T3_GENERAL_EFT1_LEADING_C5_SCALAR_THRESHOLD.md`
- [EFT1 scalar threshold: heavy-field free-Lagrangian rebase and Match audit](#source-t3-general-eft1-scalar-matching-heavy-free-rebase-md) — `docs/T3_GENERAL_EFT1_SCALAR_MATCHING_HEAVY_FREE_REBASE.md`
- [T3 EFT1 scalar-threshold: verified renormalizable reconstruction and quartic isolation](#source-t3-general-eft1-verified-reconstruction-and-quartic-isolation-md) — `docs/T3_GENERAL_EFT1_VERIFIED_RECONSTRUCTION_AND_QUARTIC_ISOLATION.md`
- [General T3 — first fermion-threshold one-loop matching probe](#source-t3-general-fermion-first-eft1-one-loop-md) — `docs/T3_GENERAL_FERMION_FIRST_EFT1_ONE_LOOP.md`
- [T3 B_-1 / C_-1: fermion-first EFT1 tree boundary](#source-t3-general-fermion-first-eft1-tree-boundary-md) — `docs/T3_GENERAL_FERMION_FIRST_EFT1_TREE_BOUNDARY.md`
- [Direct one-loop numerical C5 — experimental mixed-doublet models](#source-t3-general-t3-direct-c5-numerical-benchmark-md) — `docs/T3_GENERAL_T3_DIRECT_C5_NUMERICAL_BENCHMARK.md`
- [Fermion-first EFT1: raw hard LLHH pole/log audit](#source-t3-general-t3-eft1-hard-pole-audit-md) — `docs/T3_GENERAL_T3_EFT1_HARD_POLE_AUDIT.md`
- [Ordinary T3 `B_-1`/`C_-1`: scalar-threshold Matchete probe](#source-t3-general-t3-eft1-scalar-threshold-probe-md) — `docs/T3_GENERAL_T3_EFT1_SCALAR_THRESHOLD_PROBE.md`
- [General T3 one-scale frozen-reference Yukawa transfer (B_-1/C_-1)](#source-t3-general-t3-frozen-reference-yukawa-transfer-md) — `docs/T3_GENERAL_T3_FROZEN_REFERENCE_YUKAWA_TRANSFER.md`
- [T3 B_-1/C_-1: frozen scalar diagonal-aligned direct C5 benchmark](#source-t3-general-t3-frozen-scalar-aligned-one-scale-md) — `docs/T3_GENERAL_T3_FROZEN_SCALAR_ALIGNED_ONE_SCALE.md`
- [General T3: status and first direct one-loop Weinberg probe](#source-t3-general-t3-one-loop-status-and-probe-md) — `docs/T3_GENERAL_T3_ONE_LOOP_STATUS_AND_PROBE.md`
- [T3 B_-1/C_-1: direct one-scale C5 to low-scale SM+Weinberg endpoint](#source-t3-general-t3-one-scale-sm-weinberg-endpoint-md) — `docs/T3_GENERAL_T3_ONE_SCALE_SM_WEINBERG_ENDPOINT.md`
- [Restricted two-threshold fixed-order Weinberg projection (B_-1 and C_-1)](#source-t3-general-two-threshold-numerical-c5-md) — `docs/T3_GENERAL_TWO_THRESHOLD_NUMERICAL_C5.md`
- [T3 B_-1/C_-1: two-threshold raw Weinberg-sector pole cancellation](#source-t3-general-two-threshold-raw-pole-assembly-md) — `docs/T3_GENERAL_TWO_THRESHOLD_RAW_POLE_ASSEMBLY.md`
- [General T3: additional EFT1 scalar-pair Wilson mixing](#source-general-t3-eft1-scalar-pair-mixing-md) — `docs/GENERAL_T3_EFT1_SCALAR_PAIR_MIXING.md`
- [General T3: numerical extra-quartic EFT1 LLHH Wilson-source transport](#source-general-t3-eft1-wilson-source-transport-md) — `docs/GENERAL_T3_EFT1_WILSON_SOURCE_TRANSPORT.md`
- [General T3: grouped fermion-first threshold continuation (experimental)](#source-general-t3-fermion-first-thresholds-md) — `docs/GENERAL_T3_FERMION_FIRST_THRESHOLDS.md`
- [General T3: existing model-specific C5 and final SM+Weinberg running](#source-general-t3-final-c5-and-weinberg-running-md) — `docs/GENERAL_T3_FINAL_C5_AND_WEINBERG_RUNNING.md`
- [General T3: selected extra-quartic Weinberg projection and low-energy transport](#source-general-t3-weinberg-source-transport-md) — `docs/GENERAL_T3_WEINBERG_SOURCE_TRANSPORT.md`

---

<a id="source-t3-general-eft1-canonical-eom-handoff-fix-md"></a>

## Original document: `docs/T3_GENERAL_EFT1_CANONICAL_EOM_HANDOFF_FIX.md`

**Original SHA-256:** `d8862d371571197ef104a21c2e562f66dca8e5f3236c4e925e11aafc22035c4a`

**Status:** Historical record; check the current source and main
documentation before executing commands reproduced below.

<details>
<summary>Expand original: B_-1 / C_-1 scalar-threshold EOMEFT continuation fix</summary>

# B_-1 / C_-1 scalar-threshold EOMEFT continuation fix

## Observed live failure

The scalar-threshold probe failed at `scalar-threshold-EFT1-lagrangian-check-started`, before second-threshold matching, returning `SavedEFT1TreeHeavyScalarCheckFailed` for both B_-1 and C_-1. This does **not** invalidate the successful fermion-first tree/one-loop runs or the 68+68+68 hard-pole audit.

## Verified source-level cause / corrective design

The failed scalar-threshold Python runner supplied `*_eft1_tree.inputform.txt`, namely `matchingResult["MatchedEFT"]` from the first-stage matching. The repository's authoritative `Lagrangian/RunMatching.wl`, `CanonicalTransitionEFT`, explicitly returns `EOMEFT` for the **next threshold**, because `ReplaceEffectiveCouplings` may undo canonical mass/kinetic field redefinitions. This is a verified architectural error in the previous scalar-probe handoff; whether it was the only cause of the live failure is **not yet proven** without rerunning local Matchete.

Separately, the second-stage Wolfram code wrapped `CheckLagrangian` inside `Quiet @ Check[..., $Failed]`. `Check` returns its failure alternative on *any* message, including nonfatal messages, as noted in `Lagrangian/RunMatching.wl`. First-stage Matchete validation already uses `CheckAbort[CheckLagrangian[...],$Aborted]`; the second stage now does likewise and records diagnostics. It does **not** ignore a failed canonicalization or bypass the heavy-mass-basis check.

## Changes

- `Lagrangian/RunGeneralT3FermionFirstEFT1Tree.wl` (full replacement): additionally exports `matchingResult["EOMEFT"]` as `*_eft1_transition_eom.inputform.txt` and certifies its type. Keeps existing expanded `MatchedEFT` and LLSS diagnostics unchanged.
- `scripts/RunGeneralT3FermionFirstEFT1Tree.py` (full replacement): validates the new canonical EOMEFT certificate and artifact; cleans stale transition artifacts before a run.
- `Lagrangian/RunGeneralT3EFT1ScalarThresholdProbe.wl` (full replacement): no longer marks message-emitting successful Matchete checks as `$Failed`. Logs detailed failures and **still fails closed** unless ordinary `CheckLagrangian` and detailed `HeavyMassBasis` checks pass.
- `scripts/RunGeneralT3EFT1ScalarThresholdProbe.py` (full replacement): reads the newly exported canonical EOMEFT from the validated first-stage **tree-only** run; rejects missing/stale transition files. Prior one-loop fermion-threshold outputs remain independently validated and are not rerun.
- Both corresponding `tests/python/test_*.py` are full replacements with canonical handoff regression tests.

The root `Lagrangian/RunMatching.wl` and all Matchete CGs/couplings are untouched. No GitHub changes. The one-loop hard pole/log/finite export and the numerical C5 benchmark code are untouched.

## PowerShell, repository root

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "C:\path\to\T3_General_EFT1_Canonical_EOM_Handoff_Fix.zip" -DestinationPath . -Force

python -m pytest -q tests\python\test_general_t3_fermion_first_eft1_tree.py tests\python\test_general_t3_eft1_scalar_threshold_probe.py

# Preserve existing successful tree-EFT1 outputs before re-running only its cheap tree match:
Copy-Item -Path "output\general_t3_mass_eigen_matchete\eft1_fermion_first_tree" -Destination "output\general_t3_mass_eigen_matchete\eft1_fermion_first_tree_pre_eom_backup" -Recurse -Force

python -m scripts.RunGeneralT3FermionFirstEFT1Tree --all
python -m scripts.RunGeneralT3EFT1ScalarThresholdProbe --model B_-1
python -m scripts.RunGeneralT3EFT1ScalarThresholdProbe --model C_-1
```

The backup command assumes `eft1_fermion_first_tree_pre_eom_backup` does not already exist; give it another name if necessary. **Do not rerun** `RunGeneralT3FermionFirstEFT1OneLoop` unless its prior outputs were removed, since the loop outputs are independent of the EOMEFT export addition.

Successful first-tree re-execution should create for each model:

`output/general_t3_mass_eigen_matchete/eft1_fermion_first_tree/{B,C}_m1_final_eft1_fermion_first_tree_eft1_transition_eom.inputform.txt`

The second-stage scalar-threshold outputs remain in:

`output/general_t3_mass_eigen_matchete/eft1_scalar_threshold_probe/`

## Validation and remaining uncertainty

Offline tests verify the canonical transition contract and reject stale `MatchedEFT` inputs, incomplete EOMEFT, or failed Matchete checks; the actual native scalar-threshold matching and physical sequential C5 remain **unverified** until the user's local Wolfram run. If the second stage fails again, use the new `B_m1_final_eft1_scalar_threshold_probe.log` / `C_m1_...log`; the Matchete messages and detailed check are now visible. The next scientific steps are EFT1 tensor/RGE mapping, hard pole/RGE consistency, scalar matching and physical C5 assembly.

</details>

---

<a id="source-t3-general-eft1-hermiticity-second-threshold-md"></a>

## Original document: `docs/T3_GENERAL_EFT1_HERMITICITY_SECOND_THRESHOLD.md`

**Original SHA-256:** `f8c1c3aa9ad8475c82721b64838f4de74a5f16c5d07d3e1d104d635c81ea51ca`

**Status:** Historical record; check the current source and main
documentation before executing commands reproduced below.

<details>
<summary>Expand original: Scalar-threshold EOMEFT Hermiticity diagnostic and conservative rewrite</summary>

# Scalar-threshold EOMEFT Hermiticity diagnostic and conservative rewrite

## Live finding

After the canonical EOMEFT handoff, local Matchete passed canonical normalization,
heavy mass basis, gauge invariance and object definition in both B_-1 and C_-1.
Only `Hermiticity` was `False`. Consequently, no second-stage Match ran.

## Changes

1. Parse the saved EOMEFT after removing integrated F while X,Y are still registered Light.
   Record the **light-field** Hermiticity result in a compact JSON file.
2. Re-register X,Y Heavy. Require the original `CheckLagrangian` as before.
3. Only if Hermiticity is the sole failed check, attempt the documented identity-preserving
   Matchete simplifiers `HcSimplify` and `GreensSimplify` (with `RelabelIndices`).
4. A rewrite is accepted only if the full ordinary check becomes True, and detailed
   Hermiticity, canonical normalization, heavy mass basis, gauge invariance and
   object definition all become True. No additional HC terms are invented.
5. Export a compact `_hermiticity_diagnostics.json` even on failure. The CLI
   requires a successful, complete Hermiticity certificate for any claimed success.

If no equivalence-preserving simplification resolves the failure, the stage **still fails**
with a diagnostic JSON. In that case inspect `light_EFT1_hermiticity` and the rewrite attempts
before attempting further EFT1 integration. This patch is not a verified live repair.

## Full replacement files

- `Lagrangian/RunGeneralT3EFT1ScalarThresholdProbe.wl`
- `scripts/RunGeneralT3EFT1ScalarThresholdProbe.py`
- `tests/python/test_general_t3_eft1_scalar_threshold_probe.py`
- `docs/T3_GENERAL_EFT1_HERMITICITY_SECOND_THRESHOLD.md`

## PowerShell

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath 'C:\path\to\T3_General_EFT1_Hermiticity_Recovery.zip' -DestinationPath . -Force
python -m pytest -q tests\python\test_general_t3_eft1_scalar_threshold_probe.py
python -m scripts.RunGeneralT3EFT1ScalarThresholdProbe --model B_-1
python -m scripts.RunGeneralT3EFT1ScalarThresholdProbe --model C_-1
```

No need to rerun either fermion-first threshold; this batch consumes existing
canonical tree EOMEFT exports. No files in the production RGE path are changed.

</details>

---

<a id="source-t3-general-eft1-hermiticity-sector-isolation-md"></a>

## Original document: `docs/T3_GENERAL_EFT1_HERMITICITY_SECTOR_ISOLATION.md`

**Original SHA-256:** `b8ffa8cd1e2d5ebf316d996c0a938bc93dc18664d289e4bc1e765c63faf01a22`

**Status:** Historical record; check the current source and main
documentation before executing commands reproduced below.

<details>
<summary>Expand original: T3 EFT1 Hermiticity: exact-sector isolation</summary>

# T3 EFT1 Hermiticity: exact-sector isolation

## Verified live failure

In both `B_-1` and `C_-1`, the canonical first-threshold `EOMEFT` failed
`CheckLagrangian[...,DetailedOutput->True]["Hermiticity"]` while **X and Y
were still Light**. Reclassifying to Heavy did not cause the failure.
`HcSimplify` and `GreensSimplify` did not restore the check. Gauge invariance,
canonical normalization, heavy mass basis and object definitions passed in
the heavy-scalar check.

## This batch

The scalar-threshold probe first tests mathematically identity-preserving
rewrites of the EFT (`RelabelIndices`, `Expand` and `Expand`+`RelabelIndices`)
without changing coefficients. It accepts one only after ordinary and detailed
Matchete checks pass. Existing `HcSimplify` and `GreensSimplify` are retained,
with nonfatal simplifier messages no longer converted into `$Failed`.

If the full check still fails, the script **does not run scalar matching**.
Instead, it expands the *saved* EFT1 into disjoint operator sectors and
performs bounded detailed Matchete checks on each:

- Background (all terms other than LLSS and BSM scalar quartics)
- LLSS (all; XX, XY, YY individually)
- Heavy scalar quartics (all; purely heavy and Higgs portals individually)
- The full expanded EFT and EFT with each of the two operator classes excluded

These tests distinguish Hermiticity of each **complete operator class**;
checking a single chirality term would not be meaningful. They preserve
all original terms in the actual matching input and never add/remove HC
operators. The diagnostic output is not a certificate of pole cancellation,
physical `C5`, or a complete EFT1 RGE.

## Full replacement paths

- `Lagrangian/RunGeneralT3EFT1ScalarThresholdProbe.wl`
- `scripts/RunGeneralT3EFT1ScalarThresholdProbe.py`
- `tests/python/test_general_t3_eft1_scalar_threshold_probe.py`

New document:

- `docs/T3_GENERAL_EFT1_HERMITICITY_SECTOR_ISOLATION.md`

## Windows PowerShell

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath 'C:\path\to\T3_General_EFT1_Hermiticity_Sector_Isolation.zip' -DestinationPath . -Force
python -m pytest -q tests\python\test_general_t3_eft1_scalar_threshold_probe.py
python -m scripts.RunGeneralT3EFT1ScalarThresholdProbe --model B_-1
python -m scripts.RunGeneralT3EFT1ScalarThresholdProbe --model C_-1
```

On a genuine Hermiticity failure inspect:

`output/general_t3_mass_eigen_matchete/eft1_scalar_threshold_probe/B_m1_final_eft1_scalar_threshold_probe_hermiticity_sector_audit.json`

and its `C_m1` counterpart. In particular compare the `hermiticity` keys in
`sector_checks` for `background_no_LLSS_or_BSM_quartics`,
`background_plus_all_LLSS`, and `background_plus_all_BSM_quartics`.
The log contains a concise failure-sector summary. Upload those two small
JSON diagnostics together with the failed model's existing
`_hermiticity_diagnostics.json` for the subsequent exact physics fix.

</details>

---

<a id="source-t3-general-eft1-leading-c5-scalar-threshold-md"></a>

## Original document: `docs/T3_GENERAL_EFT1_LEADING_C5_SCALAR_THRESHOLD.md`

**Original SHA-256:** `378c4900352840beac4ed270b5bc87064853f65f68d81e0abbe0e52127853230`

**Status:** Historical record; check the current source and main
documentation before executing commands reproduced below.

<details>
<summary>Expand original: Experimental leading-C5-only scalar-threshold mode (`B_-1`, `C_-1`)</summary>

# Experimental leading-C5-only scalar-threshold mode (`B_-1`, `C_-1`)

## Why it exists

The **complete** saved first-threshold EOMEFT fails Matchete's Hermiticity
check for both models. The uploaded component diagnostics isolate the failure
in the 122 heavy-scalar-only quartic terms. Background (21 terms), the
tree-matched LLSS sector (20 / 24 terms) and all Higgs-containing quartics
(52 terms) separately pass their sector checks. The full Hermiticity problem
has **not** been fixed.

The option `--c5-relevant-only` is a separate, explicit **leading one-loop
LLHH subtheory** used to unblock a raw scalar-threshold contribution. It
retains the exact saved EFT terms in the following inventory:

- 21 background/kinetic/SM terms;
- 20 (`B_-1`) or 24 (`C_-1`) matched LLSS terms;
- all 52 quartics with two external Higgs scalar fields and two X/Y fields.

The 122 quartic terms involving **four heavy-scalar fields and no Higgs**
are excluded from this particular matching call. Such interactions cannot
contribute to the leading one-loop LLSS -> LLHH graph: with one LLSS
insertion, an additional Higgs-bearing vertex is needed to provide two
external H fields, and inserting an extra four-scalar vertex raises the loop
order. This is a power-counting restriction for the selected LLSS -> LLHH
amplitude, **not** an equality to the full EFT. Contributions to other EFT
operators are not claimed to be reproduced by the reduced matching.

No source matching coefficients or UV input files are rewritten or deleted.
The full EOMEFT remains available unchanged, and normal runs continue to
fail on its genuine Hermiticity violation. With this option, the reduced
Lagrangian must independently pass every Matchete `CheckLagrangian` check.
There is **no** automatic fallback to a restricted theory.

## PowerShell

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
python -m pytest -q tests\python\test_general_t3_eft1_c5_relevant_sector.py tests\python\test_general_t3_eft1_scalar_threshold_probe.py
python -m scripts.RunGeneralT3EFT1ScalarThresholdProbe --model B_-1 --c5-relevant-only
python -m scripts.RunGeneralT3EFT1ScalarThresholdProbe --model C_-1 --c5-relevant-only
```

Results are saved under
`output/general_t3_mass_eigen_matchete/eft1_scalar_threshold_probe/`
with `_c5_relevant` in their names so the original full-theory diagnostics are
not overwritten. The summary is `summary_final_c5_relevant.json`.

## Scientific status

**Implemented:** exact partition of the saved tree expression; check that
122 excluded terms, 20/24 LLSS terms, 52 Higgs-bearing quartics and 21
background terms reproduce the full inventory; independent validation of
the reduced subtheory; 0- and 1-loop scalar matching; raw LLHH export.

**Not implemented/verified:** full-EFT Hermiticity, general scalar-EFT
matching, independent hard pole--EFT1-RGE consistency, EFT1 Wilson running,
MSbar renormalization of the scalar-threshold contribution, or a physical
sequential three-flavour C5. The previous one-scale numerical masses are
unaffected.

**Failure handling:** an unexpected operator inventory, failed independent
`CheckLagrangian`, or unsuccessful Matchete match stops execution. The
restricted mode never marks the complete EFT as validated.

</details>

---

<a id="source-t3-general-eft1-scalar-matching-heavy-free-rebase-md"></a>

## Original document: `docs/T3_GENERAL_EFT1_SCALAR_MATCHING_HEAVY_FREE_REBASE.md`

**Original SHA-256:** `a600f5f8fe885869abd803b4096f75195d864ffe72b9fffcb6a5b38b2bbfa8e4`

**Status:** Historical record; check the current source and main
documentation before executing commands reproduced below.

<details>
<summary>Expand original: EFT1 scalar threshold: heavy-field free-Lagrangian rebase and Match audit</summary>

# EFT1 scalar threshold: heavy-field free-Lagrangian rebase and Match audit

## Scope

Experimental **one-loop LLHH-sector probe** for `B_-1`, `C_-1`. This does not repair the complete scalar potential, evolve EFT1, subtract hard UV poles, compute physical sequential C5, or reproduce the oscillation fit.

## Observed failure before this update

Both restricted-sector (`--c5-relevant-only`) inputs passed `CheckLagrangian=True`, with 20/24 saved LLSS terms, 52 Higgs-containing quartics, and 21 background terms. Matchete's subsequent `Match` emitted diagnostic messages that `NewScalarX` and `NewScalarY` were defined **Heavy** but were **Light** in the supplied Lagrangian. Although `RunT3Matching` returned `Status=Success` through its processing stages, the final tree/loop difference did not pass the explicit `hbar`/heavy-field guards and correctly failed with `ScalarThresholdLoopCorrectionInvalid`.

Matchete stage success means the stage returned a non-`$Failed` value; it does **not** by itself establish that heavy fields were integrated out.

## New handling

1. Before changing the field registry, capture the actual `FreeLag[NewScalarX,NewScalarY]` under `Mass->{Light,MSX/MSY}`.
2. After re-registering both fields with `Mass->{Heavy,MSX/MSY}`, capture their new `FreeLag` under the heavy registry.
3. In the opt-in `--c5-relevant-only` matching subtheory **only**, substitute the heavy free sector for the light one: `Expand[scalarEFT-lightEigenFree+heavyEigenFree] // RelabelIndices`. Preserve all matched LLSS and Higgs-bearing interaction terms. Do **not** alter the saved EFT1 source file.
4. Require full `CheckLagrangian=True` for this restricted subtheory. Reject any genuinely invalid Lagrangian.
5. Export a JSON audit of the Matchete `RawEFT`, `EOMEFT` and `MatchedEFT` at **both** loop orders before validating the correction. Record expression head, character count, `hbar`, scalar-field presence, `GammaCC` and a short prefix.
6. Continue only if the difference is nonzero with explicit `hbar` and the final Matchete EFT contains no `NewScalarX`, `NewScalarY` or `NewFermion`. **There is no fallback that invents a contribution if Matchete matches zero.**

This is a power-counting-aware mass/free-Lagrangian handoff attempt. It is **not an algebraic-equivalence certification** between light and heavy field power-counting conventions. The live Matchete result is still unknown until the user runs the scripts.

## Paths (relative to repository root)

**Full replacements**
- `Lagrangian/RunGeneralT3EFT1ScalarThresholdProbe.wl`
- `scripts/RunGeneralT3EFT1ScalarThresholdProbe.py`

**New files**
- `tests/python/test_general_t3_eft1_heavy_free_rebase.py`
- `docs/T3_GENERAL_EFT1_SCALAR_MATCHING_HEAVY_FREE_REBASE.md`

## PowerShell

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "C:\path\to\T3_General_EFT1_Scalar_Heavy_Free_Rebase.zip" -DestinationPath . -Force
python -m pytest -q tests\python\test_general_t3_eft1_heavy_free_rebase.py tests\python\test_general_t3_eft1_c5_relevant_sector.py tests\python\test_general_t3_eft1_scalar_threshold_probe.py
python -m scripts.RunGeneralT3EFT1ScalarThresholdProbe --model B_-1 --c5-relevant-only
python -m scripts.RunGeneralT3EFT1ScalarThresholdProbe --model C_-1 --c5-relevant-only
```

Outputs remain under `output/general_t3_mass_eigen_matchete/eft1_scalar_threshold_probe/`, with independent `_c5_relevant` stems. New diagnostics per model:

- `<stem>_free_lagrangian_rebase.json`
- `<stem>_matching_stage_audit.json`

If `ScalarThresholdLoopCorrectionInvalid` remains, upload the **two matching-stage audit JSONs** and the **two free-Lagrangian rebase JSONs**. Those identify whether the issue is the free mass power counting or Matchete's handling of dimension-five LLSS insertions. Do not rerun earlier fermion-threshold one-loop matching.

</details>

---

<a id="source-t3-general-eft1-verified-reconstruction-and-quartic-isolation-md"></a>

## Original document: `docs/T3_GENERAL_EFT1_VERIFIED_RECONSTRUCTION_AND_QUARTIC_ISOLATION.md`

**Original SHA-256:** `7189a026871f94dfa16acd51ff6c19f9c9ba2e6b1d0df2f8e025307c56e27bde`

**Status:** Historical record; check the current source and main
documentation before executing commands reproduced below.

<details>
<summary>Expand original: T3 EFT1 scalar-threshold: verified renormalizable reconstruction and quartic isolation</summary>

# T3 EFT1 scalar-threshold: verified renormalizable reconstruction and quartic isolation

## Status

This is a **fail-closed experimental repair** of `SavedEFT1TreeHeavyScalarCheckFailed`, not a demonstrated physical threshold-matching solution. It has not run against local Wolfram/Matchete in this environment. Do not infer successful two-threshold C5 until the user's Wolfram runs succeed and the subsequent renormalization/RGE checks are complete.

The local B/C logs show `Hermiticity=False` in the saved EOMEFT even when X/Y are still light, while canonical normalization, heavy-mass basis, charge/gauge checks and field inventory pass. The earlier sector audit implicates combinations containing scalar self-quartics but does not isolate an individual coupling. Neither changing the scalar classification nor `GreensSimplify` / `HcSimplify` cures it.

## Verified reconstruction attempt

`Lagrangian/RunGeneralT3EFT1ScalarThresholdProbe.wl` rebuilds the *renormalizable* SM + X,Y interactions from the same native SU(2) mass-eigenstate components as the first-stage UV Lagrangian. It keeps the **exact** saved 20 (`B_-1`) or 24 (`C_-1`) tree-matched Majorana LLSS terms, and removes the F Yukawas/propagating field. It constructs the candidate:

`L_rebuilt = L_SM + FreeLag[X,Y] + L_rotated_scalar_interactions + L_LLSS_saved`.

Crucially, **this candidate is not used unless** `(i)` `Expand[L_EOMEFT - L_rebuilt]` vanishes exactly, or strict symbolic simplification after index relabelling establishes equality to zero; `(ii)` *ordinary* `CheckLagrangian` returns `True`; `(iii)` the full detailed Matchete `Hermiticity`, `CanonicallyNormalized`, `HeavyMassBasis`, `GaugeInvariant` and `AllObjectsDefined` checks all return `True`. A further detailed check is performed after acceptance. No `PlusHc` is added, no terms are removed or made zero for matching, and no projected `(L+HC[L])/2` is accepted.

If the candidate is algebraically different from the saved EFT—even because of a legitimate yet unproven EOM identity—the runner **refuses** it. This may make the repair conservative and unsuccessful; a no-go result is diagnostic, not evidence the EFT is physically non-Hermitian.

The exact report is `<stem>_reconstruction_audit.json`; if accepted the success certificate and `<stem>_hermiticity_diagnostics.json` record `VerifiedRenormalizableReconstruction`. The Python runner verifies the proof flags before accepting the success certificate.

## Quartic coupling isolation (diagnostic only)

On failure, the runner partitions all terms exactly as before and additionally groups the 122 scalar self-quartic terms by Matchete `Coupling[lambda..., ...]` label. For each label it checks `background + coupling` and `background + all self-quartics except coupling` and writes detailed flags under `scalar_self_quartic_coupling_isolation` in `<stem>_hermiticity_sector_audit.json`. These temporary diagnostic deletions **never** feed `RunT3Matching` and cannot change a physical coefficient.

## Installation and execution (PowerShell)

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath 'C:\path\to\T3_General_EFT1_Verified_Reconstruction_Fix.zip' -DestinationPath . -Force
python -m pytest -q tests\python\test_general_t3_eft1_scalar_threshold_probe.py
python -m scripts.RunGeneralT3EFT1ScalarThresholdProbe --model B_-1
python -m scripts.RunGeneralT3EFT1ScalarThresholdProbe --model C_-1
```

Rerunning fermion-first tree/one-loop, prior hard-pole audit and earlier numerical results is unnecessary.

## If Matchete still fails

Inspect the two generated files:

- `output/general_t3_mass_eigen_matchete/eft1_scalar_threshold_probe/B_m1_final_eft1_scalar_threshold_probe_reconstruction_audit.json`
- `output/general_t3_mass_eigen_matchete/eft1_scalar_threshold_probe/B_m1_final_eft1_scalar_threshold_probe_hermiticity_sector_audit.json`

and corresponding `C_m1` files. In particular look for `algebraic_equality_proven`, `independent_matchete_check_passed` and `scalar_self_quartic_coupling_isolation`. These distinguish a true operator mismatch from a Matchete representation issue and identify the coupling to correct in the next source-level change.

## Scope

Even if the second scalar-threshold matching runs, its exported LLHH sector remains unrenormalized and cannot yet be combined with the fermion threshold without EFT1 running, matching-scale/pole checks and proper flavour normalization. The SM-only one-scale neutrino benchmark remains the currently verified numerical endpoint.

</details>

---

<a id="source-t3-general-fermion-first-eft1-one-loop-md"></a>

## Original document: `docs/T3_GENERAL_FERMION_FIRST_EFT1_ONE_LOOP.md`

**Original SHA-256:** `617b31f27b2ba1f4834d75faa6a81fa7a1633d99228b1057e943e3c24c435f9b`

**Status:** Historical record; check the current source and main
documentation before executing commands reproduced below.

<details>
<summary>Expand original: General T3 — first fermion-threshold one-loop matching probe</summary>

# General T3 — first fermion-threshold one-loop matching probe

Scope: `B_-1` and `C_-1`, native scalar eigenstates `NewScalarX`, `NewScalarY`.

## Prerequisite live results

The previously supplied **tree** first-stage runner passed locally:

| Model | Matchete tree-level EFT1 `GammaCC`/`LLSS` term count |
| --- | ---: |
| `B_-1` | 20 |
| `C_-1` | 24 |

The existing `RunGeneralT3MassEigenMatchete` detailed HeavyMassBasis check and `RunGeneralT3FermionFirstEFT1Tree` report remain prerequisites. They are never overwritten.

## What this new probe actually does

The new `RunGeneralT3FermionFirstEFT1OneLoop.wl` is a **standalone full replacement-style script** based on the live-validated eigenfield and tree-first matching source. It:

1. Rebuilds the UV Matchete model with `NewFermion` **Heavy**, and both `NewScalarX`, `NewScalarY` **Light**. Both scalar eigenstates remain in EFT1. It reuses native SU(2) contraction reconstructions for both Yukawas and all 17 scalar quartics.
2. Validates the fermion-first UV Lagrangian with ordinary and detailed `CheckLagrangian`.
3. Calls the repository's `RunT3Matching[lag, 5, 0]` and validates the tree EFT1 `LLSS` output.
4. Calls `RunT3Matching[lag, 5, 1]` **on the identical UV tree-level Lagrangian**; this uses the existing Matchete stages: `Match`, `GreensSimplify`, `EOMSimplify`, `EvaluateLoopFunctions`, `ReplaceEffectiveCouplings`.
5. Computes exactly `Expand[oneLoopMatchedEFT - matchedEFT]`. The difference is the first-threshold one-loop EFT correction, avoiding inadvertent double-counting of the tree-level EFT. Checks that the correction is nonzero, explicitly includes `hbar`, and contains no `NewFermion` field.
6. Exports the tree EFT, tree `LLSS` sector, full matched EFT through one loop, subtracted loop correction, and the **raw, unrenormalized** candidate `GammaCC`/`HHLL` sector. Candidate selection requires two explicit Higgs legs, `GammaCC`, and no remaining `X`, `Y` or `F` fields. This is a conservative textual diagnostic, **not** a final operator projection.

The candidate count may legitimately be zero. That is **not** silently promoted into a physical neutrino-mass coefficient: zero is saved and reported, and interpretation requires examining the full EFT expression.

## What this does not do

The output **does not** subtract UV poles, enforce an MSbar scheme/scale, certify the hard-pole versus EFT1 anomalous-dimension consistency, establish the complete flavour/CG normalization, run EFT1 between the heavy fermion and scalar thresholds, perform scalar-threshold matching, or calculate a physical `C5`.

The correct future hierarchical computation requires, schematically,

\[
C_5(M_S) = C_{5,\mathrm{hard}}^{\overline{\mathrm{MS}}}(M_F)
 + \Delta C_{5,\mathrm{direct}}^{\mathrm{EFT1}}(M_F\to M_S)
 + \text{scalar-threshold terms}
\]

with fixed-order loop counting and the repository's Majorana symmetry convention. The expression above is only schematic; no unjustified equality between the all-heavy one-scale result and the sequential EFT result is assumed.

## Install

Extract `T3_General_Fermion_First_EFT1_OneLoop_Probe.zip` **into the repository root**, retaining the `Lagrangian/`, `scripts/`, `tests/python/`, `docs/` directories. There are no modifications to the previously successful validators.

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "C:\path\to\T3_General_Fermion_First_EFT1_OneLoop_Probe.zip" -DestinationPath . -Force
python -m pytest -q tests\python\test_general_t3_fermion_first_eft1_one_loop.py
python -m scripts.RunGeneralT3FermionFirstEFT1OneLoop --model B_-1
python -m scripts.RunGeneralT3FermionFirstEFT1OneLoop --model C_-1
```

The run uses local tree-stage and HeavyMassBasis certificates from previous successful runs. Use `--all` only after `B_-1` has completed. The new one-loop run may take longer than tree matching; `--timeout-sec` controls the external process timeout (default 7200 seconds), without modifying an existing session.

## Outputs

`output/general_t3_mass_eigen_matchete/eft1_fermion_first_oneloop/`:

- `<stem>.json`: status, stages, source and matching flags, operator counts.
- `<stem>.progress.json` and `<stem>.log`: last checkpoint and Wolfram stdout/stderr even on failure.
- `<stem>_eft1_tree.inputform.txt`: matched tree-level EFT1.
- `<stem>_llss_tree.inputform.txt`: tree-level `GammaCC` scalar sector.
- `<stem>_eft1_tree_plus_loop.inputform.txt`: matched EFT1 through one loop.
- `<stem>_eft1_one_loop_correction.inputform.txt`: exact difference between the two matched EFTs.
- `<stem>_raw_hard_llhh_sector.inputform.txt`: raw `GammaCC`/two-Higgs candidate sector; may be `0`.
- `summary_final.json`: batch summary.

## Verification status

**Offline tested:** Python runner validation, dependency checking, output schema, nonzero/missing loop factor, zero/nonzero candidate sector, status flags, mocked subprocess; combined test collection 203/203 passed.

**Not yet verified:** Actual local one-loop Matchete matching for the F-first `B_-1` and `C_-1` models. The remote Wolfram kernel cannot run the front-end-dependent `CheckLagrangian`/`Match` stages of this project's Matchete workflow.

The upstream GitHub repository is read-only for this work. Install all changes locally; do not commit or push from this batch.

</details>

---

<a id="source-t3-general-fermion-first-eft1-tree-boundary-md"></a>

## Original document: `docs/T3_GENERAL_FERMION_FIRST_EFT1_TREE_BOUNDARY.md`

**Original SHA-256:** `feacf9bfdf174ff751500e0f6d7e88d6b9f2fd33f98816afce3728adc82a1e66`

**Status:** Historical record; check the current source and main
documentation before executing commands reproduced below.

<details>
<summary>Expand original: T3 B_-1 / C_-1: fermion-first EFT1 tree boundary</summary>

# T3 B_-1 / C_-1: fermion-first EFT1 tree boundary

## Scope / verified facts

- Prior **live** user logs: both models passed the native Matchete heavy mass-eigenfield check and exploratory all-heavy one-loop coefficient extraction. Frozen-reference scalar-aligned one-scale masses and SM+Weinberg endpoint running also passed on the user's Windows machine.
- This **new** experiment has **not** yet been run in Matchete. Offline Python tests are not a physical EFT matching certificate.
- Source checked: `Lagrangian/model/T3Fields.wl` supports `Mass -> {Heavy,mass}` / `{Light,mass}`; `Lagrangian/RunMatching.wl` exposes `RunT3Matching[lag,eftOrder,loopOrder]` and implements `Match` plus Greens/EOM/loop/effective-coupling stages.

## Physical threshold intended

The first boundary in the documented fermion-first hierarchy is

    UV(F, X, Y) -- integrate out F at mu_F --> EFT1(X, Y, LLSS)

where `NewFermion` is registered with its original heavy mass `MF`, and both scalar mass eigenstates `NewScalarX` and `NewScalarY` are registered as **Light**, with masses `MSX, MSY`. Their mixed-doublet mass matrix was already diagonalised upstream; the SU(2)-invariant source Yukawa and scalar vertices are reconstructed using native Matchete CG tensors.

The new Wolfram script builds the same 11 ordinary source operators plus 7 additional quartic source operators and verifies `CheckLagrangian`, then runs:

```wolfram
RunT3Matching[lag, 5, 0]
```

The intended tree-level EFT1 Wilson sector has two lepton fields and two **surviving** scalars, with schematic structure `LL Xdag Xdag`, `LL Xdag Ydag`, `LL Ydag Ydag` (and Hermitian conjugates), modulo Matchete's exact gamma/field conventions. No scalar is integrated out at this stage.

## Implemented fail-closed contract

- Depends on the **prior live-certified** native eigenfield result file for the same model, same numerical mass diagonalisation and same benchmark parameter export.
- Requires successful Matchete `CheckLagrangian`, including `DetailedOutput["HeavyMassBasis"]`, with only `F` declared Heavy in the new runner.
- Requires `RunT3Matching[lag,5,0]` to return `Status==Success` and a `MatchedEFT` expression.
- Exports the full matched EFT and selects a diagnostic `GammaCC` + `Field[NewScalarX|NewScalarY]` sector, following the existing `RunMatching.wl` Wilson-seed text convention. Requires at least one term and evidence of both eigenfields, with `NewFermion` field absent.
- Retains the full EFT export even if the diagnostic fails. JSON, progress file and Wolfram logs identify the failure stage; no fabricated Wilson coefficient or matching certificate is written.
- No files from earlier successful pipeline batches are modified or overwritten.

## Files in this package

- `Lagrangian/RunGeneralT3FermionFirstEFT1Tree.wl` (new full Wolfram source)
- `scripts/RunGeneralT3FermionFirstEFT1Tree.py` (new full CLI)
- `tests/python/test_general_t3_fermion_first_eft1_tree.py` (new offline tests)
- `docs/T3_GENERAL_FERMION_FIRST_EFT1_TREE_BOUNDARY.md` (this file)

## PowerShell

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "C:\path\to\T3_General_Fermion_First_EFT1_Tree_Probe.zip" -DestinationPath . -Force
python -m pytest -q tests\python\test_general_t3_fermion_first_eft1_tree.py
python -m scripts.RunGeneralT3FermionFirstEFT1Tree --model B_-1
python -m scripts.RunGeneralT3FermionFirstEFT1Tree --model C_-1
```

Saved output:

`output/general_t3_mass_eigen_matchete/eft1_fermion_first_tree/`

- `B_m1_final_eft1_fermion_first_tree.json`, `*.progress.json`, `*.log`
- `B_m1_final_eft1_fermion_first_tree_eft1_tree.inputform.txt`
- `B_m1_final_eft1_fermion_first_tree_llss_tree.inputform.txt`
- C analogues; `summary_final.json` after both succeed.

## What remains after this stage

1. Verify first-threshold one-loop hard effects and extract the full EFT1 Wilson set under fixed one-loop power counting; tree `LLSS` alone is **not** final `C5`.
2. Map the native Matchete quartics/Wilson tensor components to the RGBeta/EFT1 conventions with verified normalization; existing source counts do **not** certify physical component maps.
3. Evolve EFT1 from `mu_F` to the scalar thresholds; implement correct mass-eigenstate scalar threshold matching (including distinct mass thresholds if needed).
4. Construct symmetric final `C5`, run SM+Weinberg EFT and compare low-energy masses consistently.
5. Only then extend the matched general-model calculation beyond the two nontrivial mixed-scalar cases.

The previously reported ~0.041 / 0.051 eV low-scale numbers are **one-scale SM-only-background** endpoint results and are not replaced or reclassified by this tree-stage probe.

</details>

---

<a id="source-t3-general-t3-direct-c5-numerical-benchmark-md"></a>

## Original document: `docs/T3_GENERAL_T3_DIRECT_C5_NUMERICAL_BENCHMARK.md`

**Original SHA-256:** `08e02f73c13cb643b46a426cea7a5be051f1c268bcf3e69e0d6ae61db2ead6e3`

**Status:** Historical record; check the current source and main
documentation before executing commands reproduced below.

<details>
<summary>Expand original: Direct one-loop numerical C5 — experimental mixed-doublet models</summary>

# Direct one-loop numerical C5 — experimental mixed-doublet models

## What is computed

The evaluator reads the **actual** `B_m1` or `C_m1` 245–249 kB Matchete
`*_weinberg_ordered_coefficient.inputform.txt` generated by
`RunGeneralT3MassEigenOneLoop`. It evaluates *all* terms at arbitrary numeric
Yukawa matrices, all three quartics actually present in these two exported
coefficients, and the model's independently verified eigenvalues and `U`.
It does **not** insert a separately postulated neutrino-mass loop formula.

The source holds `hbar` as a loop-counting variable. It is evaluated **once**
at `1/(16*pi^2)` (Matchete's documented convention). The project
`Lagrangian/RunMatching.wl` defines the symmetric flavour coefficient by
`C5[i,j] = A5[i,j] + A5[j,i]` for the extracted ordered holomorphic
coefficient. The numerical mass matrix is `mnu = -(v**2 / 2) * C5`.
`numpy.linalg.svd(mnu)` returns the three nonnegative Takagi mass values
(up to unitary phase conventions), converted from GeV to eV.

The expression contains one symbolic `NFlavor` dummy label. The evaluator
explicitly sums it over the **columns** of the supplied Yukawa matrices.
If you supply multiple fermion flavours, this exploratory calculation assumes
**one common fermion mass `MF` for every column**; it cannot represent
nondegenerate fermion masses.

## Exact supplied source structure

- 118 holomorphic Weinberg terms in each model, serialized as ~245 kB (`B`)
  and ~249 kB (`C`) Mathematica InputForm files.
- Native parameter set: `MF`, `MSX`, `MSY`, `T3MEU11..22`, `y1`, `y2`,
  `lambdaT3`, `lambdaT3ExtraQ020020I1`,
  `lambdaT3ExtraQ020200I1`, and `hbar`.
- The two extra quartics are *not* separately normalized against RGBeta.
- The prepared diagonal scalar masses and U matrices are reused without edits.
- Source files are SHA256 hashed into the numeric outputs.

## Explicit illustrative benchmark

`configs/general_t3/c5_illustrative.json` chooses:

- `MF = 100000 GeV`, `v = 246 GeV`,
- `lambdaT3 = 0.1`, `lambdaT3ExtraQ020020I1 = 0`,
  `lambdaT3ExtraQ020200I1 = 0`,
- one active fermion flavour,
  `y1 = (0.005, 0.003, 0.004)^T`,
  `y2 = (0.004, 0.005, 0.003)^T`;
- scalar masses/U from the existing per-model parameter files.

Results using the **uploaded** source expressions:

| Model | Largest Takagi mass (eV) | Second (eV) | Third (eV) |
|---|---:|---:|---:|
| B_-1 | 148.4229726 | 4.6052123 | numerical zero |
| C_-1 | 86.0209794 | 2.6519176 | numerical zero |

These are *not* physical EFT1 benchmark predictions, oscillation fits, or
normalization-certified general T3 masses. Different user-supplied Yukawa
matrices and quartics will give different numbers; do not interpret the
values above as the model's viability.

## Run from repository root (PowerShell)

```powershell
python -m pytest -q tests\python\test_general_t3_numerical_c5.py
python -m scripts.RunGeneralT3NumericalC5 --all
python -m scripts.RunGeneralT3NumericalC5 --model B_-1 --config configs\general_t3\c5_illustrative.json
```

The required matched coefficient files must be located in:
`output/general_t3_mass_eigen_matchete/one_loop_probe/`;
model numeric mass parameters in:
`output/general_t3_mass_eigen_matchete/`.
These are already present on the system that generated the user's ZIP.

Output is written to:
`output/general_t3_numerical_c5_exploratory/`
with the ordered coefficient, symmetric `C5`, mass matrix, Takagi singular
values, full benchmark inputs, and SHA256 source checks. It also writes
`*_mnu_eV.csv` matrices and `summary_final.csv`.

## Not yet performed

- Physical high-scale fermion-first matching and EFT1 scalar-threshold chain;
- full 16-model numeric general quartic normalization/coupling map;
- dependence on nondegenerate heavy-fermion flavours;
- flavour fits and comparison to oscillation data;
- independent Matchete/Wolfram *numeric* substitution cross-check.

The next numerical batch should connect these direct one-scale expressions to
the actual saved EFT1 couplings and normalize the RGBeta/Matchete operator
basis; without that, this is an exploratory demonstration only.

</details>

---

<a id="source-t3-general-t3-eft1-hard-pole-audit-md"></a>

## Original document: `docs/T3_GENERAL_T3_EFT1_HARD_POLE_AUDIT.md`

**Original SHA-256:** `f726ab7167329a881228a01c18b7e42adee74728f1bd3fdbfa013b3883bbb00c`

**Status:** Historical record; check the current source and main
documentation before executing commands reproduced below.

<details>
<summary>Expand original: Fermion-first EFT1: raw hard LLHH pole/log audit</summary>

# Fermion-first EFT1: raw hard LLHH pole/log audit

This is a **diagnostic**, not a renormalized or physical matched Weinberg
coefficient. It analyses the successful local Matchete exports produced by
`RunGeneralT3FermionFirstEFT1OneLoop` for `B_-1` and `C_-1`. It reads those
files verbatim; it does **not** redo one-loop Matchete matching, evaluate
SU(2) Clebsch–Gordan tensors, alter quartic normalisations or perform EFT1
running.

## Verified directly against the uploaded user exports

| Item | B_-1 | C_-1 |
|---|---:|---:|
| Full one-loop EFT1 correction terms | 9,438 | 9,380 |
| Raw LLHH terms, signed subset of full correction | 204 | 204 |
| Pole terms containing Unicode `ϵ` | 68 | 68 |
| Log terms containing `Log[μbar2/Coupling[MF, {}, 0]^2]` | 68 | 68 |
| Finite, log-free terms | 68 | 68 |
| Independent signed tensor/flavour monomials in each sector | 68 | 68 |
| Holomorphic `Proj[-1]` terms in each sector | 34 | 34 |
| Conjugate `Proj[1]` terms in each sector | 34 | 34 |

The 68 monomials in each sector contain 28 contributions with `lambdaT3`,
20 with `lambdaT3ExtraQ020020I1`, and 20 with
`lambdaT3ExtraQ020200I1`. Their **full tensor/operator coefficients**
are identical in the pole, log and constant pieces after removal of *only*
their explicit regulator or logarithm factor. The computed raw-operator identity is

\[
\Delta\mathcal L_{LLHH}^{(1)} =
\mathcal A_{LLHH}\left(1+\frac1\epsilon+
\log\frac{\bar\mu^2}{M_F^2}\right),
\]

where \(\mathcal A_{LLHH}\) already includes Matchete's `hbar` and the
**uncontracted** operator tensors. This is a termwise identity, not a scalar
number or a physical Wilson coefficient. At the diagnostic matching point
`μbar2 = MF^2`, its **pole-discarded finite piece** is
`A_LLHH`, subject to counterterm and renormalisation consistency.

A careful distinction is essential: the identical pole and log coefficients
show that \(\partial_{\ln\mu}\log(\bar\mu^2/M_F^2)=2\) multiplies the
same tensor residue. They do **not** independently establish the EFT1
`LLSS → C5` beta function, cancellation of \(1/\epsilon\) by the
appropriate counterterm, or the physical Majorana operator normalisation.

## Source and outputs

Source directory:

`output/general_t3_mass_eigen_matchete/eft1_fermion_first_oneloop/`

The following existing source files must be present **for each model**:

- `<stem>.json` (successful Matchete matching certificate)
- `<stem>_eft1_one_loop_correction.inputform.txt` (full 9,438/9,380-term correction)
- `<stem>_raw_hard_llhh_sector.inputform.txt` (204 extracted terms)
- `<stem>_llss_tree.inputform.txt` (tree EFT1 source sector)

Run from the repository root in PowerShell:

```powershell
python -m pytest -q tests\python\test_general_t3_eft1_hard_pole_audit.py
python -m scripts.RunGeneralT3EFT1HardPoleAudit --all
```

Outputs:

`output/general_t3_mass_eigen_matchete/eft1_hard_pole_audit/`

Each model has a provenance JSON with SHA-256 hashes of the four input
artifacts, the counts and *explicit false* physics-certification fields, plus
InputForm operator expressions for (i) the raw pole/log/finite sectors,
(ii) the full-operator pole residue, log coefficient, and finite piece at
\(\bar\mu=M_F\), and (iii) the separate `Proj[-1]` and `Proj[1]` finite
branches. Each text output is an uncontracted LLHH **operator expression**,
not a flavour-space \(C_5\) matrix.

The Python parser splits only at top-level Mathematica `+`/`-` symbols and
checks exact reconstructibility. It does not use regular expressions to
rewrite CG tensors, masses, flavour indices or nested algebra; factor
removal is restricted to the two exact exported regulator/log tokens. It
asserts that every raw sector term is an identical signed term of the full
Matchete one-loop correction.

## Next genuine physics step

Use the **independent EFT1 Wilson tensor and direct Weinberg beta function**
to check whether the hard pole residual is cancelled, with the correct CG
and flavour projection. Subsequently construct the renormalized hard
threshold coefficient and direct `LLSS → C5` contribution and combine at the
scalar threshold. These are not yet performed. In particular the one-scale
neutrino masses cannot be relabelled as sequential EFT1 predictions.

</details>

---

<a id="source-t3-general-t3-eft1-scalar-threshold-probe-md"></a>

## Original document: `docs/T3_GENERAL_T3_EFT1_SCALAR_THRESHOLD_PROBE.md`

**Original SHA-256:** `6e3e6e839516c635f462ba9f0bc20e313c88f56e73ee1233f7dff18b87185001`

**Status:** Historical record; check the current source and main
documentation before executing commands reproduced below.

<details>
<summary>Expand original: Ordinary T3 `B_-1`/`C_-1`: scalar-threshold Matchete probe</summary>

# Ordinary T3 `B_-1`/`C_-1`: scalar-threshold Matchete probe

## Purpose and status

The **already live-passed** first threshold, `F` heavy with eigen-doublets `X,Y` light, produces a *tree-level* EFT1 interaction and (in a separate run) an unrenormalized one-loop hard `LLHH` sector. The hard-pole audit found 68 finite, 68 matching-log and 68 pole terms for each model. Matching those sectors coefficientwise to each other is an internal audit, **not** a hard-pole / EFT1 RGE test.

This batch tests the *second* field-elimination boundary, independently and without treating the first hard contribution as renormalized:

    EFT1(tree: X,Y light, F absent)
       -> fresh Matchete process with F absent and X,Y heavy
       -> Match(EFTOrder=5,LoopOrder=0) and Match(EFTOrder=5,LoopOrder=1)
       -> raw one-loop scalar-threshold LLHH sector

This is **not** the full fermion-first sequential EFT calculation. In particular, it does **not** propagate the one-loop hard contribution, compute EFT1 Wilson tensor running, remove UV divergences through independently checked EFT1 counterterms, split X/Y matching thresholds, or report physical `C5`.

## Actual implementation

- The one-loop first-threshold JSON and matching files are validated using the existing first-threshold certificate code; the prior heavy-mass-basis and tree matching certificates are also required.
- The **tree-only** `*_eft1_tree.inputform.txt` is taken from the successful fermion-first one-loop export (not the `tree_plus_loop` input), preventing a second one-loop counting of inherited hard terms.
- A new `wolframscript` process loads Matchete, the original model and the same native SU(2) CG construction used in the earlier live-passed scripts. It removes `F` from the field registry, preserves its mass as a real coupling, and re-registers `X,Y` as `Mass -> {Heavy,MSX/MSY}`.
- The saved tree EFT1 is parsed in the fresh process and requires ordinary and detailed `CheckLagrangian`, including `HeavyMassBasis=True`.
- Tree and one-loop matching are run on **the same tree EFT1**. Their difference is exported alongside the matched SM tree+loop EFT and its raw `LLHH` sector.
- The raw selection requires two `H` and two `l` field entries and the Matchete `GammaCC` bilinear, and excludes `X`, `Y` and `F`. A zero raw sector is reported as zero, not silently promoted to physical matching.
- Progress stages and logs survive an aborted/failed Wolfram run. Existing Lagrangian sources, saved EFT1 outputs, RGEs and numerical benchmarks are not edited.

## Project paths

New source files, installed at repository root:

- `Lagrangian/RunGeneralT3EFT1ScalarThresholdProbe.wl`
- `scripts/RunGeneralT3EFT1ScalarThresholdProbe.py`
- `tests/python/test_general_t3_eft1_scalar_threshold_probe.py`
- `docs/T3_GENERAL_T3_EFT1_SCALAR_THRESHOLD_PROBE.md`

Output directory:

`output/general_t3_mass_eigen_matchete/eft1_scalar_threshold_probe/`

Per model it contains `*.json`, `*.progress.json`, `*.log`, `*_sm_tree.inputform.txt`, `*_sm_tree_plus_loop.inputform.txt`, `*_scalar_one_loop_correction.inputform.txt`, and `*_raw_scalar_llhh.inputform.txt`.

## PowerShell

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath "C:\path\to\T3_General_EFT1_Scalar_Threshold_Probe.zip" -DestinationPath . -Force
python -m pytest -q tests\python\test_general_t3_eft1_scalar_threshold_probe.py
python -m scripts.RunGeneralT3EFT1ScalarThresholdProbe --model B_-1
python -m scripts.RunGeneralT3EFT1ScalarThresholdProbe --model C_-1
```

This needs the successful previous `RunGeneralT3FermionFirstEFT1Tree` and `RunGeneralT3FermionFirstEFT1OneLoop` results in their original output locations, the same mass-eigenstate input parameters, and your working local Matchete + `wolframscript`. A Matchete failure is expected to be diagnosed using `last_stage` in the `.progress.json` and the final `.log`; do not claim success from offline unit tests.

## Next stage if the probe passes

Build an **independent mass-eigenfield EFT1 tensor map** so that the exact Matchete `LLSS` and transformed scalar quartics populate the intermediate Wilson/renormalisable RGEs in the **same SU(2) basis**; then test whether the direct `LLSS -> LLHH` mixing and the hard-pole residue satisfy the specified MSbar consistency condition. Incorporate its RG evolution and the required counterterms before assembling the physical symmetric `C5` at the scalar threshold and running to 100 GeV.

</details>

---

<a id="source-t3-general-t3-frozen-reference-yukawa-transfer-md"></a>

## Original document: `docs/T3_GENERAL_T3_FROZEN_REFERENCE_YUKAWA_TRANSFER.md`

**Original SHA-256:** `e74c1abc5c3fa6ea65e276e5ad3926938c8e39daf0e77d5dbc3e3282b4b622bd`

**Status:** Historical record; check the current source and main
documentation before executing commands reproduced below.

<details>
<summary>Expand original: General T3 one-scale frozen-reference Yukawa transfer (B_-1/C_-1)</summary>

# General T3 one-scale frozen-reference Yukawa transfer (B_-1/C_-1)

## Purpose and limits

This is **NOT** a physical fitted-reference re-evaluation, a new fit, or
sequential fermion-first EFT1 matching. It is a controlled *hybrid* comparison:

- **Yukawas and fermion masses** are taken from the existing frozen fitted
  reference files in `Numerical/benchmarks/fitted_reference/`;
- **scalar masses and mixing U** remain those from the independently validated
  `output/general_t3_mass_eigen_matchete/{B,C}_m1_final_parameters.json`;
- **one-loop amplitude** is evaluated from the original Matchete InputForm
  coefficient, not a phenomenological substitute loop function;
- **normalisation** uses the project's ordered `A5` -> symmetric `C5`
  convention; full Matchete/RGBeta component map is *not* certified.

The frozen reference scalar parameters differ from the experimental scalar
mass-mixing inputs, even before considering the extra off-diagonal mass term.
Consequently, **the resulting neutrino masses are not the original fitted
NuFIT predictions**, and cannot be used to claim a successful fit or
reproduce the legacy `smallY_smallL` benchmark. The output includes the
original and experimental squared scalar mass values and their differences.

The two additional quartic couplings in the exported coefficients are set to
zero only because their corresponding frozen-reference input couplings are
**explicitly zero**. If either input is nonzero, the runner fails closed:
its Matchete/RGBeta operator-normalisation map is not yet certified.

## Flavour and mass treatment

The frozen `ordinary.y1`, `ordinary.y2` are 3x3 complex matrices with three
fermion-flavour columns. `ordinary.MF` is a *diagonal*, real-positive 3x3
mass matrix. The evaluator computes the scalar expression independently for
each fermion-flavour column `k`, substituting `MF = MF[k,k]` in that term:

`A5[i,j] = sum_k A5_ Matchete [y1[i,k], y2[j,k], MF[k,k], MSX, MSY, U, quartics]`.

This implementation assumes no fermion mass rotation, no off-diagonal heavy
fermion mass, and the same scalar eigenstate benchmark for all three fermion
flavours. Off-diagonal/complex `MF` fails closed rather than being treated as
a scalar number. The original one-column illustrative mode still works.

The `mnu` matrix uses `C5 = A5 + Transpose[A5]` and
`mnu = -v^2 C5 / 2`, with `hbar = 1/(16 pi^2)` once. The three Takagi masses
are the singular values of the complex-symmetric mass matrix, in eV.

`max_abs_mnu_eV` is **the largest absolute matrix entry** `max_ij |mnu_ij|`.
It is *not* the largest physical neutrino mass. For the previous illustrative
inputs, the largest Takagi mass is 148.42297 eV (B_-1) or 86.02098 eV
(C_-1); the largest matrix entries are only 61.20390 and 35.47340 eV.
The updated CLI reports both quantities by their unambiguous names.

## Installation

Extract the ZIP at the repository root. Changed file paths:

- `scripts/RunGeneralT3NumericalC5.py` — FULL replacement; supports a
  nondegenerate `fermion_masses_gev` list as well as the original common MF.
- `tests/python/test_general_t3_numerical_c5.py` — FULL replacement (retains
  existing tests).
- `scripts/RunGeneralT3ReferenceC5.py` — NEW frozen-reference transfer CLI.
- `tests/python/test_general_t3_reference_c5.py` — NEW regression tests.
- `docs/T3_GENERAL_T3_FROZEN_REFERENCE_YUKAWA_TRANSFER.md` — NEW note.

Frozen reference configurations are **read-only**. Nothing is written to their
original files, RGBeta exports, or Matchete model builder.

### Commands (PowerShell)

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
python -m pytest -q tests\python\test_general_t3_reference_c5.py tests\python\test_general_t3_numerical_c5.py
python -m scripts.RunGeneralT3ReferenceC5 --all
```

To reproduce the earlier illustrative run with clarified labels:

```powershell
python -m scripts.RunGeneralT3NumericalC5 --all
```

The transfer reads:

- `Numerical/benchmarks/fitted_reference/benchmark_manifest.json`
- `Numerical/benchmarks/fitted_reference/T3_dS1_2_dS2_2_dF_{1,3}_alpha_m1.json`
- `output/general_t3_mass_eigen_matchete/{B,C}_m1_final_parameters.json`
- `output/general_t3_mass_eigen_matchete/one_loop_probe/{B,C}_m1_final_one_loop_probe.json`
- corresponding `*_weinberg_ordered_coefficient.inputform.txt` files.

Outputs (separate from the illustrative set):

`output/general_t3_numerical_c5_reference_transfer/{B,C}_m1_final_numerical_c5.json`

`output/general_t3_numerical_c5_reference_transfer/{B,C}_m1_final_mnu_eV.csv`

`output/general_t3_numerical_c5_reference_transfer/summary_final.json`

Each model's JSON records the frozen input SHA256 and mass mismatch, with
`physical_c5_certified=false`, `sequential_EFT1_performed=false`, and
`mass_basis_compatible_with_frozen_reference=false`.

## Next physical milestone

Do not conflate this one-scale hybrid transfer with the source-consistent
fermion-first EFT1 result. Physical next steps are (i) choose/derive the
quadratic mixing input and mass basis **at the same running scale** as the
EFT1 benchmark, (ii) match the Matchete and RGBeta component-coupling basis,
(iii) run the appropriate F -> scalar thresholds, and (iv) construct the
3x3 symmetric C5 and masses at the specified comparison scale.

</details>

---

<a id="source-t3-general-t3-frozen-scalar-aligned-one-scale-md"></a>

## Original document: `docs/T3_GENERAL_T3_FROZEN_SCALAR_ALIGNED_ONE_SCALE.md`

**Original SHA-256:** `0d27429466d0c803945b97f4bd485229a1f0f8948f20f995f5f4d56e6cddbd84`

**Status:** Historical record; check the current source and main
documentation before executing commands reproduced below.

<details>
<summary>Expand original: T3 B_-1/C_-1: frozen scalar diagonal-aligned direct C5 benchmark</summary>

# T3 B_-1/C_-1: frozen scalar diagonal-aligned direct C5 benchmark

## Why

The earlier hybrid reference-Yukawa evaluation uses scalar masses from a **different**
experimental quadratic benchmark. The checked frozen fitted-reference configurations
instead specify (at their saved base scale) `mS1Sq` and `mS2Sq`, and do not specify
the new gauge-invariant scalar mixing `mu2`. Mixing **cannot** be inferred from
`lambdaT3` or the pre-existing exploratory mass-basis input.

This batch deliberately evaluates a **scenario**, not a refit:

- frozen 3x3 complex Yukawa matrices and diagonal fermion masses from
  `Numerical/benchmarks/fitted_reference/..._alpha_m1.json`;
- frozen scalar **diagonal mass squared** inputs exactly, with independent
  complex mass mixing `mu2` in GeV^2; default `mu2 = 0` explicitly assumed;
- `lambdaT3` transferred directly; the two quartics used by the exact
  one-loop coefficient are retained at their frozen **zero** values only.
  Nonzero values fail closed, since their Matchete-to-RGBeta normalization
  has not been certified;
- the **original Matchete-produced one-loop ordered coefficient** is evaluated
  with the new *verified* mass eigenvalues and mixing matrix, not replaced
  by a postulated analytic T3 loop formula.

The scalar convention is

```
Phi=(S1, epsilon Conjugate[S2]); epsilon={{0,1},{-1,0}}
A = {{mS1Sq,mu2},{Conjugate[mu2],mS2Sq}}
U^dagger A U = DiagonalMatrix[{mX^2,mY^2}], mX <= mY
Phi = U.{X,Y}
```

For `mu2 = 0` and frozen `mS1Sq > mS2Sq`, the sorted eigenbasis is
`X = epsilon Conjugate[S2]`, `Y = S1` with `U={{0,1},{1,0}}`, not
`U=IdentityMatrix[2]`.

## Installation paths

This ZIP contains **full files**, not diffs:

| Relative project path | Action |
|---|---|
| `scripts/RunGeneralT3NumericalC5.py` | Full replacement; adds strictly validated scalar-basis override |
| `scripts/RunGeneralT3FrozenScalarC5.py` | New scenario runner |
| `tests/python/test_general_t3_frozen_scalar_c5.py` | New integration and fail-closed tests |
| `docs/T3_GENERAL_T3_FROZEN_SCALAR_ALIGNED_ONE_SCALE.md` | New documentation |

It does **not** overwrite the one-loop Matchete coefficient, mass-basis
validator, or frozen fitted reference configurations.

## Usage (PowerShell, repository root)

```powershell
Expand-Archive -LiteralPath "C:\path\to\T3_General_Frozen_Scalar_Aligned_C5.zip" -DestinationPath . -Force
python -m pytest -q tests\python\test_general_t3_frozen_scalar_c5.py tests\python\test_general_t3_reference_c5.py tests\python\test_general_t3_numerical_c5.py
python -m scripts.RunGeneralT3FrozenScalarC5 --all
```

Optional *independent* mixing scenario, replacing the `mu2 = 0` default:

```powershell
python -m scripts.RunGeneralT3FrozenScalarC5 --all --mu2-real-gev2 10000 --mu2-imag-gev2 0
```

Using any specified `mu2` does **not** imply it was fitted or fixed by
existing EFT1 RGEs. Tachyonic/degenerate masses and broken input/hash
certificates are rejected.

Output directory (does not overwrite hybrid comparisons):
`output/general_t3_numerical_c5_frozen_scalar_scenario/`.

The ordered flavour coefficient is evaluated for each diagonal heavy fermion
mass and the symmetric Weinberg coefficient uses `C5 = A5 + Transpose[A5]`.
`mnu=-v^2 C5/2` and its Takagi singular values use the same conventions as
`RunGeneralT3NumericalC5.py`.

## Interpretation and remaining work

This is a **one-scale simultaneous all-heavy** result. It has correct *frozen
scalar diagonals* and flavour Yukawa input but the quadratic mixing is an
**explicit scenario**, not a frozen fit parameter. The saved Yukawa/coupling
values are taken from `base_state`, **not** evolved to each threshold in the
general EFT1 pipeline. It is **not** sequential fermion-first matching,
not a replacement for the frozen fit, and does not certify the physical
Matchete/RGBeta quartic normalization for nonzero extra couplings. In
particular, the result is **not** one-loop-matched EFT1 at the low scale.

After the scenario comparison, the next physics task remains integrating the
matched eigenbasis coefficients into the fermion-first EFT1 matching and
running sequence, with consistent mixing-quadratic running and thresholds.

</details>

---

<a id="source-t3-general-t3-one-loop-status-and-probe-md"></a>

## Original document: `docs/T3_GENERAL_T3_ONE_LOOP_STATUS_AND_PROBE.md`

**Original SHA-256:** `a3b6cc1074be39af654653d15a4f63d9591d695984e93d0640ffbc72e04ed25d`

**Status:** Historical record; check the current source and main
documentation before executing commands reproduced below.

<details>
<summary>Expand original: General T3: status and first direct one-loop Weinberg probe</summary>

# General T3: status and first direct one-loop Weinberg probe

## What is already verified

- All 16 ordinary models: general quartic candidate inventories and **experimental** quadratic+quartic EFT1 numerical runs have completed for the user's benchmarks. These are scalar-sector running results, not final neutrino-mass results.
- 14 ordinary models had diagonal Matchete heavy-basis checks; `B_-1` and `C_-1` now also pass native `HeavyMassBasis` with two mixed scalar doublet eigenfields. The user's local live runs reported the two pairs of scalar masses near (1035.58,1188.53) and (1032.61,1186.00) GeV.
- The B/C eigenbasis has 2 Yukawa source couplings and 17 scalar quartic source slots. The Matchete SU(2) tensor recoupling and numeric scalar-mass diagonalization were separately checked. Physical RGBeta-to-Matchete coupling normalization is not certified.

## What is still missing for general physical results

1. Matching to a physical Weinberg operator after the heavy states are removed, preserving the coefficient convention and flavour structure.
2. Correct numerical identification of Matchete couplings with the evolved RGBeta/EFT1 couplings and the benchmark Yukawa matrices, including relative CG signs and finite terms.
3. Consistent fermion-first and split-scalar matching at the appropriate scales, rather than comparing quantities at mismatched threshold scales.
4. Run the final SM+Weinberg RGEs and obtain symmetric 3x3 C5 matrices, neutrino masses and residual checks. Extend/check the matching and quadratic mass treatment across all 16 models; B/C do not establish this for all models.

## This experimental batch

Two **new** project files introduce a symbolic exploratory diagnostic while leaving `RunGeneralT3MassEigenMatchete.wl`, the validated heavy-basis checker, unchanged.

- `Lagrangian/RunGeneralT3MassEigenOneLoop.wl`: forks the verified eigenfield construction in a fresh Wolfram process, rechecks `HeavyMassBasis`, calls the existing `Lagrangian/RunMatching.wl` entry point `RunT3Matching[lag,5,1]`, exports the matched EFT in raw `InputForm`, and invokes `ExtractWeinbergCoefficient[matchedEFT]`. If a holomorphic dimension-five Weinberg term is found, saves the **ordered**, raw Matchete coefficient in `InputForm`. The built-in optional display canonicalizer assumes `MS1/MS2`, so it is not certified for `MSX/MSY` and is deliberately not used for a numerical C5.
- `scripts/RunGeneralT3MassEigenOneLoop.py`: requires the previous successful local native heavy-basis result, rechecks the numerical mixing input/provenance, executes the symbolic probe, and validates the status and exported artifacts. Failures report the `*.progress.json` stage and full `*.log`.
- `tests/python/test_general_t3_mass_eigen_one_loop.py`: offline interface and fail-closed checks. They do not constitute a live Matchete one-loop matching test.

The probe integrates out all heavy species **in one step at a common symbolic matching point**; it does **not** follow the project's production UV -> F -> EFT1 -> S -> Weinberg threshold sequence, does not choose a numerical matching scale and does not evaluate a numerical 3x3 C5 matrix. It therefore **does not establish physical neutrino masses or finish the general 16-model calculation.**

### Install (PowerShell, repository root)

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
Expand-Archive -LiteralPath 'C:\path\to\T3_General_One_Loop_Matching_Probe.zip' -DestinationPath . -Force
python -m pytest -q tests\python\test_general_t3_mass_eigen_one_loop.py tests\python\test_general_t3_mass_eigen_matchete.py tests\python\test_general_t3_mixed_scalar_mass_basis.py
python -m scripts.RunGeneralT3MassEigenOneLoop --model B_-1
python -m scripts.RunGeneralT3MassEigenOneLoop --model C_-1
```

### Output files

`output/general_t3_mass_eigen_matchete/one_loop_probe/` contains model-specific `*_one_loop_probe.json`, `*_one_loop_probe.progress.json`, `*_one_loop_probe.log`, and—on successful matching—`*_matched_eft.inputform.txt` and `*_weinberg_ordered_coefficient.inputform.txt`.

Next task, once the symbolic coefficient is confirmed by a **live** run: compare with the project matching convention, independently verify the loop and flavour normalization and supply benchmark inputs to evaluate C5 and `m_nu=-v^2 C5/2` at a fixed scale. Then connect to the fermion-first EFT1 evolution for final general results.

</details>

---

<a id="source-t3-general-t3-one-scale-sm-weinberg-endpoint-md"></a>

## Original document: `docs/T3_GENERAL_T3_ONE_SCALE_SM_WEINBERG_ENDPOINT.md`

**Original SHA-256:** `88c96cafe175a56a854405874f279d6e95297da26fa63799f120c005d7bfeec1`

**Status:** Historical record; check the current source and main
documentation before executing commands reproduced below.

<details>
<summary>Expand original: T3 B_-1/C_-1: direct one-scale C5 to low-scale SM+Weinberg endpoint</summary>

# T3 B_-1/C_-1: direct one-scale C5 to low-scale SM+Weinberg endpoint

## Scope

This is **not** the fermion-first EFT1 matching calculation. It continues the
**already computed** one-scale, all-heavy Matchete C5 for B_-1/C_-1 to the
frozen reference `mu_low_gev` using the project's existing
`Numerical.running.SMWeinbergEvolution.evolve_sm_weinberg` solver.

The method is deliberately a *SM-only-background diagnostic*:

1. Read immutable `Numerical/benchmarks/fitted_reference/benchmark_manifest.json`
   and the model's exact frozen reference JSON, enforcing its SHA256 hash.
2. Read the existing scalar-aligned direct-C5 JSON from
   `output/general_t3_numerical_c5_frozen_scalar_scenario/` and refuse
   anything marked physical or sequential EFT1. Its scalar diagonals come
   from the frozen input, with a separately declared quadratic mixing scenario.
3. Evolve the frozen **SM** renormalisable parameters from its stored UV scale
   (benchmark 10^7 GeV) to the stored scalar threshold (benchmark 10^3 GeV)
   with `K=0`, deliberately omitting heavy/BSM effects in this segment.
4. Insert the **unaltered** one-scale Matchete-derived symmetric C5 at the
   scalar threshold and run the existing one-loop SM+Weinberg equations to the
   saved low scale (benchmark 100 GeV).
5. Calculate the Takagi masses as singular values of `m_nu = -(v^2/2) C5`
   with the frozen benchmark's constant `vev_gev`.

The SM-only evolution above the scalar threshold is **not the actual T3
renormalisable running**. No F threshold, EFT1 LLSS, EFT1 direct Weinberg
running, scalar hard threshold or Matchete/RGBeta nonzero-quartic normalization
has been computed here. It does **not** reproduce the original fitted
oscillation result and must not be labelled as such.

## Project files

- `scripts/RunGeneralT3LowScaleOneScale.py` (new)
- `tests/python/test_general_t3_low_scale_one_scale.py` (new)
- `docs/T3_GENERAL_T3_ONE_SCALE_SM_WEINBERG_ENDPOINT.md` (new)

No existing source or previously generated results are overwritten. The runner
imports production `SMWeinbergEvolution.py` from the user's installed repo;
its high-level regression tests inject a deterministic evolution double.
Those tests are **not** independent verification of the external numerical RGE.

## Windows PowerShell

From repository root:

```powershell
Expand-Archive -LiteralPath "C:\path\to\T3_General_OneScale_SM_Weinberg_Endpoint.zip" -DestinationPath . -Force
python -m pytest -q tests\python\test_general_t3_low_scale_one_scale.py
python -m scripts.RunGeneralT3LowScaleOneScale --all
```

Prerequisite:

```powershell
python -m scripts.RunGeneralT3FrozenScalarC5 --all
```

Expected outputs:

- `output/general_t3_low_scale_one_scale/B_m1_final_low_scale_one_scale.json`
- `output/general_t3_low_scale_one_scale/C_m1_final_low_scale_one_scale.json`
- `output/general_t3_low_scale_one_scale/summary_final.json`

The JSON reports scalar-boundary and low-scale Takagi singular masses, change
ratios, C5 matrices, source/reference hashes, and explicit limitations.

## Remaining physical step

The proper fermion-first T3 calculation requires the mass-eigenfield
Lagrangian with F integrated out at `mu_F` while X,Y remain active,
the intermediate LLSS/direct Weinberg running with the appropriate tensor
basis, and the scalar matching at `mu_S` with correctly mapped couplings.
Those results must replace—not be conflated with—this diagnostic.

</details>

---

<a id="source-t3-general-two-threshold-numerical-c5-md"></a>

## Original document: `docs/T3_GENERAL_TWO_THRESHOLD_NUMERICAL_C5.md`

**Original SHA-256:** `1764293291692705ef9b7df1f2cdafd2274785bdc88b1e8ca4658951f3bd4d41`

**Status:** Historical record; check the current source and main
documentation before executing commands reproduced below.

<details>
<summary>Expand original: Restricted two-threshold fixed-order Weinberg projection (B_-1 and C_-1)</summary>

# Restricted two-threshold fixed-order Weinberg projection (B_-1 and C_-1)

## Scope and provenance

This follows the locally verified two-threshold raw operator assembly. The previous stage found 68 hard/scalar UV-pole cancellations per model; that *operator-level* cancellation is used as a prerequisite here. The scalar threshold uses the explicit `--c5-relevant-only` selection and excludes 122 four-heavy-scalar quartic terms. **The full scalar EFT has not passed Hermiticity, and the independent intermediate-EFT anomalous-dimension / RGE calculation remains outstanding.** These outputs are preliminary fixed-order estimates, not a fully renormalized sequential physical C5 prediction.

Input files, per model (`B_m1_final` / `C_m1_final`):

- `output/general_t3_mass_eigen_matchete/two_threshold_operator_assembly/<stem>_two_threshold_raw_llhh_audit.json`
- `output/general_t3_mass_eigen_matchete/two_threshold_operator_assembly/<stem>_two_threshold_raw_llhh_fixed_muF_unprojected.inputform.txt`
- Optional, for `--numeric`: `output/general_t3_numerical_c5_frozen_scalar_scenario/<stem>_numerical_c5.json`, generated by `python -m scripts.RunGeneralT3FrozenScalarC5 --all` using *explicit* zero scalar quadratic mixing unless a nonzero scenario is requested.

The projection fails closed if the prior 68/68 cancellation certificate, scoped exclusion or input operator patterns differ.

## Matchete-to-physical convention

For the 232-term fixed-muF operator, there are exactly 116 `Proj[-1]` terms with two unbarred Higgs fields and 116 conjugate `Proj[1]` terms. The exporter confirms both SU(2) epsilon tensors, the entire noncommutative Majorana chain and both lepton flavour indices for each term. For C_-1, some external-leg patterns are related to the standard orientation by *two* epsilon index exchanges, yielding no extra minus sign. The script strips these operator factors **in the same manner as** `StripWeinbergOperatorStructure` in the existing `Lagrangian/RunMatching.wl`. No arbitrary group or spinor multipliers are introduced.

The remaining ordered scalar coefficient \(A_{5,ij}\) is a full symbolic InputForm expression with `y1`, `y2`, `MF`, scalar masses, \(U_{ij}\) and quartics retained. Its `Proj[1]` counterpart is checked for the expected field/epsilon orientations, not used for an additional physical factor of two.

For numerical evaluation, the reference project's convention is

\[
C_{5,ij} = A_{5,ij}+A_{5,ji},\qquad
m_{\nu,ij}=-\frac{v^2}{2}C_{5,ij},\qquad
\hbar=\frac{1}{16\pi^2}.
\]

The 3x3 matrices sum over the fermion-generation index (three diagonal MF values in the frozen reference); neutrino Takagi masses are singular values of the complex symmetric mass matrix, converted to eV. The input Yukawas were *previously fitted* to neutrino data: output masses are not independent predictions. This does not evolve C5 to 100 GeV; it reports masses evaluated with this coefficient at the chosen common matching scale.

## PowerShell

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline

# Only needed if you do not already have the frozen scalar-aligned outputs:
python -m scripts.RunGeneralT3FrozenScalarC5 --all

# Project the already verified source operator (no Matchete run):
python -m scripts.RunGeneralT3TwoThresholdNumericalC5 --all

# Then numerically evaluate using the already computed, compatible frozen scenario:
python -m scripts.RunGeneralT3TwoThresholdNumericalC5 --all --numeric

# Focused checks:
python -m pytest -q tests\python\test_general_t3_two_threshold_numerical_c5.py
```

Main output directory:

`output/general_t3_mass_eigen_matchete/two_threshold_numerical_projection/`

For each model:

- `<stem>_restricted_fixed_muF_ordered_A5.inputform.txt`: complete 116-term ordered PL coefficient;
- `<stem>_restricted_fixed_muF_projection_audit.json`: input hashes, normalization assumptions, exact operator patterns, provenance;
- `<stem>_restricted_fixed_muF_numerical_c5.json`: optional numerical ordered A5, symmetric 3x3 C5, neutrino mass matrix, Takagi masses and comparison to **one-scale frozen** output.

No files from earlier successful stages are overwritten.

## Interpretation

The assembled pole-free expression is a fixed-order result evaluated at a common \(\bar\mu=M_F\). It retains the heavy-to-light logs \(\log(M_F^2/M_{X,Y}^2)\), giving a first sequential-scale numerical estimate. It is **not** an independent EFT1 resummation or a proof of MSbar renormalization. True two-stage EFT computation still needs the X/Y-basis Wilson-tensor normalization, EFT1 RGEs, separate threshold scales, lower-threshold matching/counterterms, and full scalar-potential consistency. The combined fixed-order mass may be compared with the simultaneous one-scale mass for a useful leading \(M_{X,Y}^2/M_F^2\) check, not interpreted as a new fitted prediction.

## Offline benchmark preview (independently assembled inputs)

Using the frozen-reference `B_-1` / `C_-1` Yukawa matrices, diagonal fermion masses 100/110/120 TeV, scalar-squared diagonals 1,306,122.448979592 and 1,000,000 GeV², zero off-diagonal mass, the corresponding unitary swap eigenbasis, and `lambdaT3=0.1`, with the two additional model quartics set to zero, the script evaluated:

| Model | Leading fixed-order largest Takagi mass | Prior aligned one-scale largest Takagi mass | Relative difference |
|---|---:|---:|---:|
| `B_-1` | 0.04359075491599 eV | 0.043679827 eV | −0.203920% |
| `C_-1` | 0.05383029966959 eV | 0.053936105 eV | −0.196168% |

These are local **reconstructed-input preview** values, not outputs from the user's Windows run. The user's own `--numeric` run reads its existing frozen-aligned one-scale output and writes provenance hashes and complete matrices. Do not treat this agreement as certification of EFT1 RG running, the complete scalar potential or independent neutrino oscillation phenomenology.

</details>

---

<a id="source-t3-general-two-threshold-raw-pole-assembly-md"></a>

## Original document: `docs/T3_GENERAL_TWO_THRESHOLD_RAW_POLE_ASSEMBLY.md`

**Original SHA-256:** `56ef48d9a96f1dce46ff0deaf9e371ff52c3bffad7eaac522e00d782f0c867af`

**Status:** Historical record; check the current source and main
documentation before executing commands reproduced below.

<details>
<summary>Expand original: T3 B_-1/C_-1: two-threshold raw Weinberg-sector pole cancellation</summary>

# T3 B_-1/C_-1: two-threshold raw Weinberg-sector pole cancellation

## Scope and verification

This batch consumes **existing**, locally generated Matchete InputForm outputs from

- `output/general_t3_mass_eigen_matchete/eft1_fermion_first_oneloop/`
- `output/general_t3_mass_eigen_matchete/eft1_hard_pole_audit/`
- `output/general_t3_mass_eigen_matchete/eft1_scalar_threshold_probe/`

It does not invoke Wolfram, rerun any matching, update GitHub, or modify source expressions. The scalar-threshold inputs are explicitly certified as `LeadingC5OneLoopSectorOnly`, with 122 four-heavy-scalar quartic terms **omitted**; the complete scalar potential still fails a separate Hermiticity check. This batch does **not** establish complete matching or the physical neutrino mass matrix.

## Exact observed term structure

| InputForm sector | B_-1 | C_-1 |
|---|---:|---:|
| Fermion threshold: UV poles | 68 | 68 |
| Fermion threshold: matching logs | 68 | 68 |
| Fermion threshold: nonlog finite terms | 68 | 68 |
| Scalar threshold: UV poles | 68 | 68 |
| Scalar threshold: scalar logs, MSX | 50 | 50 |
| Scalar threshold: scalar logs, MSY | 50 | 50 |
| Scalar threshold: nonlog finite terms | 64 | 64 |
| **Opposite UV pole-residue pairs, matched** | **68/68** | **68/68** |

On both actual datasets, the complete signed operator residue of the scalar-threshold `1/ϵ` terms is the negative of the fermion-threshold residue:

\[
\operatorname{Res}_{1/\epsilon}\!\left[\Delta\mathcal L_{LLHH}^{F,\;1\ell} + \Delta\mathcal L_{LLHH}^{S,\;1\ell}\right] = 0.
\]

Comparison retains *all* `Field`, `NCM`, `CG`, Yukawa-flavour indices, mixing spurions, and `MF` factors. It does **not** delete the operator structures to make a spurious coefficient equality.

`B_-1`: all 68 opposite signed operator residues match without changing dummy-index labels.

`C_-1`: 34 PR/conjugate terms match by SU(2) epsilon antisymmetry, while the 34 PL/holomorphic terms additionally require the change of *summed SU(2) index labels* `1 <-> 2, 3 <-> 4`. Every `Index[d$$N, SU2L[fund]]` appears exactly twice in the affected terms, so this is an allowed closed-dummy relabelling, not a change to a physical index. The two epsilon matrices are independently antisymmetric. A single epsilon argument swap therefore contributes a minus sign; two swaps cancel. The audit rejects deviations from these exact forms.

## Outputs

Run the script from the project root:

```powershell
python -m pytest -q tests\python\test_general_t3_two_threshold_pole_assembly.py
python -m scripts.RunGeneralT3TwoThresholdPoleAssembly --all
```

All output files are under

`output/general_t3_mass_eigen_matchete/two_threshold_operator_assembly/`.

For each `B_m1_final` / `C_m1_final`, the runner writes seven complete Mathematica InputForm operator expressions and a JSON report.

- `*_hard_pole_residue_canonical.inputform.txt`: first-stage `1/ϵ` residue, with `ϵ` removed but otherwise retaining the exact LLHH operator.
- `*_scalar_negative_pole_residue_canonical.inputform.txt`: negative second-stage pole residue. **Bytewise identical** to the corresponding canonical first-stage file.
- `*_hard_finite.inputform.txt`: raw first-stage terms without pole or log.
- `*_soft_finite.inputform.txt`: raw scalar-threshold terms without pole or log.
- `*_soft_logs.inputform.txt`: scalar-threshold mass-dependent log terms.
- `*_fixed_common_mu_unprojected.inputform.txt`: 300 terms from both matching stages after exact cancellation of all `1/ϵ` terms; `μbar2` remains symbolic.
- `*_fixed_muF_unprojected.inputform.txt`: 232 terms in the formal common-scale convention `μbar2 = MF^2`: the 68 first-stage hard logs vanish; second-stage logarithms become `Log[MF^2/MSX^2]` or `Log[MF^2/MSY^2]`.

**These are entire one-loop operators, not entries of the Weinberg Wilson matrix**. `hbar` remains symbolic. For a numerical Mathematica substitution, Matchete's loop bookkeeping convention in this pipeline is `hbar -> 1/(16 Pi^2)`.

## Numerical scale context (illustrative only)

The existing frozen-reference scenario uses `MSX = 1000 GeV`, `MSY = 1142.857142857143 GeV`, and reference heavy fermion masses `100`, `110`, and `120 TeV` (different fermion flavours). Therefore

| MF (GeV) | `Log[MF^2/MSX^2]` | `Log[MF^2/MSY^2]` |
|---:|---:|---:|
| 100,000 | 9.210340 | 8.943278 |
| 110,000 | 9.400961 | 9.133898 |
| 120,000 | 9.574983 | 9.307921 |

These are potentially significant logarithms, indicating why intermediate-EFT leading-log running should ultimately be included. This table is *not* a calculation of the neutrino masses. The symbolic `MF` in the files must be interpreted per fermion mass eigenstate when constructing three-generation coefficients.

## What is and is not established

**Verified from both uploaded raw expressions:**

1. Pole cancellation of 68/68 full signed LLHH operators per model, with the stated explicit SU(2) identities.
2. Every raw hard/scalar LLHH term belongs to its corresponding saved one-loop correction (accounting for a leading Mathematica unary-minus representation).
3. Both matching-stage metadata certificates and the previously generated fermion-stage pole-audit SHA256 hashes agree.
4. Outputs contain no `ϵ` UV-regulator symbol, with 300 formal common-scale or 232 formal `μbar = MF` terms.

**Not yet established:**

- An *independent* `LLSS -> C5` EFT1 anomalous dimension or full RG consistency (pole cancellation across two computations is distinct from an RGE beta-function calculation).
- The complete renormalized EFT1, as the full four-heavy-scalar quartic potential was excluded at the scalar threshold.
- SU(2), spinor, flavour and Majorana normalization of the extracted operator relative to a physical Weinberg coefficient `C5[p,q]`.
- Running from each fermion mass eigenvalue through the scalar threshold and below it, or comparison against the independently obtained direct one-scale mass eigenvalues.
- A matrix-level physical sequential neutrino mass prediction.

**Next targeted physics step:** project the pole-free full `LLHH` operator into the single consistent gauge-invariant Weinberg tensor, including chirality orientation and symmetry factors; then implement the EFT1 `LLSS` mixing and run between the two thresholds. Only after those checks is numerical matrix assembly defensible.

## Implementation provenance

New source: `scripts/RunGeneralT3TwoThresholdPoleAssembly.py`.

New tests: `tests/python/test_general_t3_two_threshold_pole_assembly.py`.

All prior files and validated matching executables remain unchanged. Do not overwrite the reference matching exports with the composite outputs. Do not present the common-scale expression as an already scheme-renormalized `C5`.

</details>

---

<a id="source-general-t3-eft1-scalar-pair-mixing-md"></a>

## Original document: `docs/GENERAL_T3_EFT1_SCALAR_PAIR_MIXING.md`

**Original SHA-256:** `e5f030bff424796a0f3bcc492f8f8fffe968660f88784fb0cd82ede220ccbf6b`

**Status:** Historical record; verify current source and core documentation
before using its old commands or provisional conclusions.

<details>
<summary>Expand original: General T3: additional EFT1 scalar-pair Wilson mixing</summary>

# General T3: additional EFT1 scalar-pair Wilson mixing

The Stage 16 B/C(-1) renormalisable UV and scalar-EFT RGE trajectories, together with the **pre-existing baseline** SM+Weinberg evolution, are already running on Windows.

Stage 17 extracts the **new quartic-dependent part** of the one-loop EFT1 psi2phi2 Wilson anomalous dimension. The existing matched F-threshold tree Wilson seed contains both PL y1*y1 and PL y2*y2 operators. The new real-scalar quartic tensor is constructed from the independent SU(2) tensors proven in Stages 9 and 13:

- `lambdaHHdagS2S2`: `V_extra = lambda2 (Hdag_1*S2_1 + Hdag_2*S2_2)^2 + h.c.`;
- `lambdaHHdagS1barS1bar`: `V_extra = lambda1 (Hdag_1*S1dag_2 - Hdag_2*S1dag_1)^2 + h.c.`;

where `lambda1`, `lambda2` are RGBeta **potential** coefficients and Stage 15 established `lambda_Matchete = -lambda_RGBeta` in the common scalar basis. `V_extra = lambda_abcd phi_a phi_b phi_c phi_d/24`; `lambda_abcd` is computed by exact real-field derivatives, with the source's identical-boson projector already incorporated.

The one-loop master tensor RGE's scalar-pair term, holding existing Wilson coefficients fixed, gives

`(16*pi^2) d C_ijab/dln(mu) |_extra = (1/2)*sum_{swap(a,b)} sum_{c,d} lambda_extra_abcd C_ijcd`.

The new script uses the existing `build_eft1_wilson_tensor` Matchete PL adapter, not an invented analytic Yukawa prefactor. It computes the nonzero real Higgs components of the **new piece only** for one lepton generation and checks that both Yukawa families occur.

**Not completed:** full 3-generation flavor Wilson ODE, inclusion of all already existing anomalous dimensions in a coupled Wilson evolution, scalar-threshold one-loop matching, pole cancellation and the final corrected C5. Do NOT add these beta components directly as the final C5; the stage is an operator-mixing audit. All production files stay unchanged. The Stage 16 legacy low-energy masses remain valid baseline values.

Requires the user's current repository `RGE/running/intermediate/ScalarOnlyTensorAdapters.py`, `output/T3_{B,C}_alpha_m1/data/eft1_after_F_wilson_seed.json`, `eft1_rgbeta_rge.json`, and Stage 15 source-convention certificate.

PowerShell from project root:

```powershell
Expand-Archive -LiteralPath .\T3_General_Scalar_Batch17.zip -DestinationPath . -Force
python -m pytest -q tests\python\test_general_t3_e_f_t1_scalar_pair_mixing.py
powershell -ExecutionPolicy Bypass -File .\scripts\RunGeneralT3EFT1ScalarPairMixing.ps1
```

Reports: `output/general_t3_stage17/{B_m1,C_m1}_extra_eft1_scalar_pair_beta.json` and `stage17_summary.json`.

</details>

---

<a id="source-general-t3-eft1-wilson-source-transport-md"></a>

## Original document: `docs/GENERAL_T3_EFT1_WILSON_SOURCE_TRANSPORT.md`

**Original SHA-256:** `c00bec8f0bf4b61cee551892e730456ffbe9f641c3e417aec69bba91d2cf4415`

**Status:** Historical record; verify current source and core documentation
before using its old commands or provisional conclusions.

<details>
<summary>Expand original: General T3: numerical extra-quartic EFT1 LLHH Wilson-source transport</summary>

# General T3: numerical extra-quartic EFT1 LLHH Wilson-source transport

## Verified preceding steps

- Stage 16: grouped two-threshold **renormalisable** UV/EFT1 one-loop RGBeta running and pre-existing baseline Weinberg evolution (both B/C(-1) run successfully on Windows).
- Stage 17: two additional scalar quartic invariants contracted into the Matchete PL EFT1 Wilson seed, producing four nonzero Higgs-real components of the **extra** scalar-pair Wilson anomalous dimension in both models (10/10 local tests and both model exports passed).
- Stage 15: potential/Lagrangian sign relation `lambda_Matchete = -lambda_RGBeta` for these couplings in the common scalar field basis. Stage18 uses the *RGBeta potential* coefficients, so it does NOT insert an additional minus sign.

## Stage18 calculation

Reuse `output/T3_{B,C}_alpha_m1/data/{uv_rgbeta_rge.json,eft1_rgbeta_rge.json}` and the genuine `output/general_t3_stage17/{B_m1,C_m1}_extra_eft1_scalar_pair_beta.json`. Run the Stage16 grouped benchmark again in memory to obtain the time-dependent EFT1 state at 17 logarithmically spaced scales. Use the full existing coupled RGEs for **all** renormalisable EFT1 couplings (including six special scalar quartics); only the Wilson insertion is new.

Stage17's scalar-pair beta components are homogeneous of degree 2 in the one-generation PL Yukawa and inverse-linear in MF. The code validates that they factor into

    beta_ijab^(1)(mu) = A_ab(mu) * conjugate(y1)^2 / MF
                       + B_ab(mu) * conjugate(y2)^2 / MF .

The only allowed free symbols in `A` and `B` are the two active complex RGBeta potential couplings and their conjugates, each to degree at most one. In the positive diagonal heavy-mass basis at the F boundary, promote the factors to symmetric 3x3 matrices:

    Y_i[p,q] = sum_r conjugate(y_i[p,r]) conjugate(y_i[q,r]) / M_r .

The Yukawas, positive heavy masses, and tree Wilson seeds are evaluated **at mu_F** and kept fixed below F, where those heavy fields are absent. Evolve the source only:

    Delta C_ijab|extra-source = (1 / (16*pi^2)) * integral_{
      ln(mu_F)}^{ln(mu_S)} beta_ijab^(1)(mu) dln(mu).

Integration uses composite Simpson quadrature over the saved EFT1 trajectory, plus an independent composite trapezoid estimate (reported as a quadrature difference). Both are downwards in logarithmic scale. Output units: `GeV^-1` for each flavour-matrix component.

### Scope (important)

This **is not the full coupled EFT1 Wilson system**: it omits homogeneous Wilson mixing (except the direct new scalar-pair source), additional operator-generated source feedback, one-loop corrections to the F threshold seed, and finite scalar-threshold matching. It exports **real** Higgs-index Wilson coefficients (`i,j` fixed to lepton SU2 indices `(1,2)`, `a,b` in Higgs real indices 1..4), not a canonically normalised Weinberg coefficient. Therefore it makes *no* new neutrino mass prediction, and does not add a contribution to the pre-existing physical baseline. This is a full **three-lepton-flavour source lift** of a tree-seed single-insertion, not a validated three-generation coupled Wilson ODE.

For radiatively generated extra quartics initialized at zero at the UV scale, this piece has nominal two-loop suppression. The computation is a selected `O(hbar^2)` effect, **not a complete two-loop prediction**.

We retain the Stage16 grouped mass hierarchy and its finite-order approximation, not individual particle thresholds. The `mu_F` tree boundary is continuous for renormalisable parameters; no new finite matching terms are assumed. The original Stage11-15 toy spectrum violates the grouped F-first ordering and is not used.

## Commands (PowerShell from project root)

```powershell
Expand-Archive -LiteralPath .\T3_General_Scalar_Batch18.zip -DestinationPath . -Force
python -m pytest -q tests\python\test_general_t3_e_f_t1_wilson_source_transport.py
powershell -ExecutionPolicy Bypass -File .\scripts\RunGeneralT3EFT1WilsonSourceTransport.ps1
```

Optional numerical consistency scan, uses more sampled points:

```powershell
powershell -ExecutionPolicy Bypass -File .\scripts\RunGeneralT3EFT1WilsonSourceTransport.ps1 -Samples 33
```

## Output files

- `output/general_t3_stage18/B_m1_extra_eft1_wilson_source_integral.json`
- `output/general_t3_stage18/C_m1_extra_eft1_wilson_source_integral.json`
- `output/general_t3_stage18/stage18_summary.json`

The script fails closed and records each model error in `stage18_summary.json`. It does not modify production physics code, the baseline outputs, or GitHub.

## Next physics task

Match these real Higgs Wilson tensors onto the correctly normalised Weinberg operator using the project's established `WeinbergTensorAdapter` and the scalar threshold, then compare with direct Matchete matching. Only after consistency checks should the extra piece be included in the final SM+Weinberg transport.

</details>

---

<a id="source-general-t3-fermion-first-thresholds-md"></a>

## Original document: `docs/GENERAL_T3_FERMION_FIRST_THRESHOLDS.md`

**Original SHA-256:** `b3dd6aff0cd59c2f8c4daf3cf48cf54d7ff8da07488efa8855c9296e900ba7d6`

**Status:** Historical record; verify current source and core documentation
before using its old commands or provisional conclusions.

<details>
<summary>Expand original: General T3: grouped fermion-first threshold continuation (experimental)</summary>

# General T3: grouped fermion-first threshold continuation (experimental)

## Implemented

1. Run the existing *complete currently coded* B/C(alpha=-1) one-loop UV RGBeta system from `mu_UV` to a common `mu_F`, including both C5-active and four other existing RGBeta complex quartics.
2. Apply `build_intermediate_scalar_boundary` at `mu_F`: remove F and both Yukawas; transport the renormalisable quartics continuously (tree-level matching; no finite threshold correction).
3. Run the existing scalar-only **EFT1 RGBeta** one-loop system to `mu_S`, including all six complex quartics.
4. Use the existing scalar-threshold projection to SM.
5. With `final_weinberg_coefficient.json` supplied, use the **previously established baseline** hierarchical C5 evaluator and run the SM+Weinberg matrix to `mu_low`. This baseline DOES NOT contain the new quartic corrections to the Wilson matching.
6. Save all running extra quartic values at UV, F boundary and S boundary, plus baseline C5/mass matrices, and a machine-readable blocked flag for extra Wilson matching.

## NOT implemented

**The generated quartics do NOT yet contribute to the physical C5 at the scalar threshold.** The missing pieces are matching of the relevant higher-dimensional LLSS Wilson tensors at F, their EFT1 anomalous dimensions with extra quartic insertions, and the associated scalar matching and pole cancellation at S. The one-step Matchete hard response is not a substitute. No unsupported extra C5 or extra low-energy neutrino mass is output. The already existing `final_weinberg_coefficient.json` is a baseline coefficient, not an extended-model result. A complete O(hbar^2) calculation also needs other two-loop effects.

**The threshold prescription is grouped**, not individual thresholds for the three different MF and two different MS. The default benchmark enforces the strict ordering `mu_UV > max(MF) > min(MF) > mu_F > max(MS) > min(MS) > mu_S > mu_low` approximately, rechecking running masses. Grouped matching is an approximation until matched finite threshold corrections and sequential heavy-generation integration are supplied. The pre-existing Stage 11–15 benchmark does **not** have this ordering (MF1=1500 GeV < MS2=1780 GeV). Its numerical results are unchanged; this batch creates a separate benchmark with MF=[5100,6000,7000] GeV.

The benchmark values are illustrative, not measured couplings or a physics fit. The source-level `lambda_Matchete=-lambda_RGBeta` relation still assumes an identical scalar-field basis. The production numerical source files are untouched.

## Run

In the repository root:

```powershell
Expand-Archive -LiteralPath .\T3_General_Scalar_Batch16.zip -DestinationPath . -Force
python -m pytest -q tests\python\test_general_t3_fermion_first_thresholds.py
powershell -ExecutionPolicy Bypass -File .\scripts\RunGeneralT3FermionFirstThresholds.ps1
```

For the two renormalisable segments only (no baseline matching):

```powershell
powershell -ExecutionPolicy Bypass -File .\scripts\RunGeneralT3FermionFirstThresholds.ps1 -RenormalisableOnly
```

Requires existing `output/T3_{B,C}_alpha_m1/data/{uv_rgbeta_rge.json,eft1_rgbeta_rge.json,final_weinberg_coefficient.json}`. Outputs `output/general_t3_stage16/{B_m1,C_m1}_hierarchical_trajectory.json` and `stage16_summary.json`.

</details>

---

<a id="source-general-t3-final-c5-and-weinberg-running-md"></a>

## Original document: `docs/GENERAL_T3_FINAL_C5_AND_WEINBERG_RUNNING.md`

**Original SHA-256:** `8f6c91ebc635829d04dd8f10477448654cb65246116e0cb7f71ee74eff2e3c3d`

**Status:** Historical record; verify current source and core documentation
before using its old commands or provisional conclusions.

<details>
<summary>Expand original: General T3: existing model-specific C5 and final SM+Weinberg running</summary>

# General T3: existing model-specific C5 and final SM+Weinberg running

## Scope

Connect the successful Stage 33 / 35 *numerical states* to existing model-specific one-loop matching artifacts. No new symbolic matching is derived. The authoritative code is reused:

- `Numerical/running/FinalC5TrajectoryAdapter.py`: heavy mass-basis conversion and construction of C5 inputs from the heavy threshold at 100 TeV and scalar masses and `lambdaT3` at 1 TeV;
- `Numerical/running/SMWeinbergStage.py`: evaluation of an **existing** symbolic C5 artifact;
- `Numerical/running/ScalarThresholdBoundary.py`: initial SM+Weinberg state;
- `Numerical/running/SMWeinbergEvolution.py`: one-loop running 1 TeV -> 100 GeV;
- `physics/NeutrinoMass.py`, `physics/NeutrinoTrajectory.py`, `physics/NeutrinoObservables.py`: low-energy neutrino mass and charged-lepton basis observables.

Full prior UV and intermediate RGEs are **not** rerun. SHA-256 validations tie Stage 33/34/35 outputs to the inputs. For each model the coefficient file is resolved from its *own* `(class,alpha)` path, never by recycling another model's matching artifact.

## Existing C5 sources

For twelve model/alpha values, the successful hierarchical fixed-order source is:

`output/full/hypercharge/T3_<family>_alpha_<pN_or_mN>/data/final_weinberg_coefficient.json`

The source must be `status=Success`, `scheme=fixed_order_one_loop_MSbar_hierarchical`, have a successful renormalized hard threshold and a validated manifestly symmetric direct-running contribution. `FinalWeinbergAdapter` combines the hard and direct terms.

Four extreme-hypercharge cases (`A_-4`, `B_-3`, `C_-3`, `E_-4`) currently do not have this hierarchical artifact in the inspected repository. They have the *different* legacy one-loop output:

`output/hypercharge/T3_<family>_alpha_<mN>/data/c5_coefficient.txt`

This uses the existing `evaluate_symbolic_c5` parser and *does not* reproduce the same hierarchical hard-plus-direct scheme. To include these four rough benchmarks the caller must explicitly supply `--allow-legacy`. The combined 16-model result is labelled **mixed matching sources**, and source kind is recorded per model. Do not treat it as a uniform precision comparison.

## Missing physics and convention caveats

- No re-evaluation of the full Matchete matching with the 38 new Stage 29 scalar quartic vertices. Those new terms have no certified EFT1 basis transfer and are absent below 100 TeV.
- The 12 reused quartic slots were evolved in Stage 34; they enter the scalar threshold state, but there is no new Wilson-operator matching.
- Existing matched C5 coefficients are evaluated on the Stage 33/35 RGE backgrounds. That is a conditional approximate prediction, **not** a new matched EFT calculation.
- Model-specific mass eigenvalues differ from the grouped matching scales (100 TeV and 1 TeV); this is an explicit common-scale approximation.
- The matching convention / coupling map between Matchete and RGBeta extra quartics remains uncertified.
- The `NO` mass eigenstate labelling follows the fitted-reference benchmark convention. No refit to oscillation data is performed.

## Files / project paths

- `scripts/FinalWeinbergRunning.py` — new complete runner.
- `tests/python/test_general_t3_final_weinberg.py` — new offline regression tests.
- `docs/GENERAL_T3_FINAL_C5_AND_WEINBERG_RUNNING.md` — this document.

All generated JSON outputs are placed under `output/general_t3_stage22`.

## Commands — PowerShell from repository root

Install ZIP (place it in repository root first):

```powershell
Expand-Archive -LiteralPath .\T3_General_Scalar_Stage36_ExistingC5_WeinbergRunning.zip -DestinationPath . -Force
```

Tests:

```powershell
python -m pytest -q tests\python\test_general_t3_final_weinberg.py
```

First model (hierarchical C5):

```powershell
python -m scripts.FinalWeinbergRunning --model A_+0
```

Run all 16 **with explicit opt-in for the four legacy approximations**:

```powershell
python -m scripts.FinalWeinbergRunning --all --allow-legacy
```

Regenerate/check the combined summary without rerunning:

```powershell
python -m scripts.FinalWeinbergRunning --summary
```

## Output

Individual model:

`output/general_t3_stage22/stage36_A_p0_low_energy_neutrinos.json`

All-model summary:

`output/general_t3_stage22/stage36_16_model_low_energy_neutrino_summary.json`

The per-model record carries provenance hashes, matching scheme, heavy/scalar numerical matching inputs, full complex C5 at 1 TeV and 100 GeV, low-energy neutrino mass matrix, `NO` neutrino mass eigenvalues, squared splittings and absolute PMNS entries. The combined summary differentiates the twelve hierarchical and four legacy artifacts. All files use the existing benchmark without refitting.

## Offline verification

The Stage 36 unit tests mock the external numerical services but validate all source routing, missing-artifact failures, alpha isolation, opt-in policy, and handoff contract. Existing Stage 30/31/33/34/35 Python regressions also pass in a local isolated staging tree. *The actual Stage 36 matching evaluations and numerical integration must run in the user's full Windows repository.*

</details>

---

<a id="source-general-t3-weinberg-source-transport-md"></a>

## Original document: `docs/GENERAL_T3_WEINBERG_SOURCE_TRANSPORT.md`

**Original SHA-256:** `f9df8e69428ae01b92dbf573cd3dad68e7b52796eb1aa48d14b58a96973ff3b4`

**Status:** Historical record; verify current source and core documentation
before using its old commands or provisional conclusions.

<details>
<summary>Expand original: General T3: selected extra-quartic Weinberg projection and low-energy transport</summary>

# General T3: selected extra-quartic Weinberg projection and low-energy transport

## Established inputs

- Stage 16: B/C(-1) grouped, fermion-first UV and scalar-only EFT RGBeta running, including all six special quartics, and **baseline** physical hierarchical C5.
- Stage 17: **new** quartic-dependent scalar-pair part of the one-loop dimension-five EFT1 Wilson-tensor anomalous dimension, based on Matchete PL LLSS tree seed.
- Stage 18: numerical integral of this **selected source**, using running EFT1 renormalisable couplings and *frozen* tree LLSS Wilson seed at mu_F. Four nonzero mixed-lepton Higgs-real components in both models (actual Windows runs).

## New verification and calculation

**Symbolic closure:** Use the project's authoritative `RGE.matching.WeinbergTensorAdapter.build_weinberg_wilson_tensor(kappa=1)` to generate the same SU(2) real-scalar template that validates the production EFT1 Wilson RGE. Independently recalculate the extra scalar-pair beta from the actual Matchete PL Wilson seed and the two SU(2) potential tensors for **all 30 independent lepton-doublet/Higgs-real index combinations** (`i<=j`, `a<=b`). Require exact identity

    delta_beta^(1)_{ijab} = beta_kappa^(1) * WeinbergTemplate_{ijab}

and require the four mixed-lepton components to agree with the saved Stage 17 expressions. If an index/sign/phase fails, **no numerical C5 is exported**. This check must be performed with the *actual local* seed; it has not been run here.

**Physical normalisation:** Independently compare the existing successful `eft1_wilson_rge.json` `weinberg_subspace_validation.beta_kappa_16pi2` with the `final_weinberg_coefficient.json` direct-running expression after inserting the canonical one-generation `C12 = conjugate(y1)*conjugate(y2)/MF` and extracting its coefficient of `log(MS/MF)`. This establishes that the RGE tensor's `kappa` uses the same one-loop direct-running C5 convention as the established physical matching adapter.

**Numerical projection:** Given the certified complex Weinberg template and the four integrated Stage 18 `C_{nu,e,a,b}`, the mixed-lepton template ratios are

    (a,b)=(1,3): -1/4
          (1,4): -i/4
          (2,3): -i/4
          (2,4): +1/4

for polynomial `kappa*(nu*H0 - e*H+)^2` and tensor `C=(1/4) fourth derivative`. Project by least squares and verify the *complete matrix* residual is smaller than 2e-8 relative, rather than using merely the four component norms.

**Scalar boundary:** Since the selected source is already LLHH with no S1 or S2 legs, the tree-level projection at the scalar threshold simply carries its coefficient to SM+Weinberg. **Finite one-loop scalar-threshold effects are not computed.**

**SM running:** Pass this **selected-source-only** symmetric 3x3 matrix through the existing `T3Trajectory.continue_with_weinberg_running`, from the grouped scalar matching scale `mu_S=1100 GeV` to `mu_low=246.22 GeV`, and convert each endpoint by `delta m_nu = -(246.22 GeV)^2/2 * delta C5`. Saved outputs remain separate from Stage 16 baseline.

## Deliberate limitations

- Full coupled EFT1 LLSS Wilson evolution and closure are not integrated (source-only, frozen LLSS seed).
- F and S matching remain **grouped**; individual heavy and scalar mass thresholds are not crossed.
- No finite extra contribution at the F or scalar threshold is included.
- With extra quartics generated from UV zero boundary, this contribution starts at **formal O(hbar^2)**; this is **not** a complete two-loop calculation, and the low-energy source-only output is **not a full physical neutrino mass prediction**.
- The running baseline and production matching results are not overwritten, nor are the source or GitHub files modified.

## Install and run (PowerShell, repository root)

```powershell
Expand-Archive -LiteralPath .\T3_General_Scalar_Batch19.zip -DestinationPath . -Force
python -m pytest -q tests\python\test_general_t3_weinberg_source_transport.py
powershell -ExecutionPolicy Bypass -File .\scripts\RunGeneralT3WeinbergSourceTransport.ps1
```

Stage 19 uses existing Stage 17/18 user-generated reports and the *production* Matchete/RGBeta exports under `output/T3_B_alpha_m1/data/` and `output/T3_C_alpha_m1/data/`; no fresh symbolic matching or Stage 18 quadrature is required. It reruns two short renormalisable RGBeta numerical segments to recover the current consistent SM matching boundary.

## New files

- `Numerical/experimental/GeneralT3WeinbergProjection.py`
- `Numerical/experimental/GeneralT3ExtraLowEnergyTransport.py`
- `scripts/RunGeneralT3WeinbergSourceTransport.py`
- `scripts/RunGeneralT3WeinbergSourceTransport.ps1`
- `tests/python/test_general_t3_weinberg_source_transport.py`
- `docs/GENERAL_T3_WEINBERG_SOURCE_TRANSPORT.md`

## Output

- `output/general_t3_stage19/B_m1_selected_extra_weinberg_transport.json`
- `output/general_t3_stage19/C_m1_selected_extra_weinberg_transport.json`
- `output/general_t3_stage19/stage19_summary.json`

A failure writes an exact error to `stage19_summary.json` and exits nonzero, without manufacturing a `delta C5` value.

</details>

---
