# General T3 matching resolution log

Purpose: document source/code changes, actual verification, and limitations as
we work from the complete symbolic scalar Lagrangians toward physically
consistent neutrino-mass predictions. **Do not convert a conditional numerical
result into a certified result without the stated independent tests.**

## Baseline before this package

- Complete symbolic scalar operator inventory exists for 16 T3 cases.
  14 have full Matchete Lagrangian checks; B(−1)/C(−1) contain off-diagonal
  scalar masses and do not pass the heavy-mass-basis check in original fields.
- All 16 candidate UV/EFT1 RG trajectories have been computed.
- B(−1)/C(−1) have numerical scalar diagonalization and a conditional direct
  one-scale mass-eigenstate `C5` procedure, but the full native interaction
  basis equivalence and sequential EFT matching are not certified.
- Six A/D/E self-bilinear cases have conditional mass-dependent anomalous
  `C5` calculations and low-energy benchmark corrections; their RGBeta–Matchete
  vertex normalization/phase map is not fully certified.
- The common complete sequential fermion-first/scalar-threshold matching is
  not certified across the full 16-model comparison.

These are historical saved outputs and audited source code, not newly generated
measurements from this package.

## Milestone BC-NATIVE-Q4-01: actual native extra-quartic tensor transport

**Scope:** B(−1) and C(−1), each with 7 saved native Matchete extra quartic
invariants. Do not conflate the 7 *Matchete* basis tensors with the source
RGBeta registry's 7 coupling slots without an explicit coupling-basis map.

**Files added:**

- `scripts/AuditT3BCNativeQuarticMassBasis.py`
- `tests/test_t3_bc_native_quartic_mass_basis.py`
- `docs/T3_BC_NATIVE_QUARTIC_MASS_BASIS_AUDIT.md`
- `docs/T3_GENERAL_SCALAR_MATCHING_RESOLUTION_LOG.md`

**Method:** Load saved Stage22 Matchete native component tensors; revalidate
ordered species, components, potential prefactors; load saved mixed-doublet
mass-basis report and recheck its mass matrix against the unitary rotation.
Transport each quartic with the existing epsilon-dependent component field map.
Check original vs rotated complex polynomial and the Hermitian-conjugate
potential at random field points. Require rank 2 for the pair of invariants
of counts `110101`. Store every transformed component tensor and SHA-256
provenance in a new audit output directory, without overwriting prior physics.

**Verification completed while authoring this package:** Python compilation
and five pytest checks using an isolated synthetic test repository populated
with the seven saved B/C tensor component patterns and a representative
mixed-doublet mass matrix. Both synthetic end-to-end audit models passed;
relative polynomial residual was about `3.6e-15` in that fixture.

**Still required from local user repository:** run pytest and the audit on the
**actual** saved Stage22 and mass-basis files. Record the local `summary.json`
produced by that run. Do **not** describe the actual project native export as
verified by this package until those commands pass locally.

**Physical certification:** NOT achieved by this milestone. The existing
Matchete/RGBeta and matching-certification flags remain false.

## Next gated milestones (not yet implemented here)

1. **BC-NATIVE-ALL-02:** transport source-native original scalar quartics and
   Yukawa Clebsch–Gordan tensors (plus relevant kinetic/gauge conventions);
   verify all interactions in the same component field map.
2. **BC-MAP-03:** obtain independent RGBeta-native invariant/potential
   coefficients; compute the exact per-block coupling-basis maps, including
   the `110101` 2×2 sector and Hermitian-conjugate orientation.
3. **BC-MATCH-04:** reassemble a real mass-eigenfield Matchete UV/EFT
   Lagrangian, obtain live `CheckLagrangian` and `HeavyMassBasis` passes and
   independently match its one-loop `C5` to the existing candidate.
4. **SELF-MAP-05:** independent component-vertex and mass-eigenstate one-loop
   matching audit for the six A/D/E self-bilinear models, including native
   RGBeta–Matchete sign, phase and coefficient normalization.
5. **EFT-THRESHOLDS-06:** test the remaining Wilson coefficient threshold
   conditions, intermediate EFT running and matching scheme alignment across
   16 models. Validate physical observables only after these checks.

For current honours-level rough benchmarks, existing candidate numbers remain
usable *as labelled approximations*. This log records only verified gains and
must be updated when each subsequent milestone passes.

## User-local confirmation: BC-NATIVE-Q4-01 (2026-10-31)

The user ran:

```powershell
python -u -X utf8 -m scripts.AuditT3BCNativeQuarticMassBasis --all --run-id bc_native_vertices_v1
```

Observed terminal output:

```text
B_-1: native tensors 7/7; 110101 rank=2; polynomial residual=3.553e-15
C_-1: native tensors 7/7; 110101 rank=2; polynomial residual=3.657e-15
```

The actual local source-based tensor-transport audit is therefore **passed**
for both models, not just the isolated synthetic fixture. The user-local
`summary.json` was not separately inspected; this record is grounded in the
reported terminal output. Native Matchete matching and physical normalization
remain uncertified.

## Milestone BC-MAP-GEOMETRY-02: native 110101 map after scalar mixing

**Added in the follow-up package** (without altering any production code):

- `scripts/AuditT3BC110101MassBasisMap.py`
- `tests/test_t3_bc_110101_mass_basis_map.py`
- `docs/T3_BC_110101_MASS_BASIS_MAP.md`
- updated `docs/T3_GENERAL_SCALAR_MATCHING_RESOLUTION_LOG.md` (this file)

Uses the two *live Wolfram-exported* registered RGBeta 110101 tensors and the
two independently saved native Matchete tensors to derive an exact 2×2
invariant-space map. The new test checks its equality **after the actual
mixed-doublet mass-basis transformation**. Inputs and output schemas include
provenance hashes and explicit false physical-certification flags.

**At package authoring:** local fixture/regression checks passed; a local
project run of `bc_110101_map_v1` has not yet been observed, so **do not mark
this milestone verified on the user's actual repository yet**.

**Remaining B/C problem:** the native RGBeta coupling function contains an
internal normalization factor, but the **physical potential/vertex prefactor**
(and absolute Lagrangian sign) has not been independently certified. The exact
basis map is **not** yet a usable map of the raw coupling constants.

**Next gated step:** establish native potential prefactors and independent
scalar/Yukawa vertex mappings. A mass-basis `CheckLagrangian`/`Match[]`
certificate must remain blocked until those succeed.

## User-local confirmation: BC-MAP-GEOMETRY-02 (2026-10-31)

The user ran the actual project command:

```powershell
python -u -X utf8 -m scripts.AuditT3BC110101MassBasisMap --all --run-id bc_110101_map_v1
```

Terminal output reported:

```text
B_-1: exact 110101 rank-2 map checked; mass-basis residual= 2.220e-16
C_-1: exact 110101 rank-2 map checked; mass-basis residual= 1.110e-16
Summary: C:\thesistemp\T3neutrinomassmodelpipeline\output\audits\t3_bc_110101_mass_basis_map\bc_110101_map_v1\summary.json
```

**Milestone now PASSED on the user's actual local repository**, based on the
reported terminal output. The `summary.json` contents were not separately
inspected. This verifies the rank-two *native invariant-tensor geometry* and
its consistency under the saved mixed-scalar mass-basis transformations.
It does **not** certify the absolute RGBeta/Matchete coupling normalization,
overall Lagrangian/potential sign, Yukawa CG mapping, native Matchete
`HeavyMassBasis` success, or sequential Wilson matching. Do not change
physical-certification flags based on this result alone.

### Follow-up: RGBeta quartic normalization source inspection (2026-10-31)

Read-only inspection of the external RGBeta source at
`aethomsen/RGBeta`, `Package/FieldsAndCouplings.m` (default branch `master`):

- `AddQuartic` (approximately lines 751--786) sets an internal quartic
  normalization `24` and multiplies it by `1/Sqrt[2]` for each complex
  scalar occurrence. Thus all four-complex-scalar 110101 operators have
  internal normalization `24(1/Sqrt[2])^4 = 6`, consistent with the saved
  native `6*Matrix[coupling][##1] &` metadata.
- `CouplingTensorMap` (approximately lines 925--941) uses
  `AvgPermutations` and optionally the conjugate coupling.
- `Lam` (approximately lines 1035--1046) includes the Hermitian-conjugate
  contribution for non-self-conjugate quartics. This is consistent with the
  source `SelfConjugate -> False` setting.
- The RGBeta publication states that the quartic normalization is defined
  explicitly by the model's `GroupInvariant`; RGBeta does not automatically
  assume conventional symmetry factors.

References:
- https://github.com/aethomsen/RGBeta/blob/master/Package/FieldsAndCouplings.m
- https://link.springer.com/article/10.1140/epjc/s10052-021-09142-4

**Interpretation:** These findings explain the internal factor `6` and
conjugate bookkeeping; they do **not** determine the absolute physical
coefficient or relative sign in the independently defined Matchete Lagrangian
without an action-level convention check. In the verified tensor basis,
`R = M A`, with `A = [[1/2,1],[-sqrt(3)/2,0]]`. A physical coupling map has
the form `c_M = sign * A * D_RGB * c_RGB`, where the potential prefactor
matrix `D_RGB` and the sign must be established independently. No assumed
unit matrix or sign should be inserted into the production matching code.

See `docs/T3_BC_RGBETA_NORMALISATION_NEXT_CHECK.md` for the next
source-level verification and its acceptance criteria.

## Milestone BC-110101-POTENTIAL-03: source-convention potential reconstruction (2026-10-31)

**Package prepared; user-local source run not yet observed.** The new files are:

- `scripts/AuditT3BC110101PotentialMap.py`
- `tests/test_t3_bc_110101_potential_map.py`
- `docs/T3_BC_110101_POTENTIAL_NORMALISATION.md`
- this complete updated `docs/T3_GENERAL_SCALAR_MATCHING_RESOLUTION_LOG.md`

Building on the user-local **PASSED** B/C native-tensor and rank-two map
reports, the script reconstructs the action-level quartic convention from
RGBeta's published generic \(\mathcal{L}_4=-\Lambda\Phi^4/4!\) and the
actual native registry metadata. The two `110101` RGBeta couplings each have
internal `6*Matrix[...][##1]`, 24 unique arrangements, and
`SelfConjugate=False`. For complex fields,
\(6(\sqrt2)^4/4!=1\); the **internal 6 is not a physical factor six**.

The Matchete source contains `PlusHc[contribution]` with unit coupling
prefactor and native quartic tensors \(M_a\). With \(R=M A\) and the previously
established exact geometry
\(A=\bigl(\begin{smallmatrix}1/2&1\\-\sqrt3/2&0\end{smallmatrix}\bigr)\),
**the source-convention candidate map** is

\[
 h_1=-g_1/2-g_2,\qquad h_2=\sqrt3\,g_1/2.
\]

A regression/fixture test of the reconstructed potential, independent complex
coupling phases, Hermitian conjugates, and both mixed mass-basis rotations
passed locally (14 tests including the previous map regressions). Source-
reconstruction polynomial residuals were O(10^-15). The fixture is NOT a
live evaluation of the installed RGBeta package, and no local project run
of this new stage was supplied at authoring time.

**Strict status:** This closes a **source-derived potential normalization
reconstruction check**. It does NOT independently certify installed RGBeta
`Lam` evaluation, exact Matchete action on a newly transformed Lagrangian,
all scalar couplings, Yukawa vertices, or physical `C5` matching. All
certification fields remain `false` and production code is unchanged.

**Next gate:** independent native RGBeta `Lam` real-field projection versus
Matchete `PlusHc` action for complex test fields, then the remaining B/C
quartic/Yukawa vertices and mass-eigenstate matching.

## User-local confirmation: BC-110101-POTENTIAL-03 (2026-10-31)

The user ran:

```powershell
python -u -X utf8 -m scripts.AuditT3BC110101PotentialMap --all --run-id bc_110101_potential_v1
```

Observed terminal output:

```text
B_-1: source reconstructed 110101 potential; sign residual=4.829e-15; mass-basis residual=1.806e-15
C_-1: source reconstructed 110101 potential; sign residual=2.220e-15; mass-basis residual=3.010e-15
Summary: C:\thesistemp\T3neutrinomassmodelpipeline\output\audits\t3_bc_110101_potential_map\bc_110101_potential_v1\summary.json
```

**Status: PASSED** for the *source-reconstructed* potential normalization and
mass-basis polynomial, on the actual user-local repository. The summary JSON
was not separately inspected. No native installed RGBeta `Lam` evaluation or
live Matchete action-level comparison occurred in this audit. The physical
coupling and matching certification flags therefore remain false.

## Follow-up package: BC-NATIVE-LAM-04 (diagnostic, not yet user-run)

Prepared new files:

- `scripts/ProbeT3BCNativeRGBetaLam.wl`
- `scripts/VerifyT3BCNativeLamProbe.py`
- `tests/test_t3_bc_native_rgbeta_lam_probe.py`
- `docs/T3_BC_NATIVE_RGBETA_LAM_PROBE.md`
- Updated complete `docs/T3_GENERAL_SCALAR_MATCHING_RESOLUTION_LOG.md`

The probe loads the installed RGBeta package in Wolfram, builds B(-1), C(-1),
loads the same registered 110101 vertices used by the RGEs, calls the actual
`Lam` tensor function, and exports structural runtime evidence. The Python
checker fails closed if the native expression or its quartic couplings cannot
be demonstrated. The runtime has not been executed by the assistant.

**Not claimed:** action-level coefficient equality, all-quartic normalization,
Yukawa vertex matching, heavy scalar `CheckLagrangian`/`HeavyMassBasis` passes,
or sequential one-loop matching. The next stage must use the actual live `Lam`
returned tensor, once its runtime representation is established, to construct
an independent quartic potential and compare it with the saved Matchete source.

## User-local issue and correction: BC-NATIVE-LAM-04 (2026-10-31)

The user ran the first Wolfram probe on their installed RGBeta v1.2.0 and
received `Cannot uniquely resolve installed RGBeta Lam symbol: {}`.
No `native_lam_v1.json` was created, so the subsequent Python verifier
raised `FileNotFoundError`. **This was a namespace-lookup failure in the
new audit script, not a physics or RGBeta tensor failure.**

Root cause verified by reading `Package/FieldsAndCouplings.m`: RGBeta declares
`PackageScope["Lam"]`. Wolfram's structured package convention places the
symbol in `RGBeta`PackageScope`Lam`. The original probe only searched the
main and `Private` contexts, so it could not find the function.

Replacement `scripts/ProbeT3BCNativeRGBetaLam.wl` uses the exact
`RGBeta`PackageScope`Lam` name, requires `Names[lamName]` to contain it,
and requires nonempty `DownValues` before calling it. The existing
`VerifyT3BCNativeLamProbe.py` is unchanged.
`tests/test_t3_bc_native_rgbeta_lam_probe.py` now contains a source regression
for this resolution, and `docs/T3_BC_NATIVE_RGBETA_LAM_PROBE.md` documents
the correction and v2 rerun commands.

**Current status:** Correction packaged and Python checks completed by the
assistant; live Wolfram validation on the user's local RGBeta is **pending**.
Do not mark native `Lam` comparison or physical matching certified based on
this namespace correction. No production Lagrangian, beta functions,
benchmark data, or GitHub repository files have been altered.

## User-local live RGBeta Lam probe PASSED: BC-NATIVE-LAM-04 (2026-10-31)

The user executed the scope-corrected v2 native Wolfram probe on RGBeta 1.2.0:

```
python -X utf8 -m pytest -q tests/test_t3_bc_native_rgbeta_lam_probe.py
# 10 passed in 0.09s
wolframscript -file scripts/ProbeT3BCNativeRGBetaLam.wl output/audits/t3_bc_native_lam/native_lam_v2.json
# Resolved installed RGBeta Lam: RGBeta`PackageScope`Lam (downvalues=1)
# B_-1: Lam tensor captured; sparse arrays=1; 110101 couplings visible=True
# C_-1: Lam tensor captured; sparse arrays=1; 110101 couplings visible=True
python -u -X utf8 -m scripts.VerifyT3BCNativeLamProbe --input output/audits/t3_bc_native_lam/native_lam_v2.json --output output/audits/t3_bc_native_lam/native_lam_v2_verified.json
# B_-1: live Lam captured; sparse arrays=1
# C_-1: live Lam captured; sparse arrays=1
```

**Passed:** live native `Lam` function resolves, evaluates, and contains the
required B/C coupling symbols and sparse array. The JSONs were not separately
inspected by the assistant. **Unverified:** native per-component coefficients,
action/potential normalization and physical matching. This is a successful
structural probe, not a physics certificate.

## New conditional gate: BC-NATIVE-LAM-05 (not yet run in user local Wolfram)

New `scripts/ProbeT3BCNativeLam110101Coefficients.wl` queries the actual
runtime sparse array at the four distinct-field ordering and its conjugate.
It tests `6/24=1/4` of the registered invariant as the coefficient at
64 entries per model (two couplings times two orientations times 16 SU(2)
components). `scripts/VerifyT3BCNativeLam110101Coefficients.py` validates the
export, with regression tests and detailed method notes in
`docs/T3_BC_NATIVE_LAM_110101_COMPONENTS.md`. The Wolfram check may uncover a
normalization discrepancy. Do not infer physical certification from a source
reconstruction or from schema tests.

**Status:** packaged for user-local execution; no actual native coefficient
result is claimed. No GitHub or production code changed.

## User-local BC-NATIVE-LAM-05 failure and extraction diagnostic (2026-10-31)

The user ran the component probe with RGBeta v1.2.0. The test suite printed
`12 passed in 0.07s`; the Wolfram run then printed:

```
Checking live Lam source coefficients for B_-1 ...
ERROR: Lam component has no native Matrix terms
Native Lam component coefficient audit failed to construct report
```

The PowerShell failure check correctly stopped before Python verification.
There is no evidence that the native component coefficient test passed.
Previous B/C native `Lam` *presence* and potential reconstruction passes
remain valid and separate; this component extraction failure does **not**
disprove their geometry, and does not certify the action normalization.

Root cause **not yet identified**: either the selected `SparseArray` entry is
not at the assumed field-index location, or the selected expression's
`Matrix[coupling]` structure is not recognized. The updated **full replacement**
`ProbeT3BCNativeLam110101Coefficients.wl` uses improved nested-head discovery,
produces JSON even on a blocked model, and records the actual `SparseArray`
entries, field index map and coupling-bearing sparse rules so the failure can
be diagnosed from the user's live runtime. The related test and instructions
have been updated. This is an audit-tool correction, not a change to any
physics formula or production RGE. Physical certification flags remain false.

## BC-NATIVE-LAM-05: local v2 diagnostics and v3 extraction fix (2026-10-31)

**Observed on the user's installed RGBeta v1.2.0:** the v2 live probe wrote
`output/audits/t3_bc_native_lam/native_lam_coeff_v2.json`, but returned
nonzero for both B(-1) and C(-1). Both requested **direct** entries were

```text
(lambdaT3GeneralQ110101I1*del[SU2L[fund], t3a, t3d]*eps[SU2L[fund], t3b, t3c])/4
+ (lambdaHHdagS1barS2barCross*del[SU2L[fund], t3a, t3c]*eps[SU2L[fund], t3b, t3d])/4
```

Each full tensor reported 48 coupling-bearing sparse-array rules. The known
`6/24=1/4` ordered coefficient is therefore visible **directly** for the two
registered quartics. The v2 extractor incorrectly required surviving
`Matrix[...]` heads even though the evaluated `Lam` sparse entry contains
**plain coupling symbols**. The v2 failure is a probe bug, not a detected
normalization mismatch. All h.c. components and the full 64/model gate remain
**unverified** at this point.

**v3 correction supplied locally (not yet run on the user machine):** full
replacement `scripts/ProbeT3BCNativeLam110101Coefficients.wl` uses native
`Coefficient` on evaluated symbols, resolves either `Bar[coupling]` or
`Conjugate[coupling]` for the conjugate entry, requires no missing/extra
contributions, and preserves fail-closed diagnostic JSON. Full replacement
Python verifier/test and instructions accompany it. The existing matching,
production RGEs, original benchmark data and remote GitHub are unchanged.
No action-level Matchete or physical matching certification is claimed.

## BC-NATIVE-LAM-05: user-local v3 scoping failure; v4 iterator correction (2026-10-31)

The user ran `ProbeT3BCNativeLam110101Coefficients.wl` v3 against installed
RGBeta v1.2.0. Both B(-1) and C(-1) produced `Part::pkspec1` /
`Extract::psl1`; `slot$...` was left unevaluated in
`positions[[slot$...]]`. The script exited nonzero before writing
`native_lam_coeff_v3.json`. No Python verifier was run and no native
component-level normalization was certified.

**Verified cause (full replacement script):** the outer Wolfram `Do` was
written as `Do[...; {slot,1,2}]`, which is **one body expression without an
iterator argument**, instead of `Do[..., {slot,1,2}]`. The v4 full replacement
uses the correct comma, adds a regression test that rejects the accidental
semicolon, and keeps the same exact native component-check logic and false
physical-certification flags. The user's previously verified native Lam
presence, source-reconstructed potential and 110101 tensor geometry are
unchanged. No production physics files or GitHub files are modified; a
successful v4 live result is not yet available.


## BC-NATIVE-LAM-05: live v4 success (user-local, 2026-10-31)

The user ran `wolframscript -file scripts/ProbeT3BCNativeLam110101Coefficients.wl
output/audits/t3_bc_native_lam/native_lam_coeff_v4.json`, then
`python -u -X utf8 -m scripts.VerifyT3BCNativeLam110101Coefficients`.
Both finished successfully. The live Mathematica probe used the installed
RGBeta v1.2.0 `RGBeta`PackageScope`Lam` function; the Python verifier printed:

```text
B_-1: live Lam component coefficients 64/64; max residual=0.000e+00
C_-1: live Lam component coefficients 64/64; max residual=0.000e+00
Native component coefficients: PASSED; full physical matching: NOT CERTIFIED
```

This establishes the actual native RGBeta direct/h.c. coefficients of the
B/C `110101` pair at the 64-component-per-model level with zero mismatch.
The source-reconstructed `6/24=1/4` ordered-tensor normalization is supported
by an **independent live RGBeta Lam evaluation**, not merely a synthetic
Python fixture. The raw JSONs were not uploaded for a separate audit of
every field; the successful verifier output is the user-local evidence.

**Limits:** No independent live Matchete `CheckLagrangian` or action-level
cross-tool identity has yet passed; no physical coupling-map or matching
certificate is claimed. Production RGE, neutrino masses, and Github remain
unchanged.

## New gate: BC-MAT-ACTION-06 (prepared; not yet run locally)

New read-only Wolfram script `scripts/ProbeT3BCMatchete110101Action.wl`
requests the existing T3 Matchete builder's seven extra scalar quartic
candidates for each B/C model, selects the two `110101` interactions,
records their live expressions, and calls `CheckLagrangian` on the two
individual terms and their sum alongside `LSM+LFree`. Its companion
`VerifyT3BCMatchete110101Action.py` rejects incorrect inventories, missing
native couplings, failed checks and any premature cross-tool or matching
certification. A detailed scope statement and run commands are in
`docs/T3_BC_LIVE_MATCHETE_110101_ACTION.md`.

This is an **independent Matchete-side action validation gate**, not yet the
cross-tool action-polynomial equality check. The scalar mixed mass is not
inserted in this gate; heavy mass-basis matching remains blocked. The new
scripts and docs do not change production or GitHub files.

## BC-MAT-ACTION-06: v1 local FrontEnd failure, v2 lifetime fix (2026-10-31)

**User-local evidence:** `ProbeT3BCMatchete110101Action.wl` v1 loaded
Matchete and attempted B(-1), then C(-1). Both printed
`FrontEndObject::notavail: A front end is not available; certain operations
require a front end.` The report construction failed and the shell guard
stopped before the Python verifier. No Matchete action validation or
cross-tool physical coupling certification was obtained.

**Source-confirmed tool bug:** v1 placed only `Needs["Matchete`"]` inside
`UsingFrontEnd`; its later model construction and `CheckLagrangian` calls
executed after that scope had ended. Existing repository exporter
`ExportGeneralT3ScalarTensorMatchete.wl` wraps the entire native tensor
calculation in `UsingFrontEnd`. The v2 full replacement applies the same
front-end lifetime to both `T3BCMOne` calls. Tests protect that scope.

**No physics change:** Original `110101` Matchete terms, CGs, mass-basis
conventions, verified RGBeta `Lam` evidence and all candidate neutrino outputs
remain unchanged. All certification flags stay false; v2 still requires
the user's live Wolfram/Matchete run. Remote GitHub was not modified.

## BC-MAT-ACTION-06: live v2 success on user machine (2026-10-31)

The user executed `wolframscript -file scripts/ProbeT3BCMatchete110101Action.wl
output/audits/t3_bc_matchete_action/matchete_110101_v2.json` and the
corresponding `scripts.VerifyT3BCMatchete110101Action` Python verifier.
Matchete loaded using the v2 full-lifetime `UsingFrontEnd` scope. Both
B(−1) and C(−1) reported **two of two individual `CheckLagrangian` passes,
and a `True` combined `CheckLagrangian` pass**; the independent Python
verifier also passed for both. Native extra-quartic inventory remains seven
per model, with two independent `110101` terms. The full scalar mixed-mass
quadratic was deliberately excluded, so this is an interaction validation,
**not** successful mixed-heavy mass-basis matching or cross-tool equality.
The actual JSON reports have not been separately inspected by the assistant;
this records the user's observed terminal output.

## BC-CROSS-ACTION-07: live-evidence-gated exported-component identity (prepared)

Added read-only `scripts/AuditT3BCLive110101ActionBridge.py`, regression
tests `tests/test_t3_bc_live_110101_action_bridge.py`, and method notes
`docs/T3_BC_LIVE_110101_ACTION_BRIDGE.md`. The new script revalidates the
Matchete v2 and RGBeta `Lam` v4 *live report schemas*, rechecks the saved
Stage22 live registered RGBeta contraction and Matchete native tensor
records, and proves the exact direct and conjugate **source-exported**
action-coefficient equality with the relative `L=-V` sign. It rejects
spurious internal `1/4`, wrong overall sign, or swapped contractions.

**Authoring validation:** 11 local regression tests passed against a
separate synthetic/live-schema fixture populated with saved native tensor
patterns. The user-local command for this new audit has not been run yet;
no successful user-local BC-CROSS-ACTION-07 report is claimed.

**Important boundary:** The previous live Matchete JSON records its
`PlusHc` interaction expressions and `CheckLagrangian` results, but does not
numerically extract the **registered** Matchete CG components. This bridge
uses the independently exported native Matchete basis tensors, not a
fresh direct live `CG` evaluation; consequently the full live cross-tool
action and physical coupling normalization remain **uncertified** until
that last independent native CG-to-action component test is done. Nothing
in the production RGEs, neutrino candidates or remote GitHub changes.

## BC-CROSS-ACTION-07: user-local PASS (2026-10-31)

User ran `pytest` on `tests/test_t3_bc_live_110101_action_bridge.py`:
**11/11 passed**, then `python -u -X utf8 -m
scripts.AuditT3BCLive110101ActionBridge --run-id bc_live_action_bridge_v1`.
B(−1) and C(−1) each passed **64/64 exact exported-basis action
coefficient identities**, with the sign and spurious extra-factor controls
rejected. The summary was written to
`output/audits/t3_bc_110101_live_action_bridge/bc_live_action_bridge_v1/summary.json`.
This records the **user's terminal output**; the raw JSON was not uploaded.
The live registered Matchete CG components were not part of that run.

## BC-LIVE-CG-08: new live registered CG gate (prepared, not run)

`ProbeT3BCMatchete110101RegisteredCG.wl` reconstructs the same two Matchete
`PlusHc` interactions using the original builder, verifies both CG symbols
occur in those interactions, then calls Matchete's **live `GetCGTensor`**
accessor (not `T3ExtendedBosonicTensorBasis`) for both B/C models inside a
single `UsingFrontEnd` scope. The saved report records 4 registered tensors,
32 components per model, and on extraction failure writes a **blocked**
diagnostic instead of falsely passing. The Python verifier revalidates the
live Matchete/RGBeta prior reports, the successful exported-basis action
bridge, source hashes and exact Stage22 native tensors; it compares all
registered CG components exactly, then repeats the signed action identity
with those **live-extracted Matchete tensors** and the saved live-RGBeta
covariance source. No production RGE or Lagrangian changes are made.

Not yet executed in live Mathematica. In particular, the GetCGTensor return
shape is a runtime check, not assumed verified by local Python fixtures.
Even on success, off-diagonal heavy mass-basis Matchete matching,
fermion/scalar threshold matching, and physical neutrino predictions are NOT
certified. No remote GitHub modifications.

## User-local confirmation: BC-CG-08 registered Matchete CGs (2026-10-31)

The user ran `ProbeT3BCMatchete110101RegisteredCG.wl` followed by
`VerifyT3BCMatchete110101RegisteredCG` with run ID
`bc_live_registered_cg_v1`.

- B(−1): native registered Matchete CG tensors **2/2**; exact components
  **32/32**; signed direct + Hermitian-conjugate coefficients **64/64**.
- C(−1): identical passes (**2/2**, **32/32**, **64/64**).
- The prior exported-component bridge also passed **64/64 per model**, with
  incorrect-sign and extra-factor controls rejected.
- This closes the B/C `110101` exported-basis quartic conventions for the
  stated checks. It does not perform direct live action subtraction or
  sequential physical matching.

## BC-MASS-EIGEN-YUKAWA-09: reuse the existing full Lagrangian and one-scale C5

Source inspection confirms that `Lagrangian/RunGeneralT3MassEigenMatchete.wl`
**already** constructs the complete B/C heavy mass-eigenstate Lagrangian,
with 17 quartics, both native Yukawa CG vertices and a detailed
`HeavyMassBasis` check. Current saved B/C reports have status
`MatcheteMassEigenHeavyBasisPassedPhysicalMatchingNotPerformed` and
`heavy_mass_basis_detailed_check_passed=true`. Thus rebuilding the same
Lagrangian would duplicate an existing implementation.

This package adds the new **read-only**
`scripts/RunT3BCMassEigenYukawaC5Handoff.py` to extract physically distinct
Y1/Y2 matrices in the X/Y basis from verified native tensor factors,
check the existing matched one-scale ordered Weinberg expression provenance,
and check/replay saved candidate B/C neutrino mass matrices.

The native Yukawa-2 CG epsilon factor differs: B(−1) `−1`, C(−1) `+1`.
Including the `−epsilon` field-replacement sign, the y2 mass-eigenstate
vertex weights are `+U[1,a]` for B and `−U[1,a]` for C. Both have y1
weights `U[0,a]`. These are coefficients in the respective native CG/spinor
bases; not a claim that the two independent Yukawa spinor vertices add.

Saved one-scale candidate Takagi masses (eV):

- B(−1): `0.027341038758356372, 0.02072521457148339, 0.003616690047615204`
- C(−1): `0.05353064456270797, 0.029391040942796477, 0.005611724480334992`

The new handoff is **prepared, not yet user-local run**. Optional
`--recompute-c5` uses the existing numerical Matchete one-loop coefficient.
No full physical EFT matching, individual X/Y thresholds or physical
prediction certification is claimed. No production files are modified.

## BC-MIXED-THRESHOLD-10: running scalar thresholds, same one-scale source (2026-10-31)

**User-local validation preceding this package:**
`scripts.RunT3BCMassEigenYukawaC5Handoff --all --run-id bc_yukawa_c5_v1`
completed for B(−1), C(−1). The repeated run
`--run-id bc_yukawa_c5_replay_v1 --recompute-c5` returned `C5 replay=True`
for both, reproducing the saved conditional one-scale Takagi spectra.
Physical sequential EFT matching remained `NOT CERTIFIED`.

**This package:** new read-only `scripts/RunT3BCSeparatedScalarThresholds.py`
reintegrates the *same* `output/general_t3_quadratic_eft1` benchmark, not the
differently initialized `continuous_v1` candidate. It checks a stored beta
source SHA, saved full scalar endpoint, native Matchete X/Y Yukawa phases,
one-scale benchmark SHA and mass spectrum. Then it solves the running
Y and X mass-eigenvalue crossing equations, exporting exact numerical mixing
matrices, scalar flavour projectors and native Yukawa weights at those scales.

**Strict limit:** Y threshold lies inside the two-scalar EFT1 and is a genuine
candidate running crossing. The X solution below Y continues that un-decoupled
EFT1, so its scale is **provisional** until a matched, X-only EFT2 is available.
No Y/X Wilson boundary was matched, no finite scalar threshold computed, no
new neutrino mass predicted, and no production/remote GitHub file modified.
The new script is prepared; the user's local numerical run is pending.

## Experimental B/C EFT2 sequential full numerical continuation (conditional, 2026-10-31)

Added local-only files `Numerical/experimental/T3BCEFT2HX.py`,
`scripts/RunT3BCEFT2Sequential.py`, and standalone regression tests.

The stage consumes existing exact Y-threshold projection outputs and replays
`continuous_v1` EFT1 scalar RG data with unchanged SHA256 sources. It rotates
and fits the complete H,X quartic tensor to the six-real-parameter Z2 inert
doublet basis, evolves the **Y-decoupled** renormalisable one-loop EFT2 and
project seven-term `psi2phi2` Wilson tensor, locates the corrected running X
threshold, and runs the conditional Weinberg coefficient to 100 GeV.

**Explicit limitations:** finite Y and X scalar matching not derived;
additional-quartic physical convention flags remain false; the initial Higgs
quadratic parameter is a tree-level proxy; external three-family flavour
indices are retained but the existing gauge-index master uses projected SM
Yukawa anomalous dimensions. The output is conditional; no physical neutrino
mass prediction is certified.

Standalone local algebra tests passed (14/14) at package creation. The
end-to-end Windows numerical run and the native Wolfram finite matching
have **not** been performed in the packaging environment. Production sources
and old outputs have not been modified.

