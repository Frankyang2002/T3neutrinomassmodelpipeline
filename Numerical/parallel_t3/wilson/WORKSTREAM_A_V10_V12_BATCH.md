# Workstream A v10–v12 — broad native quartic projection and scalar-action batch

**Status:** Code implemented and standalone Python math/tests executed. **The live
Wolfram RGBeta v10 witness has NOT been executed by the author** because this
container has no `wolframscript`. Full native results are therefore **pending
on the user's local installed RGBeta 1.2.0**, not assumed to pass.

Source reviewed read-only at GitHub HEAD
`2bc2e4b516776212645370105cc63495a3c47348` (2026-11-09),
plus RGBeta's public `Package/FieldsAndCouplings.m`,
`Package/GroupsAndIndices.m`, and `Package/TensorCalculations.m`.
The model source's reviewed Git blob is
`bc9500803f2324ad6504e71d769cb75fbbde9b07`.

## Why this batch is bigger than v9

The actual T3RGBetaModel.wl **EFT1** (`T3RGBetaBuildEFT1`, no heavy F)
registers **178 model-coupling quartic slots** across all 16 choices:

| Group | Model-coupling slots | Check |
|---|---:|---|
| Six real norm-type quartics (`lambdaH`, `lambdaS1`, `lambdaS2`, `lambdaH1`, `lambdaH2`, `lambda12`) | 96 | native v10 diagonal projector attempted |
| Higgs-adjoint portals `lambdaH{1,2}Adj` | 26 | native v10 diagonal + full 1×1/2×2 response attempted |
| Triplet self-adjoint quartics (`lambdaS1Adj`, `lambdaS2Adj`) | 14 | native v10 diagonal + v11 rank/Fierz |
| Mixed scalar adjoint `lambda12Adj` | 10 | native v10 diagonal + v11 rank |
| Mixed triplet cross `lambda12Cross` | 4 | native v10 diagonal + v11 rank |
| Topology complex quartic `lambdaT3` | 16 | v12 constructed CG/H.c. real action |
| Six complex special B(-1)/C(-1) quartics, two models | 12 | v11 constructed H.c. real actions |
| **Total** | **178** | **150 real, 28 complex** |

The 58 Higgs-portal slots of v9 are **32 norm portals (already in the 96)**
plus **26 additional adjoint portals**: do not double count.

## v10 — actual native RGBeta symbolic Lam and projectors

`scripts/RunT3ParallelANativeQuarticProjectorsV10.wl` loads the installed
RGBeta and the current repository's `T3RGBetaModel.wl`, builds each EFT1
model with all optional T3 experimental extensions disabled, and reads actual
package `$quartics` metadata and the stored `Projector` functions.

For a native quartic source coupling \(\lambda_j\) and a native projector
\(P_i\), it calculates, using the installed package functions:

\[
R_{ij}=\left.P_i[\mathrm{Lam}(\lambda_j)]\right|_{\lambda_j=1},
\qquad
P_i[\mathrm{Lam}(\lambda_j)]\big|_{\lambda_j=0}=0.
\]

The source \(\lambda_j\) is isolated using a Wolfram `Block` localizing
`$quartics` to that registered interaction **without editing RGBeta source**.
The package's own `RefineGroupStructures` reduces the remaining gauge
contractions. Linearity is checked at \(\lambda_j=0,1,2\).

- **150** real quartic native projector diagonals should give one;
- **32** Higgs-portal projector sectors produce **110** matrix elements:
  58 diagonals and 52 off-diagonals. The full matrices are tested for
  invertibility; **do not assume** different native projectors are orthogonal.
- The **28 complex couplings** are inventoried, but **their native conjugate
  channels are NOT certified** by this v10 exporter.
- `SMOKE` mode runs A(-4), B(-1), E(-2), before all16 native execution.
- A diagnostic JSON is written even if numerical/symbolic projection fails.
  The script exits nonzero on any incomplete native proof. It refuses
  overwrite and guards local Python/Wolfram source SHA256 plus the installed
  RGBeta entrypoint file SHA256.

**Native symbolic projection is not native numerical component Lam.**
The v10 Python auditor re-executes the entire v9 source-guarded chain and
builds the observed native-projector-to-Matchete conditional matrices
\(R M\), with \(M\) the 58-slot v9 conditional coupling transformation.
If native projectors are non-orthogonal, \(R^{-1}\) is recorded explicitly.
No Matchete-to-RGBeta numerical triplet frame is extracted or assumed.

## v11 — nonportal ranks and special complex quartics

`scalar_nonportal_rank_v11.py` uses a **declared** real Cartesian triplet
basis \((J^a)_{bc}=-i\epsilon_{abc}\). It checks SU(2) Lie algebra,
\(C_2(3)=2\), trace index \(T(3)=2\), independent rank of all registered
triplet self and mixed-scalar quartic sectors (14, 10, 4 respectively),
and SU(2) covariance for the real nonportal mixed-scalar bilinears.

For the triplet self-adjoint quartic,
\[
\sum_a (s^\dagger J^a s)^2
=(s^\dagger s)^2-|s^T s|^2.
\]

In a real Cartesian triplet basis, the source `lambda12Cross` tensor
produces
\[
|s_1^T s_2|^2+|s_1^\dagger s_2|^2,
\]
which is independently ranked against norm and adjoint contractions.

For each of the **six** source special B/C(-1) `SelfConjugate -> False`
quartics, the exact `GroupInvariant` delta/epsilon polynomial is transcribed
into a separate mathematical tensor, and checked at random real component
fields with complex coefficients \(c\):
\[
\frac1{4!}\,\Lambda_{abcd}x_ax_bx_cx_d
=c\,\mathcal I(z)+c^*\mathcal I(z)^*,\quad
z=(x_R+i x_I)/\sqrt2,
\]
using the source's `6*Matrix` factor and full real symmetrization. This is
**not** native RGBeta `Lam`, or a native Matchete component match.

## v12 — 16 topology `lambdaT3` complex actions

`scalar_topology_action_v12.py` derives the mathematical normalized
`Sym²(2) -> 3` CG, checks its SU(2) covariance, and constructs
`T3RGBetaS2TripleInvariant[a,b,c]` directly using the **three**
fundamental `eps` factors present in the model source. In the declared real
Cartesian triplet frame it verifies
\[
T_{abc}=-\frac{i}{\sqrt2}\epsilon_{abc}.
\]

All **16** `lambdaT3 H H S1 Bar[S2]` complex model-coupling actions are
reconstructed with `+ H.c.` and `Lam/4!` using complex random coefficients.
The B/C models' `Bar[S2]` **pseudoreal doublet index intertwiner is not
source-certified**. Their mathematical construction is therefore not a
physical SU(2) transformation-of-components certificate.

v11 and v12 refuse to reuse transcribed source CG formulas if the local
`T3RGBetaModel.wl` bytes no longer match the reviewed Git source blob
(after normalizing Windows CRLF to LF). A changed source must be reviewed
and the code regenerated, not silently accepted.

## Files and install paths

Unzip into **`C:\thesistemp\T3neutrinomassmodelpipeline`**. The ZIP provides
these **12 complete new files**, preserving each path:

1. `Numerical/parallel_t3/wilson/native_quartic_projectors_v10.py`
2. `Numerical/parallel_t3/wilson/scalar_nonportal_rank_v11.py`
3. `Numerical/parallel_t3/wilson/scalar_topology_action_v12.py`
4. `Numerical/parallel_t3/wilson/WORKSTREAM_A_V10_V12_BATCH.md`
5. `scripts/RunT3ParallelANativeQuarticProjectorsV10.wl`
6. `scripts/RunT3ParallelANativeQuarticProjectorsAuditV10.py`
7. `scripts/RunT3ParallelAScalarNonportalV11.py`
8. `scripts/RunT3ParallelATopologyScalarV12.py`
9. `scripts/RunT3ParallelAScalarBatchV10V12.ps1`
10. `tests/test_t3_parallel_a_native_quartic_projectors_v10.py`
11. `tests/test_t3_parallel_a_scalar_nonportal_v11.py`
12. `tests/test_t3_parallel_a_scalar_topology_v12.py`

(Do not delete prior v5–v9 files.)

## Single PowerShell batch

```powershell
Set-Location C:\thesistemp\T3neutrinomassmodelpipeline
$zip = "$env:USERPROFILE\Downloads\t3_workstream_A_scalar_batch_v10_v12.zip"
Expand-Archive -LiteralPath $zip -DestinationPath . -Force
powershell -ExecutionPolicy Bypass -File scripts\RunT3ParallelAScalarBatchV10V12.ps1
```

Run only the independently checked Python tasks when Wolfram is unavailable:

```powershell
powershell -ExecutionPolicy Bypass -File scripts\RunT3ParallelAScalarBatchV10V12.ps1 -SkipNative
```

The full run writes only **new** diagnostics under
`output/parallel_t3/wilson/`:

- `nonportal_scalar_invariants_conditional_v11.json`
- `topology_scalar_action_conditional_v12.json`
- `rgbeta_quartic_projectors_native_v10_smoke.json`
- `rgbeta_quartic_projectors_native_v10.json`
- `native_quartic_projectors_conditional_v10.json`

Every output refuses overwrite. To retry after a blocked diagnostic, **inspect
its recorded error first**; use fresh output paths only after identifying the
problem. No source files, production RGE, physical Wilson trajectories, or
GitHub repository are modified by this batch.

## Remaining independent blockers

1. **28 native complex-coupling `Lam` conjugate-channel projectors**: v11/v12
   reconstruct actions algebraically, but do not run their native Lam terms.
2. **Complete nonportal Matchete↔RGBeta coupling dictionary**: v11 rank is
   not a native Matchete CG extraction and v12 uses declared CGs.
3. **B(-1)/C(-1) mixed scalar mass diagonalization** and neutral self-mass
   phases; those models can be included in symbolic audits only.
4. **Signed GammaCC-to-Weyl LLSS normalization** and **complete physical
   one-loop ADM** (including crossed Yukawa correction).
5. **Certified hard fermion matching, scalar threshold at one joint boundary,
   SMEFT evolution and neutrino mass results.** Keep the required EFT order:
   `UV -> integrate out ALL fermions -> integrate out ALL scalars -> SMEFT`.

**None of the above is a physical Wilson-evolution authorization.**
