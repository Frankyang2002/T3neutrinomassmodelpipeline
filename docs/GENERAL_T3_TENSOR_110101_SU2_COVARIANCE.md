# General T3: native-source SU(2) covariance of `110101` (B_-1 and C_-1)

## Scope and evidence

Read-only test of the **actual Wolfram-exported RGBeta `Invariant` functions** for
`lambdaT3GeneralQ110101I1` and the reused
`lambdaHHdagS1barS2barCross`, compared to uploaded Matchete component tensors.
The previous 50-slot native exports are still required. No production changes;
Stage 22 physical normalisation is **UNCERTIFIED**.

Source checked: `RGE/general/generated/B_-1_extra_quartics.wl` registers the
new `lambdaT3GeneralQ110101I1`; `RGE/running/rgbeta/T3RGBetaModel.wl`
registers the reused `lambdaHHdagS1barS2barCross` (in
`T3RGBetaAddRGClosedQuartics`). B and C have the same relevant field
content and native invariant expressions. The Python test fails closed if
those exported expressions have changed, rather than assuming an old diagnosis.

## Gauge covariance calculation (exact, not fitted)

Native fields are ordered `(H, Bar[H], Bar[S1], Bar[S2])`. Under the **explicit
standard component convention**

\[
  H \longmapsto U H,\qquad
  \bar H,\bar S_1,\bar S_2 \longmapsto U^*\bar H,U^*\bar S_1,U^*\bar S_2,\qquad U\in SU(2),
\]

use `epsilon[1,2]=+1`, and generators `T_a=sigma_a/2`. For an invariant
coefficient tensor `Q_{i j k l}`, exact infinitesimal covariance requires

\[
 \sum_{p=1}^4\sum_b Q_{i_1\ldots b\ldots i_4}
              (G_p)_{b i_p}=0
 \quad\text{for every }(i_1,\ldots,i_4)\in\{1,2\}^4,
\]

where `G_1=T_a` and `G_2=G_3=G_4=-T_a^*`. The three Lie-algebra generators
are a complete test of continuous SU(2) covariance in this representation.
All calculations use exact SymPy arithmetic: no floating-point thresholds.

| Invariant | `sigma1/2` residual components | `sigma2/2` | `sigma3/2` | Covariant? |
|---|---:|---:|---:|---|
| Native `delta(4,1) epsilon(2,3)` | 0 | 0 | 0 | Yes |
| Native reused `epsilon(1,3) epsilon(2,4)` | 4 | 0 | 4 | **No** |
| Proposed `delta(3,1) epsilon(2,4)` | 0 | 0 | 0 | Yes |

Example exact violation for the reused invariant: the `sigma1/2` Lie variation
at `(i1,i2,i3,i4)=(1,1,1,2)` equals **-1**. The identity is not rescued by
an overall multiplicative factor.

**Important qualification:** RGBeta's internal convention for conjugate-field
indices has not yet been independently linked to its actual gauge-action
implementation. Thus this is an exact test of the *stated standard* SU(2)
representation convention applied to actual native invariant expressions,
not proof of an RGBeta package bug under every possible internal convention.

## Exact provisional Matchete comparison

Let `M1,M2` denote the two exported native Matchete `110101` tensors in order.
Use the source-verified native first RGBeta tensor `R1` and the **candidate,
not registered** replacement `R2* = delta(3,1) epsilon(2,4)`. For both B_-1
and C_-1:

\[
  [R_1\quad R_2^*]=[M_1\quad M_2]A,\quad
  A=\begin{pmatrix}\tfrac12&1\\-\tfrac{\sqrt3}{2}&0\end{pmatrix},\quad
  \det A=\tfrac{\sqrt3}{2},\quad
  A^{-1}=\begin{pmatrix}0&-2/\sqrt3\\1&1/\sqrt3\end{pmatrix}.
\]

Rank is 2 and the **exact component residual is zero** for each model.
With the original reused source invariant, the first RGBeta column lies in
the Matchete span, while the reused second column does not. Therefore the
proposed invariant simultaneously resolves the two earlier `110101` failures
within the stated component convention.

This is **geometry**, NOT a coupling map: the native coupling functions
`6*Matrix[coupling][##1]&` do not by themselves establish the signed scalar
potential prefactor, conjugation orientation or a globally consistent field
basis. Never send these matrices to production before those checks.

## Local installation and execution

Overlay the ZIP onto your existing repo root; do not push changes to GitHub:

```powershell
Expand-Archive -LiteralPath .\T3_General_Scalar_Stage22_110101_Covariance.zip -DestinationPath . -Force
python -m pytest -q tests\python\test_general_t3_tensor_110101_covariance.py
python -m scripts.RGBeta110101Covariance
```

The Python checker reuses the existing native evidence in
`output/general_t3_stage22/`. Output is
`output/general_t3_stage22/stage22_110101_su2_covariance.json`.

For an **additional live-Wolfram RGBeta check**, run on your Windows machine:

```powershell
wolframscript -file .\scripts\ExportGeneralT3SU2NativeCovariance.wl `
  output/general_t3_stage22/stage22_110101_wolfram_native_covariance.json
```

The Wolfram script rebuilds B_-1 and C_-1, reads their *registered* quartics
through `RGBetaNativeAccess.wl`, converts symbolic `del` and `eps`
to component evaluators **before** numerical index substitution (avoiding
RGBeta's contracted-delta trace rule), and applies all three Lie-algebra
generators. This independent
runtime script has **not** been executed here because Wolfram is unavailable.
If its results differ from the Python diagnostic, treat that as an unresolved
convention mismatch and do not modify the source until explained.

## Scope estimate after Stage 22

Approximately **3-4 substantive stages** for the *rough numerical benchmark*
goal, contingent on no further blockers:

1. Verify the source defect with live RGBeta; decide whether/how to repair the reused invariant **locally**.
2. Independently derive the signed scalar potential prefactors and 50-slot physical coupling map.
3. Integrate and smoke-test the map with the existing RGE/threshold/benchmark workflow.
4. If required, address new numerical blockers to obtain representative outputs.

The eight previously blocked complex quadratic operators may be separate
follow-up scope. This is a planning estimate, not a verified numbered roadmap
and not a forecast of effort or completion date.
