# Workstream 4 — Native source-amplitude contract (conditional batch)

## Status and scope

This is an **additive, read-only** source-to-component audit for the existing frozen T3 benchmark. It does **not** recalculate the existing 16 direct one-loop neutrino mass matrices or the tree-level Gaussian LLSS tensors. It does **not** claim physically certified Matchete/RGBeta-to-Weyl Feynman rules. Its generated report is explicitly conditional.

### Source-verified structural dictionary

The current public repository source (read-only inspection 2026-11-03) has:

| Source | Original source fields | Restricted two-Weyl construction |
|---|---|---|
| `Lagrangian/interactions/T3Topology.wl` | `y1[p,r] NCM[Bar[l],CConj[F]] S1 CG[T3Y1CG] + h.c.` | The *Hermitian-conjugate* vertex is assigned to a conjugated scalar `S1*` and the charge-conjugate heavy Weyl sector |
| same | `y2[p,r] NCM[Bar[l],F] Bar[S2] CG[T3Y2CG] + h.c.` | Its Hermitian conjugate is assigned to an unconjugated scalar `S2` and the other heavy Weyl sector |
| `RGE/running/rgbeta/T3RGBetaModel.wl` | `AddYukawa[y1,{Bar[S1],l,F}]` | RGBeta left-Weyl field `F` |
| same | `AddYukawa[y2,{S2,l,Fc}]` for Dirac, `{S2,l,F}` for Majorana | RGBeta left-Weyl field `Fc` or `F` |
| same | `AddFermionMass[MF,{F,Fc}]` for Dirac, `{F,F}` for Majorana | Unbarred doublet pair uses `epsilon`; triplet pair uses `delta`; singlet pair is scalar 1 |
| `Lagrangian/model/SU2Invariants.wl` | `T3Y1CG`: `{True,True,False}`; `T3Y2CG`: `{True,False,True}` | Distinct gauge-conjugation signatures, not interchangeable |
| same | `T3MixCG`: `{False,False,False,True}` | Symmetric two-Higgs contraction |
| `scripts/BuildT3DirectOneLoopInputs.py` | `y1_HC` recorded as `conjugate`; `y2_HC` recorded as `field` | The existing *candidate* `Y,Z` adapter, not independently native-normalised |

The `y1`, `y2`, and `lambdaT3` couplings are complex by the native Matchete source declarations. For neutral heavy singlet/triplet models `B_-1` and `C_-1`, the fermion is Majorana. The 14 other benchmarks have Dirac heavy fermions.

### Independent A_+0 doublet mass insertion

In the **unbarred RGBeta fundamental axes**, choose the explicit conventional basis

\[
\epsilon=\begin{pmatrix}0&1\\-1&0\end{pmatrix},\quad
\psi^T=(F_0,F_1,F^c_0,F^c_1),\quad
\mathcal L_m=-\frac12\psi^TM\psi+\mathrm{h.c.}
\]

for a single heavy flavour. The source-allowed invariant gives the **conditional** component matrix

\[
M=M_F\begin{pmatrix}0_2&\epsilon\\\epsilon^T&0_2\end{pmatrix},\quad
M^{-1}=\frac{M}{M_F^2},\quad
(M^{-1})_{F_0,F_1^c}=+1/M_F,\quad
(M^{-1})_{F_1,F_0^c}=-1/M_F.
\]

`epsilon` satisfies the fundamental SU(2) generator Ward identity whereas `delta` does not for two unbarred fundamentals. This is an **explicitly chosen mathematical epsilon basis**; native RGBeta component axes and the Matchete charge-conjugation phase have **not** been measured. The existing numerical input builder changes to a conjugate-index basis whose Dirac off-diagonal block is the identity, so comparing these two matrices directly without the basis transformation would be invalid.

### Independent B_-1 singlet Majorana insertion

For a canonically normalised single left-Weyl field:

\[
\mathcal L\supset -\frac12M_F FF-F(J_1+J_2)+\mathrm{h.c.},
\quad
\mathcal L_{\mathrm{eff}}=\frac{(J_1+J_2)^2}{2M_F}+\mathrm{h.c.}
\]

The coefficient of the **distinct** source product `J1 J2` is `1/MF`, not `2/MF`. The native package's signed `MF` map remains unproven, so this is not a certified physical factor for the B_-1 Lagrangian. The same distinction applies to triplet Majorana `C_-1`, with nontrivial real-representation charge-pairing metrics.

### Direct UV result and frozen completeness

The existing UV engine takes left-Weyl vertices `Y,Z`, a symmetric heavy mass matrix, and Hermitian/symmetric complex scalar Hessians. In its declared convention the one-loop kernel is `B0bar(0; MF^2, mS^2)` with full mass eigenstates, no heavy-mass expansion. Existing source inputs, matrices, decompositions, and tests are **not replaced** by this batch.

Six models allow an additional independent scalar-conjugate Yukawa: `A_-2`, `A_+0`, `D_-2`, `D_+0`, `E_-2`, `E_+0`. They are **assumed zero** in the restricted frozen benchmark, not forbidden by symmetry. Nonminimal scalar quadratic and quartic terms must likewise be listed as benchmark zeros when absent; the restricted input builder explicitly guards several such terms, including `lambdaHHdagS2S2`, `lambdaHHdagS1barS1bar`, `lambdaHHdagS1barS2barCross`, and adjoint portals. This does not constitute an exhaustive classification of all general potential invariants.

### Precise outstanding native evidence

1. A **native installed Matchete** numerical amplitude in a declared common charge-conjugation and left-Weyl basis, including `GammaCC`, `Bar`, `CConj`, `PlusHc`, both `y1` and `y2`, and the heavy Dirac/Majorana mass term.
2. A **native RGBeta** numerical fundamental/triplet tensor and epsilon orientation in the *same* component basis, plus a derivation distinguishing its real-scalar RGE Yukawa prefactor from four-component Feynman-rule coefficients.
3. Full potential-to-real-Hessian projection for `lambdaT3`, `lambdaH1`, `lambdaH2` (including any relative minus and `1/2` factors), independently of the provisional input builder.
4. Source-normalised conjugate-scalar Yukawas if the general, rather than restricted, model is requested.

Until those inputs are established, the signed UV mass matrices, their phases, and absolute normalisations remain **provisional**. Passing this contract or the previous 67 local tests is not physical certification.

### EFT boundary

At `mu0 = 1e7 GeV`, the formal all-heavy-fermion Gaussian bridge is `+J^T MF^-1 J/2`. The **first actual threshold** must integrate out **all BSM fermions simultaneously** at `muF = 1e5 GeV` while retaining *both* BSM scalar multiplets and SM fermions. The second threshold integrates out *both* BSM scalars at `muS = 1e3 GeV`, leaving SMEFT. The Gaussian bridge at `mu0` is not a scale-evolved EFT1 matching coefficient at `muF`; hard one-loop fermion matching has not been implemented in this batch.

### Installation and execution (PowerShell)

From `C:\thesistemp\T3neutrinomassmodelpipeline`, extract the supplied ZIP into the repository root **after reviewing the three new code files** (no existing engine files are replaced). Then:

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
.\.venv\Scripts\python.exe -m pytest tests\test_lo_direct_uv_native_amplitude_contract.py -q
.\.venv\Scripts\python.exe scripts\RunLODirectUVNativeAmplitudeContract.py `
    --project-root . `
    --check-frozen-inputs `
    --output-dir output\general_observable\lo_certification_direct_uv\native_amplitude_contract_01
```

The generated JSON is at

`output/general_observable/lo_certification_direct_uv/native_amplitude_contract_01/native_amplitude_source_contract.json`.

The report fails closed on changed source forms, missing frozen inputs, wrong heavy-field blocks or scalar orientations. It does **not** execute Wolfram, Matchete, or RGBeta, and does not compare UV and EFT masses.
