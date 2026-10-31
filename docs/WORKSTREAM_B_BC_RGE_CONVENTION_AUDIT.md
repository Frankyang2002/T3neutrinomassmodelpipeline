# Workstream B — B(−1)/C(−1) convention and mass-RGE audit

**Status:** scoped source audit and independent algebraic checks. **NOT a physical Matchete ↔ RGBeta certificate.**

Source inspected read-only from `Frankyang2002/T3neutrinomassmodelpipeline`, default `main`, head `305035e9df8b90e34fc74c7a3f57453ea6c6aae3c` at the start of this audit. This delivery adds only new files. The actual Windows checkout was not accessible here; the user must run the source-input audit locally. No Mathematica or Matchete invocation is claimed.

## 1. Known parameter dictionary in the **declared** source basis

For the two models, `(dS1,dS2,dF,alpha)` is B = `(2,2,1,-1)`, C = `(2,2,3,-1)`. RGBeta's complex doublets have `Y(S1)=-1/2`, `Y(S2)=+1/2`; `Y(F)=0`. `H` has `Y=+1/2`. The B/C Majorana fermion has a single left-Weyl `F` species with three flavours.

The current RGBeta scalar **candidate** potential contains

\[
V_2=m_{S1}^2 S_1^\dagger S_1+m_{S2}^2S_2^\dagger S_2+
\big[\mu_R\,S_{1i}^*\varepsilon_{ij}S_{2j}^*+\mathrm{h.c.}\big],
\qquad \varepsilon=\begin{pmatrix}0&1\\-1&0\end{pmatrix}.
\]

Here `mu_R` = `muT3GeneralS1S2` in `RunGeneralT3QuadraticRGBeta.wl`; `mS1Sq` and `mS2Sq` are the diagonal couplings. The native experimental Hessian code uses this **no-extra-1/2** convention. The Matchete extra-interaction source inserts `+PlusHc[contribution]` in the interaction **Lagrangian**. Consequently, `mu_M = -mu_R` is a *conditional same-orientation* hypothesis only; a permutation of the epsilon indices or selection of the hermitian-conjugate monomial changes the corresponding sign or conjugation. There has been no direct normalized, live Matchete-vs-RGBeta mass-vertex equality. **Do not map `muT3GeneralS1S2` onto a Matchete or numerical mass input without that test.**

The numerical mass-basis transport uses

\[
\Phi=(S_1,\varepsilon S_2^*)^T,
\quad A=\begin{pmatrix}m_{S1}^2&\mu_R\\\mu_R^*&m_{S2}^2\end{pmatrix},
\quad U^\dagger AU=\operatorname{diag}(m_-^2,m_+^2),
\quad\Phi=U\,(X,Y)^T.
\]

The analytic eigenvalues are `(mS1Sq+mS2Sq)/2 ± sqrt(((mS1Sq-mS2Sq)/2)^2+abs(mu_R)^2)`; each is fourfold degenerate in the 8-real-component gauge-unbroken Hessian. The spectrum depends on `abs(mu_R)`, so it **cannot** detect a conjugation error. A complex off-diagonal Hessian or individual vertex action is required. Scalar and Yukawa tensors require simultaneous transformation, including the `-epsilon` in `S2=-epsilon*(U[1,0] X+U[1,1] Y)^*`.

### Seven additional complex quartics: conditional Stage23/29 map

The table is a statement of current code plus saved candidate algebra, **not** a certified physical mapping. `M` denotes a Matchete coupling of `+PlusHc[lambda_M * operator]`, and `R` is the RGBeta potential coupling.

| Matchete slot | RGBeta name | Candidate relation |
|---|---|---|
| `lambdaT3ExtraQ000112I1` | `lambdaS1barS2S2bar2` | `R=-M` |
| `lambdaT3ExtraQ000202I1` | `lambdaS1bar2S2bar2` | `R=-M` |
| `lambdaT3ExtraQ001201I1` | `lambdaS1S1bar2S2bar` | `R=-M` |
| `lambdaT3ExtraQ020020I1` | `lambdaHHdagS2S2` | `R=-M` (stronger separate two-active-quartic evidence) |
| `lambdaT3ExtraQ020200I1` | `lambdaHHdagS1barS1bar` | `R=-M` (stronger separate two-active-quartic evidence) |
| `lambdaT3ExtraQ110101I1`, `...I2` | `lambdaT3GeneralQ110101I1`, `lambdaHHdagS1barS2barCross` | **Mixed 2×2 map** below |

For `110101`, ordered species are `(H,H†,S1†,S2†)`; the saved native tensor identity is `R_j = Σ_i A_ij M_i`, with

\[
A=\begin{pmatrix}1/2&1\\-\sqrt3/2&0\end{pmatrix},\quad
\begin{pmatrix}g_1\\g_2\end{pmatrix}
=-A^{-1}\begin{pmatrix}h_1\\h_2\end{pmatrix}
=\begin{pmatrix}0&2\sqrt3/3\\-1&-\sqrt3/3\end{pmatrix}
\begin{pmatrix}h_1\\h_2\end{pmatrix}.
\]

Thus `h1=-g1/2-g2`, `h2=sqrt(3)*g1/2`. **Neither individual map is a generic `-1` mapping.** The `-` sign comes from comparing RGBeta scalar potential with Matchete interaction Lagrangian, assuming the same component basis and orientation. Field rephasings, tensor normalizations, and the direct simultaneous native action still require certification. The separate `020020`/`020200` comparisons do not establish the remaining five slots.

Native RGBeta `AddQuartic` uses `6*Matrix` for four complex scalar legs in real components. The generic real-field action has `1/4!`; `6/24` per ordered real-field term is **not** an additional `1/4` multiplication of the physical complex coupling. Do not apply that erroneous rescaling. A complete matrix test must include the 24 complex-leg permutations and both conjugate terms.

## 2. Directly inspectable EFT1 one-loop quadratic results

The saved UV/EFT1 candidate exports `output/general_t3_quadratic_eft1/{B_m1,C_m1}_{uv,eft1}_quadratic_live_rgbeta.json` have all three mass beta symbols. Let `16π² dX/dln μ_RG = β_X` and abbreviate `q1=lambdaS1S1bar2S2bar`, `q2=lambdaS1barS2S2bar2`, `q3=lambdaS1bar2S2bar2`. The **EFT1** exported equations are:

\[
\begin{aligned}
\beta_{m_1^2}&=\big(-\tfrac92g_2^2-\tfrac32 g_Y^2+6\lambda_{S1}\big)m_1^2+4\lambda_{12}m_2^2+6(\mu_Rq_1^*+\mu_R^*q_1),\\
\beta_{m_2^2}&=\big(-\tfrac92g_2^2-\tfrac32 g_Y^2+6\lambda_{S2}\big)m_2^2+4\lambda_{12}m_1^2+6(\mu_Rq_2^*+\mu_R^*q_2),\\
\beta_{\mu_R}&=6q_1m_1^2+6q_2m_2^2+\big(-\tfrac92g_2^2-\tfrac32g_Y^2+2\lambda_{12}-\tfrac32\lambda_{12\mathrm{Adj}}\big)\mu_R+12q_3\mu_R^*.
\end{aligned}
\]

**Key closure consequence:** `mu_R=0` alone is **not** RG-stable when `q1*m1^2+q2*m2^2` is nonzero. Switching off only the mixed mass but retaining all extra quartics is not a closed subspace. In the UV there are additional `y1/y2` traces: coefficient `1` in B and `3/2` in C for the `mu_R` wave-function term. Other mass equations include Yukawa/MF traces. The audit only parses the simpler EFT1 mass polynomials; it checks UV **presence**, not complete UV flavour covariance.

A further verified source gap: the base `T3RGBetaModel.wl` defines the `H` scalar but does **not** register its Higgs quadratic mass `H†H` or export a corresponding beta. Thus the existing mass RGE system is not a complete renormalizable scalar-mass sector, regardless of whether its registered B/C equations evaluate. Adding this equation must be coordinated with the other workstreams rather than altering production files here.

## 3. New independent checks

`scripts/WorkstreamBBCConventionAudit.py` performs:

1. **Fail-closed source/basis inventory:** compare operator-manifest and RGBeta-registry counts, all seven candidate B/C coupling symbol names, exact Stage23/29 matrices, Matchete `PlusHc` source and RGBeta field-word declarations. It checks `A C = -I` exactly (with independent complex coefficients) and rejects silently changed certification flags.
2. **Independent complex mass Hessian:** finite differences of the potential in eight real components against a separately assembled real Hessian, full eigenvalue multiplicities, rephasing, and checks that conjugating or halving the bilinear changes the component Hessian.
3. **Actual saved EFT1 mass betas:** parse the current JSON beta strings with a restricted symbolic grammar; verify real diagonal mass beta functions, complex `U(1)×U(1)` field-rephasing covariance of the coupled system, exact nonzero `mu_R=0` source, and first-order eigenvalue evolution via independent finite differences of the `2×2` mass matrix. Compare B/C EFT1 scalar equations; UV beta keys are checked separately.
4. **Coverage gate:** report missing RGBeta Higgs quadratic beta, mark final A basis certificate absent, and explicitly set physical quartic/mass matching, sequential matching, and full 3-generation Wilson certification flags **false**.

The rephasing test uses `S1→exp(i a) S1`, `S2→exp(i b) S2`, under which `mu_R,q1,q2→exp(i(a+b))` and `q3→exp(2i(a+b))` in the compensating coupling transformation. This gives a useful independent test of the correct conjugations in the three saved beta functions. It cannot prove Matchete's inter-package field ordering.

## 4. Install and run on the local Windows checkout

The ZIP contains **complete new files**, no patches and no replacements of another workstream's files:

- `scripts/WorkstreamBBCConventionAudit.py`
- `tests/test_workstream_b_bc_convention_audit.py`
- `docs/WORKSTREAM_B_BC_RGE_CONVENTION_AUDIT.md`

```powershell
Set-Location 'C:\thesistemp\T3neutrinomassmodelpipeline'
Expand-Archive -LiteralPath 'C:\path\to\Workstream_B_BC_RGE_Audit.zip' -DestinationPath . -Force
python -X utf8 -m pytest -q tests\test_workstream_b_bc_convention_audit.py
python -X utf8 -m scripts.WorkstreamBBCConventionAudit --self-test
python -X utf8 -m scripts.WorkstreamBBCConventionAudit --output output\audits\workstream_b_bc_conventions.json
```

The final command **reads** repository inputs and writes only `output/audits/workstream_b_bc_conventions.json` as evidence; it does **not** evolve running couplings or change physical parameters. The source-input audit must be run on your exact checkout before treating the audit as passed there. If current saved inputs differ, it fails closed and must be reviewed (not bypassed).

## 5. Explicit dependency and unresolved gates

**Workstream A:** the existing `docs/general_scalar_trial_operator_manifest.json` is checked as an inventory of 16 models, but its final independently certified operator basis has not been supplied. Full operator completeness, all quartic-block normalization, quartic beta closure and the 14-model generalization are **conditional on A**. No Stage29 status label overrides this.

**Workstream B next:** compare the *installed*, normalized RGBeta mixed quadratic vertex against Matchete `PlusHc` using identical complex components, epsilon orientation, and both conjugate terms; then check `U† β_A U` with the actual evolved UV/EFT1 RG state and map the non-diagonal beta to the numerical `X,Y` mass eigenbasis. Check the missing Higgs quadratic RGE and any other missing couplings against A before a claim of complete renormalizable closure. Do not alter mass matching owned by other workstreams.

**Wilson operators:** the B/C implementation evolves six independent symmetric physical flavour pairs but uses the one-generation SM Yukawa/gauge-index projection in its seven-term Wilson beta. A full three-generation calculation requires explicit 3×3 charged-lepton Yukawa matrices and flavour-space wavefunction terms acting on both lepton indices, including transpose/conjugation and mixing between pairs, plus rederivation/verification of the complete anomalous-dimension tensor. This is *separate* from renormalizable UV/EFT1 beta closure.
