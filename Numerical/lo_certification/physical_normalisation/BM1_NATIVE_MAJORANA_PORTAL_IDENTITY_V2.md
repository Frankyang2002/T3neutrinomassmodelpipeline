# B_-1 native Majorana scalar-quadratic gate: V2 correction (2026-11-03)

## Status and provenance

- **Observed on user's installed Matchete:** `RunPhysicalBm1NativeScalarQuadratic.wl` built the restricted B_-1 Lagrangian; `CheckLagrangian` returned `True` for both Yukawas, two S1 portals, two S2 portals, both scalar self sectors, cross portals and the T3 scalar mixing. All original coefficient/mass/charge assertions passed. The script **blocked** only `portal1_identity` and `portal2_identity`; the correct non-identity Inv2 H0 projector assertions passed.
- **Inspected read-only on GitHub:** the currently indexed B_-1 `native_higgs_portal_registry.json` at `output/general_observable/lo_direct_one_loop/source_portal_cg_all16/T3_dS1_2_dS2_2_dF_1_alpha_m1/`, and unchanged `Lagrangian/interactions/ScalarPotential.wl`. Latest inspected GitHub head: `2ac58cc3de12c7af11e2e7e214415dadf0138f4f`. The public native tensor registry has the same tensor objects for both doublets. Do not treat this read as a new local extraction.
- **New independent derivation (this V2):** exact four-index decomposition of the identity in the two-dimensional B_-1 native Matchete Higgs-portal invariant basis. Verified using symbolic exact arithmetic on the full native component arrays rather than the benchmark mass matrix or provisional loop tensors. The scripts in this delivery have been syntax-checked (Python) and smoke-tested using a *source-shaped synthetic fixture*; Wolfram V2 **has not run on the user's installation**.
- **Status:** no newly source-certified neutrino masses; `0/16` remains correct. No matching, EFT or RGE.

## 1. Root cause of the B_-1 failure

The V1 native witness incorrectly expected `T3HiggsPortal1CGInv1` and `T3HiggsPortal2CGInv1` to equal the identity in the neutral Higgs projection. That was true for the **A_+0 singlet** but not for the **B_-1 doublets**. The failed checks did *not* invalidate native CG registration, Majorana mass coefficient `-1/2`, either neutral Yukawa coefficient `+1`, or the native mixing entry `-2/sqrt(3)`.

In the native ordered field axes `(Bar[H]_i,H_j,Bar[S]_a,S_b)` the two B_-1 portal invariants are, writing `I_{ijab}=delta_ij delta_ab` and `E_{ijab}=delta_ib delta_aj`,

\[
 C_{\mathrm{Inv1},ijab}=E_{ijab},\qquad
 C_{\mathrm{Inv2},ijab}=\frac{2I_{ijab}-E_{ijab}}{\sqrt3}.
\]

The components follow directly from the installed Matchete CG registry's sparse arrays. Thus the **exact full off-shell identity** is

\[
 \boxed{\ I_{ijab}=\tfrac12 C_{\mathrm{Inv1},ijab}
                      +\tfrac{\sqrt3}{2} C_{\mathrm{Inv2},ijab}.\ }
\]

No neutrino-mass fitting enters this identity. Both scalar doublets have this same source tensor basis.

For `H[2]=v/sqrt(2)` the individual Hermitian source projections are

\[
 C_{\mathrm{Inv1}}\big|_{H^0}=\begin{pmatrix}0&0\\0&1\end{pmatrix},\quad
 C_{\mathrm{Inv2}}\big|_{H^0}=\begin{pmatrix}2/\sqrt3&0\\0&1/\sqrt3\end{pmatrix}.
\]

Neither is the identity separately. Their combination exactly equals `I2`. The V1 witness checked a false physical claim; a false result was appropriate.

## 2. Correct benchmark-to-native portal assignment

The native action uses real couplings `lambdaH1Inv1`, `lambdaH1Inv2`, `lambdaH2Inv1`, `lambdaH2Inv2`. The frozen numerical benchmark instead gives two numbers `lambdaH1=0.02`, `lambdaH2=0.03`, intended as coefficients of the **full identity portal**, for the restricted benchmark point. To reproduce this *selected*, not automatically mandated, operator direction in the native basis, set

\[
\begin{aligned}
\lambda_{H1\mathrm{Inv1}}&=\tfrac12\lambda_{H1}=0.0100000000000, &
\lambda_{H1\mathrm{Inv2}}&=\tfrac{\sqrt3}{2}\lambda_{H1}=0.0173205080757,\\
\lambda_{H2\mathrm{Inv1}}&=\tfrac12\lambda_{H2}=0.0150000000000, &
\lambda_{H2\mathrm{Inv2}}&=\tfrac{\sqrt3}{2}\lambda_{H2}=0.0259807621135.
\end{aligned}
\]

A nonzero `Inv2` is **required** for an identity portal in B_-1. Explicitly setting `Inv2=0`, as the original V1 documentation asserted, does **not** reproduce the full benchmark scalar masses. For *this exact selected linear combination*, the source scalar potential has diagonal squared-mass shifts

\[
\Delta m^2_{S1}=-\frac{\lambda_{H1}v^2}{2}\,I_2,\qquad
\Delta m^2_{S2}=-\frac{\lambda_{H2}v^2}{2}\,I_2.
\]

The neutral off-diagonal scalar mixing is unaffected:

\[
 (A_{S2^0,S1^0})=-\frac{\lambda_{T3}v^2}{2}
 C_{\rm mix}[2,2,1,2]
 =+\frac{\lambda_{T3}v^2}{\sqrt3}
 =3500.1449227\;\mathrm{GeV}^2
\]

for `lambdaT3=0.1`, `v=246.22 GeV`. The external neutrino Yukawa factors are `C1=C2=+1` and the native Majorana mass term is `-MF/2`, with no extra `1/2` multiplying its one-loop self-energy.

## 3. What the two new self-contained scripts check

`RunPhysicalBm1NativeScalarQuadraticV2.wl` builds the restricted action in a fresh Matchete session, independently extracts source coefficients and registered tensors, checks both neutral CG Yukawas and the Majorana mass factor, checks **the full 16 components** of both portal identity combinations (and the four individual Higgs VEV projectors), and writes a fresh V2 native JSON. It never writes over a prior witness and does not modify source code.

`ReviewPhysicalBm1NativeLoopV2.py` refuses a wrong V1 witness; validates the V2 native field/action source hashes; **solves the full four-index identity decomposition using the measured native component tensors** with a real linear system; verifies `Inv1=1/2` and `Inv2=sqrt(3)/2` to numerical precision; and only then evaluates the independent `p=0` Feynman-parameter loop. `--portal-assignment benchmark-identity` records the four nonzero native coupling values explicitly in the JSON, with `--compare-frozen` reading the pre-existing provisional output only after the independent evaluation.

The independent loop formula remains

\[
(m_\nu)_{ab}=\frac1{16\pi^2}\sum_{r}M_r A_{21}\,{\cal D}(M_r^2;s_-,s_+)
 (g_{1,ar}g_{2,br}+g_{2,ar}g_{1,br}),
\]

where `D` is the ultraviolet-finite integral evaluated directly with `mpmath`, and no extra Majorana loop combinatorial factor is added. The identity-portal benchmark matrix and Takagi masses should reproduce the *previously derived* values (they are NOT new native evidence):

\[
 m_\nu/\mathrm{eV}\simeq\begin{pmatrix}
 -0.001231240638&-0.009878342256&+0.000945610552\\
 -0.009878342256&-0.019334391098&+0.027236878856\\
 +0.000945610552&+0.027236878856&-0.024594786699
 \end{pmatrix},
\]

with largest Takagi mass `0.05044130936 eV` for benchmark three-diagonal-mass extension. Reproducing these numbers is an **internal consistency check** after source-derived normalization; not a new physical certification. The native restricted Matchete mass is one symbol `MF`. The three distinct frozen masses (100, 110, 120 TeV) are a separate benchmark generalisation not yet a direct native mass witness. Absolute installed `GammaCC` orientation remains unresolved.

## 4. Installation and run, without overwrites

From repository root:

```powershell
$zip = Join-Path $HOME 'Downloads\T3_Bm1_Portal_Identity_V2_20261103.zip'
$paths = @(
    'scripts\RunPhysicalBm1NativeScalarQuadraticV2.wl'
    'scripts\ReviewPhysicalBm1NativeLoopV2.py'
    'Numerical\lo_certification\physical_normalisation\BM1_NATIVE_MAJORANA_PORTAL_IDENTITY_V2.md'
)
foreach ($p in $paths) { if (Test-Path -LiteralPath $p) { throw "Refuse overwrite: $p" } }
Expand-Archive -LiteralPath $zip -DestinationPath .
$stamp = Get-Date -Format 'yyyyMMdd_HHmmss'
$runDir = "output\general_observable\physical_normalisation\bm1_native_v2_$stamp"
$native = "$runDir\B_-1_native_action_v2.json"

wolframscript -file scripts\RunPhysicalBm1NativeScalarQuadraticV2.wl $native
if ($LASTEXITCODE -ne 0) { throw 'B_-1 V2 native extraction failed' }
python scripts\ReviewPhysicalBm1NativeLoopV2.py --repo-root . --native-json $native --fermion-masses native-common --output-dir "$runDir\native_common_zero"
if ($LASTEXITCODE -ne 0) { throw 'B_-1 native common mass failed' }
python scripts\ReviewPhysicalBm1NativeLoopV2.py --repo-root . --native-json $native --fermion-masses benchmark-diagonal --portal-assignment benchmark-identity --compare-frozen --output-dir "$runDir\benchmark_identity"
if ($LASTEXITCODE -ne 0) { throw 'B_-1 benchmark extension failed' }
```

## 5. Progress ledger

- **Independent physical evidence here:** exact off-shell Matchete-CG two-portal identity decomposition, derived from the actual 16 native tensor components for each scalar doublet; original V1 native run confirms which false assertions blocked. No new local installed-native success is claimed until running V2.
- **Gate movement:** G3 B_-1 blocker identified and repaired in unexecuted V2 scripts; native G3 acceptance remains pending. G1 A_+0 scalar-action subgate passed previously; G2 A_+0 remains independently numerical but not fully physically certified.
- **Physical UV source certification:** **0/16**.
- **Next decisive step:** execute the fresh V2 native Matchete extraction and then the independent source-derived B_-1 loop, keeping the native-common `MF` and benchmark-diagonal extension separate. Later resolve installed `GammaCC` and B2 coupling normalisation before certifying physical masses.
