# A_+0: native scalar quadratic witness and independent UV loop

**Scope:** restricted T3, broken phase, direct UV one-loop, B1. No EFT, RGE or GitHub write.

**Source provenance:** latest GitHub `main` inspected at `f82ee42ea2456e07e831a381fcbf2f2a3dd747d5` (2026-11-03). It is one commit after the user's reported local `1fb8f4696a5568ec80f5f6dd5c8d6c61780ce40c`. The intervening commit adds scripts/results but does not edit `T3Topology.wl`, `T3Fields.wl`, `SU2Invariants.wl`, `ScalarPotential.wl`, `LagrangianBuilder.wl`, or the native CG reports. Local working tree is **not** assumed identical; the new Wolfram witness embeds local source file SHA256 values. Read-only GitHub inspection is not the same as execution on installed Matchete.

## 1. Physics derived from actual source operators and registered CG dumps

- Native `FreeLag` (see `output/general_observable/physical_normalisation/b1_02/native/A_p0.json`): canonical `+|D S_1|^2+|D S_2|^2 -MS1^2 |S_1|^2-MS2^2 |S_2|^2` without extra factors.
- Original Yukawas in `LInt`: `+ y1 Bar[L] CConj[F] S1 CG1 + y2 Bar[L] F Bar[S2] CG2 + H.c.`. Native Dirac mass coefficient is `-1`.
- Native scalar portals are `+lambdaH1Inv1 Hdag H S1dag S1 CG` and `+(1/2) PlusHc[lambdaH2Inv{k} Hdag H S2dag S2 CG_k]`, **not** one automatically identical Matchete parameter called `lambdaH2`.
- The original mixing is `+lambdaT3 H H S1 Bar[S2] T3MixCG + H.c.`, **with no 1/2 prefactor** before substituting identical Higgs fields. The VEV contributes `v^2/2` exactly once.

With Matchete's registered original Cartesian triplet indices (not an imposed triplet weight rotation), define `c = 3^(-1/4)`. The neutral external-lepton and fermion component is `L[1]` and `F[2]` (fermion hypercharge +1/2). Registered components are:

```
T3Y1CG[1,2]                   = +1
T3Y2CG[1,2,1..3]              = c ( +1, 0, -i )
T3MixCG[2,2,1..3]             = c ( -1, 0, +i )
T3HiggsPortal1CGInv1[2,2]     = 1
Hermitian T3HiggsPortal2CGInv1[2,2] = 0_(3x3)
Hermitian T3HiggsPortal2CGInv2[2,2] = I_(3x3)
```

The first triplet Higgs portal is an independent **gauge invariant**, but its value at `H[1]=0,H[2]=v/sqrt(2)` vanishes. The second one supplies the identity mass shift. The installed-native script explicitly measures these statements again from *the assembled Matchete action* and `GetCGTensor`, not from an assumed weight basis or provisional input.

## 2. Independent physical scalar quadratic matrix

Define `Phi=(S1, S2_cartesian_1,S2_cartesian_2,S2_cartesian_3)^T`, `L_quadratic=-Phi^dag A Phi`. From the native action and `V=-L`,

\[
 A=\begin{pmatrix}
 m_{S1}^2-\frac{v^2}{2}\lambda_{H1Inv1}&
 \frac{v^2}{2}\lambda_{T3}^*\,c\,(1,0,+i)\\
 \frac{v^2}{2}\lambda_{T3}\,c\,(1,0,-i)^T&
 (m_{S2}^2-\frac{v^2}{2}\lambda_{H2Inv2})\mathbf1_3
 \end{pmatrix}.
\]

Here the native `lambdaH2Inv1` has zero H0 projection. For complex `lambdaT3`, the two off-diagonal blocks are Hermitian conjugates. The scalar mix magnitude is `sqrt(2) c |lambdaT3| v^2/2` (no additional charged-channel multiplicity). With `v=246.22 GeV`, `lambdaT3=0.1`:

```
|A_{S2_1,S1}|=|A_{S2_3,S1}| = 2303.224887177694 GeV^2
norm(A_{S2,S1})                = 3257.251872641937 GeV^2
```

If the **additional explicit parameter identification** `lambdaH1Inv1=0.02`, `lambdaH2Inv1=0`, `lambdaH2Inv2=0.03` is chosen, `A[0,0]=1305516.206095592 GeV^2`, `A[S2,S2]=(999090.635674) I_3 GeV^2`. This identifies **native** couplings by assignment; it is not an independently measured RGBeta-to-Matchete coupling dictionary.

## 3. Exact non-provisional one-loop implementation

For positive, diagonal Dirac flavour masses `M_r`, use the original native Cartesian scalar CG tensor, a direct Hermitian diagonalisation `U^dag A U=diag(m_s^2)`, and source-converted two-component couplings

\[
Y_{2;a r s}=-\overline{y_{2,ar}\,C_{2;\nu F^0 s}},\qquad
Z_{1;a r 0}=-\overline{y_{1,ar}\,C_{1;\nu F^0}}.
\]

The `Y2` vertex couples to `S2`, the `Z1` vertex couples to `S1^*`. With the all-left-handed Dirac propagator `<xi_r chi_r>`, and `Phi_s=U_{sk}(R_k+iI_k)/sqrt(2)`, the `R` and `I` contributions cancel the `YY` and `ZZ` pieces, leaving exactly one cross contraction plus its transpose:

\[
(m_\nu)_{ab}=\sum_{r,k}\frac{M_r}{16\pi^2}
B_0^{\rm fin}(0;M_r^2,m_k^2;\mu^2)
\left[(Y_{2,r}U)_{ak}(Z_{1,r}U^*)_{bk}+(a\leftrightarrow b)\right].
\]

\[
B_0^{\rm fin}(0;M^2,m^2;\mu^2)
=-\int_0^1 dx\,\ln\frac{x M^2+(1-x)m^2}{\mu^2}.
\]

The new Python reviewer implements this expression directly using only NumPy, the live witness, and the frozen *parameter* benchmark. It **does not** import `Numerical/direct_one_loop/core.py`, provisional Weyl tensors, or `RunPhysicalA0SourceLoop.py`. The UV pole coefficient is checked to vanish by summing the same contractions without `B0`.

## 4. Numerical predictions and what has actually run

The previous user terminal run measured the independent **source-defined zero-portal** matrix in eV

```
[[ 0.0026345483, -0.0013019240, 0.0107526974],
 [-0.0013019240,  0.0207211975, 0.0264372995],
 [ 0.0107526974,  0.0264372995, 0.0213538700]]
```

with descending Takagi masses `(0.048495052733, 0.010999154416, 0.007213717548) eV`, and a separate **conditional benchmark-identity** mass matrix with largest singular mass about `0.048498954012 eV`. The earlier user terminal reported conditional frozen relative residual `9.017833944886156e-12`.

In the present preparation, the new Python Cartesian-basis calculation was numerically exercised **only with an explicitly simulated native-style CG record built from the public, previously inspected native registry numbers**, not with the installed Matchete. It reproduced the zero-portal and benchmark-identity mass matrices to displayed precision and gave a UV pole norm `~4e-23 GeV`. This simulation is **NOT a new installed-native measurement**, is not included in the delivered evidence, and must not change the `0/16` certified count.

**Next new evidence:** execute `RunPhysicalA0NativeScalarQuadratic.wl` under the actual Windows Matchete installation. It fails closed if the live assembled action has different coefficients or different registered CG tensors. Feed that genuinely new JSON to `ReviewPhysicalA0NativeScalarQuadratic.py` for an independent native-Cartesian calculation. Its output directory must be new; none of the previous evidence is overwritten.

## 5. Run from the user's repository root (PowerShell)

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
$stamp = Get-Date -Format 'yyyyMMdd_HHmmss'
$native = "output\general_observable\physical_normalisation\a0_native_quadratic_$stamp\A_+0_native_action.json"

wolframscript -file scripts\RunPhysicalA0NativeScalarQuadratic.wl $native
if ($LASTEXITCODE -ne 0) { throw "Native Matchete quadratic witness blocked" }

.\.venv\Scripts\python.exe scripts\ReviewPhysicalA0NativeScalarQuadratic.py `
  --repo-root . --native-json $native `
  --output-dir "output\general_observable\physical_normalisation\a0_native_quadratic_$stamp\zero_portal"
if ($LASTEXITCODE -ne 0) { throw "Independent native-basis zero-portal calculation blocked" }

.\.venv\Scripts\python.exe scripts\ReviewPhysicalA0NativeScalarQuadratic.py `
  --repo-root . --native-json $native `
  --portal-assignment benchmark-identity --compare-frozen `
  --output-dir "output\general_observable\physical_normalisation\a0_native_quadratic_$stamp\benchmark_identity"
if ($LASTEXITCODE -ne 0) { throw "Independent native-basis benchmark-identity calculation blocked" }
```

## 6. Physical acceptance criteria

1. Installed-native `CheckLagrangian` returns Success with four experimental flags zero, canonical scalar kinetic `+1`, scalar free mass `-1`, original/HC mixing `+1`, S2 portal pairs `1/2+1/2`.
2. Registered `T3MixCG[H0,H0]` equals `(-c,0,+ic)` and H0 portal projections are `(1,0,I3)` in the **native** component basis (no independent fake weight rotations).
3. Independent Cartesian-basis loop for the zero-portal source point returns the already measured matrix to `1e-8` matrix-relative tolerance, and UV pole cancels. The optional benchmark-identity comparison should also be within `1e-8`; it remains **conditional** on coupling identification.
4. If any check fails, inspect that discrepancy as a physics-sign/normalisation issue; do not alter production or raise `physical_source_certified_models` from `0`.

**Remaining physical blockers even if this passes:** demonstrate the benchmark's actual numerical portal coupling dictionary and validate the native four-spinor-to-Weyl orientation up to consistent rephasings. B2 cross-package alignment is a separate dependency for cross-package parameters, not for an internally self-contained Matchete UV amplitude.
