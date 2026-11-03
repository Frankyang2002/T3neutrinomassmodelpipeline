# B1 native spinor convention and complex-coupling UV test — 2026-11-03

## Scope and provenance

The existing **installed-Matchete** scalar/CG witnesses are not rerun:

- A_+0: `output/general_observable/physical_normalisation/a0_native_quadratic_20261103_113729/A_+0_native_action.json`.
- B_-1: `output/general_observable/physical_normalisation/bm1_native_v2_20261103_121426/B_-1_native_action_v2.json`.
- C_-1: `output/general_observable/physical_normalisation/cm1_native_20261103_122426/C_-1_native_action.json`.

The native spinor script newly assembles only the same three **restricted** UV actions to record the actual original and Hermitian-conjugate spinor operators, canonical Dirac/Majorana mass prefactors, and what the installed Matchete version exposes about `GammaCC`. It never calls `Match`, RGBeta, RGE or EFT. Both scripts fail closed on changed source SHA256, missing native evidence, missing model, incorrect CG or altered action coefficients. They never overwrite output.

The GitHub source inspected read-only when preparing this task was `2ac58cc3de12c7af11e2e7e214415dadf0138f4f`. The local Windows `git rev-parse HEAD` is read by the native script and stored; no pull or GitHub edit occurs.

## Physics that the source fixes

The assembled Matchete `LInt` has the two distinct Yukawa monomials and **one H.c. each**, with coefficient `+1`. The complex source parameters are `y1[p,r]`, `y2[p,r]`, `lambdaT3`; the Majorana B_-1/C_-1 `LFree` has the `-MF/2` coefficient, while A_+0 has `-MF Bar[F] F`. Matchete's `GammaCC` is normally an abstract Clifford-algebra operator, not an installed numerical 4x4 matrix with an absolute phase. The script records that limitation rather than manufacturing a numerical convention.

In the *declared* two-component basis `F=(chi,bar(xi))`, choose `C=+i gamma^2 gamma^0` with `epsilon_12=+1`, and define

\[
 g_{1,ar}=-[y_{1,ar}C_1]^*,\qquad
 g_{2,ar}=-[y_{2,ar}C_2]^*,\qquad
 (m_\nu)_{ab}=\sum_r\frac{M_r A_{21}\mathcal D_r}{16\pi^2}
 (g_{1,ar}g_{2,br}+g_{2,ar}g_{1,br}).
\]

Here \(A_{21}=-\lambda_{T3}v^2 C_m/2\), with \(C_1,C_2,C_m\) projected from the previously **live** native CG tensors. The neutral coefficients are:

| Model | \(C_1\) | \(C_2\) | \(C_m\) |
|---|---|---|---|
| A_+0 | 1 | \(-i\sqrt2\,3^{-1/4}\) | \(+i\sqrt2\,3^{-1/4}\) |
| B_-1 | 1 | 1 | \(-2/\sqrt3\) |
| C_-1 | \(3^{-1/4}\) | \(-3^{-1/4}\) | \(-2/\sqrt3\) |

The loop uses the finite scalar propagator difference

\[
\mathcal D(M^2;s_-,s_+)=-\int_0^1\frac{dx}{s_+-s_-}
 \log\!\left(1+\frac{(1-x)(s_+-s_-)}{xM^2+(1-x)s_-}\right),
\]

and **independently** checks the high-precision result against the analytic finite part of \(B_0(0;M^2,s)\) at \(\mu^2=M^2\):

\[
B_0=1+\frac{s\log(s/M^2)}{M^2-s},\qquad
\mathcal D=\frac{B_0(s_+)-B_0(s_-)}{s_+-s_-}.
\]

Neither `Numerical/direct_one_loop/core.py` nor any provisional UV mass matrices or tensors are imported.

## Genuine complex-coupling point and decisive negative controls

For each model, begin with the repository's fitted-reference `y1`, `y2`, `lambdaT3` and chosen *identity-portal* parameter assignment. Deform the **actual source** couplings:

- `y1[0,1] *= exp(+0.37 i)`, `y1[1,2] *= exp(-0.19 i)`;
- `y2[2,0] *= exp(+0.41 i)`, `y2[0,2] *= exp(-0.23 i)`;
- `lambdaT3 *= exp(+0.29 i)`.

These are deliberately nonuniform phases and generally generate a genuinely complex symmetric 3x3 \(m_\nu\), not a global rephasing of a real reference matrix.

Then globally change the source couplings by \(y_1\to e^{i\theta_1}y_1\), \(y_2\to e^{i\theta_2}y_2\), \(\lambda\to e^{i\theta_\lambda}\lambda\) with \((\theta_1,\theta_2,\theta_\lambda)=(0.21,-0.17,0.43)\). The expected matrix change **in this declared convention** is

\[
 m_\nu\to e^{i(\theta_\lambda-\theta_1-\theta_2)}m_\nu
 =e^{0.39 i}m_\nu.
\]

The reviewer requires a relative phase-transport residual below \(2\times10^{-12}\). As a negative control it conjugates `y1` **without changing the native action** and requires the resulting full complex matrix to differ by more than \(10^{-4}\) in relative norm. It also exports the CP-odd invariant \(\mathrm{Im}[H_{e\mu}H_{\mu\tau}H_{\tau e}]/(\mathrm{Tr}H/3)^3\), with \(H=m_\nu m_\nu^\dagger\).

**Important boundary:** These are strong tests of the complex *chosen Weyl* implementation, but they are not an independent measurement of Matchete's numerical `GammaCC` phase. A simultaneous overall fermion charge-conjugation rephasing can change a matrix's overall sign without changing the physical Takagi spectrum. An amplitude's absolute sign therefore must not be inferred merely from comparison of Takagi masses. The source-frozen mapping to an identity portal and three independent fermion masses is still a stated choice, not a native requirement.

## Windows execution

After extracting the ZIP into `C:\thesistemp\T3neutrinomassmodelpipeline` (new paths only):

```powershell
cd C:\thesistemp\T3neutrinomassmodelpipeline
$stamp = Get-Date -Format 'yyyyMMdd_HHmmss'
$dir = "output\general_observable\physical_normalisation\b1_complex_spinor_$stamp"
$spin = "$dir\native_spinor_action.json"

wolframscript -file scripts\RunPhysicalB1SpinorConventionNative.wl $spin
if ($LASTEXITCODE -ne 0) { throw 'Native spinor audit failed' }

python scripts\ReviewPhysicalB1ComplexConvention.py `
    --repo-root . `
    --spinor-native $spin `
    --a0-native 'output\general_observable\physical_normalisation\a0_native_quadratic_20261103_113729\A_+0_native_action.json' `
    --bm1-native 'output\general_observable\physical_normalisation\bm1_native_v2_20261103_121426\B_-1_native_action_v2.json' `
    --cm1-native 'output\general_observable\physical_normalisation\cm1_native_20261103_122426\C_-1_native_action.json' `
    --output-dir "$dir\complex_loop"
if ($LASTEXITCODE -ne 0) { throw 'Independent complex UV check failed' }
```

Output path: `output/general_observable/physical_normalisation/b1_complex_spinor_<timestamp>/complex_loop/B1_complex_spinor_convention_witness.json`. No existing evidence is modified. The native spinor script is expected to report the unresolved absolute `GammaCC` phase unless a concrete installed matrix is genuinely available. The reviewer never promotes physical certification: **0/16**.

## Validation performed when preparing the package

- `ReviewPhysicalB1ComplexConvention.py` compiled with CPython and completed a **transcribed source-shaped fixture**, including three independent complex mass matrices, quadrature/analytic B0 identity and phase-covariance negative controls. Fixture values are **not live native evidence and not fitted benchmark results**.
- `RunPhysicalB1SpinorConventionNative.wl` was inspected statically against previously successful native Matchete source-script patterns; it was **not** executed in a Mathematica/Matchete kernel in this environment.

**Next acceptance criterion:** a concrete installed spinor representation or an independently extracted Feynman rule fixing the relative (not physically unobservable overall) conjugation of the Yukawa vertices to the mass insertion. Without that, mark the relative-spinor phase gate as pending even if the numerical complex test succeeds.
