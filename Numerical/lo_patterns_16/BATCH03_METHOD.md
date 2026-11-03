# Chat A — Batch 3: exact scalar-kernel explanation of the 16 frozen T3 patterns

## Scientific objective

Identify which parts of the 16 conditional one-loop neutrino-mass matrices are consequences of SU(2) charge selection and which arise from the **particular** frozen Yukawa textures and heavy-mass spectrum. This batch is additive and READ ONLY; it does not run or change production `core.evaluate`, fit Yukawas, run RGEs, integrate out heavy particles, or claim B1/B2 physical-normalisation certification.

This package targets the user's locally reported `0b210f18f6a2cfb10fff03059e2ec4a73d2195a4` checkout and checks the **actual local source/file SHA256** against Batch 1 before evaluation. GitHub may subsequently advance. Do not replace local files with newer GitHub snapshots. The 16 genuine benchmark evaluations and 288 scans were already reported as successful in the user's Windows checkout; this module is a *new analytic regression*, not a repeat of that verification.

## 1. Exact mixing matrix function in the frozen scalar sector

The source-conditional Weyl input uses

\[
\mathcal L\supset -\nu_i\psi_a\big[Y_{ias}\phi_s+Z_{ias}\phi_s^*\big]+\text{h.c.}, \qquad
V_2=\phi^\dagger A\phi +\operatorname{Re}(\phi^T B\phi).
\]

The **restricted** frozen inputs have `B=0`. The complex scalar mass matrix `A` has diagonal S1 and S2 blocks and (only) electric-charge-preserving S1–S2 mixing. Since each of the singlet, doublet and triplet multiplets has exactly one field with each electric charge, each nonzero mixing block is an **isolated 2×2 Hermitian matrix**

\[
 A_q=\begin{pmatrix}a_q&\delta_q^*\\\delta_q&b_q\end{pmatrix}, \qquad
 \delta_q=A_{S_2(t_q),S_1(s_q)},
 \quad\Lambda_{q\pm}=\frac{a_q+b_q}{2}\pm\sqrt{\frac{(a_q-b_q)^2}{4}+|\delta_q|^2}.
\]

The frozen candidate input builder sets \(\delta_q=-\lambda_{T3}v^2(T3MixCG)_{00s_qt_q}/2\). This sign/phase is **only the provisional builder convention**; B1/B2 must separately establish its physical meaning.

For a positive heavy squared mass \(m_k^2\), write the scalar-mass-dependent part of the finite \(B_0\) integral as

\[
 f_k(x)=-\frac{x}{x-m_k^2}\ln\frac{x}{m_k^2},
 \qquad f_k(m_k^2)=-1.
\]

The scalar-independent terms cancel from the radiatively generated cross-mass because the offdiagonal element of an identity matrix is zero. For any analytic function \(f\), spectral decomposition of the 2×2 scalar block gives the **exact identity**

\[
 [f(A_q)]_{S_2,S_1}
 =\delta_q \frac{f(\Lambda_{q+})-f(\Lambda_{q-})}{\Lambda_{q+}-\Lambda_{q-}}.
\]

Thus the full scalar kernel is

\[
 K_{qk}^{\rm exact}=\delta_q
 \frac{B_0(0;m_k^2,\Lambda_{q+})-B_0(0;m_k^2,\Lambda_{q-})}
 {\Lambda_{q+}-\Lambda_{q-}}.
\]

At leading small mixing,

\[
 K_{qk}^{(1)}=\delta_q\frac{B_0(0;m_k^2,a_q)-B_0(0;m_k^2,b_q)}{a_q-b_q},
 \quad \text{equal masses: } \delta_q\frac{\partial B_0}{\partial a_q}.
\]

The difference quotient is finite even when \(a_q=b_q\). The exact kernel is an odd function of \(\delta_q\) for a real common rescaling while its magnitude is affected at cubic and higher order. Therefore

\[
 m_\nu(-\lambda)=-m_\nu(\lambda),\quad m_\nu(0)=0,
 \qquad m_\nu(2\lambda)-2m_\nu(\lambda)=O(\lambda^3).
\]

These are **conditional numerical input-space identities**, not independent proofs of the physical Lagrangian normalisation.

## 2. Explicit heavy propagator, charge and flavour contractions

The code does **not** call the original Takagi diagonaliser. Instead it diagonalises the Hermitian matrix \(F^2=MM^\dagger\), using \(F^2=Q\,\mathrm{diag}(m_k^2)Q^\dagger\) and projectors \(P_k=Q_{:k}Q_{:k}^\dagger\). With the frozen vertices `y1_HC` in `Z` (S1*) and `y2_HC` in `Y` (S2), the exact mass is

\[
 m_\nu=\frac{1}{16\pi^2}\sum_{q,k}K_{qk}
 \left[Z_{s_q}\,(M^\dagger P_k)\,Y_{t_q}^{T}+
       \big(Z_{s_q}\,(M^\dagger P_k)\,Y_{t_q}^{T}\big)^T\right].
\]

Here \(Z_s,Y_t\) are 3×\(n_F\) matrices carrying the active-lepton flavour index, and the charged component index is fixed **before** any heavy diagonalisation. The two summands ensure symmetric Majorana flavour indices with no invented extra factor of two. The native conditional reference uses the same propagator but obtains scalar eigenmodes by diagonalising the **full Cartesian Hessian**; Batch 3 tests the result independently using the closed 2×2 formula.

Grouping exactly/approximately degenerate heavy squared-mass eigenvectors into complete projectors produces a basis-invariant set of heavy-singular-mass contributions \(m_{\nu,g}\). The code reports their coherence ratio

\[
 C_F=\frac{\left\|\sum_g m_{\nu,g}\right\|_F}
             {\sum_g\|m_{\nu,g}\|_F}\leq1.
\]

This diagnoses possible **interference between different heavy mass groups**, as distinct from interference between charges. The grouping is by physical singular mass, not by an assumed flavour-generation index; **individual eigenvectors inside degenerate heavy spaces are not invariant**.

Likewise, the charge-channel coherence ratio is

\[
 C_Q=\frac{\left\|\sum_q m_{\nu,q}\right\|_F}
             {\sum_q\|m_{\nu,q}\|_F}\leq1.
\]

A small singular value \(m_1\) does **not** imply \(C_Q\ll1\); it may instead indicate Yukawa flavour alignment or cancellations between heavy-mass eigenspaces. Spectra are obtained from the **total complex symmetric matrix**, not from the sums of the channel Takagi masses.

## 3. What the observed patterns support

1. **A versus D:** The source CG candidate weight is the same, \(-2/\sqrt3\), and both have one active topology channel. The apparent near equality of A_-4 versus D_+2 and analogous matched pairs may arise from swapped S1/S2 representations, charge reflection and nearly identical frozen Yukawa textures. The script tests matrix *shape* after fitting one complex overall factor rather than claiming the models are identical.
2. **B:** One channel with candidate group weight \(-2/\sqrt3\); B_-1 has a neutral self-conjugate Majorana singlet, unlike B_-3 and B_+1. The associated numerical differences must be interpreted with its different fitted Yukawas. The physical Majorana 4-spinor normalization is **pending B1**.
3. **C:** One channel with candidate weight \(+2/3\) for all its alphas. Thus the enormous differences in the lightest two eigenmasses across C are **not an extra charge multiplicity**. Their frozen Yukawa matrices are substantially different (the C_-1 product of Frobenius norms is roughly 3 times that of C_+1). The script reports singular-value hierarchy and heavy-mass-group interference to test the actual flavour origin.
4. **E:** Two channels of equal candidate weight \(+\sqrt6/3\). The source's saved component partition for the four frozen E benchmarks gives two **positive, aligned** charge-matrix contributions, each approximately half of the total. This is a property of the current restricted benchmarks, not a theorem about every parameter choice in E. In addition, E uses \(|\lambda_{T3}|=0.05\), half the A–D benchmark value. The script independently reconstructs each charge matrix and tests the constructive alignment.
5. **Hypercharge alpha:** In a fixed representation the source candidate charge set translates with alpha while the CG magnitude stays constant. Other differences across alpha can come from changed fitted Yukawa structures, so no claim of universal alpha-independent physical masses follows.

There is no generic ``m3 is larger in family X'' theorem without controlled equal-parameter inputs and the B1/B2 physical cross-package map. **The checkable result here is a decomposition and reproduction of each frozen matrix**, not a universal ranking.

## 4. Validation and restrictions

- Requires complete Batch-1 (16/16) and Batch-2 (16/16, 288 scans) outputs, all source SHA256 guards and matching frozen benchmark/component hashes.
- Compares the analytic **exact** 3×3 complex matrix against saved production mass to relative tolerance \(10^{-6}\), and checks analytical predictions of `lambdaT3` scan residuals at \(t=-1,0.5,2\) against Batch 2.
- Independent synthetic tests cover Dirac and Majorana toy representations, complex delta phase, negative sign/conjugation controls, Yukawa source contract drift, scalar B/new-operator rejection, rephasing invariance, near-degenerate squared masses, Feynman-parameter integration, all 16 toy report rows, and output non-overwrite.
- Only a **new** output tree is created: `output/general_observable/lo_pattern_study16/analysis_batch03_01/`.
- No Wolfram/Matchete execution is claimed. The symbolic group weights are source-audited candidate weights, not a B1/B2 certificate.

## Windows commands

```powershell
Set-Location 'C:\thesistemp\T3neutrinomassmodelpipeline'
.\.venv\Scripts\python.exe -m pytest -q tests\test_lo_analytic_explanation_16.py
.\.venv\Scripts\python.exe scripts\ExplainLOPatterns16.py
```

If the output directory already exists, the script refuses to overwrite it. To run again, use a **new** directory within `output/general_observable/lo_pattern_study16` with `--output-dir`.
