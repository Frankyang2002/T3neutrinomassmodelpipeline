# Chat A Batch 8 — Scalar-splitting convergence and factorization-breaking projection

## Scope

Read-only postprocessing of the **conditional** candidate Weyl UV direct-loop inputs and saved Batch-6/7 results. No fermion/scalar integration, no SMEFT, no production UV source changes, and no physically certified Matchete/RGBeta-to-Weyl matching. The split diagnostic was introduced by Batch 6 and is an **ad hoc, U(1)em-preserving mass deformation**, not a proven gauge-invariant UV interaction.

The result applies to the restricted scalar quadratic with zero anomalous scalar bilinear (`scalar_B_gev2=0`), diagonal within-multiplet masses, and disjoint equal-charge two-scalar blocks. All 16 model coefficients use the *same* source candidate tensor convention as Batches 5–7.

## Exact scalar kernel

For a fixed heavy eigenspace \(g\) and scalar charge \(q\), let

\[
A_q(t)=\begin{pmatrix}a_q+t\epsilon_{a,q}&\delta_q^*\\
\delta_q&b_q+t\epsilon_{b,q}\end{pmatrix},\qquad
s=\sqrt{(a-b)^2+4|\delta_q|^2},\quad
\lambda_\pm=\tfrac12(a+b\pm s).
\]

All derivatives below hold \(\delta_q\) fixed, and use \(f(x)=B_0(0;M_g^2,x)\):

\[
K_q=\delta_q D,\qquad D=\frac{f(\lambda_+)-f(\lambda_-)}{s},
\]

\[
D_+=\frac{f'(\lambda_+)-D}{s},\qquad D_-=\frac{D-f'(\lambda_-)}{s},
\]

\[
\partial_a K_q=\delta_q\bigg[D_+\frac{1+(a-b)/s}{2}
+D_-\frac{1-(a-b)/s}{2}\bigg],
\]

\[
\partial_b K_q=\delta_q\bigg[D_+\frac{1-(a-b)/s}{2}
+D_-\frac{1+(a-b)/s}{2}\bigg].
\]

The exact finite-mixing first derivative of the generation-space coefficient tensor is

\[
B_{g,1}^{(R)}=\sum_q
\left[\epsilon_{a,q}\partial_a K_q+\epsilon_{b,q}\partial_b K_q\right]
T_{qg}^{(R)},\qquad
T_{qg}^{(R)}=C^{(1)}_{q}(M^\dagger P_g) C^{(2)T}_{q}.
\]

Thus \(B_g(t)=B_g(0)+tB_{g,1}+O(t^2)\) as \(t\to0\). The exact derivative carries every odd order in scalar mixing. Batch 7's derivative of the *leading-in-\(\delta\)* kernel is a different approximation. Replacing the exact derivative with it leaves an \(O(t|\delta|^3)\) term, so a strict quadratic-convergence claim needs the exact derivative.

## Transverse (non-proportional) part

Stack the three 3×3 complex \(B_g\) matrices into 27-dimensional vectors: \(R(t)\) for A_-2, and \(B(t)\) for another model. The least-squares complex multiplier and residual are

\[
\kappa(t)=\frac{\langle R(t),B(t)\rangle}{\langle R(t),R(t)\rangle},
\qquad D(t)=B(t)-\kappa(t)R(t).
\]

With \(R(t)=R_0+tR_1+\cdots\), \(B(t)=B_0+tB_1+\cdots\),

\[
\kappa_1=\frac{\langle R_1,B_0\rangle+\langle R_0,B_1\rangle
-\kappa_0\big(\langle R_1,R_0\rangle+\langle R_0,R_1\rangle\big)}
{\langle R_0,R_0\rangle},
\]

\[
D_1=B_1-\kappa_0R_1-\kappa_1R_0.
\]

These are *complex* inner products. The derivative of \(\kappa\) is necessary: subtracting a fixed baseline multiplier would confound a pure rescaling with genuine nonproportionality.

For every \(t\) we save: exact and leading-mixing tangent remainders normalized by \(\|B(0)\|_F\); \(\|D(t)-D_0\|/\|B(0)\|\); its first-order prediction \(t\|D_1\|/\|B(0)\|\); the transverse remainder; and the pointwise residual \(\|D(t)\|/\|B(t)\|\). Relative transverse errors are `null` when its observable increment is beneath a 10^-9 relative signal floor. Convergence exponents are **not estimated** for noise-floor quantities.

## Exact Batch-6 scenario continuation

`split_input(t)` uses the same diagonal perturbation as Batch 6, multiplied by \(t\). The endpoints \(t=0\) and \(t=1\) are verified against the saved baseline and `charge_split_diagnostic` 3×3×3 tensors; the source inputs are untouched. The sampled values are \(0, 1/16, 1/8, 1/4, 1/2, 1\). A central finite difference at \(t=\pm10^{-4}\) independently validates the analytic derivative.

The all-model report checks the measured exact remainder follows a quadratic local power law on \(1/8\to1/4\), provided it is above the numerical floor; and the measured nonproportional increment follows a linear local power law when resolvable.

## Interpretation limits

- An \(O(t^2)\) remainder confirms the first derivative of this **particular diagnostic path**; it does not imply physical SU(2)-breaking scalar masses are realized in the actual T3 Lagrangian.
- Near-zero nonproportional response for some A models is correctly classified as an algebraically aligned/noise-limited case; it is not a fitted power law.
- The Batch-6 factorization residual compares coefficients across heavy eigenspaces, not only the final neutrino masses.
- Physical coefficient phases and overall normalization remain **0/16 source-certified** pending B1/B2.
