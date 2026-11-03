# Chat A Batch 8 — Controlled splitting convergence and transverse factorization

**Source-conditional, physically uncertified: 0/16.** This uses the saved candidate Weyl tensors, the original scalar split diagnostic and no EFT matching. The split is a numerical U(1)em-preserving mass deformation, **not** a validated gauge-invariant UV operator.

## Exact scalar-kernel differentiation

For a charge block A_q(t)=[[a_q+t Delta a_q, delta_q^*],[delta_q,b_q+t Delta b_q]], let \(\lambda_\pm=(a+b\pm s)/2\), \(s=\sqrt{(a-b)^2+4|\delta|^2}\), \(f(x)=B_0(0;M_g^2,x)\), and \(D=(f(\lambda_+)-f(\lambda_-))/s\).

\[K_q=\delta D,\qquad D_+=(f^{\prime}(\lambda_+)-D)/s,\quad D_-=(D-f^{\prime}(\lambda_-))/s.\]
\[\partial_aK_q=\delta\left[D_+{1+(a-b)/s\over2}+D_-{1-(a-b)/s\over2}\right],\]\[\partial_bK_q=\delta\left[D_+{1-(a-b)/s\over2}+D_-{1+(a-b)/s\over2}\right].\]

The exact slope is \(B_{g,1}=\sum_q(\Delta a_q\partial_aK_q+\Delta b_q\partial_bK_q)T_{qg}\). The diagonal mass shifts and the complex scalar mixing are held fixed as specified, while t varies.

Crucially, Batch 7 used the **leading-in-mixing** version of this slope. That is accurate for small mixing, but its difference from the exact slope produces an \(O(t|\delta|^3)\) term. Only the exact finite-mixing derivative gives a clean \(O(t^2)\) remainder.

## Nonproportionality: the quantity relevant to factorization

Flatten all heavy-generation coefficient matrices into complex vectors \(R(t)\) (reference) and \(B(t)\) (model). Define

\[\kappa(t)={\langle R(t),B(t)\rangle\over\langle R(t),R(t)\rangle},\qquad D(t)=B(t)-\kappa(t)R(t).\]

\[\kappa_1={\langle R_1,B_0\rangle+\langle R_0,B_1\rangle-\kappa_0(\langle R_1,R_0\rangle+\langle R_0,R_1\rangle)\over\langle R_0,R_0\rangle},\]\[D_1=B_1-\kappa_0R_1-\kappa_1R_0,\qquad D(t)=D_0+tD_1+O(t^2).\]

This is a complex least-squares projection, and the derivative includes the change in the *best-fit multiplier*. It compares the factorization-breaking signal itself rather than disguising it in the much larger total coefficient norm.

## All-16 scalar-splitting convergence

The first two errors below are normalized to the t=0 coefficient norm. They are evaluated at t=1/4. The last is the relative error of the **transverse increment** against its measured magnitude; n/a indicates a numerically vanishing transverse signal. A slope of 2 means quadratic convergence.

| Model | exact-derivative error at t=1/4 | Batch7-leading-derivative error | O(t²) convergence slope | transverse increment error | t=1 nonparallel residual |
|---|---:|---:|---:|---:|---:|
| A_-4 | 1.9814e-05 | 1.9815e-05 | 2.0059 | n/a | 6.9894e-16 |
| A_-2 | 1.9814e-05 | 1.9815e-05 | 2.0059 | n/a | 1.8138e-16 |
| A_+0 | 1.9814e-05 | 1.9815e-05 | 2.0059 | n/a | 6.6925e-16 |
| B_-3 | 2.4707e-05 | 2.4714e-05 | 2.0122 | 0.0025978 | 1.0840e-04 |
| B_-1 | 2.4707e-05 | 2.4714e-05 | 2.0122 | 0.0025978 | 1.0840e-04 |
| B_+1 | 2.4707e-05 | 2.4714e-05 | 2.0122 | 0.0025978 | 1.0840e-04 |
| C_-3 | 2.4707e-05 | 2.4714e-05 | 2.0122 | 0.0025978 | 1.0840e-04 |
| C_-1 | 2.4707e-05 | 2.4714e-05 | 2.0122 | 0.0025978 | 1.0840e-04 |
| C_+1 | 2.4707e-05 | 2.4714e-05 | 2.0122 | 0.0025978 | 1.0840e-04 |
| D_-2 | 5.6706e-05 | 5.6716e-05 | 2.0149 | 0.011471 | 2.2349e-04 |
| D_+0 | 5.6706e-05 | 5.6716e-05 | 2.0149 | 0.011471 | 2.2349e-04 |
| D_+2 | 5.6706e-05 | 5.6716e-05 | 2.0149 | 0.011471 | 2.2349e-04 |
| E_-4 | 1.8922e-05 | 1.8922e-05 | 2.0047 | 0.0022273 | 2.2412e-05 |
| E_-2 | 1.8922e-05 | 1.8922e-05 | 2.0047 | 0.0022273 | 2.2412e-05 |
| E_+0 | 1.8922e-05 | 1.8922e-05 | 2.0047 | 0.0022273 | 2.2412e-05 |
| E_+2 | 1.8922e-05 | 1.8922e-05 | 2.0047 | 0.0022273 | 2.2412e-05 |

## Summary and interpretation

- All 16 baseline and t=1 endpoint tensors match the saved Batch-6 coefficients within the declared independent reconstruction tolerance.
- Largest exact finite-mixing slope / central-difference discrepancy: 1.587e-05.
- Maximum exact-derivative t=1/4 full-tensor remainder / baseline norm: 5.671e-05.
- Largest saved Batch-6 t=1 nonproportional residual: 2.235e-04.

A tiny residual of the complete coefficient tensor is not the same as a tiny relative error in its **transverse** component. The report includes both. Near-zero transverse channels are reported as noise-limited rather than assigned spurious power-law exponents.

These results establish numerical convergence for the *specific* controlled splitting path. They are not a symbolic identity for an arbitrary T3 scalar potential. The physical 4-spinor/CG phase dictionary remains unresolved (0/16).
