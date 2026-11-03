# Chat A — Batch 3: analytical explanation of all 16 T3 patterns

**Basis:** frozen restricted direct UV, fixed `mu0=10^7 GeV`. The 2×2 scalar kernel is
derived exactly in the **provisional Weyl input conventions**. Physical cross-package normalization: **0/16 certified**.

## Analytical kernel and its observable consequences

For each equal-charge pair `s∈S1`, `t∈S2` the complex-scalar squared-mass block is
\[A_q=\begin{pmatrix}a_q&\delta_q^*\\\delta_q&b_q\end{pmatrix},\quad\lambda_{q,\pm}=\tfrac12(a_q+b_q)\pm\sqrt{\tfrac14(a_q-b_q)^2+|\delta_q|^2}.\]
The **exact** scalar insertion entering `Z(S1*) × Y(S2)` is
\[K_q(m_F^2)=\delta_q\frac{B_0(0;m_F^2,\lambda_+)-B_0(0;m_F^2,\lambda_-)}{\lambda_+-\lambda_-}.\]
At leading order in mixing, replace the divided difference at `λ±` by the one at `a_q,b_q`. Higher corrections are odd in `δ_q` for a common phase, with relative size `O(|δ_q|² / mass-gap²)` away from degeneracy.
The source-conditional amplitude is
\[(m_\nu)_{ij}=\frac1{16\pi^2}\sum_{q,k}K_q(m_k^2)\big[Z_{i,s_q}(M^\dagger P_k)Y^T_{t_q,j}+(Z_{i,s_q}(M^\dagger P_k)Y^T_{t_q,j})^T\big],\quadP_k=Q_{:k}Q_{:k}^\dagger,\quad M M^\dagger=Q\,\mathrm{diag}(m_k^2)Q^\dagger.\]
This retains the heavy-flavour matrix and its relative phases. It uses neither a Takagi rotation nor the real-scalar Hessian diagonalisation.

## All-16 per-model interpretation

| Model | charges | m1/m3 | m2/m3 | |y1||y2| | |λT3| | exact kernel error | small-mixing error | λ×2 predicted / measured |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| A_-4 | -2 | 0.010882 | 0.17391 | 3.6553e-08 | 0.1 | 8.33e-12 | 1.63e-07 | 9.81e-07 / 9.81e-07 |
| A_-2 | -1 | 0.14857 | 0.2267 | 3.4422e-08 | 0.1 | 7.35e-12 | 1.62e-07 | 9.74e-07 / 9.74e-07 |
| A_+0 | 0 | 0.14875 | 0.22681 | 3.4412e-08 | 0.1 | 8.71e-12 | 1.62e-07 | 9.74e-07 / 9.74e-07 |
| B_-3 | -1 | 0.14935 | 0.22721 | 3.4401e-08 | 0.1 | 6.42e-12 | 1.87e-07 | 1.12e-06 / 1.12e-06 |
| B_-1 | 0 | 0.089659 | 0.19435 | 5.9205e-08 | 0.1 | 1.48e-11 | 1.93e-07 | 1.16e-06 / 1.16e-06 |
| B_+1 | 1 | 0.1501 | 0.22769 | 3.4734e-08 | 0.1 | 7.06e-12 | 1.87e-07 | 1.12e-06 / 1.12e-06 |
| C_-3 | -1 | 0.48958 | 0.51221 | 6.3062e-08 | 0.1 | 5.95e-12 | 1.89e-07 | 1.14e-06 / 1.14e-06 |
| C_-1 | 0 | 0.60203 | 0.61955 | 1.9371e-07 | 0.1 | 8.08e-12 | 1.89e-07 | 1.14e-06 / 1.14e-06 |
| C_+1 | 1 | 0.0076484 | 0.17304 | 6.1457e-08 | 0.1 | 1.13e-11 | 1.86e-07 | 1.12e-06 / 1.12e-06 |
| D_-2 | 0 | 0.1488 | 0.22684 | 3.4393e-08 | 0.1 | 5.87e-12 | 1.62e-07 | 9.74e-07 / 9.74e-07 |
| D_+0 | 1 | 0.1486 | 0.22672 | 3.4421e-08 | 0.1 | 5.4e-12 | 1.62e-07 | 9.74e-07 / 9.74e-07 |
| D_+2 | 2 | 0.01089 | 0.17391 | 3.6553e-08 | 0.1 | 6.9e-12 | 1.63e-07 | 9.81e-07 / 9.81e-07 |
| E_-4 | -1, -2 | 0.086244 | 0.19271 | 5.9152e-08 | 0.05 | 1.04e-11 | 3.62e-08 | 2.17e-07 / 2.17e-07 |
| E_-2 | 0, -1 | 0.111 | 0.20435 | 7.0533e-08 | 0.05 | 1.63e-11 | 3.57e-08 | 2.14e-07 / 2.14e-07 |
| E_+0 | 1, 0 | 0.10922 | 0.20342 | 6.9895e-08 | 0.05 | 2.1e-11 | 3.57e-08 | 2.14e-07 / 2.14e-07 |
| E_+2 | 2, 1 | 0.08609 | 0.19264 | 5.9301e-08 | 0.05 | 1.5e-11 | 3.62e-08 | 2.17e-07 / 2.17e-07 |

### Channel-by-channel explanation

- **A_-4** (Dirac): Q=-2: |δ|=3257.3 GeV², |W_CG|=1.1547, ||M_Q||=0.044258 eV. m₁/m₃=0.01088; single-channel topology; no inter-charge interference.
- **A_-2** (Dirac): Q=-1: |δ|=3257.3 GeV², |W_CG|=1.1547, ||M_Q||=0.050312 eV. m₁/m₃=0.1486; single-channel topology; no inter-charge interference.
- **A_+0** (Dirac): Q=0: |δ|=3257.3 GeV², |W_CG|=1.1547, ||M_Q||=0.050251 eV. m₁/m₃=0.1488; single-channel topology; no inter-charge interference.
- **B_-3** (Dirac): Q=-1: |δ|=3500.1 GeV², |W_CG|=1.1547, ||M_Q||=0.050087 eV. m₁/m₃=0.1493; single-channel topology; no inter-charge interference.
- **B_-1** (Majorana): Q=0: |δ|=3500.1 GeV², |W_CG|=1.1547, ||M_Q||=0.051584 eV. m₁/m₃=0.08966; single-channel topology; no inter-charge interference.
- **B_+1** (Dirac): Q=1: |δ|=3500.1 GeV², |W_CG|=1.1547, ||M_Q||=0.04983 eV. m₁/m₃=0.1501; single-channel topology; no inter-charge interference.
- **C_-3** (Dirac): Q=-1: |δ|=3500.1 GeV², |W_CG|=0.66667, ||M_Q||=0.067457 eV. m₁/m₃=0.4896; single-channel topology; no inter-charge interference.
- **C_-1** (Majorana): Q=0: |δ|=3500.1 GeV², |W_CG|=0.66667, ||M_Q||=0.082308 eV. m₁/m₃=0.602; single-channel topology; no inter-charge interference.
- **C_+1** (Dirac): Q=1: |δ|=3500.1 GeV², |W_CG|=0.66667, ||M_Q||=0.048498 eV. m₁/m₃=0.007648; single-channel topology; no inter-charge interference.
- **D_-2** (Dirac): Q=0: |δ|=3257.3 GeV², |W_CG|=1.1547, ||M_Q||=0.050235 eV. m₁/m₃=0.1488; single-channel topology; no inter-charge interference.
- **D_+0** (Dirac): Q=1: |δ|=3257.3 GeV², |W_CG|=1.1547, ||M_Q||=0.050302 eV. m₁/m₃=0.1486; single-channel topology; no inter-charge interference.
- **D_+2** (Dirac): Q=2: |δ|=3257.3 GeV², |W_CG|=1.1547, ||M_Q||=0.044252 eV. m₁/m₃=0.01089; single-channel topology; no inter-charge interference.
- **E_-4** (Dirac): Q=-1: |δ|=1515.6 GeV², |W_CG|=0.8165, ||M_Q||=0.019668 eV; Q=-2: |δ|=1515.6 GeV², |W_CG|=0.8165, ||M_Q||=0.019668 eV. m₁/m₃=0.08624; two-channel overlap cosine=1.
- **E_-2** (Dirac): Q=0: |δ|=1515.6 GeV², |W_CG|=0.8165, ||M_Q||=0.021801 eV; Q=-1: |δ|=1515.6 GeV², |W_CG|=0.8165, ||M_Q||=0.021801 eV. m₁/m₃=0.111; two-channel overlap cosine=1.
- **E_+0** (Dirac): Q=1: |δ|=1515.6 GeV², |W_CG|=0.8165, ||M_Q||=0.021817 eV; Q=0: |δ|=1515.6 GeV², |W_CG|=0.8165, ||M_Q||=0.021817 eV. m₁/m₃=0.1092; two-channel overlap cosine=1.
- **E_+2** (Dirac): Q=2: |δ|=1515.6 GeV², |W_CG|=0.8165, ||M_Q||=0.019742 eV; Q=1: |δ|=1515.6 GeV², |W_CG|=0.8165, ||M_Q||=0.019742 eV. m₁/m₃=0.08609; two-channel overlap cosine=1.

## Matched-shape comparisons (NOT exact representation symmetries)

| Pair | best amplitude scale | shape error after complex rescaling |
|---|---:|---:|
| A_-4 vs D_+2 | 0.999867 | 2.805e-05 |
| A_-2 vs D_+0 | 0.999796 | 3.311e-05 |
| A_+0 vs D_-2 | 0.999671 | 5.399e-05 |
| B_-3 vs B_+1 | 0.994866 | 0.0008378 |
| C_-3 vs C_+1 | 0.109 | 0.9884 |
| E_-4 vs E_+2 | 1.00375 | 0.0001607 |
| E_-2 vs E_+0 | 1.00071 | 0.001923 |

## What these tests mean

- The *exact kernel error* is an internal numerical regression test of the source-conditional loop, not native physical certification.
- The *small-mixing error* distinguishes a leading-order analytic approximation from the exact 1-loop scalar spectrum.
- A small m₁/m₃ does not mean a forbidden charge channel: flavour misalignment, determinant suppression, or interference may cause it.
- The heavy-mass-group coherence, `||sum_g M_g||/sum_g ||M_g||`, diagnoses cancellation **between distinct heavy singular-mass eigenspaces**. Individual vectors inside a degenerate eigenspace have no invariant meaning.
- For E, the two-channel cosine measures matrix-direction interference. A value near +1 means the charge matrices reinforce; it does **not** imply an extra normalization factor of 2.
- Numerical masses are from individual frozen Yukawa fits. Changing alpha also changes the fit; a family trend is not by itself a hypercharge theorem.
- Source-CG candidate weights follow `T3Y1CG*T3Y2CG*conjugate(T3MixCG)`. The physical four-spinor/portal dictionary remains B1/B2 dependent.
- The code refuses changed scalar B, new Yukawa vertices, source hashes, partial scans, and false certification promotion.

## Open physics

Physical normalisation, spinor phases, native Majorana factors and portal sign must be resolved in B1/B2. After the coordinator supplies their certified dictionary, rerun the original all-16 numerical calculation before asserting physical masses.
