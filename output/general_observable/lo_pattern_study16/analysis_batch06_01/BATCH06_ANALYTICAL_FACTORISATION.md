# Chat A Batch 6 — CG-weighted heavy-kernel factorization

**Conditional only:** candidate Matchete/RGBeta→Weyl source dictionary; no physical normalization certificate and no EFT calculation. Do not interpret counterfactual masses as fitted or gauge-invariant models.

## Source-derived algebra

The extracted candidate Weyl vertices obey \(Z_{iIs}=\sum_r Y_{1,ir}^{*}C^{(1)}_{rIs}\), \(Y_{jIt}=\sum_sY_{2,js}^{*}C^{(2)}_{sIt}\). The two-scalar kernel and heavy-mass spectral projector give

\[B^{(R)}_{g,rs}=\sum_{q}K_q(M_g^2)\, C^{(1)}_{r,I,s_q}(M^\dagger P_g)_{IJ}C^{(2)}_{s,J,t_q},\]

\[K_q(M^2)=\delta_q\frac{B_0(0;M^2,\lambda_{q+})-B_0(0;M^2,\lambda_{q-})}{\lambda_{q+}-\lambda_{q-}},\]

\[M_\nu={1\over16\pi^2}\sum_g\left(Y_1^* B^{(R)}_gY_2^{*T} + (Y_1^*B^{(R)}_gY_2^{*T})^T\right).\]

Here each B is a **3×3 source-generation coefficient matrix** (not a neutrino mass matrix). B_g is evaluated before choosing any Yukawa matrices. The heavy group projector P_g includes the *whole degenerate mass eigenspace*.

## Exact factorization criterion

If \(B^{(R)}_g=\kappa_R B^{(0)}_g\) with the **same** complex \(\kappa_R\) for every heavy group g, then \(M_\nu^{(R)}=\kappa_R M_\nu^{(0)}\) for **arbitrary** Y1,Y2 at those scalar and heavy mass parameters. No Yukawa sampling is needed to prove this finite-dimensional statement. A tiny single-spectrum residual does not prove factorization for all mass parameters.

For small scalar mixing, each charge contribution has \(K_q=\delta_qD(M_g^2,a_q,b_q)+O(\delta_q|\delta_q|^2)\). Common unmixed scalar masses and charge-independent mixing magnitudes can make the common loop function factor out of the CG contractions; this is a **conditional sufficient scenario**, not an assertion about every possible T3 scalar potential.

## Scenarios and interpretation

| Scenario | Change | Interpretation |
|---|---|---|
| baseline | Matched |lambdaT3|=0.10, original scalar and heavy spectra | Source-conditional benchmark comparison |
| scalar_common_shift | S1 diagonal +28%, S2 diagonal −22% | Common SU(2)-multiplet masses shifted |
| heavy_generation_tilt | Three heavy Takagi mass groups ×(0.71,1.19,1.64) | Mass eigenvalues changed with invariant projectors |
| large_scalar_mixing | S1–S2 scalar mixing ×20 | Stress nonlinear exact kernel; not a fitted point |
| charge_split_diagnostic | Multiplet component diagonals deliberately split | Outside original restricted SU(2)-degenerate mass ansatz; **not a derived UV potential** |

Every scenario uses the exact two-complex-scalar kernel and leaves vertices/CG tensors unchanged. Counterfactual scenarios do not alter original outputs.

## All-16 Yukawa-independent coefficient comparisons

Reference is `A_-2`; each coefficient fit uses **all nine complex entries in all three heavy groups**, not only the final neutrino matrix. `eps` is the relative Frobenius residual, and `spread` is the maximum relative disagreement of group-specific multipliers with the common multiplier.

| Model | baseline eps | baseline spread | common-shift eps | heavy-tilt eps | large-mix eps | charge-split eps |
|---|---:|---:|---:|---:|---:|---:|
| A_-4 | 6.453e-16 | 5.551e-16 | 6.582e-16 | 7.351e-16 | 7.635e-16 | 6.989e-16 |
| A_-2 | 0.000e+00 | 0.000e+00 | 0.000e+00 | 0.000e+00 | 0.000e+00 | 1.814e-16 |
| A_+0 | 6.358e-16 | 2.220e-16 | 6.592e-16 | 5.801e-16 | 6.599e-16 | 6.692e-16 |
| B_-3 | 4.530e-10 | 5.952e-10 | 4.359e-10 | 2.365e-09 | 1.818e-07 | 1.084e-04 |
| B_-1 | 4.530e-10 | 5.952e-10 | 4.359e-10 | 2.365e-09 | 1.818e-07 | 1.084e-04 |
| B_+1 | 4.530e-10 | 5.952e-10 | 4.359e-10 | 2.365e-09 | 1.818e-07 | 1.084e-04 |
| C_-3 | 4.530e-10 | 5.952e-10 | 4.359e-10 | 2.365e-09 | 1.818e-07 | 1.084e-04 |
| C_-1 | 4.530e-10 | 5.952e-10 | 4.359e-10 | 2.365e-09 | 1.818e-07 | 1.084e-04 |
| C_+1 | 4.530e-10 | 5.952e-10 | 4.359e-10 | 2.365e-09 | 1.818e-07 | 1.084e-04 |
| D_-2 | 5.688e-16 | 3.331e-16 | 5.148e-16 | 3.411e-16 | 5.313e-16 | 2.235e-04 |
| D_+0 | 5.130e-16 | 2.220e-16 | 4.986e-16 | 3.550e-16 | 4.931e-16 | 2.235e-04 |
| D_+2 | 1.492e-15 | 6.661e-16 | 1.518e-15 | 1.252e-15 | 1.534e-15 | 2.235e-04 |
| E_-4 | 3.928e-10 | 5.142e-10 | 3.772e-10 | 2.051e-09 | 1.574e-07 | 2.241e-05 |
| E_-2 | 3.928e-10 | 5.142e-10 | 3.772e-10 | 2.051e-09 | 1.574e-07 | 2.241e-05 |
| E_+0 | 3.928e-10 | 5.142e-10 | 3.772e-10 | 2.051e-09 | 1.574e-07 | 2.241e-05 |
| E_+2 | 3.928e-10 | 5.142e-10 | 3.772e-10 | 2.051e-09 | 1.574e-07 | 2.241e-05 |

## Heavy-group multipliers at the frozen common-mixing spectrum

A single κ across all three heavy groups is the test; κ itself is a **candidate basis-dependent complex coefficient ratio**, not a physical CG normalization.

| Model | Re(κ) | Im(κ) | Re(κ1) | Re(κ2) | Re(κ3) |
|---|---:|---:|---:|---:|---:|
| A_-4 | 1 | 0 | 1 | 1 | 1 |
| A_-2 | 1 | 0 | 1 | 1 | 1 |
| A_+0 | 1 | 0 | 1 | 1 | 1 |
| B_-3 | 1.00000003 | 0 | 1.00000003 | 1.00000003 | 1.00000002 |
| B_-1 | 1.00000003 | 0 | 1.00000003 | 1.00000003 | 1.00000002 |
| B_+1 | 1.00000003 | 0 | 1.00000003 | 1.00000003 | 1.00000002 |
| C_-3 | -0.577350284 | 0 | -0.577350284 | -0.577350284 | -0.577350283 |
| C_-1 | -0.577350284 | 0 | -0.577350284 | -0.577350284 | -0.577350283 |
| C_+1 | -0.577350284 | 0 | -0.577350284 | -0.577350284 | -0.577350283 |
| D_-2 | 1 | 0 | 1 | 1 | 1 |
| D_+0 | 1 | 0 | 1 | 1 | 1 |
| D_+2 | 1 | 0 | 1 | 1 | 1 |
| E_-4 | -1.41421353 | 0 | -1.41421353 | -1.41421353 | -1.41421353 |
| E_-2 | -1.41421353 | 0 | -1.41421353 | -1.41421353 | -1.41421353 |
| E_+0 | -1.41421353 | 0 | -1.41421353 | -1.41421353 | -1.41421353 |
| E_+2 | -1.41421353 | 0 | -1.41421353 | -1.41421353 | -1.41421353 |

## Mirror/symmetry partner tests at B_g level

Unlike the Batch-5 final-matrix fit, this compares all heavy-group coefficient entries and tests the result under spectral deformations.

| Pair | baseline eps | heavy tilt eps | scalar shift eps | charge split eps |
|---|---:|---:|---:|---:|
| A_-4/D_+2 | 1.492e-15 | 1.180e-15 | 1.484e-15 | 2.235e-04 |
| A_-2/D_+0 | 5.130e-16 | 3.550e-16 | 4.986e-16 | 2.235e-04 |
| A_+0/D_-2 | 3.606e-16 | 4.718e-16 | 3.609e-16 | 2.235e-04 |
| B_-3/B_+1 | 1.819e-16 | 1.724e-18 | 1.828e-16 | 1.987e-18 |
| C_-3/C_+1 | 2.454e-16 | 2.144e-16 | 2.457e-16 | 2.662e-16 |
| E_-4/E_+2 | 4.310e-16 | 2.574e-16 | 3.360e-16 | 3.341e-16 |
| E_-2/E_+0 | 1.716e-14 | 2.108e-14 | 1.717e-14 | 1.719e-14 |

## Exact versus small-mixing group tensors

The file `batch06_linear_approximation.csv` records the relative Frobenius error of the first-order scalar-mixing approximation for each model and stress scenario. The exact kernel is used for all other results.

## Same-texture full neutrino matrix checks

For each coefficient set, two reference Yukawa pairs are inserted into the algebraic expression for Mnu. The real frozen source was also reconstructed and checked against Batch 1/3/5. Matrix alignment across representations is a consequence of B_g proportionality, not an independent UV measurement.

## E charge-channel equality test

Baseline E has two constructively aligned charge sectors (Batches 2/3). The deliberate component diagonal split changes their scalar kernel functions; the reported per-charge coefficient norms show if that equality survives. These are **coefficient** norms before choosing a Yukawa texture.

| E model | scenario | charge 1 norm | charge 2 norm | norm ratio |
|---|---|---:|---:|---:|
| E_-4 | baseline | 3.237021e-01 | 3.237021e-01 | 1 |
| E_-2 | baseline | 3.237021e-01 | 3.237021e-01 | 1 |
| E_+0 | baseline | 3.237021e-01 | 3.237021e-01 | 1 |
| E_+2 | baseline | 3.237021e-01 | 3.237021e-01 | 1 |
| E_-4 | large_scalar_mixing | 6.474408e+00 | 6.474408e+00 | 1 |
| E_-2 | large_scalar_mixing | 6.474408e+00 | 6.474408e+00 | 1 |
| E_+0 | large_scalar_mixing | 6.474408e+00 | 6.474408e+00 | 1 |
| E_+2 | large_scalar_mixing | 6.474408e+00 | 6.474408e+00 | 1 |
| E_-4 | charge_split_diagnostic | 3.266137e-01 | 3.218017e-01 | 1.0149531 |
| E_-2 | charge_split_diagnostic | 3.266137e-01 | 3.218017e-01 | 1.0149531 |
| E_+0 | charge_split_diagnostic | 3.266137e-01 | 3.218017e-01 | 1.0149531 |
| E_+2 | charge_split_diagnostic | 3.266137e-01 | 3.218017e-01 | 1.0149531 |

## Interpretation guardrails

- A small epsilon for one frozen spectrum proves only a finite-dimensional identity of the saved candidate Weyl tensors at that spectrum.
- If epsilon is small for **all** mass deformations, that is stronger evidence for structural factorization, not a proof for an arbitrary potential.
- If a deformation raises epsilon, it provides an explicit counterexample to extrapolating the baseline identity beyond its assumptions.
- Different models retain different fermion representations and Dirac/Majorana constructions. No shared physical coupling convention is implied.
- The charge-split diagnostic is deliberately outside the restricted source multiplet-degenerate scalar mass ansatz; its gauge-invariant UV completion has **not** been established.
- Numerical precision at ~1e−14 is not a physical amplitude certificate.
- The 0/16 B1/B2 normalization certificate count is unchanged.
