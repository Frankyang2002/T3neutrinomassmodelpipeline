# Step 2: Ma-model UV RGE validation

## Reference and conventions

The implementation in `MaUVRGE.py` was checked term by term against
A. Merle and M. Platscher, *Running of Radiative Neutrino Masses: The
Scotogenic Model -- REVISITED*, JHEP 11 (2015) 148,
[arXiv:1507.06314v3](https://arxiv.org/abs/1507.06314), Appendix A.

Both the paper and the code use

$$
\mathcal D=(4\pi)^2\mu\frac{d}{d\mu}
=16\pi^2\frac{d}{d\ln\mu}.
$$

The code variable `gY` is the paper's ordinary-hypercharge coupling $g_1$;
it is not GUT normalized. The scalar potential convention is

$$
V=m_H^2H^\dagger H+m_\eta^2\eta^\dagger\eta
+\frac{\lambda_1}{2}(H^\dagger H)^2
+\frac{\lambda_2}{2}(\eta^\dagger\eta)^2
+\lambda_3(H^\dagger H)(\eta^\dagger\eta)
+\lambda_4(H^\dagger\eta)(\eta^\dagger H)
+\frac{\lambda_5}{2}\left[(H^\dagger\eta)^2+\mathrm{h.c.}\right].
$$

## Equation audit

| Sector | Paper equation | Implementation result |
|---|---|---|
| Gauge couplings | (A-1), $b=(7,-3,-7)$ | Exact match |
| $Y_u,Y_d$ | Standard SM form with modified gauge running | Exact match |
| $Y_e$ | (A-2a) | Exact match |
| Scotogenic Yukawa $h$ | (A-2b) | Exact match |
| Majorana matrix $M$ | (A-3) | Exact match |
| $\lambda_1$ | (A-4a) | Exact match |
| $\lambda_2$ | (A-4b) | Exact match |
| $\lambda_3$ | (A-4c) | Exact match |
| $\lambda_4$ | (A-4d) | Exact match |
| $\lambda_5$ | (A-4e) | Exact match; multiplicative |
| $m_H^2$ | (A-5a) | Exact match |
| $m_\eta^2$ | (A-5b) | Match after the matrix-order correction below |

The trace abbreviations also match:

$$
\begin{aligned}
T&=\operatorname{Tr}(Y_e^\dagger Y_e+3Y_u^\dagger Y_u+3Y_d^\dagger Y_d),\\
T_\nu&=\operatorname{Tr}(h^\dagger h),\\
T_4&=\operatorname{Tr}[(Y_e^\dagger Y_e)^2+3(Y_u^\dagger Y_u)^2+3(Y_d^\dagger Y_d)^2],\\
T_{4\nu}&=\operatorname{Tr}[(h^\dagger h)^2],\\
T_{\nu e}&=\operatorname{Tr}(h^\dagger hY_e^\dagger Y_e).
\end{aligned}
$$

## Correction made during the audit

For a diagonal, real Majorana matrix the last term in (A-5b) is

$$
-4\sum_i M_i^2(hh^\dagger)_{ii}.
$$

The code uses the Takagi convention

$$
U^TMU=D,\qquad h_{\rm mass}=U^Th.
$$

Therefore its basis-independent matrix form is

$$
-4\operatorname{Tr}(MM^\dagger hh^\dagger),
$$

not $-4\operatorname{Tr}(M^\dagger Mhh^\dagger)$. The ordering has been
corrected. The two expressions agree when $M$ is real and diagonal, explaining
why the original example produced an apparently correct result.

## Numerical checks

The corrected implementation passes the following focused checks:

1. $\beta_M$ is symmetric and covariant under a unitary heavy-field basis
   transformation.
2. $\beta_{m_\eta^2}$ is invariant under the same transformation.
3. $\operatorname{Tr}(MM^\dagger hh^\dagger)$ reproduces
   $\sum_iM_i^2(h_{\rm mass}h_{\rm mass}^\dagger)_{ii}$ after Takagi
   diagonalisation.
4. $\lambda_5=0$ remains an RGE fixed point.
5. The full example still integrates successfully.

For the supplied nearly diagonal benchmark, the correction changes the final
result only at approximately $3\times10^{-11}$ relative level. It is
nevertheless required for general complex flavour input.

## Scope

This completes validation of the full-theory beta functions used above the
common threshold. The loop matching formula, conjugation conventions, and
relation to the Matchete T3-B coefficient belong to Step 3 and are not claimed
as validated by this audit.
