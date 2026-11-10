# Conventions used in Matchete

We employ the conventions as used in Peskin and Schroeder:

* Metric: The *mostly-minus* signature metric $g_{\mu \nu} = \mathrm{diag}(+1, -1, -1, -1)$ is used

* Levi-Civita tensor: for spacetime indices $\varepsilon^{0 1 2 3} = - \varepsilon_{0 1 2 3} = +1$

* Chiral projectors: $P_L = \frac{1 - \gamma_5}{2}$, $P_R = \frac{1 + \gamma_5}{2}$

* Covariant derivatives: $D_\mu = \partial_\mu - i T^a A^a_\mu$ with Hermitian generators $T^a$

* Trace normalization: for the fundamental representations $\mathrm{Tr}[T^a T^b] = \frac{1}{2} \delta^{ab}$

* Dimensional regularization:
	* Spacetime dimension: $d = 4 - 2 \epsilon$

	* Renormalization scheme: $\overline{\text{MS}}$ or evanescent-free version of $\overline{\text{MS}}$

	* Anti-commuting $\gamma_5$ subject to $\mathrm{tr}[ \gamma^\mu \gamma^\nu \gamma^\rho \gamma^\sigma \gamma_5 ] = - 4 i \varepsilon^{\mu \nu \rho \sigma}$ (Naive Dimensional Regularization). This results in a loss of trace cyclicity for $ \gamma_5 $-odd traces with 6 or more regular gamma matrices.

* Products of gamma matrices: $\Gamma_{\mu_1,... \mu_n} = \gamma_{[\mu_1} ... \gamma_{\mu_n]}$ , where the antisymmetrized indices include a normalization factor $1/n!$. Thus, $\Gamma_{\mu\nu} = (\gamma_\mu \gamma_\nu - \gamma_\nu \gamma_\mu)/2 = -i \sigma_{\mu\nu}$.

* For the conventions used in the UFO interface see Sec. 2.4 of arXiv:2609.XXXXX.