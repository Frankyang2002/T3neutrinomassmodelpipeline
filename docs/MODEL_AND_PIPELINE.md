# T3 models and effective field theory

This page defines the model conventions and physical calculation used throughout the project. It distinguishes the **required theoretical description** from the calculations currently demonstrated by the available software.

## Field content and electroweak conventions

The ordinary T3 radiative model extends the Standard Model (SM) by heavy fermions \(F\) and heavy scalars \(S_1,S_2\), with a \(\mathbb Z_2\) symmetry under which the new fields are odd and the SM fields are even. The project uses

\[
Q=T_3+Y,\qquad
Y(S_1)=\frac{\alpha}{2},\quad
Y(S_2)=\frac{\alpha+2}{2},\quad
Y(F)=\frac{\alpha+1}{2}.
\]

Some references instead define \(Y_{\mathrm{lit}}=2Y\). This difference must be handled explicitly when comparing charge assignments. The five ordinary electroweak representation families are

| Family | \((d_{S_1},d_{S_2},d_F)\) | Allowed \(\alpha\) |
| --- | --- | --- |
| A | (1,3,2) | −4, −2, 0 |
| B | (2,2,1) | −3, −1, +1 |
| C | (2,2,3) | −3, −1, +1 |
| D | (3,1,2) | −2, 0, +2 |
| E | (3,3,2) | −4, −2, 0, +2 |

These assignments define the ordinary 16-model catalogue. A shared-scalar/scotogenic construction has one physical scalar even when the diagram uses two formal scalar roles; it must not be counted as two independent scalar multiplets.

## General renormalisable interactions

The starting theory is schematically

\[
\mathcal L_{\rm UV}=\mathcal L_{\rm SM}
+\mathcal L_{\rm kin,mass}^{\rm heavy}
+\mathcal L_{\rm Yukawa}^{\rm allowed}
-V_{\rm scalar}^{\rm allowed}.
\]

The word *allowed* means every independent interaction compatible with Lorentz symmetry, electroweak gauge invariance and the imposed \(\mathbb Z_2\). It is not enough to include only the vertices visible in the T3 neutrino-mass diagram.

Relevant interactions can include both ordinary T3 Yukawa couplings and independent conjugate-scalar Yukawas, general allowed heavy-fermion mass terms, independent scalar quartic contractions, Higgs–scalar mixing interactions and additional complex scalar bilinears. The precise set depends on the representation and hypercharge. In particular, B(−1) and C(−1) permit additional scalar structures and a hypercharge-neutral Majorana fermion; separate flavour contributions to the Weinberg coefficient can occur.

A numerical benchmark may choose zero values for some allowed couplings. Such a choice is a **boundary condition**, not a statement that the interaction is forbidden. Any reported calculation should specify these restrictions.

## Required effective-theory sequence

The physical matching order is

\[
\underbrace{\mathcal L_{\rm UV}(\mathrm{SM},F,S_1,S_2)}_{\text{all heavy fields}}
\xrightarrow{\text{integrate out all }F}
\underbrace{\mathcal L_{\rm scalar\ EFT}(\mathrm{SM},S_1,S_2)}_{\text{heavy scalars retained}}
\xrightarrow{\text{integrate out all }S}
\underbrace{\mathcal L_{\rm SMEFT}}_{\text{SM fields}}.
\]

**Above the fermion thresholds:** Evolve all relevant UV couplings and masses. Fermion flavour eigenvalues need not be degenerate; matching must reflect the prescribed physical masses and field basis.

**After fermion decoupling:** The scalar EFT contains the remaining heavy scalar sector and the induced higher-dimensional operators. Evolve all relevant coefficients using the correct field content and perturbative order. A single Weinberg coefficient need not capture all effects in this theory.

**At scalar decoupling:** Match the remaining scalars onto SMEFT, including allowed scalar mixing, inherited operators and the loop contributions needed for the requested precision. Avoid double counting contributions that are already contained in effective operators.

**Below the scalar thresholds:** Evolve the SMEFT coefficients, including the Weinberg operator. In the documented convention,

\[
m_\nu=-\frac{v^2}{2}C_5.
\]

The resulting coefficient is symmetric in flavour indices and must be expressed in a consistent charged-lepton basis. Takagi factorisation gives non-negative Majorana neutrino masses; the PMNS mixing parameters follow after accounting for the charged-lepton diagonalisation.

The documented one-loop RGE convention is

\[
16\pi^2\frac{\mathrm dX}{\mathrm d\ln\mu}=\beta_X^{(1)}.
\]

Apply the loop factor exactly once, in agreement with the convention of each stored or newly generated RGE expression.

## Scope of the available computations

The codebase includes model construction, group factors, symbolic expressions, candidate matching calculations, numerical running and parameter fitting. It **does not establish one verified, fully general, complete-one-loop numerical calculation through both thresholds for every model**. Some available numerical results use restricted interactions or conditional effective-operator prescriptions. General symbolic expressions alone do not certify a complete general numerical computation.

A correct assessment of an individual study needs its interaction inventory, fermion and scalar mass bases, Wilson coefficient normalisations, RGE sources, matching assumptions, input parameters, fit definition and numerical checks. A small fit residual is evidence only within that configuration.

For the current source-level boundary, see [Capabilities and limitations](IMPLEMENTATION_STATUS.md). For numerical inputs and results, see [Results overview](RESULTS_OVERVIEW.md) and [Reproducibility](REPRODUCIBILITY.md).
