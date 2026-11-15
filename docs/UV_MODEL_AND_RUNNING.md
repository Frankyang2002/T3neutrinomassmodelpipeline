# T3 ultraviolet model and renormalisation-group evolution

This document consolidates the scientific conventions and validation questions underlying the ultraviolet (UV) model. The [model and EFT contract](MODEL_AND_PIPELINE.md) remains authoritative for the intended calculation, while [capabilities](IMPLEMENTATION_STATUS.md) distinguish available code from complete physical verification.

## Field content and representations

For an ordinary T3 model the heavy fields are fermions \(F\) and scalar multiplets \(S_1,S_2\), with the appropriate \(Z_2\) assignment. The electroweak convention is

\[
Q=T_3+Y,\qquad Y(S_1)=\frac{\alpha}{2},\qquad
Y(S_2)=\frac{\alpha+2}{2},\qquad Y(F)=\frac{\alpha+1}{2}.
\]

The 16-model catalogue spans five representation families:

| Family | \((\dim S_1,\dim S_2,\dim F)\) | \(\alpha\) |
|---|---|---|
| A | (1,3,2) | −4, −2, 0 |
| B | (2,2,1) | −3, −1, +1 |
| C | (2,2,3) | −3, −1, +1 |
| D | (3,1,2) | −2, 0, +2 |
| E | (3,3,2) | −4, −2, 0, +2 |

The *general* UV interaction inventory is the set of independent renormalisable Lorentz-, gauge- and \(Z_2\)-invariant terms, not merely the vertices shown in the neutrino-mass topology. This includes additional Yukawa contractions and complex scalar bilinears or quartics where permitted. B(−1) and C(−1) require particular care with scalar mixing and independent flavour contributions.

## RG evolution and conventions

At one loop, when a symbolic exporter supplies the numerator of the beta function,

\[
16\pi^2\frac{dX}{d\ln\mu}=\beta_X^{(1)}.
\]

Apply the factor \(16\pi^2\) exactly once. The coupling dictionary must identify the same operators, Clebsch–Gordan normalisations, complex conjugations and mass/flavour basis on the symbolic and numerical sides. A saved symbolic RGE is not, by itself, a validated numerical input for a general-interaction fit.

## Representation and gauge-consistency checks

The relevant checks include: tensor symmetry and independent contractions; charge and hypercharge conservation; SU(2) generator and invariant-tensor conventions; scalar-field phase/basis transformations; fermion flavour masses; and the correspondence between Matchete and RGBeta operator definitions.

An RGBeta symbolic invariance diagnostic can remain unresolved even when elementary charge sums are consistent. A failed symbolic reduction must **not** be called either a verified Ward-identity violation or a successful Ward-identity check. Registration of a candidate interaction and physical certification are different claims.

## Evidence limits

The repository has ordinary-model cataloguing, extended interaction builders, RGE generators and selected numerical evolution paths. The current documented source does not establish one complete, general-interaction, one-loop-matched neutrino calculation for all 16 models. For a result, retain the exact model, active interaction basis, coefficient source, initial couplings and scale prescription. See [verification](TEST_STRATEGY.md) and [reproducibility](REPRODUCIBILITY.md).

Historical calculation-specific source notes are retained in a reversible external archive when safely removable. Their preservation does not promote an exploratory calculation to a certified physics result.
