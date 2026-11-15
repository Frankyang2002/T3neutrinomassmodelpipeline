# Scalar effective theory, scalar matching and SMEFT

After all heavy fermions are integrated out, the intermediate effective theory contains the Standard Model fields, the heavy scalar sector and the induced effective operators. Its running and matching must use a consistent operator and scalar-mass basis.

\[
\mathrm{UV}
\xrightarrow{\;F\;}
\underbrace{\mathrm{SM}+S_1+S_2+\text{induced operators}}_{\text{scalar EFT}}
\xrightarrow{\;S_1,S_2\;}
\mathrm{SMEFT}.
\]

The arrows denote integrating out the indicated complete heavy sector. The displayed notation presumes two distinct physical scalars; shared-scalar models must not double-count the same field.

## Intermediate-theory operator evolution

Scalar self- and Higgs-portal interactions, permitted mixed bilinears, mass mixing and Wilson coefficients are representation dependent. Their numerical evolution must respect the explicit mass/basis dictionary. An effective operator containing two lepton fields and two scalar fields cannot generally be evolved correctly by running only its overall coupling: anomalous dimensions and potential operator mixing matter.

For the leading Wilson tensor, check lepton-flavour symmetry, scalar-index symmetries and the transformation of every index under basis changes. A seven-term self-class flow or an externally supplied scalar anomalous dimension is evidence for those implemented terms, **not** evidence that all independent allowed operator classes have been included. The complete basis at the chosen dimension/loop order must be established separately.

## Scalar threshold

Integrate out the remaining physical scalars, including permitted scalar mixing. The scalar matching prescription must keep finite contributions, inherited Wilson operators and appropriate subtraction terms without double counting. A direct fixed-order UV neutrino-mass estimate and a sequential EFT estimate must not be compared as the *same physical prediction* when their inputs, coefficient normalisations, matching prescriptions or scales differ.

For a conditional benchmark, report the exact restrictions (for example, omitted extra scalar interactions or a frozen input spectrum) rather than describing the result as the general-model matching.

## Final SMEFT and observables

Below the last scalar threshold, evolve the Weinberg coefficient to the low-energy scale. In the repository's documented convention,

\[
m_\nu=-\frac{v^2}{2}C_5,\qquad C_5=C_5^{\mathsf T}.
\]

Extract physical Majorana neutrino masses using Takagi factorisation and determine mixing with the charged-lepton flavour convention consistently applied. A different convention in an old diagnostic must be translated explicitly rather than silently reused.

## Verification requirements

A trustworthy scalar-to-SMEFT calculation needs agreement of the native scalar mass basis, SU(2) contractions and source coefficient normalisations; matching-scale continuity at the recorded order; numerical stability; and independent checks that the full required operator basis is represented. See [fermion matching](FERMION_MATCHING.md), [numerical verification](NUMERICAL_VALIDATION.md) and [implementation status](IMPLEMENTATION_STATUS.md).
