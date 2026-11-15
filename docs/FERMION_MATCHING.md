# Integrating out heavy fermions: matching and normalisation

This document explains the fermion matching requirements for the T3 effective-field-theory calculation. The correct physical sequence is

\[
\mathcal L_{\rm UV}(\mathrm{SM},F,S_1,S_2)
\xrightarrow{\text{integrate out all }F}
\mathcal L_{\rm scalar\ EFT}(\mathrm{SM},S_1,S_2).
\]

Both heavy scalar multiplets remain dynamical throughout the intermediate scalar EFT (unless the model has one physical shared scalar). Non-degenerate heavy fermion masses may require separate fermion thresholds; the particle *sector* ordering remains fermions before scalars.

## Matching data required at the fermion threshold

The matching calculation must specify the complete allowed UV interaction basis, physical fermion mass and flavour conventions, light-field basis, effective-operator basis and perturbative order. The effective action contains inherited operators and threshold-induced Wilson coefficients, not just a single triangle contribution.

Native symbolic output from Matchete must be passed through a verified normalisation and coupling map before numerical evaluation. For example, a correctly recorded one-loop Laurent series can still contain ultraviolet poles or unremoved effective-theory soft contributions. It is **not** automatically a finite renormalised Wilson coefficient.

For one-loop renormalised matching, distinguish:

- The full UV amplitude or effective action with the relevant field content.
- The effective-theory contribution in the same operator and mass convention.
- Hard/soft separation, pole cancellation and appropriate local counterterms.
- Coefficient inheritance and the Wilson-operator mixing required at the retained order.
- Independent numerical checks after flavour and scalar-basis conversion.

A result that verifies only a restricted scalar-quartic/fermion-triangle topology, a leading heavy-mass expansion, or a formal Laurent decomposition must be labelled accordingly. Agreement between two implementations of the same restriction does not demonstrate completeness of the general UV model.

## Coefficient and flavour conventions

The final SMEFT lepton-number-violating coefficient is complex symmetric in lepton flavour; in the documented project convention,

\[
C_5=C_5^{\mathsf T},\qquad m_\nu=-\frac{v^2}{2}\,C_5.
\]

Do not infer the normalisation of a raw ordered amplitude from this equation alone. Verify how the original source converts its ordered amplitude into the symmetric coefficient, including relative signs, conjugations, SU(2) factors and any extra flavour structures. Different intermediate operator bases require an explicit matching dictionary.

## Matching status and interpretation

The repository contains native symbolic fermion matching routes, source inventories, independent numerical triangle checks and exploratory subtraction calculations. Several retained research documents explicitly state that their calculations stop short of a complete renormalised, general-interaction numerical matching boundary. The repository's [implementation status](IMPLEMENTATION_STATUS.md) and each source's recorded input/normalisation determine whether a specific prediction is usable.

For scalar evolution following the fermion threshold, continue with [scalar EFT and SMEFT](SCALAR_EFT_AND_WEINBERG.md). Historical matching evidence, including checks that did not establish finite Wilson coefficients, must remain recoverable outside the public documentation tree.
