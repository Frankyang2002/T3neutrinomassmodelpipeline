# Results and interpretation

This page summarises the scientific questions addressed by the T3 model calculations and distinguishes recorded numerical benchmarks from general predictions. Detailed historical datasets, assumptions and output locations are recorded in [Numerical results and provenance](RESULTS.md).

## Questions investigated

The project studies how the choice of electroweak representation, running between heavy-particle mass thresholds and changes in scalar and Yukawa couplings affect the light-neutrino masses and mixing parameters in the 16 ordinary T3 models.

The November 2026 Honours thesis reports that, in many of the studied models, the neutrino flavour structure is approximately preserved while the overall mass scale changes. The B(−1) and C(−1) cases show distinctive behaviour consistent with additional interactions and independent flavour contributions. The thesis also identifies important effects from Weinberg-operator evolution on the neutrino mass scale. These conclusions concern the benchmarks and approximations investigated in that work; they do not establish a universal behaviour for every allowed interaction choice.

## What the saved numerical comparisons demonstrate

The repository records a **restricted, individually refitted 16-model comparison** with two chosen scalar-coupling scenarios. In that saved study, **32/32 individual fits** met the study's screening criterion. Reusing one reference model's fitted Yukawa matrices without refitting other representations met the same criterion for **10/32** evaluated cases.

This shows that representation-dependent effects and model-specific Yukawa choices are important in that benchmark. The fitted oscillation parameters are not parameter-free predictions: the fit used considerably more Yukawa degrees of freedom than measured observables.

Other saved comparisons examine the relative effect of running in the UV theory, the scalar effective theory and SMEFT, and how results respond to changes in mass scales or coupling inputs. Their numerical interpretation depends on the exact dataset, available source coefficients, threshold convention and whether inputs were refitted.

## Limits on interpretation

The saved results are not all generated from one uniform, fully general one-loop matching calculation. Restricted scalar potentials, conditional Wilson evolution, incomplete source coverage or numerical fitting assumptions can affect a comparison. Missing calculated entries should not be silently filled in, and a run completing successfully does not establish a physical result beyond the specified approximation.

For reproducibility, retain the model choice, fitted Yukawa matrices, scales, numerical configuration, source coefficient files, output manifest and plotting source. See [Capabilities and limitations](IMPLEMENTATION_STATUS.md), [Numerical methods](NUMERICAL.md) and [Detailed saved results](RESULTS.md).
