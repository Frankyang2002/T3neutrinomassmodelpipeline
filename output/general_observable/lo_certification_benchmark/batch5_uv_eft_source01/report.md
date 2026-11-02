# Batch 5: A_-4 fitted-input UV/EFT numerical diagnostics

**All UV amplitudes are provisional, not physical matching.**

UV frozen source was reproduced before substituting any evolved input.
Reproduction relative matrix difference: 0.000e+00.

| Scalar inputs | UV singular values (eV) | Best-scale residual | Extra portals omitted |
|---|---|---:|---:|
| UV | 0.047047863, 0.0081729993, 0.00051449454 | 0.017615 | 0 |
| F | 0.049853599, 0.0086656384, 0.00054376295 | 0.017837 | 2 |
| S | 0.051249277, 0.0089144025, 0.00055732254 | 0.018092 | 2 |

EFT conditional singular masses per |N_GammaCC| (eV): 0.0276613742, 0.00457629944, 0.000366426698

The one-complex-scale fit is descriptive only; it is NOT a physical GammaCC determination.
The F/S scenarios omit unsupported evolved adjoint portals if present;
their omitted values are preserved in summary.json. The UV scenario
uses the original UV scalar sector, not a common scalar parameter scale.
The S scenario is explicitly hybrid (muF heavy fermions, muS scalars).

UV -> integrate out ALL F at muF -> EFT1(S1,S2) -> integrate out BOTH scalars at muS -> SMEFT.
Unresolved: signed GammaCC, full hard fermion matching, complete Wilson running,
native scalar/portal coupling normalisation, full EFT/UV perturbative consistency.
