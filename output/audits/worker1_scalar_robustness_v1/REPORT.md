# Worker 1 scalar quartic robustness — local diagnostic

**Status:** diagnostic, not physical certification.

EFT order: UV → all fermions out → EFT1 with both scalars → both scalars out → SMEFT.

| Model | Raw discrepancy | Best real coefficient | Residual after fit |
|---|---:|---:|---:|
| A_-4 | 1.9604369711180984 | -0.9216565588633524 | 0.38800668487045953 |

## A_-4 blocks

| Block | Best real coefficient | Rescaled residual |
|---|---:|---:|
| S1S1 | -1.0 | 0.0 |
| S1S2 | -0.09763107293781742 | 0.9952226753832595 |
| S2S1 | -0.09763107293781742 | 0.9952226753832595 |
| S2S2 | -1.0000000000000002 | 0.0 |

### Nonzero coupling directions

| Coupling | Axis | Coefficient | Shape residual |
|---|---|---:|---:|
| lambdaH2 | real | -1.0000000000000002 | 2.220446049250313e-16 |
| lambdaH2Adj | real | -1.0000000000000002 | 1.2819751242557095e-16 |
| lambdaT3 | real | -0.09763107293781739 | 0.9952226753832581 |
| lambdaH1 | real | -0.9999999999999998 | 2.220446049250313e-16 |
| lambdaT3 | imag | -0.0976310729378175 | 0.9952226753832584 |
| D_-2 | 1.8848955969735197 | -0.7764157057450806 | 0.6302211134771419 |

## D_-2 blocks

| Block | Best real coefficient | Rescaled residual |
|---|---:|---:|
| S1S1 | -1.0 | 0.0 |
| S1S2 | 0.41178128083801047 | 0.9112827095646059 |
| S2S1 | 0.4117812808380106 | 0.9112827095646059 |
| S2S2 | -0.011865399420561747 | 0.9999296036704738 |

### Nonzero coupling directions

| Coupling | Axis | Coefficient | Shape residual |
|---|---|---:|---:|
| lambdaH1 | real | -1.0000000000000002 | 2.220446049250313e-16 |
| lambdaH1Adj | real | -1.0000000000000004 | 2.398355846665908e-16 |
| lambdaT3 | real | -0.09763107293781748 | 0.9952226753832583 |
| lambdaT3GeneralQ020101I1 | real | 1.0000000000000002 | 9.06493303673679e-17 |
| lambdaH2 | real | -0.9999999999999998 | 2.830524433501838e-16 |
| lambdaT3GeneralQ110002I1 | real | 0.9999999999999998 | 2.220446049250313e-16 |
| lambdaT3 | imag | -0.0976310729378175 | 0.9952226753832584 |
| lambdaT3GeneralQ020101I1 | imag | 1.0000000000000002 | 2.0269806481627062e-16 |
| lambdaT3GeneralQ110002I1 | imag | 0.9999999999999998 | 2.220446049250313e-16 |

## Physical assessment

A global sign/factor is ruled out as a full explanation. The lambdaT3 direction contains non-proportional structures. The native/RGBeta map and full scalar threshold are not certified; no fitted tensor was injected into the solver.
**No certified physical alternative and paired identical-UV IVP runs were supplied; low-energy C5 shape, mixing, and stage fractions cannot be attributed to the tensor discrepancy.**
