# Sixteen general-scalar representations: shared A(-2) Yukawas

Scenario: `large_uv_frozen_0p5`; frozen UV, evolved EFT1 and SM+Weinberg.
These are **conditional exploratory** calculations, not certified full general-scalar matching.

Successful numerical evaluations: **16/16**.

| Model | Status | chi2 | Δm21² [eV²] | Δm3l² [eV²] | C5 residual vs A(-2) |
|---|---|---:|---:|---:|---:|
| A_+0 | EvaluatedConditional | 0.00545308 | 7.54084e-05 | 0.00251228 | 3.36905e-16 |
| A_-2 | EvaluatedConditional | 0.000132883 | 7.5364e-05 | 0.0025108 | 0 |
| A_-4 | EvaluatedConditional | 11551.1 | 1.93402e-05 | 0.000653052 | 0.00122687 |
| B_+1 | EvaluatedConditional | 2.19297 | 7.61395e-05 | 0.00253663 | 2.16926e-16 |
| B_-1 | EvaluatedConditional | 1.87099e+06 | 0.00097169 | 0.00105499 | 0.568449 |
| B_-3 | EvaluatedConditional | 12189.6 | 1.78668e-05 | 0.000601638 | 0.000975636 |
| C_+1 | EvaluatedConditional | 9256.2 | 2.53798e-05 | 0.000845545 | 1.76652e-16 |
| C_-1 | EvaluatedConditional | 196094 | 0.00049261 | 0.00117357 | 0.346358 |
| C_-3 | EvaluatedConditional | 17823.8 | 5.95561e-06 | 0.000200546 | 0.000975636 |
| D_+0 | EvaluatedConditional | 8.75453e-06 | 7.53715e-05 | 0.00251105 | 2.49064e-16 |
| D_+2 | EvaluatedConditional | 238.673 | 8.33973e-05 | 0.00277843 | 2.3487e-16 |
| D_-2 | EvaluatedConditional | 0.00948981 | 7.54206e-05 | 0.00251269 | 2.48983e-16 |
| E_+0 | EvaluatedConditional | 4108.03 | 4.20669e-05 | 0.00140148 | 1.61104e-16 |
| E_+2 | EvaluatedConditional | 3106.77 | 4.64084e-05 | 0.00154612 | 1.99382e-16 |
| E_-2 | EvaluatedConditional | 4101.43 | 4.20936e-05 | 0.00140238 | 1.74385e-16 |
| E_-4 | EvaluatedConditional | 15511.3 | 1.05219e-05 | 0.000356941 | 0.00164909 |

## Pairwise C5 scaling at 1 TeV

- B_-3_to_C_-3: k=(-0.577350269)+i(0); residual=1.90762e-16; solar ratio=0.333333333; atmospheric ratio=0.333333333.
- B_-1_to_C_-1: k=(-0.545106875)+i(0); residual=0.818218; solar ratio=0.506961716; atmospheric ratio=1.1123929.
- A_-2_to_D_-2: k=(1.00037436)+i(0); residual=2.48983e-16; solar ratio=1.00075111; atmospheric ratio=1.00075111.

A small matrix residual establishes approximate threshold C5 proportionality within the specific implemented schemes; it does not certify representation physics.
Some fits may fail experimental screening even when their forward numerical evaluation succeeds.
Individual failures and full matrices are saved in `models/` JSON files.
