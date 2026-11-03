# Track B — archive-source pre-run nine-model batch

**Measured:** nine saved native Matchete CG tensors (with all physical charge sectors) and source benchmark couplings; original Matchete action/Dirac indexed mass still pending for these nine. **Certification:** 0/9.

| Model | Internal fermion Q | Takagi masses, descending (eV) | Relative frozen residual | Eligible |
|---|---|---|---:|---|
| A_-4 | -2 | 0.043600644, 0.0075827281, 0.00047448293 | 8.87e-12 | NO |
| A_-2 | -1 | 0.048560209, 0.011008458, 0.0072147535 | 8.53e-12 | NO |
| D_-2 | +0 | 0.048482372, 0.010997754, 0.0072141749 | 7.54e-12 | NO |
| D_+0 | +1 | 0.048549903, 0.011007039, 0.0072146703 | 7.29e-12 | NO |
| D_+2 | +2 | 0.04359482, 0.0075817378, 0.00047473317 | 7.84e-12 | NO |
| E_-4 | -1, -2 | 0.038488455, 0.0074169371, 0.0033193971 | 1.5e-11 | NO |
| E_-2 | +0, -1 | 0.042469478, 0.0086788288, 0.0047141911 | 2.12e-11 | NO |
| E_+0 | +1, +0 | 0.042515357, 0.0086485894, 0.0046435843 | 2.58e-11 | NO |
| E_+2 | +2, +1 | 0.03863356, 0.0074424848, 0.0033259668 | 2.03e-11 | NO |

Every report includes the full scalar quadratic Hessian in the source SU(2) weight basis, both zero-portal and separately labelled conditional identity-portal neutrino matrices, scalar masses, Takagi singular values, and read-only frozen comparison.

The high-precision B0/Feynman checks in these numerical reports were performed using benchmark scalar eigenvalues; adaptive quadrature is required due to a boundary layer in the loop parameter integral. Engine tests include complex and charged synthetic three-family tests, not claimed as nine separate native complex-source measurements.

**Unresolved physical gates:** installed Matchete GammaCC relative phase, individual native Yukawa/H.c./scalar action and indexed Dirac mass declarations, independently certified benchmark-to-native portal coupling dictionary. The reusable Python script runs those nine model-specific numerical checks on the local repository.
