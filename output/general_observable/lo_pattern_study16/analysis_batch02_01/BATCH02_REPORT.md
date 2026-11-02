# Chat A — Batch 2 numerical pattern analysis

Scope: the **existing** all-16 restricted direct-UV output at fixed `mu0=1e7 GeV`, not a new numerical evaluation. No EFT or RGE.

**Certification: 0/16 physically certified.** All masses and cross-family rankings are conditional on the unverified source-to-Weyl dictionary.

## Frozen benchmark comparison

| Model | alpha | m1 (eV) | m2 (eV) | m3 (eV) | channels | charge coherence | lambda t=2 residual |
|---|---:|---:|---:|---:|---:|---:|---:|
| A_-4 | -4 | 0.00047448 | 0.0075827 | 0.0436006 | 1 | 1 | 9.81e-07 |
| A_-2 | -2 | 0.0072148 | 0.011008 | 0.0485602 | 1 | 1 | 9.74e-07 |
| A_+0 | +0 | 0.0072143 | 0.011 | 0.048499 | 1 | 1 | 9.74e-07 |
| B_-3 | -3 | 0.0072182 | 0.010981 | 0.0483324 | 1 | 1 | 1.12e-06 |
| B_-1 | -1 | 0.0045225 | 0.0098034 | 0.0504413 | 1 | 1 | 1.16e-06 |
| B_+1 | +1 | 0.0072158 | 0.010946 | 0.0480743 | 1 | 1 | 1.12e-06 |
| C_-3 | -3 | 0.026947 | 0.028192 | 0.0550407 | 1 | 1 | 1.14e-06 |
| C_-1 | -1 | 0.037497 | 0.038589 | 0.062285 | 1 | 1 | 1.14e-06 |
| C_+1 | +1 | 0.00036549 | 0.0082687 | 0.0477861 | 1 | 1 | 1.12e-06 |
| D_-2 | -2 | 0.0072142 | 0.010998 | 0.0484824 | 1 | 1 | 9.74e-07 |
| D_+0 | +0 | 0.0072147 | 0.011007 | 0.0485499 | 1 | 1 | 9.74e-07 |
| D_+2 | +2 | 0.00047473 | 0.0075817 | 0.0435948 | 1 | 1 | 9.81e-07 |
| E_-4 | -4 | 0.0033194 | 0.0074169 | 0.0384885 | 2 | 1 | 2.17e-07 |
| E_-2 | -2 | 0.0047142 | 0.0086788 | 0.0424695 | 2 | 1 | 2.14e-07 |
| E_+0 | +0 | 0.0046436 | 0.0086486 | 0.0425154 | 2 | 1 | 2.14e-07 |
| E_+2 | +2 | 0.003326 | 0.0074425 | 0.0386336 | 2 | 1 | 2.17e-07 |

`m1,m2,m3` are sorted singular values, **not** experimental normal/inverted ordering labels.

## Family statistics of frozen m3

| Family | count | mean m3 (eV) | min | max |
|---|---:|---:|---:|---:|
| A | 3 | 0.0468866 | 0.0436006 | 0.0485602 |
| B | 3 | 0.0489493 | 0.0480743 | 0.0504413 |
| C | 3 | 0.0550373 | 0.0477861 | 0.062285 |
| D | 3 | 0.0468757 | 0.0435948 | 0.0485499 |
| E | 4 | 0.0405267 | 0.0384885 | 0.0425154 |

## Mechanism and robust numerical interpretation

1. **Yukawa scan:** with fixed spectra, `m(t*y1,y2)=t^2 M11+t M12+M22` and `m(y1,t*y2)=M11+t M12+t^2 M22`. The report checks all three sorted singular values against these exact quadratic identities. It separately reports whether the special topology-only linear law is numerically accurate.
2. **Scalar mixing:** in the restricted input the off-diagonal S1–S2 block is scaled. With B=0 and only topology cross-terms the full mass is odd under a mixing-sign flip, but linear scaling at finite mixing is approximate. The t=2 residual quantifies departures without interpreting them as source bugs.
3. **Heavy mass:** `M -> t M` changes both mass insertions and loop functions. The fitted t=0.5 to t=2 log slope is empirical and is not a universal -1 law.
4. **Charge channels:** `||sum M_Q||/sum ||M_Q||` is a coherent/interfering Frobenius diagnostic, NOT a probabilistic branching fraction. Any remainder outside the source-CG topology charges is included in the denominator. Negative pairwise real overlap identifies destructive matrix interference.
5. **Individually sorted masses do not add:** the mass matrix, not its Takagi spectrum, is additive across internal charge sectors.

## E-family two-channel overlap

- `E_-4`: Q=-2,-1; normalised real Frobenius overlap = +1; charge coherence = 1.
- `E_-2`: Q=-1,0; normalised real Frobenius overlap = +1; charge coherence = 1.
- `E_+0`: Q=0,1; normalised real Frobenius overlap = +1; charge coherence = 1.
- `E_+2`: Q=1,2; normalised real Frobenius overlap = +1; charge coherence = 1.

## Source and audit limits

- No cross-package Matchete/RGBeta spinor, portal, or CG normalisation is promoted.
- All old files are read-only. Original file SHA256 values are recorded in `SOURCE_MANIFEST.json`, and the existing source guard and frozen inputs are checked.
- Cross-model differences can reflect different Yukawa matrices, heavy masses, portal masses, scalar mixing, interference, and any still-unverified source mapping.
- No EFT, threshold matching, running, new fitting, or altered physics input is performed.
