# B/C candidate mass-mixed C5 and matching-convention diagnostic

Historical study: `small_uv_running_5e-9`; new run: `matching_with_mu_v1`

**Not a physical or sequential EFT1 C5 prediction.** Uses the saved Matchete one-scale ordered expression with candidate RGBeta scalar endpoint values and historical frozen Yukawas/MF.
Two convention hypotheses: prior identity m=r, and Stage23/26 conditional signed map m=-r for the two direct quartic operators. Neither physically certified.

| Model | Historical fit passes | abs(mu²)/min(diagonal m²) | Identity joint ΔC5 / C5_0 | Conditional-minus joint ΔC5 / C5_0 | Conditional-minus mixing-only |
|---|---|---:|---:|---:|---:|
| B_-1 | yes | 3.50394e-10 | 0.59893 | 0.598582 | 4.69568e-10 |
| C_-1 | yes | 2.10159e-10 | 0.375731 | 0.386374 | 1.89463e-10 |

## Per-model decomposition (all norms divided by ||C5_ref||_F)

| Model | Sign | Background diagonal | Mixing at q=0 | Direct quartics, unmixed | Quartic×mixing interference | Total |
|---|---|---:|---:|---:|---:|---:|
| B_-1 | identity_prior_preview | 0.000709874 | 4.69568e-10 | 0.598756 | 3.10251e-11 | 0.59893 |
| B_-1 | stage23_conditional_minus | 0.000709874 | 4.69568e-10 | 0.598756 | 3.10251e-11 | 0.598582 |
| C_-1 | identity_prior_preview | 0.00655668 | 1.89463e-10 | 0.381033 | 2.63272e-11 | 0.375731 |
| C_-1 | stage23_conditional_minus | 0.00655668 | 1.89463e-10 | 0.381033 | 2.63273e-11 | 0.386374 |

## Explicit quadratic-mixing scenario

Additional input-mixing run is compared to the same signed/identity quartic-only endpoint; its RG-evolved quartics and masses may also differ.
| Model | Map | abs(mu²)/min(m²) | Relative C5 shift vs quartics-only |
|---|---|---:|---:|
| B_-1 | identity_prior_preview | 1.92928e-06 | 2.46223e-06 |
| B_-1 | stage23_conditional_minus | 1.92928e-06 | 2.05699e-06 |
| C_-1 | identity_prior_preview | 5.74312e-07 | 6.31447e-07 |
| C_-1 | stage23_conditional_minus | 5.74312e-07 | 4.32589e-07 |

## Interpretation boundaries

- The reference C5 omits radiatively regenerated quartics in the zero-input branch, consistent with the Batch9 fixed-order preview. Their endpoint values remain saved upstream.
- Background term includes diagonal scalar masses and lambdaT3 changes. Yukawas/MF remain at saved historical fitted boundary, not RG endpoints.
- The mixed-mass U and eigenmasses are inserted in an **existing** one-scale Matchete coefficient, not a freshly matched general-potential EFT.
- All four decomposition components add as complex matrices, not as their displayed Frobenius norms.
- 1D signed maps are conditionally source-supported; a physical global Matchete/RGBeta field identification and nonzero mixed threshold matching remain unverified.
- Additional five scalar quartics can affect RG evolution but are not directly inserted into this saved one-scale C5 expression.
- Historical results and preceding Batch9 runs remain untouched.
