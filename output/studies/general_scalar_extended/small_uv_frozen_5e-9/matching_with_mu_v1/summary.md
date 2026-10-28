# B/C candidate mass-mixed C5 and matching-convention diagnostic

Historical study: `small_uv_frozen_5e-9`; new run: `matching_with_mu_v1`

**Not a physical or sequential EFT1 C5 prediction.** Uses the saved Matchete one-scale ordered expression with candidate RGBeta scalar endpoint values and historical frozen Yukawas/MF.
Two convention hypotheses: prior identity m=r, and Stage23/26 conditional signed map m=-r for the two direct quartic operators. Neither physically certified.

| Model | Historical fit passes | abs(mu²)/min(diagonal m²) | Identity joint ΔC5 / C5_0 | Conditional-minus joint ΔC5 / C5_0 | Conditional-minus mixing-only |
|---|---|---:|---:|---:|---:|
| B_-1 | yes | 2.11165e-10 | 0.904533 | 0.904533 | 9.97394e-11 |
| C_-1 | NO | 2.11165e-10 | 0.418277 | 0.418277 | 4.67971e-11 |

## Per-model decomposition (all norms divided by ||C5_ref||_F)

| Model | Sign | Background diagonal | Mixing at q=0 | Direct quartics, unmixed | Quartic×mixing interference | Total |
|---|---|---:|---:|---:|---:|---:|
| B_-1 | identity_prior_preview | 3.52929e-16 | 9.97394e-11 | 0.904533 | 4.50967e-12 | 0.904533 |
| B_-1 | stage23_conditional_minus | 3.52929e-16 | 9.97394e-11 | 0.904533 | 4.5098e-12 | 0.904533 |
| C_-1 | identity_prior_preview | 3.5306e-16 | 4.67971e-11 | 0.418277 | 4.44534e-12 | 0.418277 |
| C_-1 | stage23_conditional_minus | 3.5306e-16 | 4.67971e-11 | 0.418277 | 4.44547e-12 | 0.418277 |

## Explicit quadratic-mixing scenario

Additional input-mixing run is compared to the same signed/identity quartic-only endpoint; its RG-evolved quartics and masses may also differ.
| Model | Map | abs(mu²)/min(m²) | Relative C5 shift vs quartics-only |
|---|---|---:|---:|
| B_-1 | identity_prior_preview | 0.0102166 | 0.00510553 |
| B_-1 | stage23_conditional_minus | 0.0102166 | 0.00296236 |
| C_-1 | identity_prior_preview | 0.0102166 | 0.00266721 |
| C_-1 | stage23_conditional_minus | 0.0102166 | 0.0018385 |

## Interpretation boundaries

- The reference C5 omits radiatively regenerated quartics in the zero-input branch, consistent with the Batch9 fixed-order preview. Their endpoint values remain saved upstream.
- Background term includes diagonal scalar masses and lambdaT3 changes. Yukawas/MF remain at saved historical fitted boundary, not RG endpoints.
- The mixed-mass U and eigenmasses are inserted in an **existing** one-scale Matchete coefficient, not a freshly matched general-potential EFT.
- All four decomposition components add as complex matrices, not as their displayed Frobenius norms.
- 1D signed maps are conditionally source-supported; a physical global Matchete/RGBeta field identification and nonzero mixed threshold matching remain unverified.
- Additional five scalar quartics can affect RG evolution but are not directly inserted into this saved one-scale C5 expression.
- Historical results and preceding Batch9 runs remain untouched.
