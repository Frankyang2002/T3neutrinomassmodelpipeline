# B/C candidate mass-mixed C5 and matching-convention diagnostic

Historical study: `large_uv_frozen_0p5`; new run: `matching_conventions_v1`

**Not a physical or sequential EFT1 C5 prediction.** Uses the saved Matchete one-scale ordered expression with candidate RGBeta scalar endpoint values and historical frozen Yukawas/MF.
Two convention hypotheses: prior identity m=r, and Stage23/26 conditional signed map m=-r for the two direct quartic operators. Neither physically certified.

| Model | Historical fit passes | abs(mu²)/min(diagonal m²) | Identity joint ΔC5 / C5_0 | Conditional-minus joint ΔC5 / C5_0 | Conditional-minus mixing-only |
|---|---|---:|---:|---:|---:|
| B_-1 | NO | 0.0205559 | 0.770518 | 0.788217 | 0.00865642 |
| C_-1 | yes | 0.0205559 | 0.971368 | 0.993325 | 0.0107985 |

## Per-model decomposition (all norms divided by ||C5_ref||_F)

| Model | Sign | Background diagonal | Mixing at q=0 | Direct quartics, unmixed | Quartic×mixing interference | Total |
|---|---|---:|---:|---:|---:|---:|
| B_-1 | identity_prior_preview | 0.000458131 | 0.00865642 | 0.779164 | 0.000454306 | 0.770518 |
| B_-1 | stage23_conditional_minus | 0.000458131 | 0.00865642 | 0.779164 | 0.000454306 | 0.788217 |
| C_-1 | identity_prior_preview | 0.000462706 | 0.0107985 | 0.98215 | 0.000433954 | 0.971368 |
| C_-1 | stage23_conditional_minus | 0.000462706 | 0.0107985 | 0.98215 | 0.000433954 | 0.993325 |

## Interpretation boundaries

- The reference C5 omits radiatively regenerated quartics in the zero-input branch, consistent with the Batch9 fixed-order preview. Their endpoint values remain saved upstream.
- Background term includes diagonal scalar masses and lambdaT3 changes. Yukawas/MF remain at saved historical fitted boundary, not RG endpoints.
- The mixed-mass U and eigenmasses are inserted in an **existing** one-scale Matchete coefficient, not a freshly matched general-potential EFT.
- All four decomposition components add as complex matrices, not as their displayed Frobenius norms.
- 1D signed maps are conditionally source-supported; a physical global Matchete/RGBeta field identification and nonzero mixed threshold matching remain unverified.
- Additional five scalar quartics can affect RG evolution but are not directly inserted into this saved one-scale C5 expression.
- Historical results and preceding Batch9 runs remain untouched.
