# Workstream 3 — final conditional numerical closeout

**Status:** `ALL_16_BASELINES_AUDITED_SOME_UNIT_RESPONSES_BLOCKED`

This report describes **exploratory source-backed numerical responses**, not the
complete physical one-loop SMEFT Weinberg coefficient. No physical
`GammaCC` factor or extra HHSS named coefficient is assigned.

## Source-backed completion

- Baseline completed: 16/16 models.
- Missing baseline portal directions incorporated: 14/14.
- Extra named source symbols: 24/24.
- Extra unit-reference directions: 24 total; 14 valid conditional projected responses, 10 Ward-rejected, 0 other blocked.

## Definition of the sensitivity numbers

For the saved baseline matrix $C_{5,0}$ and a unit-reference invariant $r$,
the code records $A_r=\partial C_5/\partial\operatorname{Re}c_r$
and $B_r=\partial C_5/\partial\operatorname{Im}c_r$ in the declared V3 convention.
For hypothetical $|c_r|\leq\kappa$,

$$\max_{|c_r|\leq\kappa}\frac{\|\operatorname{Re}c_r A_r+\operatorname{Im}c_r B_r\|_F}{\|C_{5,0}\|_F}
=\kappa\frac{\sqrt{\lambda_{\max}(G_r)}}{\|C_{5,0}\|_F},$$

where $(G_r)_{ij}=\operatorname{Re}\operatorname{Tr}(M_i^\dagger M_j)$ for $M=(A_r,B_r)$.
The sum over reference directions is a triangle-inequality upper envelope for
independent hypothetical unit-reference coefficients of that size. It is
**not a physical upper limit** and need not use the correct Matchete/RGBeta
phase or coupling basis. A missing/rejected direction makes the full-sector
envelope unavailable rather than equal to zero.

## Model comparison

| Model | Valid/total extra refs | Largest per-unit response / baseline | Envelope $|c|\leq 0.001$ | Envelope $|c|\leq 0.01$ | Envelope $|c|\leq 0.1$ |
|---|---:|---:|---:|---:|---:|
| A_-4 | 0/0 | 0 | 0 | 0 | 0 |
| A_-2 | 1/2 | 0 | not available | not available | not available |
| A_+0 | 1/2 | 0 | not available | not available | not available |
| B_-3 | 1/1 | 0 | 0 | 0 | 0 |
| B_-1 | 2/4 | 78.34 | not available | not available | not available |
| B_+1 | 1/1 | 0 | 0 | 0 | 0 |
| C_-3 | 1/1 | 0 | 0 | 0 | 0 |
| C_-1 | 2/4 | 82.46 | not available | not available | not available |
| C_+1 | 1/1 | 0 | 0 | 0 | 0 |
| D_-2 | 1/2 | 0 | not available | not available | not available |
| D_+0 | 1/2 | 0 | not available | not available | not available |
| D_+2 | 0/0 | 0 | 0 | 0 | 0 |
| E_-4 | 0/0 | 0 | 0 | 0 | 0 |
| E_-2 | 1/2 | 0 | not available | not available | not available |
| E_+0 | 1/2 | 0 | not available | not available | not available |
| E_+2 | 0/0 | 0 | 0 | 0 | 0 |

These ratios are **dimensionless reference-basis sensitivities**;
a value of 1 does not imply a 100% physical correction.

## Largest valid single-reference sensitivities

| Model | Reference | Unassigned source names | Per-unit relative norm | Reference $|c|$ for 10% baseline response |
|---|---|---|---:|---:|
| C_-1 | reference_invariant_1 [0, 2, 0, 2, 0, 0] | lambdaHHdagS1barS1bar | 82.46 | 0.001213 |
| B_-1 | reference_invariant_1 [0, 2, 0, 0, 2, 0] | lambdaHHdagS2S2 | 78.34 | 0.001276 |
| B_-1 | reference_invariant_1 [0, 2, 0, 2, 0, 0] | lambdaHHdagS1barS1bar | 76.4 | 0.001309 |
| C_-1 | reference_invariant_1 [0, 2, 0, 0, 2, 0] | lambdaHHdagS2S2 | 75.31 | 0.001328 |
| A_-2 | reference_invariant_1 [1, 1, 0, 0, 0, 2] | lambdaT3GeneralQ110002I1 | 0 | n/a |
| A_+0 | reference_invariant_1 [1, 1, 0, 2, 0, 0] | lambdaT3GeneralQ110200I1 | 0 | n/a |
| B_-3 | reference_invariant_1 [0, 2, 0, 0, 0, 2] | lambdaT3GeneralQ020002I1 | 0 | n/a |
| B_+1 | reference_invariant_1 [0, 2, 2, 0, 0, 0] | lambdaT3GeneralQ022000I1 | 0 | n/a |

## Conclusion and limits

The all-16 completed **baseline** portal calculation is an adequate
exploratory numerical deliverable in its stated scope. Sensitivities above
describe how strongly the saved unit-reference directions could contribute
**if independently specified coefficients were supplied**. They cannot
determine the actual relative physical impact without those coefficients.

**Unresolved physical work (outside this numerical closeout):**

- GammaCC-to-Weyl signed normalization
- Matchete/RGBeta reference coefficient normalization, phases, orientation
- special scalar original-to-boundary physical basis maps
- 24 named extra HHSS coefficients and reliable full EFT1 evolution
- full one-loop hard finite/counterterm/operator closure
- earlier independent scalar projection failures in A_-2, D_+0, E_-2, E_+0

**EFT ordering:** UV -> integrate out ALL heavy fermions @muF -> EFT1(H,S1,S2) -> integrate out BOTH heavy scalars @muS -> SMEFT

**Final classification:** numerical Workstream 3 reported with explicit
limitations; neither physical $C_5$ nor complete one-loop matching certified.
