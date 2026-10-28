# T3 scalar relevance — Batch 4: quadratic RG source audit

- 16 model reports checked; **8** additional bilinears (6 self, 2 mixed).
- Saved UV and EFT1 **candidate RGBeta beta functions** were inspected.
- The eight extra quadratic beta functions have generically nonzero sources at mu2=0 when their source quartics are nonzero.
- When **all additional quartics** are also set zero, each beta-mass source vanishes algebraically at this one-loop order.
- Fixed-slope figures below use the saved EFT1 upper-scale boundary, with mu2 manually set zero; no independent zero-mu2 trajectory was integrated.

| Model | Kind | Source quartics | Estimated |mu2|/min(m²) from one EFT1 log |
|---|---|---:|---:|
| A_+0 | self_s1 | 2 | 0.38750% |
| A_-2 | self_s2 | 2 | 0.44165% |
| A_-4 | none | 0 | — |
| B_+1 | none | 0 | — |
| B_-1 | mixed | 2 | 0.43545% |
| B_-3 | none | 0 | — |
| C_+1 | none | 0 | — |
| C_-1 | mixed | 2 | 0.43402% |
| C_-3 | none | 0 | — |
| D_+0 | self_s1 | 2 | 0.59715% |
| D_+2 | none | 0 | — |
| D_-2 | self_s2 | 2 | 0.46642% |
| E_+0 | self_s1 | 3 | 0.74706% |
| E_+2 | none | 0 | — |
| E_-2 | self_s2 | 3 | 0.73685% |
| E_-4 | none | 0 | — |

## Perturbative interpretation

For mu2(mu_F)=0 and nonzero source quartics at mu_F,

    mu2(mu_S) ≈ beta_mu2^(1)(mu_F;mu2=0) ln(mu_S/mu_F)/(16*pi²).

This is an **illustrative fixed-slope leading-log term**; it neglects coupling running, feedback, and the UV segment. It is not the difference between two complete physical parameter trajectories.

If extra quartics are themselves set zero at the high-scale boundary and generated only at one loop, induced mu2 may first appear at higher iterated order. A full proof of RG closure requires considering the complete coupled beta system, not only beta_mu2.

## Scope limits

- Matchete–RGBeta coupling and scalar-field normalizations are NOT certified.
- The existing quadratic RG trajectories started from a nonzero mu2=0.01*min(m²), not from zero.
- No physical one-loop C5 with nonzero quadratic mixing is computed here.
- A small induced mass parameter is not a validated relative neutrino-mass error.
- Source quartics can affect scalar RG evolution even if absent in the direct one-vertex C5 graph.
