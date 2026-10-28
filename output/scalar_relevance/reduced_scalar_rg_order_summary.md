# T3 reduced scalar-sector RG generation orders — Batch 5

- Models: **16**; additional quartic slots: **50**; bilinear slots: **8**.
- UV quartic slots sourced immediately with all additional couplings zero: **6**.
- UV additional bilinears first reachable via two iterated one-loop beta insertions: **0**.
- These counts refer to candidate RGBeta beta functions, not complete C5 matching or exact RG closure.

| Model | UV immediate quartics | UV direct C5 quartic orders | UV bilinear order | EFT1 immediate quartics | EFT1 bilinear order |
|---|---:|---|---:|---:|---:|
| A_+0 | 0 | — | — | 0 | — |
| A_-2 | 0 | — | — | 0 | — |
| A_-4 | 0 | — | — | 0 | — |
| B_+1 | 0 | — | — | 0 | — |
| B_-1 | 3 | 1, 1 | — | 0 | — |
| B_-3 | 0 | — | — | 0 | — |
| C_+1 | 0 | — | — | 0 | — |
| C_-1 | 3 | 1, 1 | — | 0 | — |
| C_-3 | 0 | — | — | 0 | — |
| D_+0 | 0 | — | — | 0 | — |
| D_+2 | 0 | — | — | 0 | — |
| D_-2 | 0 | — | — | 0 | — |
| E_+0 | 0 | — | — | 0 | — |
| E_+2 | 0 | — | — | 0 | — |
| E_-2 | 0 | — | — | 0 | — |
| E_-4 | 0 | — | — | 0 | — |

## Interpretation

For a UV boundary where all additional couplings vanish, a first-order induced quartic has leading-log scaling
`delta_lambda ~ beta_lambda^(1) log(mu/muUV)/(16*pi^2)`.
Insertion into a one-loop neutrino-mass diagram is formally of second loop order if the boundary value is zero.
A bilinear generated through a sourced quartic appears after a second beta-function insertion and has a double-log scaling,
`delta_mu2 ~ mS^2 * (log(mu/muUV)/(16*pi^2))^2` times source coefficients.
These are formal orders, not numerical error estimates. Accidental cancellations or large logarithms require separate analysis.

For each of the UV and EFT1 segments, the loop-log factor for a 100:1 scale ratio is `0.029163`.

## Limitations

- Counts use an RGBeta candidate interaction basis, not a certified Matchete/RGBeta physical coupling map.
- Radiative generation order counts iterated one-loop RG insertions, not one-loop matching diagrams.
- A zero one-loop beta at the boundary does not imply exact RG closure or a protecting symmetry.
- Separate UV and EFT1 chains do not implement threshold matching, field redefinitions, or a two-stage integration.
- No physical neutrino mass, C5 correction, or quantitative relative error bar has been computed.
