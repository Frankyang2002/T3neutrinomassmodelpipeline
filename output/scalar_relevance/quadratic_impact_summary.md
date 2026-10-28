# T3 scalar relevance - Batch 3: quadratic mass impact

- Models: **16**; allowed additional bilinears: **8** (6 self, 2 mixed).
- Candidate initial mass coefficient: **0.01 x min(mS1Sq,mS2Sq)** for allowed cases.
- All candidate spectra have positive mass squared: **True**.
- Independent zero-mu2 RG trajectory computed: **no**.
- Matchete/RGBeta physical normalization certified: **no**.
- General C5 with nonzero mixing matched: **no**.

| Model | Kind | Final |mu2| [GeV^2] | Max relative m^2 shift | sin(2theta) |
|---|---|---:|---:|---:|
| A_+0 | self_s1 | 16112.1 | 0.0253 | - |
| A_-2 | self_s2 | 19107.7 | 0.0305 | - |
| A_-4 | none | 0 | 0 | - |
| B_+1 | none | 0 | 0 | - |
| B_-1 | mixed | 15936.3 | 0.000697 | 0.0937 |
| B_-3 | none | 0 | 0 | - |
| C_+1 | none | 0 | 0 | - |
| C_-1 | mixed | 15808.8 | 0.00069 | 0.0929 |
| C_-3 | none | 0 | 0 | - |
| D_+0 | self_s1 | 19590.5 | 0.0237 | - |
| D_+2 | none | 0 | 0 | - |
| D_-2 | self_s2 | 14871.1 | 0.0307 | - |
| E_+0 | self_s1 | 24818.9 | 0.0296 | - |
| E_+2 | none | 0 | 0 | - |
| E_-2 | self_s2 | 23991.1 | 0.0374 | - |
| E_-4 | none | 0 | 0 | - |

## Interpretation

A nonzero quadratic is a tree-level insertion in a scalar propagator; it can enter an otherwise one-loop neutrino-mass diagram at one-loop order. Its omission is NOT justified by perturbative loop counting.
A smaller eigenvalue shift is not a neutrino-mass error estimate: in the mixed B/C models the eigenvector/Yukawa rotation is linear in mu2/(m1^2-m2^2) while eigenvalue shifts begin quadratically away from degeneracy.
Setting an allowed quadratic coefficient to zero is a benchmark boundary assumption. Establish RG generation and nonzero-mixing C5 before claiming a quantitative neutrino-mass error bound.

## What was verified

The saved candidate quadratic spectra were checked against independent closed-form eigenvalues and their zero-mu2 counterparts at fixed diagonal masses. This does NOT re-run native Matchete/RGBeta or certify the physical component basis.
