# Batch 17 — constrained UV-running B/C Sobol + optimisation

Model: `B_-1`; run: `b_uv_constrained_pilot_v1`; status: `CandidateConstrainedAccepted`.
Candidate accepted if chi2 <= 5.86, rhoF <= 0.5, and mS(1 TeV) in [0.5, 5.0] TeV.

Sobol: 3 anchors × 8 points each.
Mass-safe matched candidate evaluations: 294; failed numerical evaluations: 0.

Best screened chi2: **0.0028880425**
rhoF: 0.028489
scalar masses at 1 TeV: 1.83888–2.90727 TeV
UV lambdaT3: 8.9975198e-06

| Observable | Value |
|---|---:|
| delta_m21_sq_ev2 | 7.535429979e-05 |
| delta_m3l_sq_ev2 | 0.002511380667 |
| sin2_theta12 | 0.3091107529 |
| sin2_theta13 | 0.022483698 |
| sin2_theta23 | 0.4701154621 |

## Interpretation

- All five NuFIT residuals and mass screens are recomputed at each successful trial.
- An invalid mass screen is never included as a valid chi2 point.
- The sampled UV mass ordering and 1 TeV scalar window are numerical conditions, not full threshold certification.
- The Matchete–RGBeta coupling normalization, one-scale C5 and sequential matching remain conditional.
- No historical or earlier exploratory outputs are overwritten.

