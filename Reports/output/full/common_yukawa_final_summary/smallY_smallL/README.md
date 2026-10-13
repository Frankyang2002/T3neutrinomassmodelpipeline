# Common-Yukawa T3 numerical comparison — final preliminary summary

## Inputs and scope
- 16 T3 representations, common real 3×3 Yukawa matrices Y1 and Y2.
- λT3 = 0.05; MF = diag(100,110,120) TeV; mS1² = 1.306122448979592×10⁶ GeV²; mS2² = 10⁶ GeV² (fixed, as in uploaded common-input checkpoint).
- Model-specific representation-dependent scalar couplings are retained. Thus **not every ordinary parameter is identical** across representations.
- Five NuFIT 6.1 (IC24 with SK-atm, NO) Gaussian-approximated observables; diagonal covariance, not official full Δχ² likelihood.
- Plots and table generated from user-provided least-squares outputs; no fresh model evaluations or independent rerun.

## Results
- 16/16 successful model evaluations.
- Mean χ² = 6563.2459; range 1997.2416 to 14921.8474.
- Previous stochastic best = 6584.3860; least-squares improvement = 0.321%.
- Maximum difference between reported χ² and reconstructed diagonal Gaussian χ² = 6.71e-09.

## By representation class

| Class | n | Mean χ² | Min χ² | Max χ² |
|---|---:|---:|---:|---:|
| (1,3,2) | 3 | 4472.68 | 3157.55 | 5143.80 |
| (2,2,1) | 3 | 5354.63 | 4926.13 | 6096.27 |
| (2,2,3) | 3 | 14652.17 | 14489.04 | 14921.85 |
| (3,1,2) | 3 | 4467.85 | 3154.89 | 5139.25 |
| (3,3,2) | 4 | 4542.49 | 1997.24 | 7211.56 |

## Interpretation
- Atmospheric and solar mass-squared splittings show pronounced representation-dependent differences under shared Yukawa input.
- Mixing-angle predictions are much more tightly clustered, and generally near the Gaussian NuFIT central values.
- The common fit remains poor (mean χ² ≫ 1). This does not rule out individually fitted models; the 16 individually fitted reference benchmarks are a distinct calculation.
- No additional optimisation is needed for this preliminary comparison.

## Reproducibility
- Source: uploaded optimized_predictions.csv and optimized_predictions.json from bounded least-squares run.
- Original fitted benchmarks, stochastic checkpoint, and all earlier results were not modified.
- Plot horizontal shading is ±1 marginal Gaussian σ, not a joint confidence region.
- This package belongs at `Reports/output/full/common_yukawa_final_summary/smallY_smallL/`.
