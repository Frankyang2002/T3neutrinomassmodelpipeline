# Ma full-running comparison

Input configuration: `ma_full_running_fitted.json`

This benchmark separates running above and below one common matching threshold. The fitted point is a benchmark reconstruction of the selected neutrino target; it is not an independent prediction of those target values.

## Main result

- UV-only relative shift in $C_5$: **0.865%**
- EFT-only relative shift in $C_5$: **4.294%**
- Full relative shift in $C_5$: **6.898%**
- UV–EFT interaction relative to frozen $C_5$: **3.467%**
- Full/frozen neutrino-mass norm ratio: **0.931026**

## Low-energy comparison

| Case | UV | EFT | m1 (eV) | m2 (eV) | m3 (eV) | Δm²21 (eV²) | Δm²31 (eV²) | sin²θ12 | sin²θ23 | sin²θ13 |
|---|:---:|:---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Frozen | no | no | 0.00107434 | 0.00934189 | 0.0538571 | 8.61167e-05 | 0.00289943 | 0.307952 | 0.469707 | 0.0221215 |
| UV only | yes | no | 0.00108336 | 0.00943831 | 0.0543195 | 8.7908e-05 | 0.00294943 | 0.308 | 0.470001 | 0.0221501 |
| EFT only | no | yes | 0.0010282 | 0.00894075 | 0.0515445 | 7.88797e-05 | 0.00265578 | 0.307951 | 0.469706 | 0.0221214 |
| Full | yes | yes | 0.001 | 0.00871206 | 0.0501398 | 7.49e-05 | 0.002513 | 0.308 | 0.47 | 0.02215 |

## Reproducibility and limitations

The calculation uses one-loop Ma RGEs above the threshold, small-$\lambda_5$ matching, and one-loop SM+Weinberg running below it. All heavy fields are integrated out at one common scale; split thresholds and finite threshold corrections are not included.

Generated files:

- `ma_full_running_comparison.json`
- `ma_case_comparison.csv`
- `ma_uv_trajectory.csv`
- `ma_eft_trajectory.csv`
- `ma_running_effects.png`
- `ma_running_effects.pdf`
- `ma_uv_trajectory.png`
- `ma_uv_trajectory.pdf`
- `ma_eft_trajectory.png`
- `ma_eft_trajectory.pdf`
- `ma_running_report.md`
