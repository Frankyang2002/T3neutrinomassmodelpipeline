| Model | $R$ | $k$ | $\epsilon$ | Finite drift | Predicted $\delta(\Delta m^2)$ | Actual solar | Actual atmospheric |
|:--|--:|--:|--:|--:|--:|--:|--:|
| A(-4) | 26.3 | -26.3 | 7.97e-07 | 54.3% | 6.38e+04% | 1.45e+05% | 1.45e+05% |
| B(-3) | 32.8 | -32.8 | 8.47e-07 | 36.2% | 1.01e+05% | 1.79e+05% | 1.79e+05% |
| C(-3) | 34.9 | -34.9 | 7.88e-07 | 29.2% | 1.15e+05% | 1.83e+05% | 1.83e+05% |
| D(-2) | 32.2 | -32.2+0.00058i | 1.9e-06 | 42.9% | 9.74e+04% | 1.9e+05% | 1.9e+05% |
| E(-4) | 16.9 | -16.9 | 8.07e-07 | 97.8% | 2.53e+04% | 9.46e+04% | 9.46e+04% |

$R=\|C_5^{\mathrm{ind}}\|_F/\|C_5^{\mathrm{fin,full}}\|_F$; $k=\mathrm{Tr}[(C_5^{\mathrm{fin,full}})^\dagger C_5^{\mathrm{ind}}]/\|C_5^{\mathrm{fin,full}}\|_F^2$; $\epsilon=\|C_5^{\mathrm{ind}}-kC_5^{\mathrm{fin,full}}\|_F/\|C_5^{\mathrm{ind}}\|_F$.
Predicted correction is $100(|1+k|^2-1)\%$. Actual corrections use `1111` versus `1001`.
Finite drift (%) is $100\|C_5^{\mathrm{fin,full}}-C_5^{\mathrm{fin,1001}}\|_F/\|C_5^{\mathrm{fin,1001}}\|_F$.
The prediction is justified only for small proportionality residual and small finite drift; it neglects differences in SMEFT evolution.
Results are exploratory and inherit the matching and scalar-potential approximations of the saved runs.
