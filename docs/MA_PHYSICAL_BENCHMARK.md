# Step 5: fitted physical Ma benchmark

## Definition

The benchmark holds the gauge couplings, SM Yukawas, scalar parameters, and
high-scale Majorana matrix fixed to `ma_full_running_example.json`. It fits only
the high-scale scotogenic Yukawa matrix $h(\Lambda)$.

| Input | Value |
|---|---:|
| Ordering | Normal |
| $m_1$ | $0.001\ \mathrm{eV}$ |
| Dirac phase $\delta$ | $0$ |
| Majorana phases $\alpha_{21},\alpha_{31}$ | $0,0$ |
| Oscillation inputs | Stored NuFIT 6.0 central values |

Every objective evaluation runs the full chain

$$
h(10^6\ \mathrm{GeV})
\longrightarrow \text{Ma RGEs}
\longrightarrow \text{matching at }10^3\ \mathrm{GeV}
\longrightarrow \text{SMEFT RGEs}
\longrightarrow m_\nu(M_Z).
$$

## Fitted high-scale Yukawa matrix

The fitted point is real because all selected CP phases are zero:

$$
h(10^6\ \mathrm{GeV})=
\begin{pmatrix}
-0.02205573 & 0.01310885 & -0.00778776\\
-0.04140601 & -0.04141922 & 0.04758617\\
0.02626394 & 0.11964873 & 0.12711208
\end{pmatrix}.
$$

The largest entry is $0.1271$ at the UV scale and $0.1322$ at matching, safely
below $\sqrt{4\pi}$.

## Low-energy result

| Observable | Fitted result |
|---|---:|
| $m_1$ | $0.0010000\ \mathrm{eV}$ |
| $m_2$ | $0.00871206\ \mathrm{eV}$ |
| $m_3$ | $0.0501398\ \mathrm{eV}$ |
| $\Delta m_{21}^2$ | $7.4900\times10^{-5}\ \mathrm{eV}^2$ |
| $\Delta m_{31}^2$ | $2.5130\times10^{-3}\ \mathrm{eV}^2$ |
| $\sin^2\theta_{12}$ | $0.3080$ |
| $\sin^2\theta_{23}$ | $0.4700$ |
| $\sin^2\theta_{13}$ | $0.02215$ |

All five observables are inside their stored NuFIT $3\sigma$ intervals. The
relative low-scale $C_5$ residual is $7.9\times10^{-14}$ and the neutrino
Takagi residual is $1.2\times10^{-16}$.

## Running effects at the fitted point

| Comparison with frozen parameters | Relative $C_5$ change |
|---|---:|
| UV only | $0.8649\%$ |
| EFT only | $4.2940\%$ |
| Full UV and EFT | $6.8976\%$ |

The full neutrino-mass Frobenius norm is $0.9310$ times the frozen result.

## Theory checks over the UV trajectory

The code sampled 64 logarithmically spaced scales from $10^6$ to
$10^3\ \mathrm{GeV}$.

| Check | Result |
|---|---|
| Bounded-from-below inequalities | Pass at all samples |
| Positive charged and neutral inert squared masses | Pass at all samples |
| Gauge and Yukawa couplings below $\sqrt{4\pi}$ | Pass at all samples |
| Scalar quartics below $4\pi$ | Pass at all samples |

The smallest bounded-from-below margin is $\lambda_2=0.1731$. The smallest
inert squared mass is $1.4406\times10^6\ \mathrm{GeV}^2$.

## Interpretation

This is a proof-of-existence benchmark, not a statistical prediction. The
oscillation central values are targets, so their near-zero diagnostic pulls
measure fitting accuracy. A later scan must vary masses, scalar couplings,
lightest neutrino mass, and CP phases and then apply dark-matter,
charged-lepton-flavour-violation, collider, and cosmological constraints.
