# Workstream B: coupled UV → EFT1 mass benchmark

Native-RGBeta candidate benchmark; read-only display of previously 
integrated JSON trajectories. All listed mass quantities are in GeV².
No RGEs were recomputed and no physical matching is certified.

| Quantity (GeV²) | B: 10⁷ GeV | B: 10⁵ GeV | B: 10³ GeV | C: 10⁷ GeV | C: 10⁵ GeV | C: 10³ GeV |
|---|---:|---:|---:|---:|---:|---:|
| $m_H^2$ | -7,800 | -12,048 | -12,426.8 | -7,800 | -12,084.1 | -12,774.6 |
| $m_1^2$ | 1.30612e+06 | 1.35513e+06 | 1.41193e+06 | 1.30612e+06 | 1.35549e+06 | 1.40593e+06 |
| $m_2^2$ | 1e+06 | 1.03349e+06 | 1.07328e+06 | 1e+06 | 1.03229e+06 | 1.06713e+06 |
| $\mathrm{Re}\,\mu_R$ | 8,253.36 | 13,108.1 | 18,867.7 | 8,253.36 | 13,079.4 | 18,709.6 |
| $\mathrm{Im}\,\mu_R$ | 5,646.42 | 5,937.72 | 6,265.24 | 5,646.42 | 5,925.99 | 6,225.88 |
| $m_-^2$ | 999,674 | 1.03285e+06 | 1.07212e+06 | 999,674 | 1.03166e+06 | 1.06598e+06 |
| $m_+^2$ | 1.30645e+06 | 1.35577e+06 | 1.41309e+06 | 1.30645e+06 | 1.35613e+06 | 1.40707e+06 |

## Endpoint comparison at Q = 10³ GeV

| Quantity (GeV²) | C(-1) − B(-1) |
|---|---:|
| $m_H^2$ | -347.807 |
| $m_1^2$ | -5,995.88 |
| $m_2^2$ | -6,154.4 |
| $\mathrm{Re}\,\mu_R$ | -158.074 |
| $\mathrm{Im}\,\mu_R$ | -39.3581 |
| $m_-^2$ | -6,135 |
| $m_+^2$ | -6,015.28 |

## Numerical diagnostics

- **B_-1**: UV/EFT1 nfev = 212/47; 49 distinct saved scales; tachyonic sampled points = 0; heavy-fermion singular masses at 10⁵ GeV = [100,000, 110,000, 120,000] GeV.
- **C_-1**: UV/EFT1 nfev = 212/47; 49 distinct saved scales; tachyonic sampled points = 0; heavy-fermion singular masses at 10⁵ GeV = [113,020, 124,322, 135,624] GeV.

## Interpretation and limitations

The extra scalar squared-mass eigenvalues are unbroken-electroweak-symmetry quadratic-potential eigenvalues, not physical pole masses. Positivity at sampled scales does not prove the entire trajectory is tachyon-free.
The threshold is an illustrative tree-level projection at 10⁵ GeV; it does not include finite fermion-threshold corrections or Wilson matching. The Matchete/RGBeta physical field and coupling map is uncertified, and the results are not neutrino-mass predictions.

Seven per-quantity PNG figures display B(-1) and C(-1) versus the renormalization scale, with the fermion threshold marked at 10⁵ GeV.
