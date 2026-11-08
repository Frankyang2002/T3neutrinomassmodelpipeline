# Batch 3--4: independent Standard Model Weinberg-RGE reduction

**Status:** PASS (four standalone regression tests). This verifies a necessary
limiting case and does **not** certify the general T3 A(-4) anomalous dimension.

**Independent physics reference:** S. Antusch, M. Drees, J. Kersten,
M. Lindner, M. Ratz, *Neutrino mass operator renormalization revisited*,
Phys. Lett. B **519** (2001) 238--242, Eq. (21),
https://doi.org/10.1016/S0370-2693(01)01127-3.

Their convention gives

```
16π² beta(kappa) = -(3/2)[kappa Ye†Ye + (Ye†Ye)^T kappa]
                  + lambda_Antusch kappa - 3 g2² kappa + 2 Tr(Ye†Ye) kappa
```

The test has one charged-lepton family, no quark Yukawas, and scalar
potential `V=lambda_SM*(H†H)^2`, where `lambda_Antusch=4*lambda_SM`.
Hence

```
16π² beta(kappa)/kappa = -|ye|² + 4*lambda_SM - 3*g2²
```

Numerical synthetic fixture:

| Input | Value | Contribution to 16π² beta(kappa)/kappa |
| --- | ---: | ---: |
| `ye` | `0.30 + 0.10i` | `-0.100000000000000` |
| `lambda_SM` | `0.20` | `+0.800000000000000` |
| `g2` | `0.65` | `-1.267500000000000` |
| `g1` | `0.36` | `0` (cancels in the one-loop SM result) |
| **Total** | | **`-0.567500000000000`** |

The complete beta tensor is proportional to the gauge-invariant Weinberg
input, with numerical residual at roundoff; its Yukawa coefficient is
independent of the phase of `ye`. Both crossed and conjugate-coefficient
terms vanish in this *specific* SM Weinberg limit. Therefore, agreement
**does not resolve** the general crossed-contraction `C` vs `C*` ambiguity.

Earlier, the unstarred crossed term failed a nontrivial gauge-covariance
regression. The starred, four-exchange candidate passes two gauge checks,
but differs from arXiv:2512.16890v1 Eq. (4.85), which explicitly prints
an unstarred `C` in its crossed term. The September 2026 journal publication
exists but the crossed-term equation has not been independently verified
against its rendered version in this run.

## Re-run in project checkout

```
python -m unittest discover -s tests -p "test_batch34_sm_weinberg_reference.py" -v
```

Depends on the previously installed Batch34 source, WardGate diagnostic,
and provisional crossed-candidate files. It modifies none of those files.

**Required before physical A(-4) evolution:** resolve general crossed-term
convention using an independently verified journal equation, source-level
counterterms, or a separate derivation; then test the complete gauge-invariant
operator closure with Chat A and Chat B's certified physical inputs.
