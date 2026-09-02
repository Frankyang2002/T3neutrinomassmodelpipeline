from __future__ import annotations

import numpy as np

from RGE.running.NumericalWeinbergRGE import (
    SMInitialConditions,
    evolve_weinberg,
    neutrino_mass_matrix,
)


def main() -> int:
    K0 = np.array(
        [
            [1.0e-14, 2.0e-15 + 1.0e-15j, 3.0e-15],
            [2.0e-15 + 1.0e-15j, 1.2e-14, -1.0e-15j],
            [3.0e-15, -1.0e-15j, 8.0e-15],
        ],
        dtype=complex,
    )

    initial = SMInitialConditions(
        gY=0.36,
        g2=0.64,
        g3=1.00,
        lambdaH=0.26,
        ye=np.array([2.9e-6, 6.1e-4, 1.02e-2]),
        yu=np.array([7.0e-6, 3.5e-3, 0.85]),
        yd=np.array([1.5e-5, 3.0e-4, 1.6e-2]),
        K=K0,
    )

    high = 1.0e6
    low = 1.0e2

    down = evolve_weinberg(
        initial,
        high,
        low,
    )

    if not np.allclose(
        down.K,
        down.K.T,
        rtol=1e-10,
        atol=1e-20,
    ):
        print("FAIL: K lost symmetry during running")
        return 1

    reverse_initial = SMInitialConditions(
        gY=down.gY,
        g2=down.g2,
        g3=down.g3,
        lambdaH=down.lambdaH,
        ye=down.ye,
        yu=down.yu,
        yd=down.yd,
        K=down.K,
    )

    back = evolve_weinberg(
        reverse_initial,
        low,
        high,
    )

    relative_K_error = np.max(
        np.abs(back.K - K0)
        / np.maximum(np.abs(K0), 1e-30)
    )

    mass = neutrino_mass_matrix(down.K)

    print("=" * 72)
    print("NUMERICAL WEINBERG RGE REGRESSION")
    print("=" * 72)
    print()
    print(f"High scale: {high:.3e} GeV")
    print(f"Low scale:  {low:.3e} GeV")
    print(f"Solver evaluations: {down.nfev}")
    print()
    print("Low-scale K:")
    print(down.K)
    print()
    print("Low-scale neutrino mass matrix [GeV]:")
    print(mass)
    print()
    print("Round-trip maximum relative K error:")
    print(relative_K_error)
    print()

    if relative_K_error > 1e-5:
        print("FAIL: numerical round-trip error is too large")
        return 1

    print("PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
