from __future__ import annotations

import sympy as sp

from RGE.NeutrinoMassStage import build_neutrino_mass_matrix, v


def main() -> int:
    k11, k12, k13, k22, k23, k33 = sp.symbols(
        "k11 k12 k13 k22 k23 k33"
    )

    K = sp.Matrix([
        [k11, k12, k13],
        [k12, k22, k23],
        [k13, k23, k33],
    ])

    M = build_neutrino_mass_matrix(K)

    expected = -sp.Rational(1, 2) * v**2 * K
    difference = M - expected
    symmetry_difference = M - M.T

    print("=" * 72)
    print("NEUTRINO MASS MATRIX REGRESSION")
    print("=" * 72)
    print()
    print("Difference from -(v^2/2) C5:")
    print(difference)
    print()
    print("Symmetry check M-M^T:")
    print(symmetry_difference)
    print()

    if any(sp.simplify(x) != 0 for x in difference):
        print("FAIL: mass normalization")
        return 1

    if any(sp.simplify(x) != 0 for x in symmetry_difference):
        print("FAIL: mass matrix symmetry")
        return 1

    print("PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
