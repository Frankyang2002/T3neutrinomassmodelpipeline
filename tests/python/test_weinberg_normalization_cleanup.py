from __future__ import annotations

import sympy as sp

from RGE.matching.FlavorMatchedC5 import build_flavor_matched_c5
from RGE.stages.NeutrinoMassStage import build_neutrino_mass_matrix


def main() -> int:
    A, y1, y2, v = sp.symbols("A y1 y2 v")

    kernel = A
    Y1 = sp.Matrix([[y1]])
    Y2 = sp.Matrix([[y2]])

    C5 = build_flavor_matched_c5(
        kernel,
        Y1,
        Y2,
        heavy_masses=None,
    )

    ordered = A * sp.conjugate(y1) * sp.conjugate(y2)
    c5_difference = sp.simplify(C5[0, 0] - 2 * ordered)

    m_new = build_neutrino_mass_matrix(C5, vev=v)
    m_old = -v**2 * ordered
    mass_difference = sp.simplify(m_new[0, 0] - m_old)

    print("=" * 72)
    print("WEINBERG NORMALIZATION CLEANUP REGRESSION")
    print("=" * 72)
    print()
    print("Difference C5_physical - 2*A_ordered:")
    print(c5_difference)
    print()
    print("Physical neutrino-mass difference old vs new bookkeeping:")
    print(mass_difference)
    print()

    if c5_difference != 0:
        print("FAIL: physical C5 normalization")
        return 1

    if mass_difference != 0:
        print("FAIL: physical neutrino mass changed")
        return 1

    print("PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
