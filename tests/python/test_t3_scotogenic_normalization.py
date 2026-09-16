from __future__ import annotations

import sympy as sp


def main() -> None:
    print("=" * 72)
    print("T3-B SCOTOGENIC NORMALIZATION REGRESSION")
    print("=" * 72)

    sqrt3 = sp.sqrt(3)

    # Exact T3-B T3MixCG tensor from the Matchete CG registry.
    # Index order: (H_i, H_j, S1_a, S2dag_b), with 1-based registry
    # entries converted to 0-based Python indices.
    T = sp.MutableDenseNDimArray.zeros(2, 2, 2, 2)
    # Mathematica's SparseArray InputForm encodes the first tensor index
    # through the row-pointer array {0,3,6}. The listed triples are therefore
    # (j, a, b), not full (i, j, a, b) coordinates.
    registry_entries = {
        (0, 0, 1, 0):  2 / sqrt3,
        (0, 1, 0, 0): -1 / sqrt3,
        (0, 1, 1, 1):  1 / sqrt3,
        (1, 0, 0, 0): -1 / sqrt3,
        (1, 0, 1, 1):  1 / sqrt3,
        (1, 1, 0, 1): -2 / sqrt3,
    }
    for idx, value in registry_entries.items():
        T[idx] = value

    # Project epsilon convention from the Matchete registry:
    # epsilon_12 = +1, epsilon_21 = -1.
    eps = sp.Matrix([[0, 1], [-1, 0]])

    # At the T3-B alpha=-1 scotogenic point,
    # S1_a = epsilon_ac eta_c^*, while S2dag_b = eta_b^*.
    D = sp.MutableDenseNDimArray.zeros(2, 2, 2, 2)
    for i in range(2):
        for j in range(2):
            for c in range(2):
                for b in range(2):
                    D[i, j, c, b] = sp.simplify(
                        sum(T[i, j, a, b] * eps[a, c] for a in range(2))
                    )

    expected_D = sp.MutableDenseNDimArray.zeros(2, 2, 2, 2)
    for i in range(2):
        for j in range(2):
            for c in range(2):
                for b in range(2):
                    expected_D[i, j, c, b] = -(
                        sp.KroneckerDelta(i, c) * sp.KroneckerDelta(j, b)
                        + sp.KroneckerDelta(i, b) * sp.KroneckerDelta(j, c)
                    ) / sqrt3

    tensor_residuals = [
        sp.simplify(D[i, j, c, b] - expected_D[i, j, c, b])
        for i in range(2)
        for j in range(2)
        for c in range(2)
        for b in range(2)
    ]
    if any(x != 0 for x in tensor_residuals):
        raise AssertionError(
            "T3MixCG contraction does not reduce to the expected "
            "scotogenic doublet tensor."
        )
    print("PASS: contracted T3MixCG tensor")

    # Verify the component quartic:
    #   H_i H_j eta_c^* eta_b^* D_ijcb
    #     = -(2/sqrt(3)) (H . eta^*)^2.
    H1, H2, e1, e2 = sp.symbols("H1 H2 e1 e2")
    H = [H1, H2]
    eta_star = [e1, e2]

    quartic = sp.expand(
        sum(
            H[i] * H[j] * eta_star[c] * eta_star[b] * D[i, j, c, b]
            for i in range(2)
            for j in range(2)
            for c in range(2)
            for b in range(2)
        )
    )
    expected_quartic = sp.expand(
        -(sp.Rational(2, 1) / sqrt3)
        * (H1 * e1 + H2 * e2) ** 2
    )
    quartic_difference = sp.simplify(quartic - expected_quartic)

    print("Quartic contraction difference:")
    print(quartic_difference)
    if quartic_difference != 0:
        raise AssertionError("T3-B quartic component reduction failed.")
    print("PASS: lambdaT3 quartic component normalization")

    # Literature bridge.
    # Project Lagrangian:
    #   L_mix = -(2/sqrt(3)) lambdaT3 (eta^dagger H)^2 + h.c.
    # Scotogenic literature potential:
    #   V_5 = +(lambda5/2) (H^dagger eta)^2 + h.c.
    # hence L_5 = -V_5, so for real couplings
    #   lambdaT3 = sqrt(3)/4 lambda5.
    lambda5 = sp.Symbol("lambda5", real=True)
    lambdaT3 = sqrt3 * lambda5 / 4

    project_component_coeff = sp.simplify(
        -(sp.Rational(2, 1) / sqrt3) * lambdaT3
    )
    literature_lagrangian_coeff = -lambda5 / 2

    coupling_bridge_difference = sp.simplify(
        project_component_coeff - literature_lagrangian_coeff
    )
    print("Quartic bridge difference:")
    print(coupling_bridge_difference)
    if coupling_bridge_difference != 0:
        raise AssertionError("lambdaT3 <-> lambda5 bridge failed.")
    print("PASS: lambdaT3 = sqrt(3)/4 lambda5")

    # Full C5 check in the degenerate-scalar scotogenic limit.
    #
    # Ordered T3-B Matchete coefficient:
    #   A_ij = -(2/sqrt(3)) hbar lambdaT3 M I3 h_i h_j
    #
    # Physical symmetric coefficient for y1=y2=h*:
    #   C5_ij = A_ij + A_ji = 2 A_ij
    #
    # Degenerate loop relation from the T3-B regression:
    #   M I3 = -f/M.
    #
    # hbar = 1/(16*pi^2).
    M, f = sp.symbols("M f", nonzero=True)
    hbar = 1 / (16 * sp.pi**2)
    M_I3 = -f / M

    c5_project = sp.simplify(
        2
        * (-(sp.Rational(2, 1) / sqrt3))
        * hbar
        * lambdaT3
        * M_I3
    )
    c5_literature = sp.simplify(
        lambda5 / (16 * sp.pi**2) * f / M
    )

    c5_difference = sp.simplify(c5_project - c5_literature)
    print("C5 bridge difference:")
    print(c5_difference)
    if c5_difference != 0:
        raise AssertionError("T3-B C5 does not reproduce the scotogenic formula.")
    print("PASS: T3-B physical C5 matches scotogenic literature")

    # Finally verify the VEV-convention bridge:
    #   v174 = v246/sqrt(2)
    # and project mass convention:
    #   Mnu = -(v246^2/2) C5.
    v246 = sp.Symbol("v246", positive=True)
    v174 = v246 / sp.sqrt(2)

    mnu_project = sp.simplify(
        -(v246**2 / 2) * c5_project
    )
    mnu_literature = sp.simplify(
        -(v174**2 * lambda5 / (16 * sp.pi**2)) * f / M
    )

    mass_difference = sp.simplify(mnu_project - mnu_literature)
    print("Neutrino-mass bridge difference:")
    print(mass_difference)
    if mass_difference != 0:
        raise AssertionError(
            "Project neutrino-mass normalization does not reproduce "
            "the scotogenic literature convention."
        )

    print("PASS: v246/sqrt(2) = v174 convention bridge")
    print("PASS: full T3-B scotogenic normalization")


if __name__ == "__main__":
    main()
