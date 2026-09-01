from __future__ import annotations

import numpy as np

from RGE.NeutrinoDataComparison import compare_to_nufit6_no
from RGE.T3NeutrinoTarget import build_normal_ordering_target
from RGE.T3YukawaFit import fit_three_heavy_balanced, reconstruct_c5


def takagi_observables_real_symmetric(mass_matrix_ev: np.ndarray) -> dict:
    """Recover observables for the CP-conserving target regression.

    The target used here is real symmetric and positive semidefinite, so an
    ordinary eigh is sufficient for this regression only.
    """
    eigenvalues, vectors = np.linalg.eigh(np.real_if_close(mass_matrix_ev))
    order = np.argsort(eigenvalues)
    masses = eigenvalues[order]
    u = vectors[:, order]

    # Column signs are unphysical for |U|.
    pmns_abs = np.abs(u)

    return {
        "PMNSAbs": pmns_abs.tolist(),
        "DeltaM21SqEV2": float(masses[1] ** 2 - masses[0] ** 2),
        "DeltaM31SqEV2": float(masses[2] ** 2 - masses[0] ** 2),
    }


def main() -> None:
    target = build_normal_ordering_target(lightest_mass_ev=0.01)

    assert np.allclose(
        target.mass_matrix_ev,
        target.mass_matrix_ev.T,
        atol=1e-14,
    )

    observables = takagi_observables_real_symmetric(target.mass_matrix_ev)
    comparison = compare_to_nufit6_no(observables)

    assert comparison["AllFiveInside3Sigma"] is True
    assert comparison["Chi2Diagnostic"] < 1e-16

    # Arbitrary non-degenerate smoke-test loop factors in GeV^-1.
    loop_factors = np.array([1.0e-7, 1.7e-7, 2.4e-7], dtype=complex)

    fit = fit_three_heavy_balanced(
        target.c5_matrix_gev_inv,
        loop_factors,
    )
    reconstructed = reconstruct_c5(fit.y1, fit.y2, loop_factors)

    assert fit.relative_residual < 1e-12
    assert np.allclose(
        reconstructed,
        target.c5_matrix_gev_inv,
        rtol=1e-11,
        atol=1e-28,
    )

    print("PASS: NuFIT target matrix reproduces five best-fit observables")
    print("PASS: constructive T3 Yukawa solution reconstructs target C5")
    print(f"PASS: relative C5 residual = {fit.relative_residual:.3e}")
    print(f"balanced max |y1| = {fit.max_abs_y1:.3e}")
    print(f"balanced max |y2| = {fit.max_abs_y2:.3e}")


if __name__ == "__main__":
    main()
