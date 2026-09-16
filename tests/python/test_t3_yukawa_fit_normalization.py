from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
PHENO = ROOT / "tests" / "non_pipeline_rge" / "phenomenology"
sys.path.insert(0, str(PHENO))

from T3NeutrinoTarget import build_normal_ordering_target
from T3YukawaFit import fit_three_heavy_balanced, reconstruct_c5


def main() -> None:
    print("=" * 72)
    print("T3 YUKAWA-FIT NORMALIZATION REGRESSION")
    print("=" * 72)

    target = build_normal_ordering_target(
        lightest_mass_ev=0.01,
        vev_gev=246.22,
    )

    loop_factors = np.array(
        [1.0e-7, 1.7e-7, 2.4e-7],
        dtype=complex,
    )

    fit = fit_three_heavy_balanced(
        target.c5_matrix_gev_inv,
        loop_factors,
    )
    reconstructed = reconstruct_c5(
        fit.y1,
        fit.y2,
        loop_factors,
    )

    diff = reconstructed - target.c5_matrix_gev_inv
    max_abs = float(np.max(np.abs(diff)))

    print(f"relative C5 residual = {fit.relative_residual:.3e}")
    print(f"max |C5_reconstructed - C5_target| = {max_abs:.3e} GeV^-1")

    if fit.relative_residual >= 1e-12:
        raise AssertionError("Yukawa fit residual is too large.")

    if not np.allclose(
        reconstructed,
        target.c5_matrix_gev_inv,
        rtol=1e-11,
        atol=1e-28,
    ):
        raise AssertionError("Reconstructed C5 does not equal the physical target C5.")

    a = fit.y1.conj()
    b = fit.y2.conj()
    d = np.diag(loop_factors)
    ordered = a @ d @ b.T
    physical = ordered + ordered.T

    if not np.allclose(
        reconstructed,
        physical,
        rtol=1e-12,
        atol=1e-28,
    ):
        raise AssertionError("reconstruct_c5 is not implementing C5 = A + A^T.")

    print("PASS: C5 = A + A^T")
    print("PASS: fitted Yukawas reconstruct the physical target C5")


if __name__ == "__main__":
    main()
