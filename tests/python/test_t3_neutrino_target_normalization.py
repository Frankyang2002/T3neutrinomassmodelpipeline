from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
PHENO_DIR = ROOT / "tests" / "non_pipeline_rge" / "phenomenology"
sys.path.insert(0, str(PHENO_DIR))

from T3NeutrinoTarget import EV_TO_GEV, build_neutrino_target


def check(ordering: str) -> None:
    vev = 246.22
    target = build_neutrino_target(
        ordering=ordering,
        lightest_mass_ev=0.01,
        vev_gev=vev,
    )

    reconstructed_mass_gev = (
        -0.5 * vev * vev * target.c5_matrix_gev_inv
    )
    target_mass_gev = target.mass_matrix_ev * EV_TO_GEV

    difference = reconstructed_mass_gev - target_mass_gev

    print(
        f"{ordering}: max |M_reconstructed - M_target| [GeV] = "
        f"{np.max(np.abs(difference)):.3e}"
    )

    if not np.allclose(
        reconstructed_mass_gev,
        target_mass_gev,
        rtol=1e-13,
        atol=1e-30,
    ):
        raise AssertionError(
            f"{ordering}: C5 target does not satisfy "
            "M_nu = -(v^2/2) C5"
        )


def main() -> None:
    print("=" * 72)
    print("T3 NEUTRINO TARGET NORMALIZATION REGRESSION")
    print("=" * 72)

    check("NO")
    check("IO")

    print("PASS")


if __name__ == "__main__":
    main()
