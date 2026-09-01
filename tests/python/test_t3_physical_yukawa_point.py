from __future__ import annotations

import json
import tempfile
from pathlib import Path

import numpy as np
import sympy as sp

import T3PhysicalYukawaPoint as point
from RGE.T3NeutrinoTarget import build_normal_ordering_target
from RGE.T3YukawaFit import reconstruct_c5


def main() -> None:
    # This regression isolates the integration logic from Matchete parsing by
    # replacing the loop-factor evaluator with three known model factors.
    loop_factors = np.array(
        [1.1e-7, 1.8e-7, 2.6e-7],
        dtype=complex,
    )

    original = point.evaluate_t3_loop_factors
    point.evaluate_t3_loop_factors = lambda c5_path, config: loop_factors

    try:
        config = {
            "vev_gev": 246.22,
            "t3": {
                "MS1": 1200.0,
                "MS2": 1500.0,
                "MF": [1000.0, 1800.0, 2500.0],
                "lambdaT3": 0.1,
                "hbar": 1.0 / (16.0 * np.pi**2),
                "y1_real": [[0.0] * 3 for _ in range(3)],
                "y1_imag": [[0.0] * 3 for _ in range(3)],
                "y2_real": [[0.0] * 3 for _ in range(3)],
                "y2_imag": [[0.0] * 3 for _ in range(3)],
            },
        }

        updated, diagnostics = point.build_matched_yukawa_config(
            Path("unused.txt"),
            config,
            lightest_mass_ev=0.01,
        )

        y1 = (
            np.asarray(updated["t3"]["y1_real"])
            + 1j * np.asarray(updated["t3"]["y1_imag"])
        )
        y2 = (
            np.asarray(updated["t3"]["y2_real"])
            + 1j * np.asarray(updated["t3"]["y2_imag"])
        )

        reconstructed = reconstruct_c5(
            y1,
            y2,
            loop_factors,
        )

        target = build_normal_ordering_target(
            lightest_mass_ev=0.01,
            vev_gev=246.22,
        )

        assert np.allclose(
            reconstructed,
            target.c5_matrix_gev_inv,
            rtol=1e-11,
            atol=1e-28,
        )
        assert diagnostics["RelativeC5Residual"] < 1e-12
        assert diagnostics["PerturbativeBySqrt4Pi"] is True

        # Original config must not be mutated.
        assert config["t3"]["y1_real"] == [[0.0] * 3 for _ in range(3)]

        print("PASS: physical-point integration reconstructs target C5")
        print("PASS: numerical config receives fitted y1/y2")
        print("PASS: input config remains unchanged")
        print(
            f"PASS: relative residual = "
            f"{diagnostics['RelativeC5Residual']:.3e}"
        )
    finally:
        point.evaluate_t3_loop_factors = original


if __name__ == "__main__":
    main()
