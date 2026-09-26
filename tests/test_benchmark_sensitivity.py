from __future__ import annotations

import numpy as np

from Numerical.BenchmarkSensitivity import (
    leading_prefactor_scale,
    scaled_config_payload,
)


def _payload() -> dict:
    return {
        "base_state": {
            "ordinary": {
                "lambdaT3": {"real": 2.0, "imag": -1.0},
                "y1": {
                    "real": np.ones((3, 3)).tolist(),
                    "imag": (2.0 * np.ones((3, 3))).tolist(),
                },
                "y2": {
                    "real": (3.0 * np.ones((3, 3))).tolist(),
                    "imag": (4.0 * np.ones((3, 3))).tolist(),
                },
            }
        }
    }


def test_all_scaling_directions_preserve_leading_prefactor() -> None:
    for direction in ("T1", "T2", "T3"):
        for gamma in (0.95, 1.0, 1.05):
            assert np.isclose(
                leading_prefactor_scale(direction, gamma),
                1.0,
                rtol=1.0e-14,
                atol=1.0e-14,
            )


def test_t1_scales_lambda_and_both_yukawas() -> None:
    gamma = 1.21
    result = scaled_config_payload(_payload(), "T1", gamma)
    ordinary = result["base_state"]["ordinary"]

    assert np.isclose(ordinary["lambdaT3"]["real"], 2.0 * gamma)
    assert np.isclose(ordinary["y1"]["real"][0][0], gamma ** -0.5)
    assert np.isclose(
        ordinary["y2"]["real"][0][0],
        3.0 * gamma ** -0.5,
    )


def test_t2_scales_lambda_and_y2() -> None:
    gamma = 1.1
    result = scaled_config_payload(_payload(), "T2", gamma)
    ordinary = result["base_state"]["ordinary"]

    assert np.isclose(ordinary["lambdaT3"]["real"], 2.0 * gamma)
    assert np.isclose(ordinary["y1"]["real"][0][0], 1.0)
    assert np.isclose(ordinary["y2"]["real"][0][0], 3.0 / gamma)


def test_t3_oppositely_scales_y1_and_y2() -> None:
    gamma = 1.1
    result = scaled_config_payload(_payload(), "T3", gamma)
    ordinary = result["base_state"]["ordinary"]

    assert np.isclose(ordinary["lambdaT3"]["real"], 2.0)
    assert np.isclose(ordinary["y1"]["real"][0][0], gamma)
    assert np.isclose(ordinary["y2"]["real"][0][0], 3.0 / gamma)
