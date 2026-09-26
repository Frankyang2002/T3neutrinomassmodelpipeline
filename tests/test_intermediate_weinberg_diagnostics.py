from __future__ import annotations

import numpy as np
import pytest

from Numerical.IntermediateWeinbergDiagnostics import (
    IntermediateWeinbergRunning,
)


def test_intermediate_weinberg_running_absolute_value() -> None:
    values = np.asarray(
        [
            [[0.0, 1.0j], [1.0j, 2.0]],
            [[0.0, 2.0j], [2.0j, 4.0]],
        ],
        dtype=complex,
    )

    result = IntermediateWeinbergRunning(
        mu_gev=np.asarray([1.0e10, 1.0e9]),
        delta_c5=values,
    )

    np.testing.assert_allclose(
        result.delta_c5_abs,
        np.abs(values),
    )
