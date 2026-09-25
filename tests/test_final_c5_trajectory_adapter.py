from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

import numpy as np

from Numerical.FinalC5TrajectoryAdapter import (
    FinalC5TrajectoryEvaluator,
    evaluate_final_c5_on_trajectory,
)


def test_canonical_callable_delegates_to_trajectory_evaluator() -> None:
    evaluator = FinalC5TrajectoryEvaluator(Path("final_weinberg_coefficient.json"))

    trajectory = object()
    expected = np.eye(3, dtype=complex) * 2.0e-14

    with patch(
        "Numerical.FinalC5TrajectoryAdapter.evaluate_final_c5_on_trajectory",
        return_value=expected,
    ) as mocked:
        result = evaluator({"lambdaT3": 0.01}, trajectory)

    np.testing.assert_allclose(result, expected)
    mocked.assert_called_once_with(
        Path("final_weinberg_coefficient.json"),
        trajectory,
    )
