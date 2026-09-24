from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from Numerical.OscillationFit import (
    OscillationFitResult,
)
from Numerical.ParameterScan import (
    ScanParameter,
    grid_parameter_points,
    random_parameter_points,
    run_parameter_scan,
    scan_result_to_json_dict,
    write_scan_result,
)


def _fit(
    parameters,
) -> OscillationFitResult:
    x = float(parameters["x"])
    y = float(parameters.get("y", 0.0))

    chi2 = (x - 2.0) ** 2 + y**2

    return OscillationFitResult(
        ordering="NO",
        chi2=chi2,
        observable_names=(
            "sin2_theta12",
            "sin2_theta13",
        ),
        prediction=np.array([x, y], dtype=float),
        target=np.array([2.0, 0.0], dtype=float),
        residual=np.array(
            [x - 2.0, y],
            dtype=float,
        ),
        whitened_residual=np.array(
            [x - 2.0, y],
            dtype=float,
        ),
        source="synthetic",
    )


class ParameterScanTests(unittest.TestCase):
    def test_grid_cartesian_product(self) -> None:
        points = grid_parameter_points(
            {
                "x": [1.0, 2.0],
                "y": [0.0, 1.0, 2.0],
            }
        )

        self.assertEqual(len(points), 6)
        self.assertEqual(
            points[0],
            {"x": 1.0, "y": 0.0},
        )
        self.assertEqual(
            points[-1],
            {"x": 2.0, "y": 2.0},
        )

    def test_random_scan_is_reproducible(self) -> None:
        parameters = (
            ScanParameter(
                "linear",
                0.0,
                1.0,
                "linear",
            ),
            ScanParameter(
                "log",
                1.0e-3,
                1.0e3,
                "log",
            ),
        )

        first = random_parameter_points(
            parameters,
            5,
            seed=12345,
        )
        second = random_parameter_points(
            parameters,
            5,
            seed=12345,
        )

        self.assertEqual(first, second)

        for point in first:
            self.assertGreaterEqual(
                point["linear"],
                0.0,
            )
            self.assertLess(
                point["linear"],
                1.0,
            )
            self.assertGreaterEqual(
                point["log"],
                1.0e-3,
            )
            self.assertLess(
                point["log"],
                1.0e3,
            )

    def test_best_and_accepted_points(self) -> None:
        points = grid_parameter_points(
            {
                "x": [1.0, 2.0, 3.0],
                "y": [0.0],
            }
        )

        result = run_parameter_scan(
            points,
            _fit,
            accepted_chi2_max=0.5,
        )

        self.assertEqual(
            result.best_point.parameters["x"],
            2.0,
        )
        self.assertEqual(
            result.best_point.chi2,
            0.0,
        )
        self.assertEqual(
            len(result.accepted_points),
            1,
        )

    def test_failed_point_is_retained(self) -> None:
        def evaluator(parameters):
            if parameters["x"] < 0.0:
                raise RuntimeError(
                    "unphysical point"
                )
            return _fit(parameters)

        result = run_parameter_scan(
            [
                {"x": 2.0},
                {"x": -1.0},
            ],
            evaluator,
        )

        self.assertEqual(
            len(result.successful_points),
            1,
        )
        self.assertEqual(
            len(result.failed_points),
            1,
        )
        self.assertIn(
            "unphysical point",
            result.failed_points[0].error,
        )

    def test_fail_fast_propagates_error(self) -> None:
        def evaluator(_parameters):
            raise RuntimeError("stop")

        with self.assertRaisesRegex(
            RuntimeError,
            "stop",
        ):
            run_parameter_scan(
                [{"x": 1.0}],
                evaluator,
                fail_fast=True,
            )

    def test_progress_and_checkpoint_callbacks_run_each_point(self) -> None:
        progress = []
        checkpoints = []

        result = run_parameter_scan(
            [{"x": 1.0}, {"x": 2.0}, {"x": 3.0}],
            _fit,
            progress_callback=lambda completed, total, point: progress.append(
                (completed, total, point.index, point.status)
            ),
            checkpoint_callback=lambda partial: checkpoints.append(
                len(partial.points)
            ),
        )

        self.assertEqual(len(result.points), 3)
        self.assertEqual(checkpoints, [1, 2, 3])
        self.assertEqual(
            progress,
            [
                (1, 3, 0, "Success"),
                (2, 3, 1, "Success"),
                (3, 3, 2, "Success"),
            ],
        )

    def test_json_output_contains_best_point(self) -> None:
        result = run_parameter_scan(
            [
                {"x": 1.0},
                {"x": 2.0},
                {"x": 3.0},
            ],
            _fit,
        )

        payload = scan_result_to_json_dict(
            result
        )

        self.assertEqual(
            payload["best_point_index"],
            1,
        )
        self.assertEqual(
            payload["best_chi2"],
            0.0,
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            path = (
                Path(tmpdir)
                / "scan"
                / "result.json"
            )

            written = write_scan_result(
                result,
                path,
            )

            loaded = json.loads(
                written.read_text(
                    encoding="utf-8"
                )
            )

        self.assertEqual(
            loaded["successful_point_count"],
            3,
        )
        self.assertEqual(
            loaded["best_point_index"],
            1,
        )


if __name__ == "__main__":
    unittest.main()
