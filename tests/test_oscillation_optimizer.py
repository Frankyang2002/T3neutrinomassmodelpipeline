import math
from pathlib import Path
import tempfile
import unittest

import numpy as np

from Numerical.OscillationOptimizer import (
    EvaluationRecorder,
    OptimizerParameter,
    estimate_local_sensitivities,
    load_optimizer_best,
    optimise_reduced_least_squares,
    unit_vector_from_parameters,
)


class _ToyFit:
    def __init__(self, residual):
        residual = np.asarray(residual, dtype=float)
        self.whitened_residual = residual
        self.observable_names = ("r1", "r2")
        self.prediction = residual
        self.ordering = "NO"
        self.chi2 = float(np.dot(residual, residual))


class OptimizerTests(unittest.TestCase):
    def test_linear_and_log_unit_roundtrip(self):
        parameters = (
            OptimizerParameter("x", -2.0, 2.0, "linear"),
            OptimizerParameter("y", 1.0e-3, 1.0e1, "log"),
        )
        values = {"x": 0.75, "y": 0.2}
        z = unit_vector_from_parameters(parameters, values)
        rebuilt = {
            parameter.name: parameter.from_unit(value)
            for parameter, value in zip(parameters, z)
        }
        self.assertAlmostEqual(rebuilt["x"], values["x"])
        self.assertAlmostEqual(
            math.log(rebuilt["y"]),
            math.log(values["y"]),
            places=12,
        )

    def test_sensitivity_finds_influential_parameters(self):
        parameters = (
            OptimizerParameter("x", -2.0, 2.0),
            OptimizerParameter("y", -2.0, 2.0),
            OptimizerParameter("z", -2.0, 2.0),
        )

        def evaluator(values):
            return _ToyFit(
                [
                    10.0 * (values["x"] - 0.4),
                    0.5 * (values["y"] + 0.2),
                ]
            )

        recorder = EvaluationRecorder(
            evaluator=evaluator,
            parameters=parameters,
            observable_names=("r1", "r2"),
        )
        seed = unit_vector_from_parameters(
            parameters,
            {"x": 0.0, "y": 0.0, "z": 0.0},
        )
        _, ranked = estimate_local_sensitivities(
            recorder,
            seed,
            step=0.02,
        )

        self.assertEqual(ranked[0].name, "x")
        self.assertGreater(ranked[0].score, ranked[1].score)
        self.assertEqual(ranked[-1].name, "z")
        self.assertAlmostEqual(ranked[-1].score, 0.0, places=10)

    def test_previous_optimizer_best_can_be_warm_started(self):
        payload = {
            "status": "Success",
            "best": {
                "chi2": 12.5,
                "parameters": {
                    "x": 0.25,
                    "y": -0.5,
                },
                "ordering": "NO",
                "prediction": {"dummy": 1.0},
            },
        }

        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "optimizer.json"
            path.write_text(
                __import__("json").dumps(payload),
                encoding="utf-8",
            )
            seed = load_optimizer_best(path)

        self.assertEqual(seed["source_kind"], "optimizer_result")
        self.assertIsNone(seed["index"])
        self.assertAlmostEqual(seed["chi2"], 12.5)
        self.assertAlmostEqual(seed["parameters"]["x"], 0.25)
        self.assertAlmostEqual(seed["parameters"]["y"], -0.5)

    def test_reduced_least_squares_improves_toy_fit(self):
        parameters = (
            OptimizerParameter("x", -2.0, 2.0),
            OptimizerParameter("y", -2.0, 2.0),
            OptimizerParameter("z", -2.0, 2.0),
        )

        def evaluator(values):
            return _ToyFit(
                [
                    values["x"] - 0.4,
                    values["y"] + 0.2,
                ]
            )

        with tempfile.TemporaryDirectory() as temporary:
            recorder = EvaluationRecorder(
                evaluator=evaluator,
                parameters=parameters,
                observable_names=("r1", "r2"),
                checkpoint_path=Path(temporary) / "checkpoint.json",
            )
            seed = unit_vector_from_parameters(
                parameters,
                {"x": 1.5, "y": 1.4, "z": 0.3},
            )

            baseline, ranked = estimate_local_sensitivities(
                recorder,
                seed,
                step=0.02,
            )
            baseline_chi2 = float(np.dot(baseline, baseline))

            solution, active = optimise_reduced_least_squares(
                recorder,
                seed,
                ranked,
                active_count=2,
                max_nfev=40,
            )

            self.assertTrue(solution.success)
            self.assertEqual(set(active), {"x", "y"})
            self.assertIsNotNone(recorder.best)
            self.assertLess(recorder.best["chi2"], baseline_chi2)
            self.assertLess(recorder.best["chi2"], 1.0e-10)
            self.assertTrue(
                (Path(temporary) / "checkpoint.json").is_file()
            )


if __name__ == "__main__":
    unittest.main()
