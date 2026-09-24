import json
from pathlib import Path
import tempfile
import unittest

from Numerical.PlotOptimizerResults import generate_plots


class PlotOptimizerResultsTests(unittest.TestCase):
    def test_generates_expected_files(self):
        payload = {
            "status": "Success",
            "seed_source_kind": "optimizer_result",
            "seed_stored_chi2": 10.0,
            "seed_recomputed_chi2": 10.0,
            "active_parameter_count": 2,
            "active_parameters": ["x", "y"],
            "total_pipeline_evaluations": 3,
            "successful_evaluations": 3,
            "failed_evaluations": 0,
            "sensitivity_ranking": [
                {"name": "x", "score": 3.0},
                {"name": "y", "score": 2.0},
            ],
            "best": {
                "chi2": 0.01,
                "parameters": {"x": 0.1, "y": 0.2},
                "prediction": {
                    "delta_m21_sq_ev2": 1.0,
                    "delta_m3l_sq_ev2": 2.0,
                    "sin2_theta12": 0.3,
                    "sin2_theta13": 0.02,
                    "sin2_theta23": 0.5,
                },
                "observable_names": [
                    "delta_m21_sq_ev2",
                    "delta_m3l_sq_ev2",
                    "sin2_theta12",
                    "sin2_theta13",
                    "sin2_theta23",
                ],
                "whitened_residual": [0.1, -0.1, 0.0, 0.0, 0.0],
            },
            "history": [
                {
                    "evaluation": 0,
                    "phase": "least_squares",
                    "status": "Success",
                    "chi2": 10.0,
                    "parameters": {"x": 0.5, "y": 0.5},
                },
                {
                    "evaluation": 1,
                    "phase": "least_squares",
                    "status": "Success",
                    "chi2": 1.0,
                    "parameters": {"x": 0.2, "y": 0.3},
                },
                {
                    "evaluation": 2,
                    "phase": "least_squares",
                    "status": "Success",
                    "chi2": 0.01,
                    "parameters": {"x": 0.1, "y": 0.2},
                },
            ],
        }

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "optimizer.json"
            source.write_text(json.dumps(payload), encoding="utf-8")
            output = root / "plots"

            generate_plots(source, output)

            expected = [
                "chi2_trajectory.png",
                "best_chi2_trajectory.png",
                "final_pulls.png",
                "sensitivity_ranking.png",
                "active_parameter_trajectory.png",
                "summary.txt",
            ]
            for name in expected:
                self.assertTrue((output / name).is_file(), name)


if __name__ == "__main__":
    unittest.main()
