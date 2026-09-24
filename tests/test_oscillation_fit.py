from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from Numerical.OscillationFit import (
    OSCILLATION_OBSERVABLES,
    OscillationFitTarget,
    evaluate_mass_matrix_fit,
    evaluate_oscillation_fit,
    prediction_from_mass_matrix,
)


def _mass_matrix() -> np.ndarray:
    # Generic complex symmetric matrix with non-trivial mixing.
    scale = 1.0e-11  # GeV ~ 0.01 eV
    return scale * np.array(
        [
            [2.5, 0.4 + 0.1j, 0.2 - 0.05j],
            [0.4 + 0.1j, 3.5, 0.6 + 0.2j],
            [0.2 - 0.05j, 0.6 + 0.2j, 5.0],
        ],
        dtype=complex,
    )


class OscillationFitTests(unittest.TestCase):
    def test_prediction_contains_canonical_observables(self) -> None:
        prediction = prediction_from_mass_matrix(
            _mass_matrix(),
            ordering="NO",
        )

        values = prediction.vector()

        self.assertEqual(values.shape, (5,))
        self.assertTrue(np.all(np.isfinite(values)))
        self.assertGreaterEqual(prediction.sin2_theta12, 0.0)
        self.assertLessEqual(prediction.sin2_theta12, 1.0)
        self.assertGreaterEqual(prediction.sin2_theta13, 0.0)
        self.assertLessEqual(prediction.sin2_theta13, 1.0)
        self.assertGreaterEqual(prediction.sin2_theta23, 0.0)
        self.assertLessEqual(prediction.sin2_theta23, 1.0)

    def test_target_equal_to_prediction_has_zero_chi2(self) -> None:
        prediction = prediction_from_mass_matrix(
            _mass_matrix(),
            ordering="NO",
        )
        values = prediction.as_dict()

        target = OscillationFitTarget.from_independent_errors(
            ordering="NO",
            central_values=values,
            one_sigma_errors={
                name: max(abs(values[name]) * 0.1, 1.0e-12)
                for name in values
            },
            source="unit-test",
        )

        result = evaluate_oscillation_fit(
            prediction,
            target,
        )

        self.assertAlmostEqual(result.chi2, 0.0)

    def test_one_sigma_single_pull_gives_chi2_one(self) -> None:
        prediction = prediction_from_mass_matrix(
            _mass_matrix(),
            ordering="NO",
        )

        central = prediction.as_dict()
        sigma = {
            name: max(abs(value) * 0.1, 1.0e-12)
            for name, value in central.items()
        }

        # Shift the target by exactly one sigma in one observable.
        central = dict(central)
        central["sin2_theta13"] += sigma["sin2_theta13"]

        target = OscillationFitTarget.from_independent_errors(
            ordering="NO",
            central_values=central,
            one_sigma_errors=sigma,
        )

        result = evaluate_oscillation_fit(
            prediction,
            target,
        )

        self.assertAlmostEqual(
            result.chi2,
            1.0,
            places=12,
        )

    def test_correlated_covariance_is_used(self) -> None:
        prediction = prediction_from_mass_matrix(
            _mass_matrix(),
            ordering="NO",
        )

        names = ("sin2_theta12", "sin2_theta13")
        pred = prediction.vector(names)

        covariance = np.array(
            [
                [0.01**2, 0.5 * 0.01 * 0.02],
                [0.5 * 0.01 * 0.02, 0.02**2],
            ],
            dtype=float,
        )

        target = OscillationFitTarget(
            ordering="NO",
            observable_names=names,
            central_values=pred + np.array([0.01, 0.02]),
            covariance=covariance,
            source="correlated-test",
        ).validated()

        result = evaluate_oscillation_fit(
            prediction,
            target,
        )

        expected_residual = np.array([-0.01, -0.02])
        expected = float(
            expected_residual
            @ np.linalg.inv(covariance)
            @ expected_residual
        )

        self.assertAlmostEqual(
            result.chi2,
            expected,
            places=12,
        )

    def test_json_target_loader(self) -> None:
        prediction = prediction_from_mass_matrix(
            _mass_matrix(),
            ordering="NO",
        )
        central = prediction.as_dict()

        payload = {
            "source": "synthetic-test-release",
            "ordering": "NO",
            "central_values": central,
            "one_sigma_errors": {
                name: max(abs(value) * 0.1, 1.0e-12)
                for name, value in central.items()
            },
        }

        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "target.json"
            path.write_text(
                json.dumps(payload, indent=2),
                encoding="utf-8",
            )

            target = OscillationFitTarget.from_json(
                path
            )

        self.assertEqual(
            target.observable_names,
            OSCILLATION_OBSERVABLES,
        )
        self.assertEqual(
            target.source,
            "synthetic-test-release",
        )

        result = evaluate_mass_matrix_fit(
            _mass_matrix(),
            target,
        )
        self.assertAlmostEqual(result.chi2, 0.0)

    def test_ordering_mismatch_is_rejected(self) -> None:
        prediction = prediction_from_mass_matrix(
            _mass_matrix(),
            ordering="NO",
        )

        target = OscillationFitTarget.from_independent_errors(
            ordering="IO",
            central_values=prediction.as_dict(),
            one_sigma_errors={
                name: max(abs(value) * 0.1, 1.0e-12)
                for name, value in prediction.as_dict().items()
            },
        )

        with self.assertRaisesRegex(
            ValueError,
            "mass orderings differ",
        ):
            evaluate_oscillation_fit(
                prediction,
                target,
            )


if __name__ == "__main__":
    unittest.main()
