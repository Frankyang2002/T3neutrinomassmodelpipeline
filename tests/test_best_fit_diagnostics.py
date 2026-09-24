import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from Numerical.BestFitDiagnostics import (
    _complex_matrix_payload,
    _log_save_scales,
    load_optimizer_best_parameters,
)


class BestFitDiagnosticsTests(unittest.TestCase):
    def test_load_optimizer_best_parameters(self):
        payload = {
            "status": "Success",
            "best": {
                "chi2": 1.25,
                "ordering": "NO",
                "parameters": {
                    "lambdaT3_real": 0.12,
                    "y1_11": -0.04,
                },
            },
        }

        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "optimizer.json"
            path.write_text(
                json.dumps(payload),
                encoding="utf-8",
            )
            parameters, chi2, ordering = (
                load_optimizer_best_parameters(path)
            )

        self.assertAlmostEqual(chi2, 1.25)
        self.assertEqual(ordering, "NO")
        self.assertAlmostEqual(
            parameters["lambdaT3_real"],
            0.12,
        )

    def test_complex_matrix_payload_roundtrip(self):
        matrix = np.asarray(
            [
                [1.0 + 2.0j, 3.0 - 4.0j],
                [5.0 + 0.5j, -2.0j],
            ],
            dtype=complex,
        )

        payload = _complex_matrix_payload(matrix)
        rebuilt = (
            np.asarray(payload["real"])
            + 1j * np.asarray(payload["imag"])
        )

        np.testing.assert_allclose(rebuilt, matrix)
        np.testing.assert_allclose(
            np.asarray(payload["abs"]),
            np.abs(matrix),
        )

    def test_log_save_scales_include_exact_endpoints(self):
        scales = _log_save_scales(
            7.0e9,
            1.0e3,
            12,
        )

        self.assertEqual(scales.size, 12)
        self.assertEqual(scales[0], 7.0e9)
        self.assertEqual(scales[-1], 1.0e3)
        self.assertTrue(np.all(np.diff(scales) < 0.0))


if __name__ == "__main__":
    unittest.main()
