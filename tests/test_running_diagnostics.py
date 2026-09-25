import unittest

import numpy as np

from Numerical.RunningDiagnostics import (
    _log_save_scales,
    mixing_angles_from_pmns_abs,
)


class RunningDiagnosticsTests(unittest.TestCase):
    def test_log_save_scales_are_descending_with_exact_endpoints(self):
        scales = _log_save_scales(1.0e9, 1.0e6, 10)
        self.assertEqual(scales[0], 1.0e9)
        self.assertEqual(scales[-1], 1.0e6)
        self.assertTrue(np.all(np.diff(scales) < 0.0))

    def test_pmns_magnitudes_give_standard_sin2_angles(self):
        pmns = np.asarray(
            [
                [0.82, 0.55, 0.15],
                [0.32, 0.66, 0.68],
                [0.47, 0.51, 0.72],
            ],
            dtype=float,
        )
        s12, s13, s23 = mixing_angles_from_pmns_abs(pmns)
        expected_s13 = pmns[0, 2] ** 2
        denom = 1.0 - expected_s13
        self.assertAlmostEqual(s13, expected_s13)
        self.assertAlmostEqual(s12, pmns[0, 1] ** 2 / denom)
        self.assertAlmostEqual(s23, pmns[1, 2] ** 2 / denom)


if __name__ == "__main__":
    unittest.main()
