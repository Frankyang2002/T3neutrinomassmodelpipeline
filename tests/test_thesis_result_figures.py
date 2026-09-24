import json
from pathlib import Path
import tempfile
import unittest

from Numerical.ThesisResultFigures import (
    _table_rows,
    generate_thesis_results,
)


class ThesisResultFiguresTests(unittest.TestCase):
    def _payload(self):
        matrix = {
            "real": [[1, 0, 0], [0, 2, 0], [0, 0, 3]],
            "imag": [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
            "abs": [[1, 0, 0], [0, 2, 0], [0, 0, 3]],
        }
        return {
            "status": "Success",
            "stored_best_chi2": 1.0e-6,
            "ordering": "NO",
            "scales_gev": {
                "mu_scalar_threshold": 1.0e6,
                "mu_low": 1.0e3,
            },
            "best_parameters": {
                "lambdaT3_real": 0.1,
                **{
                    f"y1_{i}{j}": 0.01 * (i + j)
                    for i in range(1, 4)
                    for j in range(1, 4)
                },
                **{
                    f"y2_{i}{j}": -0.01 * (i + j)
                    for i in range(1, 4)
                    for j in range(1, 4)
                },
            },
            "threshold_diagnostics": {
                "fermion_singular_masses_gev": [1.0e5, 1.1e5, 1.2e5],
                "scalar_running_masses_gev": [8.0e4, 8.5e4],
            },
            "low_energy": {
                "masses_ev": [0.004, 0.009, 0.05],
                "sum_masses_ev": 0.063,
                "delta_m21_sq_ev2": 7.5e-5,
                "delta_m31_sq_ev2": 2.5e-3,
                "delta_m32_sq_ev2": 2.4e-3,
                "pmns_abs": [
                    [0.82, 0.55, 0.15],
                    [0.32, 0.66, 0.68],
                    [0.47, 0.51, 0.72],
                ],
                "takagi_residual": 1.0e-15,
                "mnu_charged_lepton_basis_gev": matrix,
            },
            "running": {
                "mu_gev": [1.0e6, 1.0e4, 1.0e3],
                "masses_ev": [
                    [0.006, 0.014, 0.075],
                    [0.005, 0.011, 0.06],
                    [0.004, 0.009, 0.05],
                ],
                "c5_abs": [
                    [[[1.0e-16, 2.0e-16, 3.0e-16],
                      [2.0e-16, 4.0e-16, 5.0e-16],
                      [3.0e-16, 5.0e-16, 6.0e-16]]][0],
                    [[[9.0e-17, 1.8e-16, 2.7e-16],
                      [1.8e-16, 3.6e-16, 4.5e-16],
                      [2.7e-16, 4.5e-16, 5.4e-16]]][0],
                    [[[8.0e-17, 1.6e-16, 2.4e-16],
                      [1.6e-16, 3.2e-16, 4.0e-16],
                      [2.4e-16, 4.0e-16, 4.8e-16]]][0],
                ],
            },
        }

    def test_table_contains_expected_sections(self):
        rows = _table_rows(self._payload())
        sections = {row[0] for row in rows}
        self.assertEqual(
            sections,
            {"fit", "parameter", "neutrino", "threshold"},
        )

    def test_generate_expected_outputs(self):
        payload = self._payload()
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "best_fit.json"
            source.write_text(json.dumps(payload), encoding="utf-8")
            output = root / "results"

            generate_thesis_results(source, output)

            expected = [
                "neutrino_mass_running.png",
                "c5_running.png",
                "mnu_matrix_abs_ev.png",
                "pmns_abs.png",
                "neutrino_mass_spectrum.png",
                "best_fit_results.csv",
                "best_fit_results_table.tex",
                "best_fit_results_table.md",
                "figure_captions.txt",
                "summary.txt",
            ]
            for name in expected:
                self.assertTrue((output / name).is_file(), name)


if __name__ == "__main__":
    unittest.main()
