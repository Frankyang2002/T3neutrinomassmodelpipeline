import json
from pathlib import Path
import tempfile
import unittest

from Numerical.RunningResultFigures import generate_running_result_figures


class RunningResultFiguresTests(unittest.TestCase):
    def _payload(self):
        matrix = {
            "real": [
                [1e-11, 2e-11, 3e-11],
                [2e-11, 4e-11, 5e-11],
                [3e-11, 5e-11, 6e-11],
            ],
            "imag": [
                [0, 0, 0],
                [0, 0, 0],
                [0, 0, 0],
            ],
            "abs": [
                [1e-11, 2e-11, 3e-11],
                [2e-11, 4e-11, 5e-11],
                [3e-11, 5e-11, 6e-11],
            ],
        }
        return {
            "status": "Success",
            "scales_gev": {
                "mu_uv": 1e9,
                "mu_fermion_threshold": 1e7,
                "mu_scalar_threshold": 1e6,
                "mu_low": 1e3,
            },
            "uv_running": {
                "mu_gev": [1e9, 1e8, 1e7],
                "y1_frobenius_norm": [0.5, 0.48, 0.46],
                "y2_frobenius_norm": [0.4, 0.39, 0.38],
                "lambdaT3_abs": [0.12, 0.115, 0.11],
            },
            "intermediate_running": {
                "mu_gev": [1e7, 3e6, 1e6],
                "lambdaH1": [0.20, 0.19, 0.18],
                "lambdaH2": [0.22, 0.21, 0.20],
                "lambda12": [0.10, 0.095, 0.09],
                "lambdaT3_abs": [0.11, 0.105, 0.10],
            },
            "intermediate_direct_weinberg": {
                "mu_gev": [1e7, 3e6, 1e6],
                "delta_c5_abs": [
                    [
                        [0.0, 1e-18, 2e-18],
                        [1e-18, 3e-18, 4e-18],
                        [2e-18, 4e-18, 5e-18],
                    ],
                    [
                        [0.0, 2e-18, 4e-18],
                        [2e-18, 6e-18, 8e-18],
                        [4e-18, 8e-18, 1e-17],
                    ],
                    [
                        [0.0, 3e-18, 6e-18],
                        [3e-18, 9e-18, 1.2e-17],
                        [6e-18, 1.2e-17, 1.5e-17],
                    ],
                ],
            },
            "final_running": {
                "mu_gev": [1e6, 1e4, 1e3],
                "masses_ev": [
                    [0.006, 0.014, 0.075],
                    [0.005, 0.011, 0.06],
                    [0.004, 0.009, 0.05],
                ],
                "delta_m21_sq_ev2": [1.6e-4, 9.6e-5, 7.5e-5],
                "delta_m3l_sq_ev2": [5.5e-3, 3.5e-3, 2.5e-3],
                "sin2_theta12": [0.31, 0.305, 0.30],
                "sin2_theta13": [0.024, 0.023, 0.022],
                "sin2_theta23": [0.58, 0.57, 0.56],
                "c5_abs": [
                    [
                        [1e-16, 2e-16, 3e-16],
                        [2e-16, 4e-16, 5e-16],
                        [3e-16, 5e-16, 6e-16],
                    ],
                    [
                        [9e-17, 1.8e-16, 2.7e-16],
                        [1.8e-16, 3.6e-16, 4.5e-16],
                        [2.7e-16, 4.5e-16, 5.4e-16],
                    ],
                    [
                        [8e-17, 1.6e-16, 2.4e-16],
                        [1.6e-16, 3.2e-16, 4e-16],
                        [2.4e-16, 4e-16, 4.8e-16],
                    ],
                ],
                "mnu_charged_lepton_abs_ev": [
                    [
                        [0.01, 0.02, 0.03],
                        [0.02, 0.04, 0.05],
                        [0.03, 0.05, 0.06],
                    ],
                    [
                        [0.009, 0.018, 0.027],
                        [0.018, 0.036, 0.045],
                        [0.027, 0.045, 0.054],
                    ],
                    [
                        [0.008, 0.016, 0.024],
                        [0.016, 0.032, 0.040],
                        [0.024, 0.040, 0.048],
                    ],
                ],
            },
            "low_energy": {
                "pmns_abs": [
                    [0.82, 0.55, 0.15],
                    [0.32, 0.66, 0.68],
                    [0.47, 0.51, 0.72],
                ],
                "mnu_charged_lepton_basis_gev": matrix,
            },
        }

    def test_generates_primary_figure_set(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "diagnostics.json"
            source.write_text(
                json.dumps(self._payload()),
                encoding="utf-8",
            )
            output = root / "figures"

            generate_running_result_figures(source, output)

            expected = {
                "uv_coupling_running.png",
                "intermediate_direct_weinberg_running.png",
                "c5_running.png",
                "neutrino_mass_splitting_running.png",
                "neutrino_mixing_running.png",
                "README.txt",
            }

            generated = {path.name for path in output.iterdir()}
            self.assertEqual(generated, expected)


if __name__ == "__main__":
    unittest.main()
