from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from Numerical.ScanCLI import (
    _set_bound_value,
    _validate_rgbeta_payload_compatibility,
    build_scan_points,
    build_uv_state_from_config,
    load_scan_config,
)
from Numerical.State import T3Representation


def _write_fixture_tree(root: Path) -> Path:
    uv = root / "uv.json"
    eft1 = root / "eft1.json"
    final = root / "final.json"
    target = root / "target.json"

    uv.write_text(
        json.dumps(
            {
                "status": "Success",
                "metadata": {},
                "report_betas": {},
            }
        ),
        encoding="utf-8",
    )
    eft1.write_text(
        json.dumps(
            {
                "status": "Success",
                "metadata": {},
                "report_betas": {},
            }
        ),
        encoding="utf-8",
    )
    final.write_text(
        json.dumps(
            {
                "status": "Success",
                "combined": {
                    "ready_for_physical_majorana_numerics": True
                },
            }
        ),
        encoding="utf-8",
    )
    target.write_text(
        json.dumps(
            {
                "source": "synthetic",
                "ordering": "NO",
                "central_values": {
                    "delta_m21_sq_ev2": 1.0,
                    "delta_m3l_sq_ev2": 2.0,
                    "sin2_theta12": 0.3,
                    "sin2_theta13": 0.02,
                    "sin2_theta23": 0.5,
                },
                "one_sigma_errors": {
                    "delta_m21_sq_ev2": 0.1,
                    "delta_m3l_sq_ev2": 0.1,
                    "sin2_theta12": 0.1,
                    "sin2_theta13": 0.01,
                    "sin2_theta23": 0.1,
                },
            }
        ),
        encoding="utf-8",
    )

    zero = [[0.0, 0.0, 0.0] for _ in range(3)]

    config = {
        "representation": {
            "d_s1": 1,
            "d_s2": 3,
            "d_f": 2,
            "alpha": 0,
            "shared_scalar": False,
        },
        "scales": {
            "mu_fermion_threshold_gev": 1.0e10,
            "mu_scalar_threshold_gev": 7.0e9,
            "mu_low_gev": 1.0e3,
        },
        "files": {
            "uv_rgbeta": "uv.json",
            "eft1_rgbeta": "eft1.json",
            "final_weinberg": "final.json",
            "oscillation_target": "target.json",
        },
        "output": {
            "result_json": "out/scan.json"
        },
        "base_state": {
            "mu_gev": 1.0e12,
            "sm": {
                "gY": 0.36,
                "g2": 0.65,
                "g3": 0.9,
                "lambdaH": 0.25,
                "yu": {
                    "real": [
                        [1e-5, 0.0, 0.0],
                        [0.0, 7e-3, 0.0],
                        [0.0, 0.0, 0.75],
                    ],
                    "imag": zero,
                },
                "yd": {
                    "real": [
                        [2e-5, 0.0, 0.0],
                        [0.0, 4e-4, 0.0],
                        [0.0, 0.0, 1.8e-2],
                    ],
                    "imag": zero,
                },
                "ye": {
                    "real": [
                        [3e-6, 0.0, 0.0],
                        [0.0, 6e-4, 0.0],
                        [0.0, 0.0, 1e-2],
                    ],
                    "imag": zero,
                },
            },
            "ordinary": {
                "y1": {
                    "real": [
                        [0.05, 0.0, 0.0],
                        [0.0, 0.07, 0.0],
                        [0.0, 0.0, 0.09],
                    ],
                    "imag": zero,
                },
                "y2": {
                    "real": [
                        [0.04, 0.0, 0.0],
                        [0.0, 0.06, 0.0],
                        [0.0, 0.0, 0.08],
                    ],
                    "imag": zero,
                },
                "MF": {
                    "real": [
                        [1e10, 0.0, 0.0],
                        [0.0, 1.1e10, 0.0],
                        [0.0, 0.0, 1.2e10],
                    ],
                    "imag": zero,
                },
                "mS1Sq": 6.4e19,
                "mS2Sq": 4.9e19,
                "lambdaS1": 0.1,
                "lambdaS2": 0.12,
                "lambdaH1": 0.02,
                "lambdaH2": 0.03,
                "lambda12": 0.01,
                "lambdaT3": {
                    "real": 0.005,
                    "imag": 0.002,
                },
                "lambdaH2Adj": 0.004,
                "lambdaS2Adj": 0.006,
            },
        },
        "scan": {
            "mode": "grid",
            "bindings": {
                "lambda_re": "ordinary.lambdaT3.real",
                "y11": "ordinary.y1.real.0.0",
            },
            "values": {
                "lambda_re": [0.001, 0.002],
                "y11": [0.03, 0.04],
            },
        },
    }

    config_path = root / "scan.json"
    config_path.write_text(
        json.dumps(config, indent=2),
        encoding="utf-8",
    )
    return config_path


class ScanCLITests(unittest.TestCase):
    def test_nested_binding_supports_matrix_entries(self) -> None:
        payload = {
            "ordinary": {
                "lambdaT3": {"real": 0.1, "imag": 0.0},
                "y1": {"real": [[1.0, 2.0], [3.0, 4.0]]},
            }
        }

        _set_bound_value(
            payload,
            "ordinary.lambdaT3.real",
            0.25,
        )
        _set_bound_value(
            payload,
            "ordinary.y1.real.1.0",
            9.0,
        )

        self.assertEqual(
            payload["ordinary"]["lambdaT3"]["real"],
            0.25,
        )
        self.assertEqual(
            payload["ordinary"]["y1"]["real"][1][0],
            9.0,
        )

    def test_relative_paths_resolve_from_config_directory(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            path = _write_fixture_tree(root)
            config = load_scan_config(path)

            self.assertEqual(
                config.uv_rgbeta_path,
                (root / "uv.json").resolve(),
            )
            self.assertEqual(
                config.output_path,
                (root / "out" / "scan.json").resolve(),
            )
            self.assertEqual(
                config.checkpoint_path,
                (root / "out" / "scan.checkpoint.json").resolve(),
            )

    def test_grid_points_follow_bindings(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            config = load_scan_config(
                _write_fixture_tree(Path(tmpdir))
            )
            points = build_scan_points(config)

        self.assertEqual(len(points), 4)
        self.assertEqual(
            points[0],
            {"lambda_re": 0.001, "y11": 0.03},
        )
        self.assertEqual(
            points[-1],
            {"lambda_re": 0.002, "y11": 0.04},
        )

    def test_bound_parameters_construct_valid_uv_state(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            config = load_scan_config(
                _write_fixture_tree(Path(tmpdir))
            )
            state = build_uv_state_from_config(
                config,
                {
                    "lambda_re": 0.009,
                    "y11": 0.11,
                },
            )

        self.assertAlmostEqual(
            state.lambdaT3.real,
            0.009,
        )
        self.assertAlmostEqual(
            state.lambdaT3.imag,
            0.002,
        )
        self.assertAlmostEqual(
            state.y1[0, 0].real,
            0.11,
        )
        self.assertEqual(
            state.representation.d_s2,
            3,
        )

    def test_shared_scalar_cli_config_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            path = _write_fixture_tree(root)
            payload = json.loads(
                path.read_text(encoding="utf-8")
            )
            payload["representation"]["shared_scalar"] = True
            path.write_text(
                json.dumps(payload),
                encoding="utf-8",
            )

            with self.assertRaisesRegex(
                ValueError,
                "ordinary split-scalar",
            ):
                load_scan_config(path)


    def test_uv_payload_representation_mismatch_is_rejected(self) -> None:
        representation = T3Representation(
            d_s1=1,
            d_s2=3,
            d_f=2,
            alpha=0,
            shared_scalar=False,
        ).validated()

        payload = {
            "status": "Success",
            "metadata": {
                "dS1": 2,
                "dS2": 2,
                "dF": 1,
                "alpha": -1,
                "ReportBetaConvention": (
                    "16*pi^2*dX/dln(mu); gauge BetaTerm divided by 2*g"
                ),
            },
            "report_betas": {},
        }

        with self.assertRaisesRegex(
            ValueError,
            "representation mismatch",
        ):
            _validate_rgbeta_payload_compatibility(
                payload=payload,
                representation=representation,
                payload_name="UV RGBeta payload",
                eft1=False,
            )

    def test_eft1_stage_metadata_is_checked(self) -> None:
        representation = T3Representation(
            d_s1=2,
            d_s2=2,
            d_f=1,
            alpha=-1,
            shared_scalar=False,
        ).validated()

        payload = {
            "status": "Success",
            "metadata": {
                "dS1": 2,
                "dS2": 2,
                "dF": 1,
                "alpha": -1,
                "IntegratedField": "S1",
                "ActiveBSMFields": ["F", "S2"],
                "ReportBetaConvention": (
                    "16*pi^2*dX/dln(mu); gauge BetaTerm divided by 2*g"
                ),
            },
            "report_betas": {},
        }

        with self.assertRaisesRegex(
            ValueError,
            "post-F EFT1 stage",
        ):
            _validate_rgbeta_payload_compatibility(
                payload=payload,
                representation=representation,
                payload_name="EFT1 RGBeta payload",
                eft1=True,
            )



if __name__ == "__main__":
    unittest.main()
