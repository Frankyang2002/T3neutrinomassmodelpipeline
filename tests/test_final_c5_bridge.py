from __future__ import annotations

from pathlib import Path
import unittest
from unittest.mock import patch

import numpy as np

from Numerical.FinalC5Bridge import (
    FinalC5TrajectoryBuilder,
    _heavy_mass_basis_data,
    evaluate_final_c5_from_trajectory,
    final_c5_inputs_from_trajectory,
)
from Numerical.State import (
    SMNumericalState,
    T3Representation,
    T3UVState,
)
from Numerical.T3Trajectory import (
    run_renormalisable_t3_trajectory,
)


def _uv_state() -> T3UVState:
    return T3UVState(
        mu_gev=1.0e12,
        representation=T3Representation(1, 3, 2, 0),
        sm=SMNumericalState(
            gY=0.36,
            g2=0.65,
            g3=0.90,
            lambdaH=0.25,
            yu=np.diag([1.0e-5, 7.0e-3, 0.75]).astype(complex),
            yd=np.diag([2.0e-5, 4.0e-4, 1.8e-2]).astype(complex),
            ye=np.diag([3.0e-6, 6.0e-4, 1.0e-2]).astype(complex),
        ),
        y1=np.diag([0.05, 0.07, 0.09]).astype(complex),
        y2=np.diag([0.04, 0.06, 0.08]).astype(complex),
        MF=np.diag([1.0e10, 1.1e10, 1.2e10]).astype(complex),
        mS1Sq=(8.0e9) ** 2,
        mS2Sq=(7.0e9) ** 2,
        lambdaS1=0.10,
        lambdaS2=0.12,
        lambdaH1=0.02,
        lambdaH2=0.03,
        lambda12=0.01,
        lambdaT3=0.005 + 0.002j,
        lambdaH2Adj=0.004,
        lambdaS2Adj=0.006,
    ).validated()


def _uv_payload(state: T3UVState) -> dict:
    # Keep y1, y2 and MF exactly constant in this bridge regression so their
    # scale assignment can be checked independently of beta-function details.
    betas = {
        "gY": "0*gY",
        "g2": "0*g2",
        "g3": "0*g3",
        "yu": "0*Matrix[yu][gen[$i], gen[$j]]",
        "yd": "0*Matrix[yd][gen[$i], gen[$j]]",
        "ye": "0*Matrix[ye][gen[$i], gen[$j]]",
        "y1": "0*Matrix[y1][gen[$i], heavy[$j]]",
        "y2": "0*Matrix[y2][gen[$i], heavy[$j]]",
        "MF": "0*Matrix[MF][heavy[$i], heavy[$j]]",
        "mS1Sq": "0*mS1Sq",
        "mS2Sq": "0*mS2Sq",
        "lambdaH": "0*lambdaH",
        "lambdaS1": "0*lambdaS1",
        "lambdaS2": "0*lambdaS2",
        "lambdaH1": "0*lambdaH1",
        "lambdaH2": "0*lambdaH2",
        "lambda12": "0*lambda12",
        "lambdaT3": "0*lambdaT3",
        "lambdaH2Adj": "0*lambdaH2Adj",
        "lambdaS2Adj": "0*lambdaS2Adj",
    }

    return {
        "status": "Success",
        "metadata": {
            "SharedScalar": False,
            "dS1": 1,
            "dS2": 3,
            "dF": 2,
            "alpha": 0,
            "ReportBetaConvention": (
                "16*pi^2*dX/dln(mu); gauge BetaTerm divided by 2*g"
            ),
        },
        "report_betas": betas,
    }


def _eft1_payload() -> dict:
    betas = {
        "gY": "0*gY",
        "g2": "0*g2",
        "g3": "0*g3",
        "yu": "0*Matrix[yu][gen[$i], gen[$j]]",
        "yd": "0*Matrix[yd][gen[$i], gen[$j]]",
        "ye": "0*Matrix[ye][gen[$i], gen[$j]]",
        "mS1Sq": "0*mS1Sq",
        "mS2Sq": "0*mS2Sq",
        "lambdaH": "0*lambdaH",
        "lambdaS1": "0*lambdaS1",
        "lambdaS2": "0*lambdaS2",
        "lambdaH1": "0*lambdaH1",
        "lambdaH2": "0*lambdaH2",
        "lambda12": "0*lambda12",
        "lambdaT3": "0*lambdaT3",
        "lambdaH2Adj": "0*lambdaH2Adj",
        "lambdaS2Adj": "0*lambdaS2Adj",
    }

    return {
        "status": "Success",
        "metadata": {
            "SharedScalar": False,
            "dS": None,
            "dS1": 1,
            "dS2": 3,
            "dF": 2,
            "alpha": 0,
            "IntegratedField": "F",
            "ActiveBSMFields": ["S1", "S2"],
            "ReportBetaConvention": (
                "16*pi^2*dX/dln(mu); gauge BetaTerm divided by 2*g"
            ),
        },
        "report_betas": betas,
    }


def _trajectory():
    state = _uv_state()
    return run_renormalisable_t3_trajectory(
        state,
        _uv_payload(state),
        _eft1_payload(),
        mu_fermion_threshold_gev=1.0e10,
        mu_scalar_threshold_gev=7.0e9,
    )


class FinalC5BridgeTests(unittest.TestCase):
    def test_extracts_values_at_correct_thresholds(self) -> None:
        trajectory = _trajectory()
        inputs = final_c5_inputs_from_trajectory(
            trajectory
        )

        np.testing.assert_allclose(
            inputs.heavy_masses_gev,
            np.array([1.0e10, 1.1e10, 1.2e10]),
        )
        np.testing.assert_allclose(
            inputs.y1,
            trajectory.uv_threshold_state.y1,
        )
        np.testing.assert_allclose(
            inputs.y2,
            trajectory.uv_threshold_state.y2,
        )

        self.assertEqual(
            inputs.scalar_mass1_gev,
            8.0e9,
        )
        self.assertEqual(
            inputs.scalar_mass2_gev,
            7.0e9,
        )
        self.assertEqual(
            inputs.grouped_scalar_matching_scale_gev,
            7.0e9,
        )
        self.assertEqual(
            inputs.lambda_t3,
            trajectory.eft1_threshold_state.lambdaT3,
        )

    def test_config_matches_existing_evaluator_schema(self) -> None:
        config = final_c5_inputs_from_trajectory(
            _trajectory()
        ).as_config()

        model = config["t3"]

        self.assertEqual(len(model["MF"]), 3)
        self.assertEqual(
            np.asarray(model["y1_real"]).shape,
            (3, 3),
        )
        self.assertEqual(
            np.asarray(model["y2_imag"]).shape,
            (3, 3),
        )
        self.assertEqual(model["MS"], 7.0e9)
        self.assertIn("lambdaT3", model)
        self.assertIn("hbar", model)

    def test_offdiagonal_nonmajorana_MF_is_rejected_without_silent_rotation(self) -> None:
        trajectory = _trajectory()
        state = trajectory.uv_threshold_state

        # Mutate only through construction of a replacement validated state.
        mf = state.MF.copy()
        mf[0, 1] = 1.0e6

        broken_state = T3UVState(
            mu_gev=state.mu_gev,
            representation=state.representation,
            sm=state.sm,
            y1=state.y1,
            y2=state.y2,
            MF=mf,
            mS1Sq=state.mS1Sq,
            mS2Sq=state.mS2Sq,
            lambdaS1=state.lambdaS1,
            lambdaS2=state.lambdaS2,
            lambdaH1=state.lambdaH1,
            lambdaH2=state.lambdaH2,
            lambda12=state.lambda12,
            lambdaT3=state.lambdaT3,
            lambdaH2Adj=state.lambdaH2Adj,
            lambdaS2Adj=state.lambdaS2Adj,
        ).validated()

        # Reuse the immutable trajectory container with a minimal UV-result
        # proxy whose final_state property is the broken state.
        class UVProxy:
            @property
            def final_state(self):
                return broken_state

        from dataclasses import replace
        broken_trajectory = replace(
            trajectory,
            uv=UVProxy(),
        )

        with self.assertRaisesRegex(
            ValueError,
            "diagonal positive heavy-fermion mass basis",
        ):
            final_c5_inputs_from_trajectory(
                broken_trajectory
            )


    def test_majorana_takagi_rotation_diagonalizes_MF_and_rotates_yukawas(self) -> None:
        rng = np.random.default_rng(12345)

        z = (
            rng.normal(size=(3, 3))
            + 1j * rng.normal(size=(3, 3))
        )
        q, r = np.linalg.qr(z)
        q = q @ np.diag(
            np.exp(-1j * np.angle(np.diag(r)))
        )

        masses = np.array(
            [1.0e10, 1.1e10, 1.2e10],
            dtype=float,
        )

        # Construct MF so q^T MF q = diag(masses).
        mf = q.conj() @ np.diag(masses) @ q.conj().T

        y1 = np.array(
            [
                [0.10 + 0.02j, -0.03 + 0.04j, 0.08 - 0.01j],
                [0.02 - 0.01j, 0.07 + 0.03j, -0.05 + 0.02j],
                [-0.04 + 0.06j, 0.01 - 0.02j, 0.09 + 0.05j],
            ],
            dtype=complex,
        )
        y2 = np.array(
            [
                [0.04 - 0.01j, 0.05 + 0.02j, -0.03 + 0.01j],
                [-0.06 + 0.03j, 0.08 - 0.04j, 0.02 + 0.01j],
                [0.03 + 0.02j, -0.02 + 0.05j, 0.07 - 0.03j],
            ],
            dtype=complex,
        )

        zero = 0.0 + 0.0j
        state = T3UVState(
            mu_gev=1.0e10,
            representation=T3Representation(
                2,
                2,
                1,
                -1,
            ),
            sm=_uv_state().sm,
            y1=y1,
            y2=y2,
            MF=mf,
            mS1Sq=(8.0e9) ** 2,
            mS2Sq=(7.0e9) ** 2,
            lambdaS1=0.10,
            lambdaS2=0.12,
            lambdaH1=0.02,
            lambdaH2=0.03,
            lambda12=0.01,
            lambdaT3=0.005 + 0.002j,
            lambdaH1Adj=0.0,
            lambdaH2Adj=0.0,
            lambda12Adj=0.0,
            lambdaHHdagS2S2=zero,
            lambdaHHdagS1barS1bar=zero,
            lambdaS1bar2S2bar2=zero,
            lambdaS1barS2S2bar2=zero,
            lambdaS1S1bar2S2bar=zero,
            lambdaHHdagS1barS2barCross=zero,
        ).validated()

        basis = _heavy_mass_basis_data(
            state,
            offdiagonal_atol=1.0e-10,
            imaginary_atol=1.0e-10,
            takagi_residual_rtol=1.0e-10,
            takagi_residual_atol=1.0e-10,
        )

        diagonalized = (
            basis.rotation.T
            @ state.MF
            @ basis.rotation
        )

        np.testing.assert_allclose(
            diagonalized,
            np.diag(basis.masses_gev),
            rtol=1.0e-10,
            atol=10.0,
        )
        np.testing.assert_allclose(
            basis.masses_gev,
            masses,
            rtol=1.0e-10,
            atol=10.0,
        )
        np.testing.assert_allclose(
            basis.y1,
            state.y1 @ basis.rotation,
            rtol=1.0e-12,
            atol=1.0e-12,
        )
        np.testing.assert_allclose(
            basis.y2,
            state.y2 @ basis.rotation,
            rtol=1.0e-12,
            atol=1.0e-12,
        )

    def test_diagonal_positive_majorana_basis_is_unchanged(self) -> None:
        zero = 0.0 + 0.0j
        y1 = np.array(
            [
                [0.10 + 0.02j, 0.03, 0.01],
                [0.02, 0.07 + 0.03j, 0.04],
                [0.05, 0.01, 0.09 + 0.05j],
            ],
            dtype=complex,
        )
        y2 = np.array(
            [
                [0.04 - 0.01j, 0.02, 0.03],
                [0.01, 0.08 - 0.04j, 0.02],
                [0.03, 0.02, 0.07 - 0.03j],
            ],
            dtype=complex,
        )
        state = T3UVState(
            mu_gev=1.0e10,
            representation=T3Representation(
                2,
                2,
                1,
                -1,
            ),
            sm=_uv_state().sm,
            y1=y1,
            y2=y2,
            MF=np.diag(
                [1.0e10, 1.1e10, 1.2e10]
            ).astype(complex),
            mS1Sq=(8.0e9) ** 2,
            mS2Sq=(7.0e9) ** 2,
            lambdaS1=0.10,
            lambdaS2=0.12,
            lambdaH1=0.02,
            lambdaH2=0.03,
            lambda12=0.01,
            lambdaT3=0.005 + 0.002j,
            lambdaH1Adj=0.0,
            lambdaH2Adj=0.0,
            lambda12Adj=0.0,
            lambdaHHdagS2S2=zero,
            lambdaHHdagS1barS1bar=zero,
            lambdaS1bar2S2bar2=zero,
            lambdaS1barS2S2bar2=zero,
            lambdaS1S1bar2S2bar=zero,
            lambdaHHdagS1barS2barCross=zero,
        ).validated()

        basis = _heavy_mass_basis_data(
            state,
            offdiagonal_atol=1.0e-10,
            imaginary_atol=1.0e-10,
            takagi_residual_rtol=1.0e-10,
            takagi_residual_atol=1.0e-10,
        )

        np.testing.assert_allclose(
            basis.rotation,
            np.eye(3),
        )
        np.testing.assert_allclose(
            basis.y1,
            y1,
        )
        np.testing.assert_allclose(
            basis.y2,
            y2,
        )
        np.testing.assert_allclose(
            basis.masses_gev,
            [1.0e10, 1.1e10, 1.2e10],
        )

    def test_existing_final_c5_evaluator_is_called(self) -> None:
        trajectory = _trajectory()
        expected = np.eye(3, dtype=complex) * 1.0e-14

        with patch(
            "Numerical.FinalC5Bridge.evaluate_final_weinberg_json",
            return_value=expected,
        ) as mocked:
            result = evaluate_final_c5_from_trajectory(
                Path("final_weinberg_coefficient.json"),
                trajectory,
            )

        np.testing.assert_allclose(result, expected)
        mocked.assert_called_once()

        _, config = mocked.call_args.args
        self.assertEqual(
            config["t3"]["MS"],
            7.0e9,
        )

    def test_callable_builder_matches_parameter_scan_interface(self) -> None:
        trajectory = _trajectory()
        expected = np.eye(3, dtype=complex) * 2.0e-14

        builder = FinalC5TrajectoryBuilder(
            Path("final_weinberg_coefficient.json")
        )

        with patch(
            "Numerical.FinalC5Bridge.evaluate_final_c5_from_trajectory",
            return_value=expected,
        ):
            result = builder(
                {"lambdaT3": 0.01},
                trajectory,
            )

        np.testing.assert_allclose(result, expected)


if __name__ == "__main__":
    unittest.main()
