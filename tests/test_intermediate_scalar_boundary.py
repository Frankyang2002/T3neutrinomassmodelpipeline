from __future__ import annotations

import unittest

import numpy as np

from Numerical.State import (
    SMNumericalState,
    T3Representation,
    T3UVState,
)
from Numerical.T3Trajectory import (
    continue_with_weinberg_running,
    run_renormalisable_t3_trajectory,
)


def _uv_state() -> T3UVState:
    return T3UVState(
        mu_gev=1.0e12,
        representation=T3Representation(1, 3, 2, 0),
        sm=SMNumericalState(
            gY=0.36,
            g2=0.65,
            g3=1.05,
            lambdaH=0.25,
            yu=np.diag([1.0e-5, 7.0e-3, 0.85]).astype(complex),
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
    betas = {
        "gY": "(59*gY^3)/6",
        "g2": "-1/2*g2^3",
        "g3": "-7*g3^3",
        "yu": "Matrix[yu][gen[$i], gen[$j]]",
        "yd": "Matrix[yd][gen[$i], gen[$j]]",
        "ye": "Matrix[ye][gen[$i], gen[$j]]",
        "y1": "Matrix[y1][gen[$i], heavy[$j]]",
        "y2": "Matrix[y2][gen[$i], heavy[$j]]",
        "MF": "Matrix[MF][heavy[$i], heavy[$j]]",
        "mS1Sq": "2*mS1Sq",
        "mS2Sq": "2*mS2Sq",
        "lambdaH": "lambdaH^2",
        "lambdaS1": "lambdaS1^2",
        "lambdaS2": "lambdaS2^2 + lambdaS2Adj^2",
        "lambdaH1": "lambdaH1*lambdaH",
        "lambdaH2": "lambdaH2*lambdaH",
        "lambda12": "lambda12*(lambdaS1 + lambdaS2)",
        "lambdaT3": "lambdaT3*(lambdaH + lambdaH2Adj)",
        "lambdaH2Adj": "lambdaH2Adj*lambdaH",
        "lambdaS2Adj": "lambdaS2Adj*lambdaS2",
    }
    return {
        "status": "Success",
        "metadata": {
            "SharedScalar": False,
            "dS1": state.representation.d_s1,
            "dS2": state.representation.d_s2,
            "dF": state.representation.d_f,
            "alpha": state.representation.alpha,
            "ReportBetaConvention": (
                "16*pi^2*dX/dln(mu); gauge BetaTerm divided by 2*g"
            ),
        },
        "report_betas": betas,
    }


def _intermediate_payload(state: T3UVState) -> dict:
    betas = {
        "gY": "(43*gY^3)/6",
        "g2": "-2*g2^3",
        "g3": "-7*g3^3",
        "yu": "Matrix[yu][gen[$i], gen[$j]]",
        "yd": "Matrix[yd][gen[$i], gen[$j]]",
        "ye": "Matrix[ye][gen[$i], gen[$j]]",
        "mS1Sq": "2*mS1Sq + lambda12*mS2Sq",
        "mS2Sq": "2*mS2Sq + lambda12*mS1Sq",
        "lambdaH": "lambdaH^2 + lambdaH1^2 + lambdaH2^2",
        "lambdaS1": "lambdaS1^2",
        "lambdaS2": "lambdaS2^2 + lambdaS2Adj^2",
        "lambdaH1": "lambdaH1*lambdaH",
        "lambdaH2": "lambdaH2*lambdaH",
        "lambda12": "lambda12*(lambdaS1 + lambdaS2)",
        "lambdaT3": "lambdaT3*(lambdaH + lambdaH2Adj)",
        "lambdaH2Adj": "lambdaH2Adj*lambdaH",
        "lambdaS2Adj": "lambdaS2Adj*lambdaS2",
    }
    return {
        "status": "Success",
        "metadata": {
            "SharedScalar": False,
            "dS": None,
            "dS1": state.representation.d_s1,
            "dS2": state.representation.d_s2,
            "dF": state.representation.d_f,
            "alpha": state.representation.alpha,
            "IntegratedField": "F",
            "ActiveBSMFields": ["S1", "S2"],
            "ReportBetaConvention": (
                "16*pi^2*dX/dln(mu); gauge BetaTerm divided by 2*g"
            ),
        },
        "report_betas": betas,
    }


class T3TrajectoryTests(unittest.TestCase):
    def test_complete_renormalisable_path(self) -> None:
        state = _uv_state()
        trajectory = run_renormalisable_t3_trajectory(
            state,
            _uv_payload(state),
            _intermediate_payload(state),
            mu_fermion_threshold_gev=1.0e10,
            mu_scalar_threshold_gev=7.0e9,
        )
        self.assertEqual(trajectory.uv.mu_initial_gev, 1.0e12)
        self.assertEqual(trajectory.uv.mu_final_gev, 1.0e10)
        self.assertEqual(trajectory.eft1.mu_initial_gev, 1.0e10)
        self.assertEqual(trajectory.eft1.mu_final_gev, 7.0e9)
        self.assertEqual(trajectory.final_sm_boundary.mu_gev, 7.0e9)
        self.assertFalse(hasattr(trajectory.eft1_initial_state, "MF"))

    def test_threshold_diagnostics_are_retained(self) -> None:
        state = _uv_state()
        trajectory = run_renormalisable_t3_trajectory(
            state,
            _uv_payload(state),
            _intermediate_payload(state),
            mu_fermion_threshold_gev=1.0e10,
            mu_scalar_threshold_gev=7.0e9,
        )
        self.assertEqual(trajectory.fermion_threshold.masses_gev.shape, (3,))
        self.assertEqual(len(trajectory.scalar_threshold.masses_gev), 2)

    def test_explicit_save_scales_survive_both_segments(self) -> None:
        state = _uv_state()
        trajectory = run_renormalisable_t3_trajectory(
            state,
            _uv_payload(state),
            _intermediate_payload(state),
            mu_fermion_threshold_gev=1.0e10,
            mu_scalar_threshold_gev=7.0e9,
            uv_save_scales_gev=[1.0e12, 1.0e11, 1.0e10],
            eft1_save_scales_gev=[1.0e10, 8.0e9, 7.0e9],
        )
        np.testing.assert_allclose(
            trajectory.uv.mu_gev,
            np.array([1.0e12, 1.0e11, 1.0e10]),
        )
        np.testing.assert_allclose(
            trajectory.eft1.mu_gev,
            np.array([1.0e10, 8.0e9, 7.0e9]),
        )

    def test_non_descending_thresholds_are_rejected(self) -> None:
        state = _uv_state()
        with self.assertRaisesRegex(ValueError, "mu_UV > mu_F > mu_S"):
            run_renormalisable_t3_trajectory(
                state,
                _uv_payload(state),
                _intermediate_payload(state),
                mu_fermion_threshold_gev=5.0e9,
                mu_scalar_threshold_gev=7.0e9,
            )

    def test_weinberg_continuation_uses_scalar_boundary(self) -> None:
        state = _uv_state()
        trajectory = run_renormalisable_t3_trajectory(
            state,
            _uv_payload(state),
            _intermediate_payload(state),
            mu_fermion_threshold_gev=1.0e10,
            mu_scalar_threshold_gev=7.0e9,
        )
        c5 = np.array(
            [
                [1.0e-14, 2.0e-15, 0.0],
                [2.0e-15, 3.0e-14, 1.0e-15],
                [0.0, 1.0e-15, 5.0e-14],
            ],
            dtype=complex,
        )
        full = continue_with_weinberg_running(trajectory, c5, 1.0e9)
        self.assertEqual(full.weinberg.mu_initial, 7.0e9)
        self.assertEqual(full.weinberg.mu_final, 1.0e9)
        np.testing.assert_allclose(full.c5_at_scalar_threshold, c5)

    def test_low_scale_must_be_below_scalar_threshold(self) -> None:
        state = _uv_state()
        trajectory = run_renormalisable_t3_trajectory(
            state,
            _uv_payload(state),
            _intermediate_payload(state),
            mu_fermion_threshold_gev=1.0e10,
            mu_scalar_threshold_gev=7.0e9,
        )
        with self.assertRaisesRegex(ValueError, "mu_low < mu_S"):
            continue_with_weinberg_running(
                trajectory,
                np.eye(3, dtype=complex) * 1.0e-14,
                8.0e9,
            )


if __name__ == "__main__":
    unittest.main()
