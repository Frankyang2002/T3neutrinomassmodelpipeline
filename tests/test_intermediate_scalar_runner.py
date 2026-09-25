from __future__ import annotations

import unittest

import numpy as np

from Numerical.IntermediateScalarRunner import run_intermediate_scalar_segment
from Numerical.IntermediateScalarState import T3IntermediateScalarState
from Numerical.State import (
    SMNumericalState,
    T3Representation,
)


def _state() -> T3IntermediateScalarState:
    return T3IntermediateScalarState(
        mu_gev=1.0e10,
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


def _payload(state: T3IntermediateScalarState) -> dict:
    # Structural payload matching the exact current scalar-only RGBeta keys.
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
            "YS1": "0",
            "YS2": "1",
            "YS": "Null",
            "IntegratedField": "F",
            "ActiveBSMFields": ["S1", "S2"],
            "ReportBetaConvention": (
                "16*pi^2*dX/dln(mu); gauge BetaTerm divided by 2*g"
            ),
        },
        "report_betas": betas,
    }


class IntermediateScalarRunnerTests(unittest.TestCase):
    def test_short_downward_run_reaches_endpoint(self) -> None:
        state = _state()
        result = run_intermediate_scalar_segment(state, _payload(state), 9.0e9)

        self.assertTrue(result.solver_success)
        self.assertGreaterEqual(result.n_points, 2)
        self.assertEqual(result.mu_gev[0], 1.0e10)
        self.assertEqual(result.mu_gev[-1], 9.0e9)

    def test_requested_save_scales_are_retained(self) -> None:
        state = _state()
        scales = [1.0e10, 9.8e9, 9.5e9, 9.0e9]

        result = run_intermediate_scalar_segment(
            state,
            _payload(state),
            9.0e9,
            save_scales_gev=scales,
        )

        np.testing.assert_allclose(
            result.mu_gev,
            np.asarray(scales),
            rtol=0.0,
            atol=0.0,
        )

    def test_final_state_reconstructs(self) -> None:
        state = _state()
        result = run_intermediate_scalar_segment(state, _payload(state), 9.0e9)

        final_state = result.final_state

        self.assertIsInstance(final_state, T3IntermediateScalarState)
        self.assertEqual(final_state.mu_gev, 9.0e9)
        self.assertTrue(np.isfinite(final_state.sm.gY))
        self.assertTrue(np.isfinite(final_state.lambdaH2Adj))
        self.assertFalse(hasattr(final_state, "MF"))

    def test_running_changes_couplings(self) -> None:
        state = _state()
        result = run_intermediate_scalar_segment(state, _payload(state), 9.0e9)

        final_state = result.final_state

        self.assertNotEqual(final_state.sm.gY, state.sm.gY)
        self.assertNotEqual(final_state.lambdaS1, state.lambdaS1)

    def test_upward_running_is_supported(self) -> None:
        state = _state()

        result = run_intermediate_scalar_segment(
            state,
            _payload(state),
            1.1e10,
            save_scales_gev=[1.0e10, 1.05e10, 1.1e10],
        )

        self.assertEqual(result.mu_gev[-1], 1.1e10)
        self.assertGreater(result.mu_gev[-1], result.mu_gev[0])

    def test_equal_scales_are_rejected(self) -> None:
        state = _state()

        with self.assertRaisesRegex(ValueError, "must be different"):
            run_intermediate_scalar_segment(state, _payload(state), state.mu_gev)


if __name__ == "__main__":
    unittest.main()
