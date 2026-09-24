from __future__ import annotations

from pathlib import Path
import unittest

import numpy as np

from Numerical.State import (
    SMNumericalState,
    T3Representation,
    T3UVState,
)
from Numerical.RGBetaEvaluator import load_rgbeta_payload
from Numerical.UVRunner import run_uv_segment


FIXTURE = (
    Path(__file__).resolve().parent
    / "data"
    / "rgbeta_class_a_sample.json"
)


def _class_a_state() -> T3UVState:
    return T3UVState(
        mu_gev=1.0e12,
        representation=T3Representation(
            d_s1=1,
            d_s2=3,
            d_f=2,
            alpha=0,
        ),
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


class UVRunnerTests(unittest.TestCase):
    def test_short_downward_run_reaches_requested_endpoint(self) -> None:
        payload = load_rgbeta_payload(FIXTURE)
        state = _class_a_state()

        result = run_uv_segment(
            state,
            payload,
            9.0e11,
        )

        self.assertTrue(result.solver_success)
        self.assertGreaterEqual(result.n_points, 2)
        self.assertAlmostEqual(result.mu_gev[0], state.mu_gev)
        self.assertAlmostEqual(result.mu_gev[-1], 9.0e11)
        self.assertEqual(result.y.shape[0], result.layout.size)

    def test_requested_save_scales_are_retained(self) -> None:
        payload = load_rgbeta_payload(FIXTURE)
        state = _class_a_state()

        requested = [1.0e12, 9.8e11, 9.5e11, 9.0e11]

        result = run_uv_segment(
            state,
            payload,
            9.0e11,
            save_scales_gev=requested,
        )

        np.testing.assert_allclose(
            result.mu_gev,
            np.asarray(requested),
            rtol=1.0e-13,
            atol=0.0,
        )

    def test_final_state_is_reconstructible_and_finite(self) -> None:
        payload = load_rgbeta_payload(FIXTURE)
        state = _class_a_state()

        result = run_uv_segment(
            state,
            payload,
            9.0e11,
        )

        final_state = result.final_state

        self.assertAlmostEqual(final_state.mu_gev, 9.0e11)
        self.assertTrue(np.isfinite(final_state.sm.gY))
        self.assertTrue(np.all(np.isfinite(final_state.y1.real)))
        self.assertTrue(np.all(np.isfinite(final_state.y1.imag)))
        self.assertEqual(final_state.y1.shape, (3, 3))

    def test_running_changes_a_coupling(self) -> None:
        payload = load_rgbeta_payload(FIXTURE)
        state = _class_a_state()

        result = run_uv_segment(
            state,
            payload,
            9.0e11,
        )

        final_state = result.final_state

        self.assertNotEqual(final_state.sm.gY, state.sm.gY)
        self.assertNotEqual(final_state.sm.g2, state.sm.g2)

    def test_upward_running_is_supported(self) -> None:
        payload = load_rgbeta_payload(FIXTURE)
        state = _class_a_state()

        result = run_uv_segment(
            state,
            payload,
            1.1e12,
            save_scales_gev=[1.0e12, 1.05e12, 1.1e12],
        )

        self.assertAlmostEqual(result.mu_gev[-1], 1.1e12)
        self.assertGreater(result.mu_gev[-1], result.mu_gev[0])

    def test_equal_scales_are_rejected(self) -> None:
        payload = load_rgbeta_payload(FIXTURE)
        state = _class_a_state()

        with self.assertRaisesRegex(
            ValueError,
            "must be different",
        ):
            run_uv_segment(
                state,
                payload,
                state.mu_gev,
            )


if __name__ == "__main__":
    unittest.main()
