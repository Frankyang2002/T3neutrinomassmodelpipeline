from __future__ import annotations

from pathlib import Path
import unittest

import numpy as np

from Numerical.BetaVector import LOOP_FACTOR
from Numerical.RGBetaEvaluator import (
    derivative_from_rgbeta_payload,
    evaluate_inputform_expression,
    evaluate_rgbeta_payload,
    load_rgbeta_payload,
)
from Numerical.State import (
    SMNumericalState,
    T3Representation,
    T3UVState,
)
from Numerical.StateVector import pack_uv_state


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


class RGBetaEvaluatorTests(unittest.TestCase):
    def test_matrix_bar_trans_trace_syntax(self) -> None:
        y = np.array(
            [
                [1.0 + 2.0j, 2.0 - 1.0j],
                [0.5 + 0.25j, -3.0 + 0.0j],
            ],
            dtype=complex,
        )

        env = {"y": y}

        trace_value = evaluate_inputform_expression(
            "Tr[y . Trans[Bar[y]]]",
            env,
        )
        expected_trace = np.trace(y @ y.conj().T)

        self.assertAlmostEqual(trace_value, expected_trace)

        matrix_value = evaluate_inputform_expression(
            "Matrix[y, Trans[Bar[y]], y][gen[$i], gen[$j]]",
            env,
        )
        expected_matrix = y @ y.conj().T @ y

        np.testing.assert_allclose(
            matrix_value,
            expected_matrix,
            rtol=0.0,
            atol=0.0,
        )

    def test_real_class_a_export_evaluates_every_beta(self) -> None:
        payload = load_rgbeta_payload(FIXTURE)
        state = _class_a_state()

        evaluated = evaluate_rgbeta_payload(payload, state)
        _, layout = pack_uv_state(state)

        self.assertEqual(set(evaluated), set(layout.names))

        for name in ("yu", "yd", "ye", "y1", "y2", "MF"):
            self.assertEqual(evaluated[name].shape, (3, 3))

        for name in (
            "gY",
            "g2",
            "g3",
            "mS1Sq",
            "mS2Sq",
            "lambdaH",
            "lambdaS1",
            "lambdaS2",
            "lambdaH1",
            "lambdaH2",
            "lambda12",
            "lambdaT3",
            "lambdaH2Adj",
            "lambdaS2Adj",
        ):
            self.assertEqual(np.ndim(evaluated[name]), 0)

    def test_report_gauge_convention_is_used(self) -> None:
        payload = load_rgbeta_payload(FIXTURE)
        state = _class_a_state()
        evaluated = evaluate_rgbeta_payload(payload, state)

        self.assertAlmostEqual(
            evaluated["gY"],
            (59.0 / 6.0) * state.sm.gY**3,
        )
        self.assertAlmostEqual(
            evaluated["g2"],
            -0.5 * state.sm.g2**3,
        )
        self.assertAlmostEqual(
            evaluated["g3"],
            -7.0 * state.sm.g3**3,
        )

    def test_derivative_matches_state_vector_length_and_loop_factor(self) -> None:
        payload = load_rgbeta_payload(FIXTURE)
        state = _class_a_state()

        derivative, layout = derivative_from_rgbeta_payload(
            payload,
            state,
        )

        vector, state_layout = pack_uv_state(state)

        self.assertEqual(layout, state_layout)
        self.assertEqual(derivative.shape, vector.shape)

        gY_block = layout.block("gY")
        expected = (
            (59.0 / 6.0)
            * state.sm.gY**3
            / LOOP_FACTOR
        )
        self.assertAlmostEqual(
            derivative[gY_block.start],
            expected,
        )

    def test_payload_representation_mismatch_is_rejected(self) -> None:
        payload = load_rgbeta_payload(FIXTURE)
        state = _class_a_state()

        broken = dict(payload)
        broken["metadata"] = dict(payload["metadata"])
        broken["metadata"]["alpha"] = 1

        with self.assertRaisesRegex(
            ValueError,
            "does not match the numerical state",
        ):
            evaluate_rgbeta_payload(broken, state)


if __name__ == "__main__":
    unittest.main()
