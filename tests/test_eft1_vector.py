from __future__ import annotations

import unittest

import numpy as np

from Numerical.EFT1RGBetaEvaluator import (
    eft1_derivative_from_payload,
    evaluate_eft1_rgbeta_payload,
)
from Numerical.EFT1State import T3EFT1State
from Numerical.EFT1StateVector import (
    pack_eft1_state,
    unpack_eft1_state,
)
from Numerical.State import (
    SMNumericalState,
    T3Representation,
)


def _state() -> T3EFT1State:
    return T3EFT1State(
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


def _synthetic_payload(state: T3EFT1State) -> dict:
    # This payload is structural only: it checks the real parser/evaluator
    # against the exact coupling set exported by current T3RGBetaEFT1OneLoopBetas.
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


class EFT1VectorTests(unittest.TestCase):
    def test_round_trip(self) -> None:
        state = _state()
        vector, layout = pack_eft1_state(state)
        rebuilt = unpack_eft1_state(vector, layout)

        self.assertEqual(rebuilt.representation, state.representation)
        self.assertEqual(rebuilt.lambdaT3, state.lambdaT3)
        self.assertEqual(rebuilt.lambdaH2Adj, state.lambdaH2Adj)
        np.testing.assert_allclose(rebuilt.sm.yu, state.sm.yu)

        self.assertNotIn("MF", layout.names)
        self.assertNotIn("y1", layout.names)
        self.assertNotIn("y2", layout.names)

    def test_layout_matches_current_eft1_export_keys(self) -> None:
        state = _state()
        _, layout = pack_eft1_state(state)
        payload = _synthetic_payload(state)

        self.assertEqual(
            set(layout.names),
            set(payload["report_betas"]),
        )

    def test_evaluator_maps_all_eft1_betas(self) -> None:
        state = _state()
        payload = _synthetic_payload(state)

        evaluated = evaluate_eft1_rgbeta_payload(payload, state)
        _, layout = pack_eft1_state(state)

        self.assertEqual(set(evaluated), set(layout.names))
        np.testing.assert_allclose(evaluated["yu"], state.sm.yu)
        self.assertEqual(evaluated["lambdaT3"], state.lambdaT3 * (
            state.sm.lambdaH + state.lambdaH2Adj
        ))

    def test_derivative_has_state_vector_shape(self) -> None:
        state = _state()
        payload = _synthetic_payload(state)

        derivative, layout = eft1_derivative_from_payload(
            payload,
            state,
        )
        vector, expected_layout = pack_eft1_state(state)

        self.assertEqual(layout, expected_layout)
        self.assertEqual(derivative.shape, vector.shape)
        self.assertTrue(np.all(np.isfinite(derivative)))

    def test_wrong_stage_metadata_is_rejected(self) -> None:
        state = _state()
        payload = _synthetic_payload(state)
        payload["metadata"]["IntegratedField"] = "S1"

        with self.assertRaisesRegex(
            ValueError,
            "does not match",
        ):
            evaluate_eft1_rgbeta_payload(payload, state)


if __name__ == "__main__":
    unittest.main()
