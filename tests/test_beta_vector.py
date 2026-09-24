from __future__ import annotations

import unittest

import numpy as np

from Numerical.BetaVector import (
    LOOP_FACTOR,
    beta_values_to_derivative,
    validate_beta_keys,
)
from Numerical.State import (
    SMNumericalState,
    SharedT3UVState,
    T3Representation,
)
from Numerical.StateVector import (
    pack_uv_state,
)


def _sm_state() -> SMNumericalState:
    return SMNumericalState(
        gY=0.36,
        g2=0.65,
        g3=1.1,
        lambdaH=0.26,
        yu=np.diag([1.0e-5, 7.0e-3, 0.9]).astype(complex),
        yd=np.diag([2.0e-5, 4.0e-4, 2.0e-2]).astype(complex),
        ye=np.diag([3.0e-6, 6.0e-4, 1.0e-2]).astype(complex),
    )


def _shared_state() -> SharedT3UVState:
    return SharedT3UVState(
        mu_gev=1.0e12,
        representation=T3Representation(
            2,
            2,
            1,
            -1,
            shared_scalar=True,
        ),
        sm=_sm_state(),
        h=np.eye(3, dtype=complex) * (0.1 + 0.02j),
        MF=np.diag([1.0e10, 1.1e10, 1.2e10]).astype(complex),
        mSSq=(8.0e9) ** 2,
        lambdaS=0.1,
        lambda3=0.01,
        lambda4=0.02,
        lambda5=0.03,
    ).validated()


def _mock_betas(layout) -> dict[str, object]:
    """Create shape-correct beta values with easy-to-recognise entries."""
    result: dict[str, object] = {}

    for index, block in enumerate(layout.blocks, start=1):
        scale = float(index)

        if block.kind == "real_scalar":
            result[block.name] = scale
        elif block.kind == "complex_scalar":
            result[block.name] = scale + 1j * (scale + 0.5)
        elif block.kind == "complex_matrix":
            real = np.full(block.shape, scale, dtype=float)
            imag = np.full(block.shape, scale + 0.5, dtype=float)
            result[block.name] = real + 1j * imag
        else:
            raise AssertionError(block.kind)

    return result


class BetaVectorTests(unittest.TestCase):
    def test_derivative_has_same_real_layout_as_state(self) -> None:
        vector, layout = pack_uv_state(_shared_state())
        betas = _mock_betas(layout)

        derivative = beta_values_to_derivative(
            betas,
            layout,
        )

        self.assertEqual(derivative.shape, vector.shape)
        self.assertTrue(np.all(np.isfinite(derivative)))

    def test_loop_factor_is_applied_exactly_once(self) -> None:
        _, layout = pack_uv_state(_shared_state())
        betas = _mock_betas(layout)

        raw = beta_values_to_derivative(
            betas,
            layout,
            divide_by_loop_factor=False,
        )
        physical = beta_values_to_derivative(
            betas,
            layout,
            divide_by_loop_factor=True,
        )

        np.testing.assert_allclose(
            physical,
            raw / LOOP_FACTOR,
            rtol=0.0,
            atol=0.0,
        )

    def test_complex_matrix_packing_matches_state_vector_convention(self) -> None:
        _, layout = pack_uv_state(_shared_state())
        betas = _mock_betas(layout)

        raw = beta_values_to_derivative(
            betas,
            layout,
            divide_by_loop_factor=False,
        )

        h_block = layout.block("h")
        h_beta = np.asarray(betas["h"], dtype=complex)
        n = h_beta.size
        block_values = raw[h_block.start:h_block.stop]

        np.testing.assert_allclose(
            block_values[:n],
            h_beta.real.reshape(-1),
            rtol=0.0,
            atol=0.0,
        )
        np.testing.assert_allclose(
            block_values[n:],
            h_beta.imag.reshape(-1),
            rtol=0.0,
            atol=0.0,
        )

    def test_missing_beta_is_rejected(self) -> None:
        _, layout = pack_uv_state(_shared_state())
        betas = _mock_betas(layout)
        del betas["lambda5"]

        with self.assertRaisesRegex(
            ValueError,
            "missing state variables: lambda5",
        ):
            validate_beta_keys(betas, layout)

    def test_extra_beta_is_rejected_by_default(self) -> None:
        _, layout = pack_uv_state(_shared_state())
        betas = _mock_betas(layout)
        betas["notAStateVariable"] = 1.0

        with self.assertRaisesRegex(
            ValueError,
            "variables absent from the state layout",
        ):
            validate_beta_keys(betas, layout)

    def test_real_state_variable_rejects_complex_beta(self) -> None:
        _, layout = pack_uv_state(_shared_state())
        betas = _mock_betas(layout)
        betas["g2"] = 1.0 + 0.2j

        with self.assertRaisesRegex(
            ValueError,
            "must be real",
        ):
            beta_values_to_derivative(betas, layout)


if __name__ == "__main__":
    unittest.main()
