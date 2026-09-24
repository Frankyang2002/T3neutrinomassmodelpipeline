"""Map evaluated T3 beta functions onto the real numerical state vector.

This module deliberately does not parse RGBeta's Wolfram ``InputForm`` text.
It handles the next unambiguous layer:

    evaluated beta values -> real ODE derivative vector

The convention expected here is the repository's ``report_betas`` convention,

    16*pi^2 dX/dln(mu) = beta_X^(1),

including the gauge-coupling conversion already performed by
``T3RGBetaConventionalReportBetas`` in the Wolfram runner.

A later RGBeta evaluator only needs to return a dictionary whose keys are the
state-layout coupling names and whose values have the corresponding scalar or
matrix shapes.  This module then performs all real/complex packing and divides
by 16*pi^2 exactly once.
"""

from __future__ import annotations

from typing import Any, Mapping

import numpy as np

from Numerical.StateVector import StateVectorLayout, VectorBlock


LOOP_FACTOR = 16.0 * np.pi**2


def _finite_real_scalar(name: str, value: Any) -> float:
    """Return one finite real beta value."""
    if np.ndim(value) != 0:
        raise ValueError(f"Beta {name} must be a scalar.")

    result = complex(value)

    if not (np.isfinite(result.real) and np.isfinite(result.imag)):
        raise ValueError(f"Beta {name} must be finite.")

    if not np.isclose(result.imag, 0.0, rtol=0.0, atol=1.0e-13):
        raise ValueError(
            f"Beta {name} must be real for a real state variable; "
            f"got imaginary part {result.imag}."
        )

    return float(result.real)


def _finite_complex_scalar(name: str, value: Any) -> complex:
    """Return one finite complex beta value."""
    if np.ndim(value) != 0:
        raise ValueError(f"Beta {name} must be a scalar.")

    result = complex(value)

    if not (np.isfinite(result.real) and np.isfinite(result.imag)):
        raise ValueError(f"Beta {name} must be finite.")

    return result


def _finite_complex_matrix(
    name: str,
    value: Any,
    shape: tuple[int, ...],
) -> np.ndarray:
    """Return one finite matrix-valued beta with exactly the layout shape."""
    matrix = np.asarray(value, dtype=complex)

    if matrix.shape != shape:
        raise ValueError(
            f"Beta {name} must have shape {shape}, got {matrix.shape}."
        )

    if not np.all(np.isfinite(matrix.real)) or not np.all(
        np.isfinite(matrix.imag)
    ):
        raise ValueError(f"Beta {name} must contain only finite entries.")

    return matrix


def _pack_beta_block(
    block: VectorBlock,
    value: Any,
) -> np.ndarray:
    """Pack one evaluated beta into the matching real-vector block."""

    if block.kind == "real_scalar":
        return np.asarray(
            [_finite_real_scalar(block.name, value)],
            dtype=float,
        )

    if block.kind == "complex_scalar":
        z = _finite_complex_scalar(block.name, value)
        return np.asarray([z.real, z.imag], dtype=float)

    if block.kind == "complex_matrix":
        matrix = _finite_complex_matrix(
            block.name,
            value,
            block.shape,
        )
        return np.concatenate(
            [
                matrix.real.reshape(-1),
                matrix.imag.reshape(-1),
            ]
        )

    raise ValueError(
        f"Unsupported vector block kind {block.kind!r} for {block.name}."
    )


def validate_beta_keys(
    beta_values: Mapping[str, Any],
    layout: StateVectorLayout,
    *,
    allow_extra: bool = False,
) -> None:
    """Require complete beta coverage of the numerical state layout.

    The current UV ODE state is defined by ``StateVectorLayout``.  Every block
    must therefore have exactly one evaluated beta before an ODE derivative can
    be built.  Extra keys are rejected by default so a mismatch between RGBeta
    export and numerical state cannot be silently ignored.
    """

    expected = set(layout.names)
    supplied = set(beta_values)

    missing = sorted(expected.difference(supplied))
    extra = sorted(supplied.difference(expected))

    if missing:
        raise ValueError(
            "Evaluated beta dictionary is missing state variables: "
            + ", ".join(missing)
        )

    if extra and not allow_extra:
        raise ValueError(
            "Evaluated beta dictionary contains variables absent from the "
            "state layout: "
            + ", ".join(extra)
        )


def beta_values_to_derivative(
    beta_values: Mapping[str, Any],
    layout: StateVectorLayout,
    *,
    divide_by_loop_factor: bool = True,
    allow_extra: bool = False,
) -> np.ndarray:
    """Return ``dy/dln(mu)`` in the exact layout used by ``pack_uv_state``.

    Parameters
    ----------
    beta_values:
        Evaluated one-loop beta functions keyed by RGBeta coupling name.

        By default these are assumed to obey

            16*pi^2 dX/dln(mu) = beta_X^(1).

    layout:
        State-vector layout returned by ``pack_uv_state``.

    divide_by_loop_factor:
        If True (default), divide the packed one-loop beta vector by
        ``16*pi^2``.  Set False only for an evaluator that already returns the
        physical derivative ``dX/dln(mu)``.

    allow_extra:
        Whether keys not present in the current state layout may be ignored.

    Returns
    -------
    numpy.ndarray
        A finite one-dimensional real derivative vector with length
        ``layout.size``.
    """

    validate_beta_keys(
        beta_values,
        layout,
        allow_extra=allow_extra,
    )

    derivative = np.empty(layout.size, dtype=float)

    for block in layout.blocks:
        packed = _pack_beta_block(
            block,
            beta_values[block.name],
        )

        if packed.size != block.size:
            raise RuntimeError(
                f"Internal beta packing error for {block.name}: "
                f"packed {packed.size} entries, expected {block.size}."
            )

        derivative[block.start:block.stop] = packed

    if divide_by_loop_factor:
        derivative = derivative / LOOP_FACTOR

    if not np.all(np.isfinite(derivative)):
        raise ValueError("Packed ODE derivative contains non-finite entries.")

    return derivative
