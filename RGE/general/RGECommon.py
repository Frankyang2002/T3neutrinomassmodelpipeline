from __future__ import annotations

from collections.abc import Callable
import sympy as sp

from RGE.general.RGEModel import RGEModel

# Shared symbolic definitions and validation helpers for the general RGE.

HALF = sp.Rational(1, 2)


C = sp.IndexedBase("C")


QuarticComponent = Callable[[int, int, int, int], sp.Expr]


def _validate_output_component(
    model: RGEModel,
    output_component: tuple[int, int, int, int],
) -> None:
    if len(output_component) != 4:
        raise ValueError("output_component must have the form (i, j, a, b).")

    _, _, a, b = output_component
    n = model.total_real_scalar_dimension

    if not 1 <= a <= n:
        raise IndexError(f"Scalar index a={a} is outside 1,...,{n}.")
    if not 1 <= b <= n:
        raise IndexError(f"Scalar index b={b} is outside 1,...,{n}.")
