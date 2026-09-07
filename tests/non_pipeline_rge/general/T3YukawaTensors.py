from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import sympy as sp

from tests.non_pipeline_rge.general import FermionBasis

# Conversion between complex T3 Yukawa components and the real-scalar
# Yukawa tensor used by the master Weinberg RGE.

@dataclass(frozen=True)
class ComplexYukawaComponent:
    """One component coefficient multiplying psi_i psi_j times a complex scalar."""

    fermion_i: int
    fermion_j: int
    scalar_name: str
    scalar_component: int
    scalar_conjugated: bool
    coefficient: sp.Expr


def build_real_yukawa_tensor(
    scalar_model: RGEModel,
    fermion_basis: FermionBasis,
    components: Iterable[ComplexYukawaComponent],
    symmetrize_fermions: bool = False,
) -> dict[tuple[int, int, int], sp.Expr]:
    """Convert complex-scalar Yukawa components into y_ija in the real basis.

    For z=(R+iI)/sqrt(2):
      y_R = Y/sqrt(2),
      y_I = +i Y/sqrt(2).

    For z^*=(R-iI)/sqrt(2):
      y_R = Y/sqrt(2),
      y_I = -i Y/sqrt(2).

    The function does not add any extra 1/2 or sign associated with a chosen
    two-component-fermion Lagrangian convention.  The exported complex
    coefficient must already match Eq. (4.85)'s y convention.
    """

    result: dict[tuple[int, int, int], sp.Expr] = {}

    for item in components:
        if not 1 <= item.fermion_i <= fermion_basis.dimension:
            raise IndexError("fermion_i is outside the global fermion basis.")
        if not 1 <= item.fermion_j <= fermion_basis.dimension:
            raise IndexError("fermion_j is outside the global fermion basis.")

        block = scalar_model.block(item.scalar_name)
        if not 1 <= item.scalar_component <= block.scalar.su2_dimension:
            raise IndexError("scalar_component is outside its SU(2) multiplet.")

        local_r = 2 * item.scalar_component - 1
        local_i = 2 * item.scalar_component
        r = block.local_to_global(local_r)
        im = block.local_to_global(local_i)

        coefficient = sp.sympify(item.coefficient)
        imaginary_sign = -1 if item.scalar_conjugated else 1

        contributions = {
            (item.fermion_i, item.fermion_j, r): coefficient / SQRT2,
            (
                item.fermion_i,
                item.fermion_j,
                im,
            ): imaginary_sign * sp.I * coefficient / SQRT2,
        }

        if symmetrize_fermions and item.fermion_i != item.fermion_j:
            contributions.update(
                {
                    (item.fermion_j, item.fermion_i, r): coefficient / SQRT2,
                    (
                        item.fermion_j,
                        item.fermion_i,
                        im,
                    ): imaginary_sign * sp.I * coefficient / SQRT2,
                }
            )

        for key, value in contributions.items():
            result[key] = sp.simplify(result.get(key, sp.S.Zero) + value)

    return result


def yukawa_component_function(
    components: Mapping[tuple[int, int, int], sp.Expr],
):
    """Return the y(i,j,a) callable expected by MasterRGEInputs."""

    def y(i: int, j: int, a: int) -> sp.Expr:
        return sp.sympify(components.get((i, j, a), sp.S.Zero))

    return y
