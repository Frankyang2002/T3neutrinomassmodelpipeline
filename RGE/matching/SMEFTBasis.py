from __future__ import annotations

"""One-generation SM basis and Yukawa tensor used by production RGE stages."""

import sympy as sp

from RGE.general.FermionBasis import FermionBasis, WeylFermion
from RGE.general.RGEModel import ComplexScalar, RGEModel


ye, yu, yd = sp.symbols("ye yu yd")


def _add_symmetric_yukawa(
    entries: dict[tuple[int, int, int], sp.Expr],
    i: int,
    j: int,
    a: int,
    value: sp.Expr,
) -> None:
    entries[(i, j, a)] = value
    entries[(j, i, a)] = value


def build_sm_yukawa(
    scalar_model: RGEModel,
    fermion_basis: FermionBasis,
):
    """Build the one-generation SM y_ija tensor including colour multiplicity."""

    entries: dict[tuple[int, int, int], sp.Expr] = {}
    h = scalar_model.block("H")
    hp_r = h.local_to_global(1)
    hp_i = h.local_to_global(2)
    h0_r = h.local_to_global(3)
    h0_i = h.local_to_global(4)
    root2 = sp.sqrt(2)

    nu = fermion_basis.global_index("L", 1)
    e = fermion_basis.global_index("L", 2)
    ec = fermion_basis.global_index("eC", 1)

    _add_symmetric_yukawa(
        entries,
        nu,
        ec,
        hp_r,
        ye / root2,
    )
    _add_symmetric_yukawa(
        entries,
        nu,
        ec,
        hp_i,
        -sp.I * ye / root2,
    )
    _add_symmetric_yukawa(
        entries,
        e,
        ec,
        h0_r,
        ye / root2,
    )
    _add_symmetric_yukawa(
        entries,
        e,
        ec,
        h0_i,
        -sp.I * ye / root2,
    )

    for colour in range(1, 4):
        q_up = fermion_basis.global_index(
            f"Q{colour}",
            1,
        )
        q_down = fermion_basis.global_index(
            f"Q{colour}",
            2,
        )
        uc = fermion_basis.global_index(
            f"uC{colour}",
            1,
        )
        dc = fermion_basis.global_index(
            f"dC{colour}",
            1,
        )

        _add_symmetric_yukawa(
            entries,
            q_up,
            uc,
            h0_r,
            yu / root2,
        )
        _add_symmetric_yukawa(
            entries,
            q_up,
            uc,
            h0_i,
            sp.I * yu / root2,
        )
        _add_symmetric_yukawa(
            entries,
            q_down,
            uc,
            hp_r,
            -yu / root2,
        )
        _add_symmetric_yukawa(
            entries,
            q_down,
            uc,
            hp_i,
            -sp.I * yu / root2,
        )

        _add_symmetric_yukawa(
            entries,
            q_up,
            dc,
            hp_r,
            yd / root2,
        )
        _add_symmetric_yukawa(
            entries,
            q_up,
            dc,
            hp_i,
            -sp.I * yd / root2,
        )
        _add_symmetric_yukawa(
            entries,
            q_down,
            dc,
            h0_r,
            yd / root2,
        )
        _add_symmetric_yukawa(
            entries,
            q_down,
            dc,
            h0_i,
            -sp.I * yd / root2,
        )

    def yukawa(i: int, j: int, a: int) -> sp.Expr:
        return sp.sympify(
            entries.get((i, j, a), sp.S.Zero)
        )

    return yukawa


def build_sm_eft() -> tuple[RGEModel, FermionBasis]:
    """Construct the active one-generation SM EFT basis below the T3 threshold."""

    scalar_model = RGEModel(
        (
            ComplexScalar(
                name="H",
                su2_dimension=2,
                hypercharge=sp.Rational(1, 2),
            ),
        )
    )

    fermions = [
        WeylFermion(
            name="L",
            su2_dimension=2,
            hypercharge=-sp.Rational(1, 2),
        ),
        WeylFermion(
            name="eC",
            su2_dimension=1,
            hypercharge=sp.Integer(1),
        ),
    ]

    for colour in range(1, 4):
        fermions.extend(
            [
                WeylFermion(
                    name=f"Q{colour}",
                    su2_dimension=2,
                    hypercharge=sp.Rational(1, 6),
                ),
                WeylFermion(
                    name=f"uC{colour}",
                    su2_dimension=1,
                    hypercharge=-sp.Rational(2, 3),
                ),
                WeylFermion(
                    name=f"dC{colour}",
                    su2_dimension=1,
                    hypercharge=sp.Rational(1, 3),
                ),
            ]
        )

    return scalar_model, FermionBasis(
        fermions=tuple(fermions)
    )
