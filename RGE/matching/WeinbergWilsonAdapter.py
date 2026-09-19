from __future__ import annotations

"""Construct the Weinberg-operator C_ijab tensor in the production RGE basis."""

from dataclasses import dataclass
from itertools import product
from typing import Mapping

import sympy as sp

from RGE.general.FermionBasis import FermionBasis, WeylFermion
from RGE.general.RGEModel import ComplexScalar, RGEModel


# ---------------------------------------------------------------------------
# Sparse Wilson-coefficient lookup
# Consolidated from the former standalone helper module.
# ---------------------------------------------------------------------------

class SparseWilsonLookup:
    """Sparse C_ijab lookup with zero for components absent from the mapping."""

    def __init__(
        self,
        components: Mapping[tuple[int, int, int, int], sp.Expr],
    ):
        self.components = components

    def __getitem__(
        self,
        key: tuple[int, int, int, int],
    ) -> sp.Expr:
        return sp.sympify(
            self.components.get(key, sp.S.Zero)
        )

    def __call__(
        self,
        i: int,
        j: int,
        a: int,
        b: int,
    ) -> sp.Expr:
        return self[i, j, a, b]


def wilson_component_function(
    components: Mapping[tuple[int, int, int, int], sp.Expr],
) -> SparseWilsonLookup:
    """Return a sparse Wilson-coefficient lookup callable."""

    return SparseWilsonLookup(components)


# ---------------------------------------------------------------------------
# One-generation SM EFT basis and Yukawa tensor
# Consolidated from the former SMEFTBasis helper.
# ---------------------------------------------------------------------------

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

@dataclass(frozen=True)
class WeinbergEmbedding:
    """Names of the SM blocks inside the global RGE bases."""

    lepton_name: str = "L"
    higgs_name: str = "H"


def _complex_scalar_from_real_pair(
    real_symbol: sp.Symbol,
    imag_symbol: sp.Symbol,
) -> sp.Expr:
    return (real_symbol + sp.I * imag_symbol) / sp.sqrt(2)


def weinberg_polynomial(
    scalar_model: RGEModel,
    fermion_basis: FermionBasis,
    kappa: sp.Expr,
    embedding: WeinbergEmbedding = WeinbergEmbedding(),
) -> tuple[sp.Expr, dict[int, sp.Symbol], dict[int, sp.Symbol]]:
    """Return the expanded Weinberg polynomial and global placeholders."""

    lepton_block = fermion_basis.block(embedding.lepton_name)

    if lepton_block.fermion.su2_dimension != 2:
        raise ValueError(
            "The Weinberg operator requires an SU(2) lepton doublet."
        )

    higgs_block = scalar_model.block(embedding.higgs_name)

    if higgs_block.scalar.su2_dimension != 2:
        raise ValueError(
            "The Weinberg operator requires an SU(2) Higgs doublet."
        )

    fermion_symbols = {
        index: sp.Symbol(f"psi{index}", commutative=True)
        for index in range(1, fermion_basis.dimension + 1)
    }

    scalar_symbols = {
        index: sp.Symbol(f"phi{index}", real=True)
        for index in range(
            1,
            scalar_model.total_real_scalar_dimension + 1,
        )
    }

    nu_index = fermion_basis.global_index(
        embedding.lepton_name,
        1,
    )
    e_index = fermion_basis.global_index(
        embedding.lepton_name,
        2,
    )

    hp_r = higgs_block.local_to_global(1)
    hp_i = higgs_block.local_to_global(2)
    h0_r = higgs_block.local_to_global(3)
    h0_i = higgs_block.local_to_global(4)

    h_plus = _complex_scalar_from_real_pair(
        scalar_symbols[hp_r],
        scalar_symbols[hp_i],
    )
    h_zero = _complex_scalar_from_real_pair(
        scalar_symbols[h0_r],
        scalar_symbols[h0_i],
    )

    nu = fermion_symbols[nu_index]
    electron = fermion_symbols[e_index]

    polynomial = sp.expand(
        sp.sympify(kappa)
        * (nu * h_zero - electron * h_plus) ** 2
    )

    return polynomial, fermion_symbols, scalar_symbols


def build_weinberg_wilson_tensor(
    scalar_model: RGEModel,
    fermion_basis: FermionBasis,
    kappa: sp.Expr,
    embedding: WeinbergEmbedding = WeinbergEmbedding(),
    simplify: bool = True,
) -> dict[tuple[int, int, int, int], sp.Expr]:
    """Construct all nonzero C_ijab components in the global real basis."""

    polynomial, fermions, scalars = weinberg_polynomial(
        scalar_model,
        fermion_basis,
        kappa,
        embedding,
    )

    lepton_indices = tuple(
        fermion_basis.global_index(
            embedding.lepton_name,
            component,
        )
        for component in (1, 2)
    )

    higgs_block = scalar_model.block(embedding.higgs_name)
    higgs_indices = tuple(
        higgs_block.local_to_global(local)
        for local in range(1, 5)
    )

    result: dict[tuple[int, int, int, int], sp.Expr] = {}

    for i, j, a, b in product(
        lepton_indices,
        lepton_indices,
        higgs_indices,
        higgs_indices,
    ):
        derivative = sp.diff(
            polynomial,
            fermions[i],
            fermions[j],
            scalars[a],
            scalars[b],
        ) / 4

        value = (
            sp.simplify(derivative)
            if simplify
            else derivative
        )

        if value != 0:
            result[(i, j, a, b)] = value

    return result


def validate_weinberg_tensor_symmetry(
    tensor: Mapping[tuple[int, int, int, int], sp.Expr],
) -> None:
    """Require symmetry under i<->j and a<->b."""

    for (i, j, a, b), value in tensor.items():
        if sp.simplify(
            tensor.get((j, i, a, b), 0) - value
        ) != 0:
            raise AssertionError(
                f"Fermion symmetry failed at {(i, j, a, b)}."
            )

        if sp.simplify(
            tensor.get((i, j, b, a), 0) - value
        ) != 0:
            raise AssertionError(
                f"Scalar symmetry failed at {(i, j, a, b)}."
            )
