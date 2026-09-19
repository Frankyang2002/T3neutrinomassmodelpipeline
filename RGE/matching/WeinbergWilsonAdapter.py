from __future__ import annotations

"""Construct the Weinberg-operator C_ijab tensor in the production RGE basis."""

from dataclasses import dataclass
from itertools import product
from typing import Mapping

import sympy as sp

from RGE.general.FermionBasis import FermionBasis
from RGE.general.RGEModel import RGEModel


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
