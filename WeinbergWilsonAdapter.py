from __future__ import annotations

"""
Construct the Weinberg-operator Wilson tensor C_ijab used by the general
psi^2 phi^2 RGE engine.

The electroweak Weinberg structure is

    kappa (L . H)^2
      = kappa (nu H0 - e H+)^2.

The RGE engine uses real scalar components and a symmetric coefficient tensor
C_ijab.  We therefore introduce commuting placeholders for the two Weyl
components of L and the real components of H, expand the operator, and define

    C_ijab = (1/4)
             d^4 P /
             (d psi_i d psi_j d phi_a d phi_b),

where

    P = kappa (nu H0 - e H+)^2.

The factor 1/4 automatically accounts for the two symmetric fermion slots and
the two symmetric scalar slots.  It reproduces the legacy neutral conversion

    C_RR = kappa/2
    C_RI = C_IR = i kappa/2
    C_II = -kappa/2.
"""

from dataclasses import dataclass
from itertools import product
from typing import Mapping

import sympy as sp

from GeneralWeinbergRGEGenerator import RGEModel
from T3RGETensors import FermionBasis


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
    """Return the expanded Weinberg polynomial and global placeholder symbols."""

    lepton_block = fermion_basis.block(embedding.lepton_name)

    if lepton_block.fermion.su2_dimension != 2:
        raise ValueError("The Weinberg operator requires an SU(2) lepton doublet.")

    higgs_block = scalar_model.block(embedding.higgs_name)

    if higgs_block.scalar.su2_dimension != 2:
        raise ValueError("The Weinberg operator requires an SU(2) Higgs doublet.")

    fermion_symbols = {
        index: sp.Symbol(f"psi{index}", commutative=True)
        for index in range(1, fermion_basis.dimension + 1)
    }

    scalar_symbols = {
        index: sp.Symbol(f"phi{index}", real=True)
        for index in range(1, scalar_model.total_real_scalar_dimension + 1)
    }

    nu_index = fermion_basis.global_index(embedding.lepton_name, 1)
    e_index = fermion_basis.global_index(embedding.lepton_name, 2)

    # H is ordered in the usual SU(2) weight basis:
    #
    #   component 1 = H+
    #   component 2 = H0.
    #
    # Each complex component occupies the local real pair (R,I).
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
        sp.sympify(kappa) * (nu * h_zero - electron * h_plus) ** 2
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

    lepton_block = fermion_basis.block(embedding.lepton_name)
    lepton_indices = tuple(
        fermion_basis.global_index(embedding.lepton_name, component)
        for component in (1, 2)
    )

    higgs_block = scalar_model.block(embedding.higgs_name)
    higgs_indices = tuple(
        higgs_block.local_to_global(local)
        for local in range(1, 5)
    )

    result: dict[tuple[int, int, int, int], sp.Expr] = {}

    # Only L,L,H,H can be nonzero at the matching scale, so there is no reason
    # to scan all BSM fermion/scalar indices.
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

        value = sp.simplify(derivative) if simplify else derivative

        if value != 0:
            result[(i, j, a, b)] = value

    return result


def validate_weinberg_tensor_symmetry(
    tensor: Mapping[tuple[int, int, int, int], sp.Expr],
) -> None:
    """Require symmetry under i<->j and a<->b."""

    for (i, j, a, b), value in tensor.items():
        if sp.simplify(tensor.get((j, i, a, b), 0) - value) != 0:
            raise AssertionError(f"Fermion symmetry failed at {(i,j,a,b)}.")

        if sp.simplify(tensor.get((i, j, b, a), 0) - value) != 0:
            raise AssertionError(f"Scalar symmetry failed at {(i,j,a,b)}.")


def _self_check() -> None:
    """Regression against the legacy neutral C <-> kappa mapping."""

    from GeneralWeinbergRGEGenerator_complete import RGEModel
    from T3RGETensors import WeylFermion, FermionBasis

    kappa = sp.Symbol("kappa")

    scalar_model = RGEModel.t3(
        d_s1=2,
        y_s1=-sp.Rational(1, 2),
        d_s2=2,
        y_s2=sp.Rational(1, 2),
    )

    fermion_basis = FermionBasis(
        fermions=(
            WeylFermion(
                name="L",
                su2_dimension=2,
                hypercharge=-sp.Rational(1, 2),
            ),
            WeylFermion(
                name="F",
                su2_dimension=1,
                hypercharge=0,
            ),
        )
    )

    tensor = build_weinberg_wilson_tensor(
        scalar_model,
        fermion_basis,
        kappa,
    )

    validate_weinberg_tensor_symmetry(tensor)

    nu = fermion_basis.global_index("L", 1)
    e = fermion_basis.global_index("L", 2)

    H = scalar_model.block("H")

    hp_r = H.local_to_global(1)
    hp_i = H.local_to_global(2)
    h0_r = H.local_to_global(3)
    h0_i = H.local_to_global(4)

    # Neutral-neutrino block: legacy mapping.
    assert sp.simplify(tensor[(nu, nu, h0_r, h0_r)] - kappa / 2) == 0
    assert sp.simplify(tensor[(nu, nu, h0_r, h0_i)] - sp.I * kappa / 2) == 0
    assert sp.simplify(tensor[(nu, nu, h0_i, h0_r)] - sp.I * kappa / 2) == 0
    assert sp.simplify(tensor[(nu, nu, h0_i, h0_i)] + kappa / 2) == 0

    # Charged block from + kappa e e H+ H+.
    assert sp.simplify(tensor[(e, e, hp_r, hp_r)] - kappa / 2) == 0
    assert sp.simplify(tensor[(e, e, hp_r, hp_i)] - sp.I * kappa / 2) == 0
    assert sp.simplify(tensor[(e, e, hp_i, hp_i)] + kappa / 2) == 0

    # Mixed block from -2 kappa nu e H0 H+.
    assert sp.simplify(tensor[(nu, e, h0_r, hp_r)] + kappa / 4) == 0
    assert sp.simplify(tensor[(e, nu, h0_r, hp_r)] + kappa / 4) == 0

    # Two imaginary factors give the opposite sign relative to RR.
    assert sp.simplify(tensor[(nu, e, h0_i, hp_i)] - kappa / 4) == 0

    # One imaginary factor carries -i kappa/4.
    assert sp.simplify(tensor[(nu, e, h0_i, hp_r)] + sp.I * kappa / 4) == 0

    print(
        "WeinbergWilsonAdapter self-check passed "
        f"({len(tensor)} nonzero C_ijab components)."
    )


if __name__ == "__main__":
    _self_check()
