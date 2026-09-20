# Defines the fermion basis for the ijkl part of our RGEs
# We assign the indices to the correct fields and whichever component it is

from __future__ import annotations
from dataclasses import dataclass

import sympy as sp

from RGE.general.GaugeGenerators import (
    GaugeSector,
    g1,
    g2,
    global_su2_generators,
    global_u1_generator,
    su2_complex_generators,
)
from RGE.general.ScalarBasis import ScalarBasis


@dataclass(frozen=True)
class WeylFermion:
    """Weyl multiplet used by the y_ija/C_ijab convention."""

    name: str
    su2_dimension: int
    hypercharge: sp.Expr
    multiplicity: int = 1 # How many flavours, where each flavour gets a block
    conjugated_representation: bool = False

    def __post_init__(self) -> None:
        if self.su2_dimension < 1:
            raise ValueError("SU(2) representation dimension must be positive.")
        if self.multiplicity < 1:
            raise ValueError("Fermion multiplicity must be positive.")


@dataclass(frozen=True)
class FermionBasisBlock:
    """One fermion multiplet/flavour copy inside the global Weyl basis."""

    fermion: WeylFermion
    copy: int
    first: int
    last: int

    @property
    def indices(self) -> range:
        return range(self.first, self.last + 1)


@dataclass
class FermionBasis:
    """Global basis and gauge generators for all left-handed Weyl fermions.
    It will have the dimension of 
    sum of all fermions with dim*multiplicity"""

    fermions: tuple[WeylFermion, ...]

    
    def __post_init__(self) -> None:
        start = 1
        blocks: list[FermionBasisBlock] = []

        for fermion in self.fermions:
            for copy in range(1, fermion.multiplicity + 1):
                stop = start + fermion.su2_dimension - 1
                blocks.append(
                    FermionBasisBlock(
                        fermion=fermion,
                        copy=copy,
                        first=start,
                        last=stop,
                    )
                )
                start = stop + 1

        self.blocks = tuple(blocks)
        self.dimension = start - 1


    def block(self, name: str, copy: int = 1) -> FermionBasisBlock:
        '''We get the fermion basis block from the name'''
        for block in self.blocks:
            if block.fermion.name == name and block.copy == copy:
                return block
        raise KeyError(f"No fermion block {name!r}, copy={copy}.")

    def global_index(
        self,
        name: str,
        component: int,
        copy: int = 1,
    ) -> int:
        '''We get the global index from the 
        species + flavour copy + su2 component given'''
        block = self.block(name, copy)
        if not 1 <= component <= block.fermion.su2_dimension:
            raise IndexError(
                f"{name} component must lie in "
                f"1,...,{block.fermion.su2_dimension}."
            )
        return block.first + component - 1


def _fermion_local_su2_generators(
    fermion: WeylFermion,
) -> tuple[sp.Matrix, sp.Matrix, sp.Matrix]:
    """Return generators for the requested Weyl representation."""

    generators = su2_complex_generators(fermion.su2_dimension)

    if fermion.conjugated_representation:
        generators = tuple(
            -generator.conjugate()
            for generator in generators
        )

    return generators


def fermion_global_su2_generators(
    basis: FermionBasis,
) -> tuple[sp.Matrix, sp.Matrix, sp.Matrix]:
    """Embed SU(2)_L generators into the full Weyl-fermion basis.
    We get a block diagonal of generators for each field"""

    result = [sp.zeros(basis.dimension) for _ in range(3)]

    for block in basis.blocks:
        local = _fermion_local_su2_generators(block.fermion)
        offset = block.first - 1

        for generator_number in range(3):
            for row in range(block.fermion.su2_dimension):
                for column in range(block.fermion.su2_dimension):
                    result[generator_number][
                        offset + row,
                        offset + column,
                    ] = local[generator_number][row, column]

    return tuple(sp.simplify(generator) for generator in result)


def fermion_global_u1_generator(basis: FermionBasis) -> sp.Matrix:
    """Return the U(1)_Y generator in the global Weyl-fermion basis."""

    result = sp.zeros(basis.dimension)

    for block in basis.blocks:
        hypercharge = sp.sympify(block.fermion.hypercharge)
        if block.fermion.conjugated_representation:
            hypercharge = -hypercharge

        for index in block.indices:
            result[index - 1, index - 1] = hypercharge

    return result


def build_gauge_sectors(
    scalar_model: ScalarBasis,
    fermion_basis: FermionBasis,
) -> tuple[GaugeSector, GaugeSector]:
    """Build the SU(2)_L and U(1)_Y sectors required by the master RGE."""

    return (
        GaugeSector(
            coupling=g2,
            scalar_generators=global_su2_generators(scalar_model),
            fermion_generators=fermion_global_su2_generators(
                fermion_basis
            ),
        ),
        GaugeSector(
            coupling=g1,
            scalar_generators=(global_u1_generator(scalar_model),),
            fermion_generators=(
                fermion_global_u1_generator(fermion_basis),
            ),
        ),
    )
