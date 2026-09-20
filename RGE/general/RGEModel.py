# Tells us all representations and indices that exists, mostly for scalars
# Our RGE formulae requires scalar indices and real scalar components
# So for d dimensional SU(2) multiplet we have 2d components
# We assign which indices goes to higgs and which scalar fields 
# And how many real scalar fields exist etc
# When we talk about scalar basis we are for example saying:
# For SU(2) H = 2, S1 = 1, S2 = 3, 
# 1-4 -> H, 5-6 -> S1, 7-12 -> S2 
# Same fields share indices, like if we have shared scalar we do not split in the interaction
# Important for interactions and making sure couplings are assigned correctly as well
from __future__ import annotations

from dataclasses import dataclass
import sympy as sp

@dataclass(frozen=True)
class ComplexScalar:
    """A class for a complex scalar multiplet"""

    name: str
    su2_dimension: int
    hypercharge: sp.Rational

    def __post_init__(self) -> None:
        if self.su2_dimension < 1:
            raise ValueError("SU(2) representation dimension must be positive.")


@dataclass(frozen=True)
class ScalarBasisBlock:
    """Location of one complex multiplet inside the global real-scalar basis."""

    scalar: ComplexScalar
    first: int
    last: int

    @property
    def real_dimension(self) -> int:
        return 2 * self.scalar.su2_dimension # Normal complex to real conversion dimensions

    @property
    def indices(self) -> range:
        return range(self.first, self.last + 1) # Get indices

    # We map real component index of a multiplet to the complete scalar tensor
    def local_to_global(self, local_index: int) -> int:
        if not 1 <= local_index <= self.real_dimension:
            raise IndexError(
                f"{self.scalar.name} local real index must lie in "
                f"1,...,{self.real_dimension}."
            )
        # Global = first + local coord - 1
        return self.first + local_index - 1


@dataclass
class RGEModel:
    """All scalar data to be used by RGE system."""

    scalars: tuple[ComplexScalar, ...]

    # This does 2 things
    # 1. Enforce unique scalar names
    # 2. We create our scalar basis based on what scalar fields we have
    def __post_init__(self) -> None:
        names = [scalar.name for scalar in self.scalars]
        if len(names) != len(set(names)):
            raise ValueError("Scalar names must be unique.")

        start = 1
        blocks: dict[str, ScalarBasisBlock] = {}

        for scalar in self.scalars:
            stop = start + 2 * scalar.su2_dimension - 1
            blocks[scalar.name] = ScalarBasisBlock(
                scalar=scalar,
                first=start,
                last=stop,
            )
            start = stop + 1

        self.blocks = blocks
        self.total_real_scalar_dimension = start - 1

    
    def block(self, scalar_name: str) -> ScalarBasisBlock:
        '''We get the scalar basis block'''
        try:
            return self.blocks[scalar_name]
        except KeyError as exc:
            raise KeyError(
                f"Unknown scalar {scalar_name!r}. "
                f"Available scalars: {tuple(self.blocks)}"
            ) from exc

    @classmethod
    def t3(
        cls,
        d_s1: int,
        y_s1,
        d_s2: int,
        y_s2,
        include_higgs: bool = True,
    ) -> "RGEModel":
        """Construct the scalar sector of a generalized T3 model."""

        scalars: list[ComplexScalar] = []

        if include_higgs:
            scalars.append(
                ComplexScalar(
                    name="H",
                    su2_dimension=2,
                    hypercharge=sp.Rational(1, 2),
                )
            )

        scalars.extend(
            [
                ComplexScalar(
                    name="S1",
                    su2_dimension=int(d_s1),
                    hypercharge=sp.Rational(y_s1),
                ),
                ComplexScalar(
                    name="S2",
                    su2_dimension=int(d_s2),
                    hypercharge=sp.Rational(y_s2),
                ),
            ]
        )

        return cls(tuple(scalars))

    @classmethod
    def t3_shared(
        cls,
        d_s: int,
        y_s=sp.Rational(1, 2),
        include_higgs: bool = True,
    ) -> "RGEModel":
        """Construct the physical one-scalar scotogenic-style scalar sector."""
        scalars: list[ComplexScalar] = []
        if include_higgs:
            scalars.append(
                ComplexScalar(
                    name="H",
                    su2_dimension=2,
                    hypercharge=sp.Rational(1, 2),
                )
            )
        scalars.append(
            ComplexScalar(
                name="S",
                su2_dimension=int(d_s),
                hypercharge=sp.Rational(y_s),
            )
        )
        return cls(tuple(scalars))

