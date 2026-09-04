from __future__ import annotations

from dataclasses import dataclass
import sympy as sp

@dataclass(frozen=True)
class ComplexScalar:
    """One complex colour-singlet scalar multiplet."""

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
        return 2 * self.scalar.su2_dimension

    @property
    def indices(self) -> range:
        return range(self.first, self.last + 1)

    def local_to_global(self, local_index: int) -> int:
        if not 1 <= local_index <= self.real_dimension:
            raise IndexError(
                f"{self.scalar.name} local real index must lie in "
                f"1,...,{self.real_dimension}."
            )
        return self.first + local_index - 1


@dataclass
class RGEModel:
    """Minimal model data required by the generic scalar-side RGE machinery."""

    scalars: tuple[ComplexScalar, ...]

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
