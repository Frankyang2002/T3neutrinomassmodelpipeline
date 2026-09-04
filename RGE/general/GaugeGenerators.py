from __future__ import annotations

from dataclasses import dataclass
import sympy as sp

from RGE.general.RGEModel import RGEModel

# Gauge generators, gauge couplings, and gauge-sector helpers used by
# the general Weinberg RGE.

g1, g2, g3 = sp.symbols("g1 g2 g3")


def su2_complex_generators(dimension: int) -> tuple[sp.Matrix, sp.Matrix, sp.Matrix]:
    """Return Hermitian SU(2) generators (T1, T2, T3) for dimension 2j+1."""

    if not isinstance(dimension, int) or dimension < 1:
        raise ValueError("SU(2) representation dimension must be a positive integer.")

    j = sp.Rational(dimension - 1, 2)
    m_values = [j - position for position in range(dimension)]
    t_plus = sp.zeros(dimension)

    for column, m in enumerate(m_values):
        target_m = m + 1
        if target_m in m_values:
            row = m_values.index(target_m)
            t_plus[row, column] = sp.sqrt((j - m) * (j + m + 1))

    t_minus = t_plus.conjugate().T
    t1 = (t_plus + t_minus) / 2
    t2 = (t_plus - t_minus) / (2 * sp.I)
    t3 = sp.diag(*m_values)

    return tuple(sp.simplify(generator) for generator in (t1, t2, t3))


def real_scalar_generator(generator: sp.Matrix) -> sp.Matrix:
    """Convert a complex generator to the interleaved (R1,I1,R2,I2,...) basis."""

    dimension = generator.rows
    theta = sp.zeros(2 * dimension)

    for p in range(dimension):
        for q in range(dimension):
            entry = sp.expand_complex(generator[p, q])
            real_part = sp.re(entry)
            imaginary_part = sp.im(entry)

            r_p, i_p = 2 * p, 2 * p + 1
            r_q, i_q = 2 * q, 2 * q + 1

            theta[r_p, r_q] += sp.I * imaginary_part
            theta[r_p, i_q] += sp.I * real_part
            theta[i_p, r_q] += -sp.I * real_part
            theta[i_p, i_q] += sp.I * imaginary_part

    return sp.simplify(theta)


def su2_real_scalar_generators(dimension: int) -> tuple[sp.Matrix, sp.Matrix, sp.Matrix]:
    """Return the three SU(2) generators in the real-scalar basis."""

    return tuple(
        real_scalar_generator(generator)
        for generator in su2_complex_generators(dimension)
    )


def u1_real_scalar_generator(dimension: int, hypercharge) -> sp.Matrix:
    """Return the U(1)_Y generator in the real-scalar basis."""

    return real_scalar_generator(sp.Rational(hypercharge) * sp.eye(dimension))


def embed_local_generator(
    model: RGEModel,
    scalar_name: str,
    local_generator: sp.Matrix,
) -> sp.Matrix:
    """Embed one multiplet generator into the model's full real-scalar basis."""

    block = model.block(scalar_name)

    if local_generator.rows != block.real_dimension:
        raise ValueError(
            f"Generator for {scalar_name} has dimension {local_generator.rows}, "
            f"expected {block.real_dimension}."
        )

    result = sp.zeros(model.total_real_scalar_dimension)
    start = block.first - 1

    for row in range(block.real_dimension):
        for column in range(block.real_dimension):
            result[start + row, start + column] = local_generator[row, column]

    return result


def global_su2_generators(model: RGEModel) -> tuple[sp.Matrix, sp.Matrix, sp.Matrix]:
    """Return SU(2) generators acting on every complex scalar in the model."""

    generators = [sp.zeros(model.total_real_scalar_dimension) for _ in range(3)]

    for scalar in model.scalars:
        local = su2_real_scalar_generators(scalar.su2_dimension)

        for generator_index in range(3):
            generators[generator_index] += embed_local_generator(
                model,
                scalar.name,
                local[generator_index],
            )

    return tuple(sp.simplify(generator) for generator in generators)


def global_u1_generator(model: RGEModel) -> sp.Matrix:
    """Return the U(1)_Y generator acting on the full real-scalar basis."""

    result = sp.zeros(model.total_real_scalar_dimension)

    for scalar in model.scalars:
        local = u1_real_scalar_generator(
            scalar.su2_dimension,
            scalar.hypercharge,
        )
        result += embed_local_generator(model, scalar.name, local)

    return sp.simplify(result)


@dataclass(frozen=True)
class GaugeSector:
    """One gauge-group contribution to Eq. (4.85).

    ``scalar_generators`` and ``fermion_generators`` contain matching matrices
    for every adjoint generator A of the gauge group.
    """

    coupling: sp.Expr
    scalar_generators: tuple[sp.Matrix, ...]
    fermion_generators: tuple[sp.Matrix, ...]

    def __post_init__(self) -> None:
        if len(self.scalar_generators) != len(self.fermion_generators):
            raise ValueError(
                "Scalar and fermion generator lists must contain the same "
                "number of gauge generators."
            )


def gauge_sectors_from_generators(
    model: RGEModel,
    fermion_su2_generators: tuple[sp.Matrix, sp.Matrix, sp.Matrix],
    fermion_u1_generator: sp.Matrix,
) -> tuple[GaugeSector, GaugeSector]:
    """Build the SU(2)_L and U(1)_Y sectors used by Eq. (4.85).

    This helper assumes the fermion basis supplied by the caller already
    contains every fermion component relevant to C_ijab and y_ija.
    """

    su2_sector = GaugeSector(
        coupling=g2,
        scalar_generators=global_su2_generators(model),
        fermion_generators=fermion_su2_generators,
    )

    u1_sector = GaugeSector(
        coupling=g1,
        scalar_generators=(global_u1_generator(model),),
        fermion_generators=(fermion_u1_generator,),
    )

    return su2_sector, u1_sector


def quadratic_casimir_matrix(
    generators: tuple[sp.Matrix, ...],
) -> sp.Matrix:
    r"""Return C2 = sum_A T^A T^A for one representation.

    This implements the matrix definitions in Eqs. (A.7) and (A.9).
    """

    if not generators:
        raise ValueError("At least one generator is required.")

    dimension = generators[0].rows
    result = sp.zeros(dimension)

    for generator in generators:
        if generator.shape != (dimension, dimension):
            raise ValueError("All generators must have the same square dimension.")
        result += generator * generator

    return sp.simplify(result)
