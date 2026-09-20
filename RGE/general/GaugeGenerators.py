# We get gauge generators from our scalar blocks 
# We get the 3 SU(2) generators

from __future__ import annotations

from dataclasses import dataclass
import sympy as sp

from RGE.general.ScalarBasis import ScalarBasis

g1, g2, g3 = sp.symbols("g1 g2 g3") # Symbolic gauge couplings

def su2_complex_generators(dimension: int) -> tuple[sp.Matrix, sp.Matrix, sp.Matrix]:
    """Return Hermitian SU(2) generators (T1, T2, T3) for dimension 2j+1."""

    if not isinstance(dimension, int) or dimension < 1:
        raise ValueError("SU(2) representation dimension must be a positive integer.")

    j = sp.Rational(dimension - 1, 2) # j = (d - 1)/2
    m_values = [j - position for position in range(dimension)] # m = j,...-j
    t_plus = sp.zeros(dimension)

    # We get a t_plus matrix, where for each value of m and j
    # we have the raising operator
    # T+|j,m>=\sqrt{(j-m)(j+m+1)}|j.m+1> 
    for column, m in enumerate(m_values):
        target_m = m + 1
        if target_m in m_values:
            row = m_values.index(target_m)
            t_plus[row, column] = sp.sqrt((j - m) * (j + m + 1))

    t_minus = t_plus.conjugate().T # T- is just the conjugate of T+
    # From our defined diagonal basis for T3, we have the following
    t1 = (t_plus + t_minus) / 2 
    t2 = (t_plus - t_minus) / (2 * sp.I)
    t3 = sp.diag(*m_values)

    # Now we have all our generators
    return tuple(sp.simplify(generator) for generator in (t1, t2, t3))


def real_scalar_generator(generator: sp.Matrix) -> sp.Matrix:
    """Convert a complex generator to the interleaved (R1,I1,R2,I2,...) basis.
    So just complex generator -> real basis
    Gives us a 2d x 2d matrix"""

    dimension = generator.rows
    theta = sp.zeros(2 * dimension)

    # Make our hermitian matrix and split our real and imaginary components into separate real fields
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
    """Return the three SU(2) generators in the real-scalar basis"""

    return tuple(
        real_scalar_generator(generator)
        for generator in su2_complex_generators(dimension)
    )


def u1_real_scalar_generator(dimension: int, hypercharge) -> sp.Matrix:
    """Return the U(1)_Y generator in the real-scalar basis."""

    return real_scalar_generator(sp.Rational(hypercharge) * sp.eye(dimension))


def embed_local_generator(
    model: ScalarBasis,
    scalar_name: str,
    local_generator: sp.Matrix,
) -> sp.Matrix:
    """Embed one multiplet generator into the model's full real-scalar basis.
    We make sure the generator is in the right form for a field.
    Eg: if H: 1-4, S:5-8. 
    Our block diagonal generators will be 0 for the ones it is not applying on
    So applying a generator on scalar will only have values in the 5-8 indices"""

    block = model.block(scalar_name)

    if local_generator.rows != block.real_dimension:
        raise ValueError(
            f"Generator for {scalar_name} has dimension {local_generator.rows}, "
            f"expected {block.real_dimension}."
        )

    result = sp.zeros(model.total_real_scalar_dimension)
    start = block.first - 1

    # We just add the local generator to the right indices
    for row in range(block.real_dimension):
        for column in range(block.real_dimension):
            result[start + row, start + column] = local_generator[row, column]

    return result


def global_su2_generators(model: ScalarBasis) -> tuple[sp.Matrix, sp.Matrix, sp.Matrix]:
    """Gets full model-wide SU(2) generators. We construct local generators
    and place them into the right block in a matrix."""

    # Initialise our generators to be 0
    generators = [sp.zeros(model.total_real_scalar_dimension) for _ in range(3)]

    # For each scalar, we fill each generator spaces with its local generator
    for scalar in model.scalars:
        local = su2_real_scalar_generators(scalar.su2_dimension)

        for generator_index in range(3):
            generators[generator_index] += embed_local_generator(
                model,
                scalar.name,
                local[generator_index],
            )

    return tuple(sp.simplify(generator) for generator in generators)


def global_u1_generator(model: ScalarBasis) -> sp.Matrix:
    """Return the U(1)_Y generator acting on the full real-scalar basis."""

    result = sp.zeros(model.total_real_scalar_dimension)
    # For each scalar fill in the single generator with its local generator
    for scalar in model.scalars:
        local = u1_real_scalar_generator(
            scalar.su2_dimension,
            scalar.hypercharge,
        )
        result += embed_local_generator(model, scalar.name, local)

    return sp.simplify(result)


@dataclass(frozen=True)
class GaugeSector:
    """
    Gets our gauge sector
    Couplings, scalar generators and fermion generators
    We get one for U1 and SU2
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


def quadratic_casimir_matrix(
    generators: tuple[sp.Matrix, ...],
) -> sp.Matrix:
    r"""Return C2 = sum_A T^A T^A for one representation.
    It will be block diagonal  
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
