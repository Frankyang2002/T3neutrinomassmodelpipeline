from __future__ import annotations

"""
General symbolic helpers for psi^2 phi^2 / Weinberg-operator RGE work.

This file is the generalized successor to the scotogenic-specific
WeinbergRGEGenerator.py.  The legacy file is intentionally left untouched and
can be used as a reference for conventions and previously checked special-case
results.

Current scope:
- arbitrary complex SU(2) scalar multiplets;
- dynamic global real-scalar basis;
- arbitrary SU(2) and U(1)_Y real-basis generators;
- generic lambda_abcd and C_ijab tensors;
- generic gauge and scalar-quartic contributions already implemented in the
  legacy prototype;
- automatic neutral-component identification and C -> K substitutions.

The tensor structure of the complete general psi^2 phi^2 beta function,
Eq. (4.85), and the collinear anomalous dimensions from Eqs. (A.2) and
(A.3) are implemented here.
"""

from dataclasses import dataclass
from functools import lru_cache
from fractions import Fraction
from typing import Callable, Iterable, Mapping

import sympy as sp


HALF = sp.Rational(1, 2)

g1, g2, g3 = sp.symbols("g1 g2 g3")
C = sp.IndexedBase("C")
LAMBDA = sp.IndexedBase("lambda")


# -----------------------------------------------------------------------------
# Model and scalar-basis definitions
# -----------------------------------------------------------------------------


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


# -----------------------------------------------------------------------------
# Gauge generators
# -----------------------------------------------------------------------------


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


# -----------------------------------------------------------------------------
# Global scalar generators
# -----------------------------------------------------------------------------


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


# -----------------------------------------------------------------------------
# Generic quartic tensor lambda_abcd
# -----------------------------------------------------------------------------


QuarticComponent = Callable[[int, int, int, int], sp.Expr]


def symbolic_quartic_component(
    a: int,
    b: int,
    c: int,
    d: int,
    lambda_tensor=LAMBDA,
) -> sp.Expr:
    """Return a formal lambda_abcd component without assuming a UV coupling basis."""

    return lambda_tensor[a, b, c, d]


def mapping_quartic_component(
    components: Mapping[tuple[int, int, int, int], sp.Expr],
) -> QuarticComponent:
    """Build a symmetric lambda_abcd lookup from an explicit component mapping."""

    normalized = {
        tuple(sorted(key)): sp.sympify(value)
        for key, value in components.items()
    }

    def component(a: int, b: int, c: int, d: int) -> sp.Expr:
        return normalized.get(tuple(sorted((a, b, c, d))), sp.S.Zero)

    return component


# -----------------------------------------------------------------------------
# Generic RGE contributions currently inherited from the legacy prototype
# -----------------------------------------------------------------------------


def _validate_output_component(
    model: RGEModel,
    output_component: tuple[int, int, int, int],
) -> None:
    if len(output_component) != 4:
        raise ValueError("output_component must have the form (i, j, a, b).")

    _, _, a, b = output_component
    n = model.total_real_scalar_dimension

    if not 1 <= a <= n:
        raise IndexError(f"Scalar index a={a} is outside 1,...,{n}.")
    if not 1 <= b <= n:
        raise IndexError(f"Scalar index b={b} is outside 1,...,{n}.")


def _gauge_generator_contribution(
    theta: sp.Matrix,
    gauge_coupling: sp.Symbol,
    i: int,
    j: int,
    a_index: int,
    b_index: int,
    real_dimension: int,
    coefficient,
    numerical_coefficient,
) -> sp.Expr:
    """Compute a generic symmetrized scalar-generator contribution."""

    result = sp.S.Zero

    for c_index in range(real_dimension):
        for d_index in range(real_dimension):
            c = c_index + 1
            d = d_index + 1

            symmetrized_generators = (
                theta[a_index, c_index] * theta[b_index, d_index]
                + theta[b_index, c_index] * theta[a_index, d_index]
            )

            result += (
                numerical_coefficient
                * gauge_coupling**2
                * symmetrized_generators
                * coefficient[i, j, c, d]
            )

    return sp.simplify(result)


def first_gauge_component(
    model: RGEModel,
    output_component: tuple[int, int, int, int],
    coefficient=C,
) -> sp.Expr:
    r"""Return the gauge-generator contribution already checked in the legacy code.

    Implements

      - sum_alpha g_alpha^2 sum_cd
        (theta^alpha_ac theta^alpha_bd + theta^alpha_bc theta^alpha_ad) C_ijcd.

    SU(3) is absent because the scalar multiplets considered here are colour singlets.
    """

    _validate_output_component(model, output_component)
    i, j, a, b = output_component
    a_index, b_index = a - 1, b - 1
    n = model.total_real_scalar_dimension

    result = sp.S.Zero

    for theta in global_su2_generators(model):
        result += _gauge_generator_contribution(
            theta=theta,
            gauge_coupling=g2,
            i=i,
            j=j,
            a_index=a_index,
            b_index=b_index,
            real_dimension=n,
            coefficient=coefficient,
            numerical_coefficient=-1,
        )

    result += _gauge_generator_contribution(
        theta=global_u1_generator(model),
        gauge_coupling=g1,
        i=i,
        j=j,
        a_index=a_index,
        b_index=b_index,
        real_dimension=n,
        coefficient=coefficient,
        numerical_coefficient=-1,
    )

    return sp.simplify(result)


def first_scalar_component(
    model: RGEModel,
    output_component: tuple[int, int, int, int],
    quartic_component: QuarticComponent = symbolic_quartic_component,
    coefficient=C,
) -> sp.Expr:
    r"""Return sum_cd lambda_abcd C_ijcd in the global real-scalar basis."""

    _validate_output_component(model, output_component)

    i, j, a, b = output_component
    n = model.total_real_scalar_dimension
    result = sp.S.Zero

    for c in range(1, n + 1):
        for d in range(1, n + 1):
            quartic = quartic_component(a, b, c, d)

            if quartic != 0:
                result += quartic * coefficient[i, j, c, d]

    return sp.simplify(result)



# -----------------------------------------------------------------------------
# Full psi^2 phi^2 master RGE, Eq. (4.85)
# -----------------------------------------------------------------------------


YukawaComponent = Callable[[int, int, int], sp.Expr]
ScalarAnomalousDimension = Callable[[int, int], sp.Expr]
FermionAnomalousDimension = Callable[[int, int], sp.Expr]


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


@dataclass(frozen=True)
class MasterRGEInputs:
    """Tensor inputs required to evaluate the general psi^2 phi^2 beta function."""

    fermion_dimension: int
    yukawa: YukawaComponent
    quartic: QuarticComponent
    gauge_sectors: tuple[GaugeSector, ...] = ()
    gamma_scalar: ScalarAnomalousDimension | None = None
    gamma_fermion: FermionAnomalousDimension | None = None

    def __post_init__(self) -> None:
        if self.fermion_dimension < 1:
            raise ValueError("fermion_dimension must be positive.")


def symbolic_yukawa_component(
    i: int,
    j: int,
    a: int,
    yukawa_tensor=sp.IndexedBase("y"),
) -> sp.Expr:
    """Return a formal Yukawa component y_ija."""

    return yukawa_tensor[i, j, a]


def symbolic_scalar_anomalous_dimension(
    b: int,
    d: int,
    gamma_tensor=sp.IndexedBase("gamma_s"),
) -> sp.Expr:
    """Return a formal scalar collinear anomalous dimension gamma_s^bd."""

    return gamma_tensor[b, d]


def symbolic_fermion_anomalous_dimension(
    j: int,
    k: int,
    gamma_tensor=sp.IndexedBase("gamma_f"),
) -> sp.Expr:
    """Return a formal fermion collinear anomalous dimension gamma_f^jk."""

    return gamma_tensor[j, k]


def _swap_pair(first, second):
    """Return the identity and swap permutations of a two-element pair."""

    return ((first, second), (second, first))


def _pair_product_permutations(a, b, i, j):
    """Return the four independent swaps in sigma({a,b} x {i,j})."""

    return tuple(
        (a_perm, b_perm, i_perm, j_perm)
        for a_perm, b_perm in _swap_pair(a, b)
        for i_perm, j_perm in _swap_pair(i, j)
    )


def _simultaneous_pair_permutations(a, b, i, j):
    """Return the two simultaneous pair exchanges used by Eq. (4.85)."""

    return (
        (a, b, i, j),
        (b, a, j, i),
    )


def _validate_master_dimensions(
    model: RGEModel,
    inputs: MasterRGEInputs,
    output_component: tuple[int, int, int, int],
) -> None:
    """Validate external indices and tensor dimensions for Eq. (4.85)."""

    _validate_output_component(model, output_component)
    i, j, _, _ = output_component

    if not 1 <= i <= inputs.fermion_dimension:
        raise IndexError(
            f"Fermion index i={i} is outside 1,...,{inputs.fermion_dimension}."
        )
    if not 1 <= j <= inputs.fermion_dimension:
        raise IndexError(
            f"Fermion index j={j} is outside 1,...,{inputs.fermion_dimension}."
        )

    scalar_dimension = model.total_real_scalar_dimension

    for sector in inputs.gauge_sectors:
        for theta in sector.scalar_generators:
            if theta.shape != (scalar_dimension, scalar_dimension):
                raise ValueError(
                    "A scalar gauge generator does not match the model's "
                    f"{scalar_dimension}-dimensional real-scalar basis."
                )

        for generator in sector.fermion_generators:
            if generator.shape != (
                inputs.fermion_dimension,
                inputs.fermion_dimension,
            ):
                raise ValueError(
                    "A fermion gauge generator does not match fermion_dimension."
                )


def yukawa_wavefunction_term(
    model: RGEModel,
    inputs: MasterRGEInputs,
    output_component: tuple[int, int, int, int],
    coefficient=C,
) -> sp.Expr:
    r"""Return y_ijd y^*_{kld} C_klab from Eq. (4.85)."""

    _validate_master_dimensions(model, inputs, output_component)
    i, j, a, b = output_component
    nf = inputs.fermion_dimension
    ns = model.total_real_scalar_dimension
    y = inputs.yukawa

    result = sp.S.Zero

    for k in range(1, nf + 1):
        for l in range(1, nf + 1):
            for d in range(1, ns + 1):
                result += (
                    y(i, j, d)
                    * sp.conjugate(y(k, l, d))
                    * coefficient[k, l, a, b]
                )

    return sp.simplify(result)


def scalar_pair_term(
    model: RGEModel,
    inputs: MasterRGEInputs,
    output_component: tuple[int, int, int, int],
    coefficient=C,
) -> sp.Expr:
    r"""Return the -1/2! sigma({a,b}) gauge/quartic term in Eq. (4.85).

    This contains both

      - sum_alpha g_alpha^2
        (theta_ac theta_bd + theta_bc theta_ad) C_ijcd

    and

      + sum_cd lambda_abcd C_ijcd,

    when the scalar quartic tensor is fully symmetric.
    """

    _validate_master_dimensions(model, inputs, output_component)
    i, j, a, b = output_component
    ns = model.total_real_scalar_dimension
    quartic = inputs.quartic
    result = sp.S.Zero

    for a_perm, b_perm in _swap_pair(a, b):
        for c in range(1, ns + 1):
            for d in range(1, ns + 1):
                gauge_piece = sp.S.Zero

                for sector in inputs.gauge_sectors:
                    for theta in sector.scalar_generators:
                        gauge_piece += (
                            sector.coupling**2
                            * theta[a_perm - 1, c - 1]
                            * theta[b_perm - 1, d - 1]
                        )

                result += -HALF * (
                    2 * gauge_piece - quartic(a_perm, b_perm, c, d)
                ) * coefficient[i, j, c, d]

    return sp.simplify(result)


def mixed_yukawa_gauge_term(
    model: RGEModel,
    inputs: MasterRGEInputs,
    output_component: tuple[int, int, int, int],
    coefficient=C,
) -> sp.Expr:
    r"""Return the third line's first sigma({a,b} x {i,j}) term in Eq. (4.85).

    Implements

      sum_sigma [
        2 y_jkc y^*_{klb}
        + y_jkb y^*_{klc}
        + 4 sum_alpha g_alpha^2 t^A_lj theta^A_bc
      ] C_ilac.
    """

    _validate_master_dimensions(model, inputs, output_component)
    i, j, a, b = output_component
    nf = inputs.fermion_dimension
    ns = model.total_real_scalar_dimension
    y = inputs.yukawa
    result = sp.S.Zero

    # The Yukawa contractions require the four independent external-pair swaps.
    for a_perm, b_perm, i_perm, j_perm in _pair_product_permutations(
        a, b, i, j
    ):
        for k in range(1, nf + 1):
            for l in range(1, nf + 1):
                for c in range(1, ns + 1):
                    yukawa_bracket = (
                        2
                        * y(j_perm, k, c)
                        * sp.conjugate(y(k, l, b_perm))
                        + y(j_perm, k, b_perm)
                        * sp.conjugate(y(k, l, c))
                    )

                    result += (
                        yukawa_bracket
                        * coefficient[i_perm, l, a_perm, c]
                    )

    # The gauge contraction carries the same four independent external-pair
    # swaps.  Unlike the Yukawa contraction above, it has no summed k index.
    for a_perm, b_perm, i_perm, j_perm in _pair_product_permutations(
        a, b, i, j
    ):
        for l in range(1, nf + 1):
            for c in range(1, ns + 1):
                gauge_bracket = sp.S.Zero

                for sector in inputs.gauge_sectors:
                    for t_generator, theta in zip(
                        sector.fermion_generators,
                        sector.scalar_generators,
                    ):
                        gauge_bracket += (
                            4
                            * sector.coupling**2
                            * t_generator[l - 1, j_perm - 1]
                            * theta[b_perm - 1, c - 1]
                        )

                result += (
                    gauge_bracket
                    * coefficient[i_perm, l, a_perm, c]
                )

    return sp.simplify(result)


def crossed_yukawa_term(
    model: RGEModel,
    inputs: MasterRGEInputs,
    output_component: tuple[int, int, int, int],
    coefficient=C,
) -> sp.Expr:
    r"""Return sum_sigma y_jla y_ikc C_lkbc from Eq. (4.85)."""

    _validate_master_dimensions(model, inputs, output_component)
    i, j, a, b = output_component
    nf = inputs.fermion_dimension
    ns = model.total_real_scalar_dimension
    y = inputs.yukawa
    result = sp.S.Zero

    for a_perm, b_perm, i_perm, j_perm in _simultaneous_pair_permutations(
        a, b, i, j
    ):
        for k in range(1, nf + 1):
            for l in range(1, nf + 1):
                for c in range(1, ns + 1):
                    result += (
                        y(j_perm, l, a_perm)
                        * y(i_perm, k, c)
                        * coefficient[l, k, b_perm, c]
                    )

    return sp.simplify(result)


def conjugate_coefficient_yukawa_term(
    model: RGEModel,
    inputs: MasterRGEInputs,
    output_component: tuple[int, int, int, int],
    coefficient=C,
) -> sp.Expr:
    r"""Return (2 y_ild y_jkd + y_ijd y_kld) C^*_{klab} from Eq. (4.85)."""

    _validate_master_dimensions(model, inputs, output_component)
    i, j, a, b = output_component
    nf = inputs.fermion_dimension
    ns = model.total_real_scalar_dimension
    y = inputs.yukawa
    result = sp.S.Zero

    for k in range(1, nf + 1):
        for l in range(1, nf + 1):
            for d in range(1, ns + 1):
                result += (
                    2 * y(i, l, d) * y(j, k, d)
                    + y(i, j, d) * y(k, l, d)
                ) * sp.conjugate(coefficient[k, l, a, b])

    return sp.simplify(result)


def scalar_anomalous_dimension_term(
    model: RGEModel,
    inputs: MasterRGEInputs,
    output_component: tuple[int, int, int, int],
    coefficient=C,
) -> sp.Expr:
    r"""Return sum_sigma({a,b}) gamma_{c,s}^{bd} C_ijad from Eq. (4.85).

    The ``c`` in gamma_{c,s} labels the *collinear* anomalous dimension; it is
    not a summed scalar index.
    """

    _validate_master_dimensions(model, inputs, output_component)

    if inputs.gamma_scalar is None:
        return sp.S.Zero

    i, j, a, b = output_component
    ns = model.total_real_scalar_dimension
    gamma_s = inputs.gamma_scalar
    result = sp.S.Zero

    for a_perm, b_perm in _swap_pair(a, b):
        for d in range(1, ns + 1):
            result += (
                gamma_s(b_perm, d)
                * coefficient[i, j, a_perm, d]
            )

    return sp.simplify(result)


def fermion_anomalous_dimension_term(
    model: RGEModel,
    inputs: MasterRGEInputs,
    output_component: tuple[int, int, int, int],
    coefficient=C,
) -> sp.Expr:
    r"""Return sum_sigma({i,j}) gamma_{c,f}^{jk} C_ikab from Eq. (4.85).

    The ``c`` in gamma_{c,f} labels the *collinear* anomalous dimension; it is
    not a summed fermion index.
    """

    _validate_master_dimensions(model, inputs, output_component)

    if inputs.gamma_fermion is None:
        return sp.S.Zero

    i, j, a, b = output_component
    nf = inputs.fermion_dimension
    gamma_f = inputs.gamma_fermion
    result = sp.S.Zero

    for i_perm, j_perm in _swap_pair(i, j):
        for k in range(1, nf + 1):
            result += (
                gamma_f(j_perm, k)
                * coefficient[i_perm, k, a, b]
            )

    return sp.simplify(result)


def calculate_master_rge(
    model: RGEModel,
    inputs: MasterRGEInputs,
    output_component: tuple[int, int, int, int],
    coefficient=C,
    simplify_each: bool = True,
) -> dict[str, sp.Expr]:
    """Evaluate every term displayed in Eq. (4.85).

    The anomalous-dimension terms are included whenever ``gamma_scalar`` and
    ``gamma_fermion`` are supplied.  Their explicit formulae should be taken
    from Eqs. (A.2) and (A.3), rather than guessed here.
    """

    functions = (
        ("yukawa_wavefunction", yukawa_wavefunction_term),
        ("scalar_pair", scalar_pair_term),
        ("mixed_yukawa_gauge", mixed_yukawa_gauge_term),
        ("crossed_yukawa", crossed_yukawa_term),
        ("conjugate_coefficient_yukawa", conjugate_coefficient_yukawa_term),
        ("scalar_anomalous_dimension", scalar_anomalous_dimension_term),
        ("fermion_anomalous_dimension", fermion_anomalous_dimension_term),
    )

    contributions: dict[str, sp.Expr] = {}

    for name, function in functions:
        value = function(
            model=model,
            inputs=inputs,
            output_component=output_component,
            coefficient=coefficient,
        )
        contributions[name] = sp.simplify(value) if simplify_each else value

    total = sum(contributions.values(), sp.S.Zero)
    contributions["total"] = sp.simplify(total)

    return contributions


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




# -----------------------------------------------------------------------------
# Collinear anomalous dimensions, Eqs. (A.2) and (A.3)
# -----------------------------------------------------------------------------


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


def scalar_collinear_anomalous_dimension(
    model: RGEModel,
    inputs: MasterRGEInputs,
    a: int,
    b: int,
) -> sp.Expr:
    r"""Return gamma_{c,s}^{ab} from Eq. (A.2).

    Implements

      gamma_{c,s}^{ab}
        = -4 sum_alpha g_alpha^2 [C2(S_alpha)]_{ab}
          + 1/2 sum_{i,j}
            (y_{jia} y^*_{jib} + y_{jib} y^*_{jia}).
    """

    ns = model.total_real_scalar_dimension
    nf = inputs.fermion_dimension

    if not 1 <= a <= ns or not 1 <= b <= ns:
        raise IndexError(f"Scalar indices must lie in 1,...,{ns}.")

    result = sp.S.Zero

    for sector in inputs.gauge_sectors:
        casimir = quadratic_casimir_matrix(sector.scalar_generators)
        result += -4 * sector.coupling**2 * casimir[a - 1, b - 1]

    y = inputs.yukawa

    yukawa_piece = sp.S.Zero

    for i in range(1, nf + 1):
        for j in range(1, nf + 1):
            yukawa_piece += (
                y(j, i, a) * sp.conjugate(y(j, i, b))
                + y(j, i, b) * sp.conjugate(y(j, i, a))
            )

    result += HALF * yukawa_piece

    return sp.simplify(result)


def fermion_collinear_anomalous_dimension(
    model: RGEModel,
    inputs: MasterRGEInputs,
    i: int,
    j: int,
) -> sp.Expr:
    r"""Return gamma_{c,f}^{ij} from Eq. (A.3).

    The source equation is implemented exactly as printed:

      gamma_{c,f}^{ij}
        = -3 sum_alpha [C2(F_alpha)]_{ij}
          + 1/2 sum_{k,a} y_{ika} y^*_{jka}.

    In particular, unlike Eq. (A.2), the displayed Eq. (A.3) does not show an
    explicit factor g_alpha^2 multiplying C2(F_alpha).
    """

    nf = inputs.fermion_dimension
    ns = model.total_real_scalar_dimension

    if not 1 <= i <= nf or not 1 <= j <= nf:
        raise IndexError(f"Fermion indices must lie in 1,...,{nf}.")

    result = sp.S.Zero

    for sector in inputs.gauge_sectors:
        casimir = quadratic_casimir_matrix(sector.fermion_generators)
        result += (
            -3
            * sector.coupling**2
            * casimir[i - 1, j - 1]
        )

    y = inputs.yukawa

    yukawa_piece = sp.S.Zero

    for k in range(1, nf + 1):
        for a in range(1, ns + 1):
            yukawa_piece += (
                y(i, k, a)
                * sp.conjugate(y(j, k, a))
            )

    result += HALF * yukawa_piece

    return sp.simplify(result)


def with_collinear_anomalous_dimensions(
    model: RGEModel,
    inputs: MasterRGEInputs,
) -> MasterRGEInputs:
    """Return a copy of the master-RGE inputs with A.2 and A.3 wired in."""

    @lru_cache(maxsize=None)
    def gamma_scalar(a: int, b: int) -> sp.Expr:
        return scalar_collinear_anomalous_dimension(
            model=model,
            inputs=inputs,
            a=a,
            b=b,
        )

    @lru_cache(maxsize=None)
    def gamma_fermion(i: int, j: int) -> sp.Expr:
        return fermion_collinear_anomalous_dimension(
            model=model,
            inputs=inputs,
            i=i,
            j=j,
        )

    return MasterRGEInputs(
        fermion_dimension=inputs.fermion_dimension,
        yukawa=inputs.yukawa,
        quartic=inputs.quartic,
        gauge_sectors=inputs.gauge_sectors,
        gamma_scalar=gamma_scalar,
        gamma_fermion=gamma_fermion,
    )


def calculate_complete_master_rge(
    model: RGEModel,
    inputs: MasterRGEInputs,
    output_component: tuple[int, int, int, int],
    coefficient=C,
    simplify_each: bool = True,
) -> dict[str, sp.Expr]:
    """Evaluate Eq. (4.85) using the collinear anomalous dimensions A.2 and A.3."""

    complete_inputs = with_collinear_anomalous_dimensions(
        model=model,
        inputs=inputs,
    )

    return calculate_master_rge(
        model=model,
        inputs=complete_inputs,
        output_component=output_component,
        coefficient=coefficient,
        simplify_each=simplify_each,
    )



# -----------------------------------------------------------------------------
# Neutral components and C -> K conversion
# -----------------------------------------------------------------------------


def su2_weights(dimension: int) -> tuple[sp.Rational, ...]:
    """Return T3 weights ordered consistently with su2_complex_generators."""

    j = sp.Rational(dimension - 1, 2)
    return tuple(j - position for position in range(dimension))


def neutral_complex_component_position(scalar: ComplexScalar) -> int | None:
    """Return the 1-based complex-component position with Q=T3+Y=0."""

    target_t3 = -scalar.hypercharge
    weights = su2_weights(scalar.su2_dimension)

    try:
        return weights.index(target_t3) + 1
    except ValueError:
        return None


def neutral_real_pair(model: RGEModel, scalar_name: str) -> tuple[int, int] | None:
    """Return the global (R,I) indices of a scalar's neutral complex component."""

    block = model.block(scalar_name)
    complex_position = neutral_complex_component_position(block.scalar)

    if complex_position is None:
        return None

    local_r = 2 * complex_position - 1
    local_i = 2 * complex_position

    return (
        block.local_to_global(local_r),
        block.local_to_global(local_i),
    )


def neutral_k_symbols(model: RGEModel) -> dict[str, sp.Symbol]:
    """Create one symbolic K coefficient for every scalar with a neutral component."""

    result: dict[str, sp.Symbol] = {}

    for scalar in model.scalars:
        if neutral_real_pair(model, scalar.name) is not None:
            result[scalar.name] = sp.Symbol(f"K_{scalar.name}")

    return result


def c_to_k_substitutions(
    model: RGEModel,
    fermion_indices: tuple[int, int],
    coefficient=C,
    k_symbols: Mapping[str, sp.Symbol] | None = None,
) -> dict[sp.Expr, sp.Expr]:
    """Return neutral-neutral C -> K substitutions for every scalar multiplet."""

    i, j = fermion_indices
    symbols = dict(k_symbols or neutral_k_symbols(model))
    substitutions: dict[sp.Expr, sp.Expr] = {}

    n = model.total_real_scalar_dimension

    for a in range(1, n + 1):
        for b in range(1, n + 1):
            substitutions[coefficient[i, j, a, b]] = sp.S.Zero

    for scalar_name, k_symbol in symbols.items():
        pair = neutral_real_pair(model, scalar_name)

        if pair is None:
            continue

        r, im = pair

        substitutions[coefficient[i, j, r, r]] = k_symbol / 2
        substitutions[coefficient[i, j, r, im]] = sp.I * k_symbol / 2
        substitutions[coefficient[i, j, im, r]] = sp.I * k_symbol / 2
        substitutions[coefficient[i, j, im, im]] = -k_symbol / 2

    return substitutions


def convert_c_expression_to_k(
    expression: sp.Expr,
    model: RGEModel,
    fermion_indices: tuple[int, int],
    coefficient=C,
    k_symbols: Mapping[str, sp.Symbol] | None = None,
) -> sp.Expr:
    """Convert supported neutral C components into model-generated K symbols."""

    substitutions = c_to_k_substitutions(
        model=model,
        fermion_indices=fermion_indices,
        coefficient=coefficient,
        k_symbols=k_symbols,
    )

    return sp.simplify(sp.expand(expression.xreplace(substitutions)))


# -----------------------------------------------------------------------------
# RGE aggregation
# -----------------------------------------------------------------------------


RGEContribution = Callable[..., sp.Expr]


def calculate_partial_rge(
    model: RGEModel,
    output_component: tuple[int, int, int, int],
    quartic_component: QuarticComponent = symbolic_quartic_component,
    coefficient=C,
) -> dict[str, sp.Expr]:
    """Evaluate the generalized contributions implemented so far.

    This is deliberately called *partial* because the full master RGE has not yet
    been implemented in this new file.
    """

    contributions = {
        "first_gauge_component": first_gauge_component(
            model=model,
            output_component=output_component,
            coefficient=coefficient,
        ),
        "first_scalar_component": first_scalar_component(
            model=model,
            output_component=output_component,
            quartic_component=quartic_component,
            coefficient=coefficient,
        ),
    }

    contributions["total"] = sp.simplify(
        sum(contributions.values(), sp.S.Zero)
    )

    return contributions


# -----------------------------------------------------------------------------
# Small self-checks
# -----------------------------------------------------------------------------


def _self_check() -> None:
    """Check dimensions and neutral-component bookkeeping on representative T3 models."""

    model = RGEModel.t3(
        d_s1=3,
        y_s1=0,
        d_s2=3,
        y_s2=1,
    )

    assert model.total_real_scalar_dimension == 16
    assert model.block("H").indices == range(1, 5)
    assert model.block("S1").indices == range(5, 11)
    assert model.block("S2").indices == range(11, 17)

    assert neutral_real_pair(model, "H") is not None
    assert neutral_real_pair(model, "S1") is not None
    assert neutral_real_pair(model, "S2") is not None


if __name__ == "__main__":
    _self_check()
    print("GeneralWeinbergRGEGenerator self-check passed.")
