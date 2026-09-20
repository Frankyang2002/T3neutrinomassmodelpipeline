# Does 1 loop RGE for our wilson coefficient

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
import sympy as sp

from RGE.general.GaugeGenerators import GaugeSector
from RGE.general.ScalarBasis import ScalarBasis

# Shared symbolic definitions and validation helpers for the general RGE.
HALF = sp.Rational(1, 2)
C = sp.IndexedBase("C")
# Quartic(a,b,c,d), yukawa(i,j,a)
QuarticComponent = Callable[[int, int, int, int], sp.Expr]
YukawaComponent = Callable[[int, int, int], sp.Expr]
# Scalar(a,b), fermion(i,j)
ScalarAnomalousDimension = Callable[[int, int], sp.Expr]
FermionAnomalousDimension = Callable[[int, int], sp.Expr]



@dataclass(frozen=True)
class WilsonRGEInputs:
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

# Some helpers
def _swap_pair(first, second):
    """Get all permutation of a pair (itself and 1 swap)."""

    return ((first, second), (second, first))


def _pair_product_permutations(a, b, i, j):
    """Return the four independent swaps in sigma({a,b} x {i,j})."""

    return tuple(
        (a_perm, b_perm, i_perm, j_perm)
        for a_perm, b_perm in _swap_pair(a, b)
        for i_perm, j_perm in _swap_pair(i, j)
    )


def _simultaneous_pair_permutations(a, b, i, j):
    """Pair exchange, not all combinations
    Switching scalar indices, and also switching fermion indices once"""

    return (
        (a, b, i, j),
        (b, a, j, i),
    )




def _validate_rge_dimensions(
    model: ScalarBasis,
    inputs: WilsonRGEInputs,
    output_component: tuple[int, int, int, int],
) -> None:
    """Make sure everything is valid
    We first validate scalar indices and inputs
    Then we validate fermion indices and inputs
    Then we check if the generators are in the right dimensions for our field."""

    _, _, a, b = output_component
    n = model.total_real_scalar_dimension

    if len(output_component) != 4:
        raise ValueError("output_component must have the form (i, j, a, b).")
    
    if not 1 <= a <= n:
        raise IndexError(f"Scalar index a={a} is outside 1,...,{n}.")
    if not 1 <= b <= n:
        raise IndexError(f"Scalar index b={b} is outside 1,...,{n}.")

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
    model: ScalarBasis,
    inputs: WilsonRGEInputs,
    output_component: tuple[int, int, int, int],
    coefficient=C,
) -> sp.Expr:
    r"""Return y_ijd y^*_{kld} C_klab
    We sum over k,l,d"""

    _validate_rge_dimensions(model, inputs, output_component)
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
    model: ScalarBasis,
    inputs: WilsonRGEInputs,
    output_component: tuple[int, int, int, int],
    coefficient=C,
) -> sp.Expr:
    r"""Return the -1/2! sigma({a,b})sigma(cd)[2sigma(\alpha,A)
    g^2_alpha \theta^A_ac \theta^A_bd
    -\lambda_abcd] C_ijcd 
    term
    """

    _validate_rge_dimensions(model, inputs, output_component)
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
    model: ScalarBasis,
    inputs: WilsonRGEInputs,
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

    _validate_rge_dimensions(model, inputs, output_component)
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
    model: ScalarBasis,
    inputs: WilsonRGEInputs,
    output_component: tuple[int, int, int, int],
    coefficient=C,
) -> sp.Expr:
    r"""Return sum_sigma y_jla y_ikc C_lkbc from Eq. (4.85)."""

    _validate_rge_dimensions(model, inputs, output_component)
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
    model: ScalarBasis,
    inputs: WilsonRGEInputs,
    output_component: tuple[int, int, int, int],
    coefficient=C,
) -> sp.Expr:
    r"""Return (2 y_ild y_jkd + y_ijd y_kld) C^*_{klab} from Eq. (4.85)."""

    _validate_rge_dimensions(model, inputs, output_component)
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
    model: ScalarBasis,
    inputs: WilsonRGEInputs,
    output_component: tuple[int, int, int, int],
    coefficient=C,
) -> sp.Expr:
    r"""Return sum_sigma({a,b}) gamma_{c,s}^{bd} C_ijad from Eq. (4.85).

    The ``c`` in gamma_{c,s} labels the *collinear* anomalous dimension; it is
    not a summed scalar index.
    """

    _validate_rge_dimensions(model, inputs, output_component)

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
    model: ScalarBasis,
    inputs: WilsonRGEInputs,
    output_component: tuple[int, int, int, int],
    coefficient=C,
) -> sp.Expr:
    r"""Return sum_sigma({i,j}) gamma_{c,f}^{jk} C_ikab from Eq. (4.85).

    The ``c`` in gamma_{c,f} labels the *collinear* anomalous dimension; it is
    not a summed fermion index.
    """

    _validate_rge_dimensions(model, inputs, output_component)

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


def calculate_wilson_tensor_rge(
    model: ScalarBasis,
    inputs: WilsonRGEInputs,
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
