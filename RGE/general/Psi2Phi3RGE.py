"""One-loop self-running of dimension-six ``psi^2 phi^3`` operators.

This module implements the ``psi^2 phi^3 -> psi^2 phi^3`` anomalous-dimension
sector used for scalar-first T3 studies.  The coefficient convention is

    C6[i,j,a,b,c] psi_i psi_j phi_a phi_b phi_c,

with separate symmetry in the two fermion indices and the three real-scalar
indices.  The tensor contractions follow Eq. (4.25) of Aebischer, Bresciani
and Selimovic, *Anomalous Dimension of a General Effective Gauge Theory II:
Fermionic Sector* (arXiv:2512.16890).

Scope
-----
This is the explicit ``psi^2 phi^3 <- psi^2 phi^3`` kernel in Eq. (4.25).
The source explicitly states that the ``psi^2 phi^3`` RGEs must also be
supplemented by the EOM/redundant-operator expressions in Eq. (2.36), built
from the ``D^2 phi^4`` sector.  Those contributions are not implemented here.
The T3 scalar-first production backend must therefore remain blocked until
the exact Matchete boundary and every required contribution for that boundary
have been verified.
"""

from __future__ import annotations

from itertools import permutations

import sympy as sp

from RGE.general.ScalarBasis import ScalarBasis
from RGE.general.WilsonTensorRGE import HALF, WilsonRGEInputs


C6 = sp.IndexedBase("C6")


def _scalar_permutations(a: int, b: int, c: int):
    """Return the six terms in the scalar symmetrisation sum.

    Repeated external component labels deliberately produce repeated tuples.
    The source formula sums permutations, not unique numerical tuples, and the
    explicit ``1/2!`` prefactor accounts for the symmetry of the remaining
    scalar pair.
    """

    return tuple(permutations((a, b, c), 3))


def _fermion_permutations(i: int, j: int):
    """Return the two terms in the fermion symmetrisation sum."""

    return ((i, j), (j, i))


def _mixed_permutations(a: int, b: int, c: int, i: int, j: int):
    """Return the product of scalar and fermion external permutations."""

    return tuple(
        (*scalar_perm, *fermion_perm)
        for scalar_perm in _scalar_permutations(a, b, c)
        for fermion_perm in _fermion_permutations(i, j)
    )


def _validate_dimensions(
    model: ScalarBasis,
    inputs: WilsonRGEInputs,
    output_component: tuple[int, int, int, int, int],
) -> None:
    """Validate one ``(i,j,a,b,c)`` component and all gauge generators."""

    if len(output_component) != 5:
        raise ValueError(
            "output_component must have the form (i, j, a, b, c)."
        )

    i, j, a, b, c = output_component
    nf = inputs.fermion_dimension
    ns = model.total_real_scalar_dimension

    for name, index in (("i", i), ("j", j)):
        if not 1 <= index <= nf:
            raise IndexError(
                f"Fermion index {name}={index} is outside 1,...,{nf}."
            )

    for name, index in (("a", a), ("b", b), ("c", c)):
        if not 1 <= index <= ns:
            raise IndexError(
                f"Scalar index {name}={index} is outside 1,...,{ns}."
            )

    for sector in inputs.gauge_sectors:
        for theta in sector.scalar_generators:
            if theta.shape != (ns, ns):
                raise ValueError(
                    "A scalar gauge generator does not match the model's "
                    f"{ns}-dimensional real-scalar basis."
                )
        for generator in sector.fermion_generators:
            if generator.shape != (nf, nf):
                raise ValueError(
                    "A fermion gauge generator does not match "
                    "fermion_dimension."
                )


def yukawa_wavefunction_term(
    model: ScalarBasis,
    inputs: WilsonRGEInputs,
    output_component: tuple[int, int, int, int, int],
    coefficient=C6,
) -> sp.Expr:
    """Return ``y_ijd y^*_{kld} C6_klabc``."""

    _validate_dimensions(model, inputs, output_component)
    i, j, a, b, c = output_component
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
                    * coefficient[k, l, a, b, c]
                )

    return sp.simplify(result)


def scalar_pair_term(
    model: ScalarBasis,
    inputs: WilsonRGEInputs,
    output_component: tuple[int, int, int, int, int],
    coefficient=C6,
) -> sp.Expr:
    r"""Return the scalar-pair gauge/quartic contraction.

    Implements

    ``-1/2! sum_sigma(a,b,c) [2 g^2 theta_bd theta_ce - lambda_bcde]
    C6_ijade``

    with an implicit sum over gauge sectors/generators and internal ``d,e``.
    """

    _validate_dimensions(model, inputs, output_component)
    i, j, a, b, c = output_component
    ns = model.total_real_scalar_dimension
    quartic = inputs.quartic
    result = sp.S.Zero

    for a_perm, b_perm, c_perm in _scalar_permutations(a, b, c):
        for d in range(1, ns + 1):
            for e in range(1, ns + 1):
                gauge_piece = sp.S.Zero

                for sector in inputs.gauge_sectors:
                    for theta in sector.scalar_generators:
                        gauge_piece += (
                            sector.coupling**2
                            * theta[b_perm - 1, d - 1]
                            * theta[c_perm - 1, e - 1]
                        )

                result += -HALF * (
                    2 * gauge_piece
                    - quartic(b_perm, c_perm, d, e)
                ) * coefficient[i, j, a_perm, d, e]

    return sp.simplify(result)


def mixed_yukawa_gauge_term(
    model: ScalarBasis,
    inputs: WilsonRGEInputs,
    output_component: tuple[int, int, int, int, int],
    coefficient=C6,
) -> sp.Expr:
    r"""Return the mixed Yukawa/gauge contraction of Eq. (4.25)."""

    _validate_dimensions(model, inputs, output_component)
    i, j, a, b, c = output_component
    nf = inputs.fermion_dimension
    ns = model.total_real_scalar_dimension
    y = inputs.yukawa
    result = sp.S.Zero

    for (
        a_perm,
        b_perm,
        c_perm,
        i_perm,
        j_perm,
    ) in _mixed_permutations(a, b, c, i, j):
        for k in range(1, nf + 1):
            for l in range(1, nf + 1):
                for d in range(1, ns + 1):
                    bracket = (
                        2
                        * y(j_perm, k, d)
                        * sp.conjugate(y(k, l, c_perm))
                        + y(j_perm, k, c_perm)
                        * sp.conjugate(y(k, l, d))
                    )

                    gauge_piece = sp.S.Zero
                    for sector in inputs.gauge_sectors:
                        for t_generator, theta in zip(
                            sector.fermion_generators,
                            sector.scalar_generators,
                            strict=True,
                        ):
                            gauge_piece += (
                                4
                                * sector.coupling**2
                                * t_generator[l - 1, j_perm - 1]
                                * theta[c_perm - 1, d - 1]
                            )

                    result += HALF * (
                        bracket + gauge_piece
                    ) * coefficient[
                        i_perm,
                        l,
                        a_perm,
                        b_perm,
                        d,
                    ]

    return sp.simplify(result)


def conjugate_coefficient_yukawa_term(
    model: ScalarBasis,
    inputs: WilsonRGEInputs,
    output_component: tuple[int, int, int, int, int],
    coefficient=C6,
) -> sp.Expr:
    """Return ``(2 y_ild y_jkd + y_ijd y_kld) C6^*_klabc``."""

    _validate_dimensions(model, inputs, output_component)
    i, j, a, b, c = output_component
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
                ) * sp.conjugate(coefficient[k, l, a, b, c])

    return sp.simplify(result)


def scalar_anomalous_dimension_term(
    model: ScalarBasis,
    inputs: WilsonRGEInputs,
    output_component: tuple[int, int, int, int, int],
    coefficient=C6,
) -> sp.Expr:
    """Return the scalar collinear anomalous-dimension contribution."""

    _validate_dimensions(model, inputs, output_component)

    if inputs.gamma_scalar is None:
        return sp.S.Zero

    i, j, a, b, c = output_component
    ns = model.total_real_scalar_dimension
    gamma_s = inputs.gamma_scalar
    result = sp.S.Zero

    for a_perm, b_perm, c_perm in _scalar_permutations(a, b, c):
        for d in range(1, ns + 1):
            result += (
                HALF
                * gamma_s(c_perm, d)
                * coefficient[i, j, a_perm, b_perm, d]
            )

    return sp.simplify(result)


def fermion_anomalous_dimension_term(
    model: ScalarBasis,
    inputs: WilsonRGEInputs,
    output_component: tuple[int, int, int, int, int],
    coefficient=C6,
) -> sp.Expr:
    """Return the fermion collinear anomalous-dimension contribution."""

    _validate_dimensions(model, inputs, output_component)

    if inputs.gamma_fermion is None:
        return sp.S.Zero

    i, j, a, b, c = output_component
    nf = inputs.fermion_dimension
    gamma_f = inputs.gamma_fermion
    result = sp.S.Zero

    for i_perm, j_perm in _fermion_permutations(i, j):
        for k in range(1, nf + 1):
            result += (
                gamma_f(j_perm, k)
                * coefficient[i_perm, k, a, b, c]
            )

    return sp.simplify(result)


def calculate_psi2phi3_rge(
    model: ScalarBasis,
    inputs: WilsonRGEInputs,
    output_component: tuple[int, int, int, int, int],
    coefficient=C6,
    simplify_each: bool = True,
) -> dict[str, sp.Expr]:
    """Evaluate the implemented one-loop ``psi^2 phi^3`` self-running terms."""

    functions = (
        ("yukawa_wavefunction", yukawa_wavefunction_term),
        ("scalar_pair", scalar_pair_term),
        ("mixed_yukawa_gauge", mixed_yukawa_gauge_term),
        (
            "conjugate_coefficient_yukawa",
            conjugate_coefficient_yukawa_term,
        ),
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

    contributions["total"] = sp.simplify(
        sum(contributions.values(), sp.S.Zero)
    )
    return contributions


def calculate_complete_psi2phi3_rge(
    model: ScalarBasis,
    inputs: WilsonRGEInputs,
    output_component: tuple[int, int, int, int, int],
    coefficient=C6,
    simplify_each: bool = True,
) -> dict[str, sp.Expr]:
    """Evaluate the self-running kernel with collinear anomalous dimensions.

    ``complete`` here means Eq. (4.25) with its collinear anomalous
    dimensions included.  It does not add the Eq. (2.36) EOM contribution or
    mixing from other dimension-six operator classes.
    """

    from RGE.general.AnomalousDimensions import (
        with_collinear_anomalous_dimensions,
    )

    complete_inputs = with_collinear_anomalous_dimensions(
        model=model,
        inputs=inputs,
    )

    return calculate_psi2phi3_rge(
        model=model,
        inputs=complete_inputs,
        output_component=output_component,
        coefficient=coefficient,
        simplify_each=simplify_each,
    )
