from __future__ import annotations

from functools import lru_cache
import sympy as sp

from RGE.general.GaugeGenerators import quadratic_casimir_matrix
from RGE.general.MasterWeinbergRGE import (
    C,
    HALF,
    MasterRGEInputs,
    calculate_master_rge,
)
from RGE.general.RGEModel import RGEModel

# Collinear anomalous dimensions and the complete one-loop Weinberg RGE.

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
