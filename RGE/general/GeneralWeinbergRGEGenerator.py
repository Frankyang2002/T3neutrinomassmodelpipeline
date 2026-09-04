from __future__ import annotations

from RGE.general.RGEModel import ComplexScalar, RGEModel, ScalarBasisBlock

from RGE.general.GaugeGenerators import (
    GaugeSector,
    embed_local_generator,
    g1,
    g2,
    gauge_sectors_from_generators,
    global_su2_generators,
    global_u1_generator,
    quadratic_casimir_matrix,
    real_scalar_generator,
    su2_complex_generators,
    su2_real_scalar_generators,
    u1_real_scalar_generator,
)

from RGE.general.RGECommon import (
    C,
    HALF,
    QuarticComponent,
    _validate_output_component,
)

from RGE.general.MasterWeinbergRGE import (
    FermionAnomalousDimension,
    MasterRGEInputs,
    ScalarAnomalousDimension,
    YukawaComponent,
    _pair_product_permutations,
    _simultaneous_pair_permutations,
    _swap_pair,
    _validate_master_dimensions,
    calculate_master_rge,
    conjugate_coefficient_yukawa_term,
    crossed_yukawa_term,
    fermion_anomalous_dimension_term,
    mixed_yukawa_gauge_term,
    scalar_anomalous_dimension_term,
    scalar_pair_term,
    yukawa_wavefunction_term,
)

from RGE.general.AnomalousDimensions import (
    calculate_complete_master_rge,
    fermion_collinear_anomalous_dimension,
    scalar_collinear_anomalous_dimension,
    symbolic_fermion_anomalous_dimension,
    symbolic_scalar_anomalous_dimension,
    with_collinear_anomalous_dimensions,
)

from RGE.general.PartialWeinbergRGE import (
    LAMBDA,
    RGEContribution,
    _gauge_generator_contribution,
    calculate_partial_rge,
    first_gauge_component,
    first_scalar_component,
    mapping_quartic_component,
    symbolic_quartic_component,
)







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





# -----------------------------------------------------------------------------
# Model and scalar-basis definitions
# -----------------------------------------------------------------------------







from RGE.general.WeinbergKappaConversion import (
    c_to_k_substitutions,
    convert_c_expression_to_k,
    neutral_complex_component_position,
    neutral_k_symbols,
    neutral_real_pair,
    su2_weights,
)



# -----------------------------------------------------------------------------
# Gauge generators
# -----------------------------------------------------------------------------










# -----------------------------------------------------------------------------
# Global scalar generators
# -----------------------------------------------------------------------------








# -----------------------------------------------------------------------------
# Generic quartic tensor lambda_abcd
# -----------------------------------------------------------------------------








# -----------------------------------------------------------------------------
# Generic RGE contributions currently inherited from the legacy prototype
# -----------------------------------------------------------------------------











# -----------------------------------------------------------------------------
# Full psi^2 phi^2 master RGE, Eq. (4.85)
# -----------------------------------------------------------------------------








def symbolic_yukawa_component(
    i: int,
    j: int,
    a: int,
    yukawa_tensor=sp.IndexedBase("y"),
) -> sp.Expr:
    """Return a formal Yukawa component y_ija."""

    return yukawa_tensor[i, j, a]


































# -----------------------------------------------------------------------------
# Collinear anomalous dimensions, Eqs. (A.2) and (A.3)
# -----------------------------------------------------------------------------













# -----------------------------------------------------------------------------
# Neutral components and C -> K conversion
# -----------------------------------------------------------------------------














# -----------------------------------------------------------------------------
# RGE aggregation
# -----------------------------------------------------------------------------






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
