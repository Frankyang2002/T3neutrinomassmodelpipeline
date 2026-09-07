from __future__ import annotations

from tests.non_pipeline_rge.general.FermionBasis import (
    FermionBasis,
    FermionBasisBlock,
    WeylFermion,
    _fermion_local_su2_generators,
    build_gauge_sectors,
    fermion_global_su2_generators,
    fermion_global_u1_generator,
)

from tests.non_pipeline_rge.general.T3YukawaTensors import (
    ComplexYukawaComponent,
    build_real_yukawa_tensor,
    yukawa_component_function,
)

from tests.non_pipeline_rge.general.T3QuarticTensors import (
    ComplexQuarticComponent,
    ComplexScalarFactor,
    complex_scalar_component_expression,
    quartic_component_function,
    quartic_components_from_exchange,
    quartic_polynomial_from_components,
    quartic_tensor_from_components,
    scalar_real_symbols,
)




"""
Adapter between generalized T3 model data and the tensors used by
GeneralWeinbergRGEGenerator_complete.py.

This file deliberately does not parse pretty-printed Matchete Lagrangians.
Instead, the Wolfram/model-building side should export machine-readable
component tensors for the accepted UV interactions.  This module then:

- constructs the global real-scalar basis;
- constructs a global left-handed Weyl-fermion basis;
- embeds scalar and fermion gauge generators;
- converts complex-scalar Yukawa components into y_ija;
- expands complex quartic interactions into the fully symmetric real tensor
  lambda_abcd;
- stores/makes available the matched Wilson tensor C_ijab.

Important normalization convention
----------------------------------
The quartic conversion assumes the real-scalar convention used by Eq. (4.85),

    V_4 = (1/4!) lambda_abcd phi_a phi_b phi_c phi_d,

so lambda_abcd is obtained by four derivatives of the quartic polynomial.

For Yukawas, this adapter only converts the complex scalar into its real
(R,I) components.  The supplied complex Yukawa coefficient must already use
the same fermion ordering and overall normalization as y_ija in Eq. (4.85).
"""

from dataclasses import dataclass
import json
from pathlib import Path
import re
from typing import Iterable, Mapping, Sequence

import sympy as sp
from sympy.parsing.mathematica import parse_mathematica

from RGE.general.GaugeGenerators import (
    GaugeSector,
    g1,
    g2,
    global_su2_generators,
    global_u1_generator,
    su2_complex_generators,
)
from RGE.general.MasterWeinbergRGE import MasterRGEInputs
from RGE.general.RGECommon import C
from RGE.general.RGEModel import (
    ComplexScalar,
    RGEModel,
)


SQRT2 = sp.sqrt(2)


# Mathematica sometimes writes exact fourth roots as algebraic ``Root``
# objects instead of radicals.  For x^4 = a/b its canonical ordering is
#
#   1: -r,  2: r,  3: -I r,  4: I r,    r = (a/b)^(1/4).
#
# SymPy's Mathematica parser currently leaves these as an unevaluated
# ``Root(Lambda(...), k, 0)`` function.  In particular, it cannot then prove
# that roots 3 and 4 are complex conjugates.  The T3 tensor contractions rely
# on precisely those identities, so canonicalize this exporter-generated
# subset before parsing the rest of the expression.
_PURE_QUARTIC_ROOT_RE = re.compile(
    r"Root\[\s*-\s*(?P<numerator>\d+)\s*\+\s*"
    r"(?P<denominator>\d+)\s*\*\s*#1\^4\s*&\s*,\s*"
    r"(?P<index>[1-4])\s*,\s*0\s*\]"
)


def _canonicalize_wolfram_pure_quartic_roots(value: str) -> str:
    """Rewrite the exact ``Root[-a+b #1^4&,k,0]`` subset as radicals."""

    def replacement(match: re.Match[str]) -> str:
        numerator = match.group("numerator")
        denominator = match.group("denominator")
        root = f"(({numerator})/({denominator}))^(1/4)"
        return {
            "1": f"-({root})",
            "2": root,
            "3": f"-I*({root})",
            "4": f"I*({root})",
        }[match.group("index")]

    return _PURE_QUARTIC_ROOT_RE.sub(replacement, value)


# -----------------------------------------------------------------------------
# Exact-expression parsing
# -----------------------------------------------------------------------------


def parse_exact_expression(value) -> sp.Expr:
    """Convert a JSON number/string into a SymPy expression.

    The exporter should prefer simple exact strings such as:
      "1/2", "-sqrt(3)/2", "I/sqrt(2)", "y1[1,2]".
    """

    if isinstance(value, (int, float)):
        return sp.sympify(value)

    if not isinstance(value, str):
        raise TypeError(f"Expected an exact-expression string, got {type(value)!r}.")

    # The exchange uses Mathematica InputForm.  SymPy's Mathematica parser
    # correctly preserves exact rationals, nested Sqrt[...] expressions and
    # Conjugate[...] instead of treating them as unrelated symbol names.
    try:
        canonical_value = _canonicalize_wolfram_pure_quartic_roots(value)
        return sp.sympify(parse_mathematica(canonical_value))
    except Exception as exc:
        raise ValueError(f"Could not parse exact Wolfram expression {value!r}.") from exc


# -----------------------------------------------------------------------------
# Fermion basis
# -----------------------------------------------------------------------------
















# -----------------------------------------------------------------------------
# Scalar real-basis variables
# -----------------------------------------------------------------------------






# -----------------------------------------------------------------------------
# Yukawa tensor y_ija
# -----------------------------------------------------------------------------








# -----------------------------------------------------------------------------
# Quartic tensor lambda_abcd
# -----------------------------------------------------------------------------














# -----------------------------------------------------------------------------
# Wilson tensor C_ijab
# -----------------------------------------------------------------------------


class SparseWilsonLookup:
    def __init__(
        self,
        components: Mapping[tuple[int, int, int, int], sp.Expr],
    ):
        self.components = components

    def __getitem__(
        self,
        key: tuple[int, int, int, int],
    ) -> sp.Expr:
        return sp.sympify(
            self.components.get(key, sp.S.Zero)
        )

    def __call__(
        self,
        i: int,
        j: int,
        a: int,
        b: int,
    ) -> sp.Expr:
        return self[i, j, a, b]


def wilson_component_function(
    components: Mapping[tuple[int, int, int, int], sp.Expr],
):
    """Return a sparse Wilson-coefficient lookup callable."""

    return SparseWilsonLookup(components)


# -----------------------------------------------------------------------------
# Final assembly
# -----------------------------------------------------------------------------


@dataclass
class T3RGETensors:
    """All tensors/bases needed to call the generalized Eq. (4.85) engine."""

    scalar_model: RGEModel
    fermion_basis: FermionBasis
    yukawa_components: dict[tuple[int, int, int], sp.Expr]
    quartic_components: dict[tuple[int, int, int, int], sp.Expr]
    wilson_components: dict[tuple[int, int, int, int], sp.Expr]

    def master_inputs(self) -> MasterRGEInputs:
        """Build MasterRGEInputs from the converted UV tensors."""

        return MasterRGEInputs(
            fermion_dimension=self.fermion_basis.dimension,
            yukawa=yukawa_component_function(self.yukawa_components),
            quartic=quartic_component_function(self.quartic_components),
            gauge_sectors=build_gauge_sectors(
                self.scalar_model,
                self.fermion_basis,
            ),
        )


# -----------------------------------------------------------------------------
# Machine-readable exchange format
# -----------------------------------------------------------------------------


def load_component_exchange(path: str | Path) -> dict:
    """Load the JSON component exchange produced by the Wolfram adapter."""

    return json.loads(Path(path).read_text(encoding="utf-8"))


def scalar_model_from_exchange(data: Mapping) -> RGEModel:
    """Construct the scalar RGE model from exchange JSON."""

    scalars = tuple(
        ComplexScalar(
            name=item["name"],
            su2_dimension=int(item["su2_dimension"]),
            hypercharge=parse_exact_expression(item["hypercharge"]),
        )
        for item in data["scalars"]
    )

    return RGEModel(scalars=scalars)


def fermion_basis_from_exchange(data: Mapping) -> FermionBasis:
    """Construct the global Weyl basis from exchange JSON."""

    fermions = tuple(
        WeylFermion(
            name=item["name"],
            su2_dimension=int(item["su2_dimension"]),
            hypercharge=parse_exact_expression(item["hypercharge"]),
            multiplicity=int(item.get("multiplicity", 1)),
            conjugated_representation=bool(
                item.get("conjugated_representation", False)
            ),
        )
        for item in data["fermions"]
    )

    return FermionBasis(fermions=fermions)


# -----------------------------------------------------------------------------
# Regression/self-checks
# -----------------------------------------------------------------------------


def _self_check_quartic_normalization() -> None:
    """Recover the legacy doublet self-quartic normalization.

    For V=(lambda/2)(Phi^dagger Phi)^2, the real tensor must satisfy
      lambda_aaaa = 3 lambda,
      lambda_aabb = lambda   (a != b).
    """

    lam = sp.Symbol("lam", real=True)
    model = RGEModel(
        scalars=(
            ComplexScalar(
                name="Phi",
                su2_dimension=2,
                hypercharge=sp.Rational(1, 2),
            ),
        )
    )

    terms: list[ComplexQuarticComponent] = []

    # (lambda/2) sum_mn Phi_m^* Phi_m Phi_n^* Phi_n
    for m in range(1, 3):
        for n in range(1, 3):
            terms.append(
                ComplexQuarticComponent(
                    coefficient=lam / 2,
                    factors=(
                        ComplexScalarFactor("Phi", m, True),
                        ComplexScalarFactor("Phi", m, False),
                        ComplexScalarFactor("Phi", n, True),
                        ComplexScalarFactor("Phi", n, False),
                    ),
                )
            )

    tensor = quartic_tensor_from_components(model, terms)

    assert sp.simplify(tensor[(1, 1, 1, 1)] - 3 * lam) == 0
    assert sp.simplify(tensor[(1, 1, 2, 2)] - lam) == 0
    assert sp.simplify(tensor[(1, 1, 3, 3)] - lam) == 0


def _self_check_yukawa_conversion() -> None:
    """Check z=(R+iI)/sqrt(2) conversion used by y_ija."""

    y0 = sp.Symbol("y0")
    scalar_model = RGEModel(
        scalars=(
            ComplexScalar("Phi", 1, 0),
        )
    )
    fermion_basis = FermionBasis(
        fermions=(
            WeylFermion("f", 1, 0, multiplicity=2),
        )
    )

    tensor = build_real_yukawa_tensor(
        scalar_model,
        fermion_basis,
        (
            ComplexYukawaComponent(
                fermion_i=1,
                fermion_j=2,
                scalar_name="Phi",
                scalar_component=1,
                scalar_conjugated=False,
                coefficient=y0,
            ),
        ),
    )

    assert sp.simplify(tensor[(1, 2, 1)] - y0 / SQRT2) == 0
    assert sp.simplify(tensor[(1, 2, 2)] - sp.I * y0 / SQRT2) == 0


def _self_check() -> None:
    _self_check_quartic_normalization()
    _self_check_yukawa_conversion()


if __name__ == "__main__":
    _self_check()
    print("T3RGETensors self-check passed.")
