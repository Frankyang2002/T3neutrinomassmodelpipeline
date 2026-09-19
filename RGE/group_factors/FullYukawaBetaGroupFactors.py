from __future__ import annotations

"""Complete generic one-loop group factors for the T3 Yukawa beta functions.

Conventions
-----------
Q = T3 + Y
Y(L) = -1/2
Y(F) = (alpha + 1)/2
d_Si = d_F +/- 1

For a Yukawa invariant involving L, F and S_i define
    d_i^> = max(d_F, d_Si)

and the canonically normalised SU(2) leg factors
    G_Si = d_i^> / d_Si
    G_Li = d_i^> / 2
    G_Fi = d_i^> / d_F.

The vector-like branch is

16 pi^2 beta_y1 =
    G_S1 Tr(y1 y1^dagger) y1
  + 1/2 (G_L1 + G_F1) y1 y1^dagger y1
  + 1/2 G_L2 y2 y2^dagger y1
  + 1/2 Ye Ye^dagger y1
  - 3 g2^2 [C2(L) + C2(F)] y1
  - 3 gY^2 [Y(L)^2 + Y(F)^2] y1,

with beta_y2 obtained by 1 <-> 2 in the Yukawa pieces.

For the physical self-conjugate branch,
    alpha == -1 and dF odd,
RGBeta exposes two additional nongauge beta_y2 structures:
    + G_S2 Tr(y2 y1^dagger) y1
    + 1/2 G_F2 y2 y1^dagger y1.

Even-dF alpha=-1 RGBeta points use a different neutral-field convention and
are not exact validations of the current physical Matchete T3 model.

The SU(2) and U(1) gauge terms are the standard one-loop fermion-Casimir
contribution.  The scalar representation does not enter the gauge
coefficients.
"""

import argparse
import json
from dataclasses import asdict, dataclass
from fractions import Fraction

import sympy as sp

from RGE.group_factors.RepresentationFactors import (
    canonical_yukawa_leg_factors,
    su2_quadratic_casimir_from_dimension,
)


@dataclass(frozen=True)
class YukawaBetaCoefficients:
    trace: str
    self_matrix: str
    cross_matrix: str
    charged_lepton_matrix: str
    su2_gauge: str
    u1_gauge: str


@dataclass(frozen=True)
class T3YukawaBetaGroupFactors:
    dS1: int
    dS2: int
    dF: int
    alpha: int
    C2L: str
    C2F: str
    YL: str
    YF: str
    self_conjugate_F: bool
    y1: YukawaBetaCoefficients
    y2: YukawaBetaCoefficients
    y2_self_conjugate_extra: dict[str, str]


def _txt(value: Fraction | sp.Rational) -> str:
    value = sp.Rational(value)
    return (
        str(int(value))
        if value.q == 1
        else f"{int(value.p)}/{int(value.q)}"
    )


def physical_self_conjugate_f(dF: int, alpha: int) -> bool:
    """Current Matchete-side physical self-conjugate-F criterion."""
    return int(alpha) == -1 and int(dF) % 2 == 1


def _validate(dS1: int, dS2: int, dF: int) -> None:
    if any(d not in (1, 2, 3) for d in (dS1, dS2, dF)):
        raise ValueError("Current T3 implementation supports d in {1,2,3}.")
    if abs(dS1 - dF) != 1 or abs(dS2 - dF) != 1:
        raise ValueError("T3 Yukawa invariance requires dSi=dF+/-1.")


def complete_yukawa_group_factors(
    *,
    dS1: int,
    dS2: int,
    dF: int,
    alpha: int,
) -> T3YukawaBetaGroupFactors:
    dS1, dS2, dF, alpha = map(int, (dS1, dS2, dF, alpha))
    _validate(dS1, dS2, dF)

    y1legs = canonical_yukawa_leg_factors(dF, dS1)
    y2legs = canonical_yukawa_leg_factors(dF, dS2)

    GL1 = sp.Rational(
        y1legs.G_lepton.numerator,
        y1legs.G_lepton.denominator,
    )
    GF1 = sp.Rational(
        y1legs.G_heavy.numerator,
        y1legs.G_heavy.denominator,
    )
    GS1 = sp.Rational(
        y1legs.G_scalar.numerator,
        y1legs.G_scalar.denominator,
    )

    GL2 = sp.Rational(
        y2legs.G_lepton.numerator,
        y2legs.G_lepton.denominator,
    )
    GF2 = sp.Rational(
        y2legs.G_heavy.numerator,
        y2legs.G_heavy.denominator,
    )
    GS2 = sp.Rational(
        y2legs.G_scalar.numerator,
        y2legs.G_scalar.denominator,
    )

    C2L = sp.Rational(3, 4)
    c2f = su2_quadratic_casimir_from_dimension(dF)
    C2F = sp.Rational(c2f.numerator, c2f.denominator)

    YL = -sp.Rational(1, 2)
    YF = sp.Rational(alpha + 1, 2)

    su2 = sp.simplify(-3 * (C2L + C2F))
    u1 = sp.simplify(-3 * (YL**2 + YF**2))

    y1 = YukawaBetaCoefficients(
        trace=_txt(GS1),
        self_matrix=_txt(sp.simplify((GL1 + GF1) / 2)),
        cross_matrix=_txt(sp.simplify(GL2 / 2)),
        charged_lepton_matrix="1/2",
        su2_gauge=_txt(su2),
        u1_gauge=_txt(u1),
    )

    y2 = YukawaBetaCoefficients(
        trace=_txt(GS2),
        self_matrix=_txt(sp.simplify((GL2 + GF2) / 2)),
        cross_matrix=_txt(sp.simplify(GL1 / 2)),
        charged_lepton_matrix="1/2",
        su2_gauge=_txt(su2),
        u1_gauge=_txt(u1),
    )

    self_conjugate = physical_self_conjugate_f(dF, alpha)
    y2_self_conjugate_extra: dict[str, str] = {}
    if self_conjugate:
        y2_self_conjugate_extra = {
            "Tr_y2_y1dag*y1": _txt(GS2),
            "y2_y1dag_y1": _txt(sp.simplify(GF2 / 2)),
        }

    return T3YukawaBetaGroupFactors(
        dS1=dS1,
        dS2=dS2,
        dF=dF,
        alpha=alpha,
        C2L=_txt(C2L),
        C2F=_txt(C2F),
        YL=_txt(YL),
        YF=_txt(YF),
        self_conjugate_F=self_conjugate,
        y1=y1,
        y2=y2,
        y2_self_conjugate_extra=y2_self_conjugate_extra,
    )


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Complete generic T3 y1/y2 one-loop group-factor coefficients, "
            "including the physical self-conjugate beta_y2 extras."
        )
    )
    parser.add_argument("dS1", type=int)
    parser.add_argument("dS2", type=int)
    parser.add_argument("dF", type=int)
    parser.add_argument("alpha", type=int)
    args = parser.parse_args()

    result = complete_yukawa_group_factors(
        dS1=args.dS1,
        dS2=args.dS2,
        dF=args.dF,
        alpha=args.alpha,
    )
    print(json.dumps(asdict(result), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
