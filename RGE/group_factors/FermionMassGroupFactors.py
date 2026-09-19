from __future__ import annotations

"""Generic one-loop group factors for the T3 heavy-fermion mass beta function.

For the two T3 Yukawa vertices define

    d_i^> = max(d_F, d_Si),
    G_Fi  = d_i^> / d_F.

The heavy-fermion leg contraction therefore gives one-half G_Fi in the
fermion anomalous dimension.  With the project flavor conventions,

16 pi^2 beta_MF =
    1/2 G_F1 (y1^T y1^*) MF
  + 1/2 G_F2 MF (y2^dagger y2)
  - 6 g2^2 C2(F) MF
  - 6 gY^2 YF^2 MF.

Here
    C2(F) = (d_F^2 - 1)/4,
    YF    = (alpha + 1)/2.
"""

import argparse
from dataclasses import asdict, dataclass
import json

from RGE.group_factors.RepresentationFactors import (
    su2_quadratic_casimir_from_dimension,
)


@dataclass(frozen=True)
class FermionMassGroupFactors:
    dS1: int
    dS2: int
    dF: int
    alpha: int
    GF1: str
    GF2: str
    y1_left: str
    y2_right: str
    C2F: str
    YF: str
    su2_gauge: str
    u1_gauge: str


def _rat(num: int, den: int = 1) -> str:
    from fractions import Fraction
    value = Fraction(num, den)
    return (
        str(value.numerator)
        if value.denominator == 1
        else f"{value.numerator}/{value.denominator}"
    )


def fermion_mass_group_factors(
    dS1: int,
    dS2: int,
    dF: int,
    alpha: int,
) -> FermionMassGroupFactors:
    from fractions import Fraction

    dS1, dS2, dF, alpha = map(int, (dS1, dS2, dF, alpha))

    if abs(dS1 - dF) != 1 or abs(dS2 - dF) != 1:
        raise ValueError("T3 Yukawa invariance requires dSi=dF+/-1.")

    GF1 = Fraction(max(dF, dS1), dF)
    GF2 = Fraction(max(dF, dS2), dF)

    C2F = su2_quadratic_casimir_from_dimension(dF)
    YF = Fraction(alpha + 1, 2)

    su2 = -6 * C2F
    u1 = -6 * YF * YF

    def txt(x: Fraction) -> str:
        return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"

    return FermionMassGroupFactors(
        dS1=dS1,
        dS2=dS2,
        dF=dF,
        alpha=alpha,
        GF1=txt(GF1),
        GF2=txt(GF2),
        y1_left=txt(GF1 / 2),
        y2_right=txt(GF2 / 2),
        C2F=txt(C2F),
        YF=txt(YF),
        su2_gauge=txt(su2),
        u1_gauge=txt(u1),
    )


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Generic T3 heavy-fermion mass beta-function group factors."
    )
    parser.add_argument("dS1", type=int)
    parser.add_argument("dS2", type=int)
    parser.add_argument("dF", type=int)
    parser.add_argument("alpha", type=int)
    args = parser.parse_args()

    result = fermion_mass_group_factors(
        args.dS1, args.dS2, args.dF, args.alpha
    )
    print(json.dumps(asdict(result), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
