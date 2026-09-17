from __future__ import annotations

"""Physical one-loop group factors for the T3 heavy-scalar mass beta functions.

The supported mass parameters are mS1Sq=m_{S1}^2 and mS2Sq=m_{S2}^2.

There are two fermion-mass insertion branches:

1. Vector-like/non-self-conjugate F:
       coefficient = -4 G_Si

2. Physical self-conjugate F:
       coefficient = -16 G_Si

In the current T3 field definition, the self-conjugate branch is selected only
when Y_F=0 and the SU(2) representation is real, i.e. dF is odd.  Since
Y_F=(alpha+1)/2, this means

       self_conjugate_F <=> alpha == -1 and dF odd.

This distinction is necessary because the RGBeta model currently switches to a
neutral/self-conjugate construction whenever Y_F=0, including even dF=2.  That
RGBeta alpha=-1, dF=2 point is therefore not the same field-content branch as
the current Matchete-side T3 model.

Conventions:
    Q = T3 + Y
    Y(S1) = alpha/2
    Y(S2) = (alpha+2)/2
    C2(d) = (d^2-1)/4
    G_Si = max(dF,dSi)/dSi
"""

from dataclasses import asdict, dataclass
from fractions import Fraction
import argparse
import json


def _txt(x: Fraction) -> str:
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def _c2(d: int) -> Fraction:
    return Fraction(d * d - 1, 4)


def _validate(d1: int, d2: int, dF: int) -> None:
    if min(d1, d2, dF) < 1:
        raise ValueError("Representation dimensions must be positive.")
    if d1 not in (1, 2, 3) or d2 not in (1, 2, 3) or dF not in (1, 2, 3):
        raise ValueError("Current T3 implementation supports d in {1,2,3}.")
    if abs(d1 - dF) != 1 or abs(d2 - dF) != 1:
        raise ValueError("T3 Yukawa invariance requires dSi=dF+/-1.")


def physical_self_conjugate_f(dF: int, alpha: int) -> bool:
    """Current Matchete-side physical branch criterion."""
    return int(alpha) == -1 and int(dF) % 2 == 1


@dataclass(frozen=True)
class ScalarMassGroupFactors:
    dS1: int
    dS2: int
    dF: int
    alpha: int
    self_conjugate_F: bool
    GS1: str
    GS2: str
    C2S1: str
    C2S2: str
    YS1: str
    YS2: str
    beta_mS1Sq: dict[str, str]
    beta_mS2Sq: dict[str, str]


def scalar_mass_group_factors(
    dS1: int,
    dS2: int,
    dF: int,
    alpha: int,
    *,
    self_conjugate_F: bool | None = None,
) -> ScalarMassGroupFactors:
    d1, d2, dF, alpha = map(int, (dS1, dS2, dF, alpha))
    _validate(d1, d2, dF)

    physical_sc = physical_self_conjugate_f(dF, alpha)
    if self_conjugate_F is None:
        self_conjugate_F = physical_sc
    else:
        self_conjugate_F = bool(self_conjugate_F)

    if self_conjugate_F and alpha != -1:
        raise ValueError("A self-conjugate F requires Y_F=0, hence alpha=-1.")
    if self_conjugate_F and dF % 2 == 0:
        raise ValueError(
            "Current physical T3 branch does not treat even-dimensional "
            "SU(2) F as self-conjugate."
        )

    C1 = _c2(d1)
    C2 = _c2(d2)
    Y1 = Fraction(alpha, 2)
    Y2 = Fraction(alpha + 2, 2)
    GS1 = Fraction(max(dF, d1), d1)
    GS2 = Fraction(max(dF, d2), d2)

    heavy_factor = Fraction(-16 if self_conjugate_F else -4)

    m1 = {
        "mS1Sq*Tr_y1": _txt(2 * GS1),
        "Tr_MF_y1": _txt(heavy_factor * GS1),
        "lambdaS1*mS1Sq": _txt(Fraction(2 * (d1 + 1))),
        "g2_sq*mS1Sq": _txt(-6 * C1),
        "gY_sq*mS1Sq": _txt(-6 * Y1 * Y1),
        "lambda12*mS2Sq": _txt(Fraction(2 * d2)),
    }

    if d1 == 3:
        m1["lambdaS1Adj*mS1Sq"] = _txt(2 * C1)

    if d1 == d2 == 3:
        m1["lambda12Cross*mS2Sq"] = "4"

    m2 = {
        "mS2Sq*Tr_y2": _txt(2 * GS2),
        "Tr_MF_y2": _txt(heavy_factor * GS2),
        "lambdaS2*mS2Sq": _txt(Fraction(2 * (d2 + 1))),
        "g2_sq*mS2Sq": _txt(-6 * C2),
        "gY_sq*mS2Sq": _txt(-6 * Y2 * Y2),
        "lambda12*mS1Sq": _txt(Fraction(2 * d1)),
    }

    if d2 == 3:
        m2["lambdaS2Adj*mS2Sq"] = _txt(2 * C2)

    if d1 == d2 == 3:
        m2["lambda12Cross*mS1Sq"] = "4"

    return ScalarMassGroupFactors(
        dS1=d1,
        dS2=d2,
        dF=dF,
        alpha=alpha,
        self_conjugate_F=self_conjugate_F,
        GS1=_txt(GS1),
        GS2=_txt(GS2),
        C2S1=_txt(C1),
        C2S2=_txt(C2),
        YS1=_txt(Y1),
        YS2=_txt(Y2),
        beta_mS1Sq=m1,
        beta_mS2Sq=m2,
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("dS1", type=int)
    parser.add_argument("dS2", type=int)
    parser.add_argument("dF", type=int)
    parser.add_argument("alpha", type=int)
    args = parser.parse_args()

    result = scalar_mass_group_factors(
        args.dS1, args.dS2, args.dF, args.alpha
    )
    print(json.dumps(asdict(result), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
