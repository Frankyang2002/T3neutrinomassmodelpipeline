from __future__ import annotations

"""Generic one-loop group factors for the T3 loop-closing quartic lambdaT3.

Interaction
-----------
    lambdaT3 * H H S1 S2† + h.c.

The two Higgs doublets are in the symmetric J=1 channel.

All coefficients in this module are analytic group-factor expressions.  The
non-singlet scalar recouplings are the gauge-covariant component-basis results

    R_H1Adj = C2(S2) - C2(S1) - 2
    R_H2Adj = 2 - C2(S1) + C2(S2)
    R_12Adj = C2(S1) + C2(S2) - 2

when the corresponding operator exists.

For the T3-E-only triplet-triplet Cross operator,

    R_12Cross = -2.

The Cross operator exists only for dS1=dS2=3 in the current T3 model basis.

These are physical component-basis recouplings.  They should not be replaced
by the raw RGBeta formal-index coefficients in the pseudoreal/real cases.
"""

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
import json


@dataclass(frozen=True)
class MixingQuarticGroupFactors:
    dS1: int
    dS2: int
    dF: int
    alpha: int

    GS1: str
    GS2: str

    C2H: str
    C2S1: str
    C2S2: str

    YH: str
    YS1: str
    YS2: str

    y1_trace: str
    y2_trace: str
    yd_trace: str
    ye_trace: str
    yu_trace: str

    lambdaH: str
    lambdaH1: str
    lambdaH2: str
    lambda12: str

    lambdaH1Adj: str
    lambdaH2Adj: str
    lambda12Adj: str
    lambda12Cross: str

    has_lambdaH1Adj: bool
    has_lambdaH2Adj: bool
    has_lambda12Adj: bool
    has_lambda12Cross: bool

    su2_gauge: str
    u1_gauge: str


def _txt(value: Fraction) -> str:
    return (
        str(value.numerator)
        if value.denominator == 1
        else f"{value.numerator}/{value.denominator}"
    )


def _c2(dimension: int) -> Fraction:
    """Quadratic Casimir for the SU(2) irrep of dimension d."""
    return Fraction(dimension * dimension - 1, 4)


def _validate_t3_dimensions(dS1: int, dS2: int, dF: int) -> None:
    if min(dS1, dS2, dF) < 1:
        raise ValueError("SU(2) representation dimensions must be positive.")

    if abs(dS1 - dF) != 1 or abs(dS2 - dF) != 1:
        raise ValueError("T3 Yukawa invariance requires dSi=dF+/-1.")


def mixing_quartic_group_factors(
    dS1: int,
    dS2: int,
    dF: int,
    alpha: int,
) -> MixingQuarticGroupFactors:
    dS1, dS2, dF, alpha = map(int, (dS1, dS2, dF, alpha))
    _validate_t3_dimensions(dS1, dS2, dF)

    GS1 = Fraction(max(dF, dS1), dS1)
    GS2 = Fraction(max(dF, dS2), dS2)

    C2H = Fraction(3, 4)
    C2S1 = _c2(dS1)
    C2S2 = _c2(dS2)

    YH = Fraction(1, 2)
    YS1 = Fraction(alpha, 2)
    YS2 = Fraction(alpha + 2, 2)

    # Gauge wave-function contribution:
    #   -3 g^2 sum_external C2
    # and analogously for U(1) charges.
    su2 = -3 * (2 * C2H + C2S1 + C2S2)
    u1 = -3 * (2 * YH * YH + YS1 * YS1 + YS2 * YS2)

    has_h1_adj = dS1 > 1
    has_h2_adj = dS2 > 1
    has_12_adj = dS1 > 1 and dS2 > 1
    has_cross = dS1 == 3 and dS2 == 3

    # Gauge-covariant scalar recouplings in the physical component basis.
    r_h1 = C2S2 - C2S1 - 2 if has_h1_adj else Fraction(0)
    r_h2 = 2 - C2S1 + C2S2 if has_h2_adj else Fraction(0)
    r_12 = C2S1 + C2S2 - 2 if has_12_adj else Fraction(0)
    r_cross = Fraction(-2) if has_cross else Fraction(0)

    return MixingQuarticGroupFactors(
        dS1=dS1,
        dS2=dS2,
        dF=dF,
        alpha=alpha,
        GS1=_txt(GS1),
        GS2=_txt(GS2),
        C2H=_txt(C2H),
        C2S1=_txt(C2S1),
        C2S2=_txt(C2S2),
        YH=_txt(YH),
        YS1=_txt(YS1),
        YS2=_txt(YS2),
        y1_trace=_txt(GS1),
        y2_trace=_txt(GS2),
        yd_trace="6",
        ye_trace="2",
        yu_trace="6",
        lambdaH="2",
        lambdaH1="4",
        lambdaH2="4",
        lambda12="2",
        lambdaH1Adj=_txt(r_h1),
        lambdaH2Adj=_txt(r_h2),
        lambda12Adj=_txt(r_12),
        lambda12Cross=_txt(r_cross),
        has_lambdaH1Adj=has_h1_adj,
        has_lambdaH2Adj=has_h2_adj,
        has_lambda12Adj=has_12_adj,
        has_lambda12Cross=has_cross,
        su2_gauge=_txt(su2),
        u1_gauge=_txt(u1),
    )


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Analytic group factors in the one-loop beta_lambdaT3."
    )
    parser.add_argument("dS1", type=int)
    parser.add_argument("dS2", type=int)
    parser.add_argument("dF", type=int)
    parser.add_argument("alpha", type=int)
    args = parser.parse_args()

    result = mixing_quartic_group_factors(
        args.dS1, args.dS2, args.dF, args.alpha
    )
    print(json.dumps(asdict(result), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
