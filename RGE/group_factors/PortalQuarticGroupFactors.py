from __future__ import annotations

"""Analytic backbone group factors for beta_lambdaH1, beta_lambdaH2, beta_lambda12.

This module deliberately contains only structures whose A--E pattern has
already been reduced to generic SU(2)/hypercharge factors.  Pure-gauge g^4
terms, adjoint-quartic recouplings, and T3-E Cross mixing are kept for the
next recoupling stage rather than fitted model-by-model.

Conventions
-----------
Q = T3 + Y
Y(H)  = 1/2
Y(S1) = alpha/2
Y(S2) = (alpha+2)/2

C2(d) = (d^2-1)/4
G_Si  = max(dF,dSi)/dSi

The resolved coefficients multiply the following monomials in 16 pi^2 beta.
"""

from dataclasses import asdict, dataclass
from fractions import Fraction
import argparse
import json

from RGE.group_factors.RepresentationFactors import (
    su2_quadratic_casimir_from_dimension,
)


def _txt(x: Fraction) -> str:
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def _validate_t3(dS1: int, dS2: int, dF: int) -> None:
    if min(dS1, dS2, dF) < 1:
        raise ValueError("Representation dimensions must be positive.")
    if abs(dS1 - dF) != 1 or abs(dS2 - dF) != 1:
        raise ValueError("T3 Yukawa invariance requires dSi=dF+/-1.")


@dataclass(frozen=True)
class PortalQuarticBackbone:
    dS1: int
    dS2: int
    dF: int
    alpha: int

    GS1: str
    GS2: str
    C2S1: str
    C2S2: str
    YS1: str
    YS2: str

    # beta_lambdaH1
    H1_y1_trace: str
    H1_Yd_trace: str
    H1_Ye_trace: str
    H1_Yu_trace: str
    H1_y1Ye_mixed: str
    H1_lambdaT3_sq: str
    H1_g2_linear: str
    H1_gY_linear: str
    H1_lambdaH_lambdaH1: str
    H1_lambdaH1_sq: str
    H1_lambdaH1_lambdaS1: str
    H1_lambda12_lambdaH2: str

    # beta_lambdaH2
    H2_y2_trace: str
    H2_Yd_trace: str
    H2_Ye_trace: str
    H2_Yu_trace: str
    H2_y2Ye_mixed: str
    H2_lambdaT3_sq: str
    H2_g2_linear: str
    H2_gY_linear: str
    H2_lambdaH_lambdaH2: str
    H2_lambdaH2_sq: str
    H2_lambdaH2_lambdaS2: str
    H2_lambda12_lambdaH1: str

    # beta_lambda12
    L12_y1_trace: str
    L12_y2_trace: str
    L12_y1y2_mixed: str
    L12_lambdaT3_sq: str
    L12_g2_linear: str
    L12_gY_linear: str
    L12_lambda12_sq: str
    L12_lambda12_lambdaS1: str
    L12_lambda12_lambdaS2: str
    L12_lambdaH1_lambdaH2: str


def portal_quartic_backbone(
    dS1: int,
    dS2: int,
    dF: int,
    alpha: int,
) -> PortalQuarticBackbone:
    dS1, dS2, dF, alpha = map(int, (dS1, dS2, dF, alpha))
    _validate_t3(dS1, dS2, dF)

    GS1 = Fraction(max(dF, dS1), dS1)
    GS2 = Fraction(max(dF, dS2), dS2)
    C1 = su2_quadratic_casimir_from_dimension(dS1)
    C2 = su2_quadratic_casimir_from_dimension(dS2)
    CH = Fraction(3, 4)
    YH = Fraction(1, 2)
    Y1 = Fraction(alpha, 2)
    Y2 = Fraction(alpha + 2, 2)

    return PortalQuarticBackbone(
        dS1=dS1,
        dS2=dS2,
        dF=dF,
        alpha=alpha,
        GS1=_txt(GS1),
        GS2=_txt(GS2),
        C2S1=_txt(C1),
        C2S2=_txt(C2),
        YS1=_txt(Y1),
        YS2=_txt(Y2),

        H1_y1_trace=_txt(2 * GS1),
        H1_Yd_trace="6",
        H1_Ye_trace="2",
        H1_Yu_trace="6",
        H1_y1Ye_mixed=_txt(-2 * GS1),
        H1_lambdaT3_sq=_txt(Fraction(12, dS1)),
        H1_g2_linear=_txt(-6 * (CH + C1)),
        H1_gY_linear=_txt(-6 * (YH * YH + Y1 * Y1)),
        H1_lambdaH_lambdaH1="6",
        H1_lambdaH1_sq="4",
        H1_lambdaH1_lambdaS1=_txt(Fraction(2 * (dS1 + 1))),
        H1_lambda12_lambdaH2=_txt(Fraction(2 * dS2)),

        H2_y2_trace=_txt(2 * GS2),
        H2_Yd_trace="6",
        H2_Ye_trace="2",
        H2_Yu_trace="6",
        H2_y2Ye_mixed=_txt(-2 * GS2),
        H2_lambdaT3_sq=_txt(Fraction(12, dS2)),
        H2_g2_linear=_txt(-6 * (CH + C2)),
        H2_gY_linear=_txt(-6 * (YH * YH + Y2 * Y2)),
        H2_lambdaH_lambdaH2="6",
        H2_lambdaH2_sq="4",
        H2_lambdaH2_lambdaS2=_txt(Fraction(2 * (dS2 + 1))),
        H2_lambda12_lambdaH1=_txt(Fraction(2 * dS1)),

        L12_y1_trace=_txt(2 * GS1),
        L12_y2_trace=_txt(2 * GS2),
        L12_y1y2_mixed=_txt(-2 * GS1 * GS2),
        L12_lambdaT3_sq=_txt(Fraction(12, dS1 * dS2)),
        L12_g2_linear=_txt(-6 * (C1 + C2)),
        L12_gY_linear=_txt(-6 * (Y1 * Y1 + Y2 * Y2)),
        L12_lambda12_sq="4",
        L12_lambda12_lambdaS1=_txt(Fraction(2 * (dS1 + 1))),
        L12_lambda12_lambdaS2=_txt(Fraction(2 * (dS2 + 1))),
        L12_lambdaH1_lambdaH2="4",
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("dS1", type=int)
    parser.add_argument("dS2", type=int)
    parser.add_argument("dF", type=int)
    parser.add_argument("alpha", type=int)
    args = parser.parse_args()

    result = portal_quartic_backbone(args.dS1, args.dS2, args.dF, args.alpha)
    print(json.dumps(asdict(result), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
