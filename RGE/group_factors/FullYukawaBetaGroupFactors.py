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

Then

16 pi^2 beta_y1 =
    G_S1 Tr(y1 y1^dagger) y1
  + 1/2 (G_L1 + G_F1) y1 y1^dagger y1
  + 1/2 G_L2 y2 y2^dagger y1
  + 1/2 Ye Ye^dagger y1
  - 3 g2^2 [C2(L) + C2(F)] y1
  - 3 gY^2 [Y(L)^2 + Y(F)^2] y1,

and beta_y2 is obtained by 1 <-> 2 in the Yukawa pieces, with the same
gauge coefficients because both vertices contain the same two fermion
multiplets L and F.

The SU(2) and U(1) gauge terms are the standard one-loop fermion-Casimir
contribution to a Yukawa beta function.  The scalar representation does not
enter these gauge coefficients.
"""

import argparse
import json
from dataclasses import asdict, dataclass
from pathlib import Path
import sys

import sympy as sp

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from RGE.group_factors.YukawaBetaGroupFactors import _leg_factors


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
    y1: YukawaBetaCoefficients
    y2: YukawaBetaCoefficients


def _c2_from_dimension(d: int) -> sp.Expr:
    d = int(d)
    j = sp.Rational(d - 1, 2)
    return sp.simplify(j * (j + 1))


def complete_yukawa_group_factors(
    *,
    dS1: int,
    dS2: int,
    dF: int,
    alpha: int,
) -> T3YukawaBetaGroupFactors:
    dS1, dS2, dF, alpha = map(int, (dS1, dS2, dF, alpha))

    y1legs = _leg_factors(dF, dS1)
    y2legs = _leg_factors(dF, dS2)

    GL1 = sp.Rational(max(dF, dS1), 2)
    GF1 = sp.Rational(max(dF, dS1), dF)
    GS1 = sp.Rational(max(dF, dS1), dS1)

    GL2 = sp.Rational(max(dF, dS2), 2)
    GF2 = sp.Rational(max(dF, dS2), dF)
    GS2 = sp.Rational(max(dF, dS2), dS2)

    C2L = sp.Rational(3, 4)
    C2F = _c2_from_dimension(dF)

    YL = -sp.Rational(1, 2)
    YF = sp.Rational(alpha + 1, 2)

    su2 = sp.simplify(-3 * (C2L + C2F))
    u1 = sp.simplify(-3 * (YL**2 + YF**2))

    y1 = YukawaBetaCoefficients(
        trace=str(GS1),
        self_matrix=str(sp.simplify((GL1 + GF1) / 2)),
        cross_matrix=str(sp.simplify(GL2 / 2)),
        charged_lepton_matrix="1/2",
        su2_gauge=str(su2),
        u1_gauge=str(u1),
    )

    y2 = YukawaBetaCoefficients(
        trace=str(GS2),
        self_matrix=str(sp.simplify((GL2 + GF2) / 2)),
        cross_matrix=str(sp.simplify(GL1 / 2)),
        charged_lepton_matrix="1/2",
        su2_gauge=str(su2),
        u1_gauge=str(u1),
    )

    return T3YukawaBetaGroupFactors(
        dS1=dS1,
        dS2=dS2,
        dF=dF,
        alpha=alpha,
        C2L=str(C2L),
        C2F=str(C2F),
        YL=str(YL),
        YF=str(YF),
        y1=y1,
        y2=y2,
    )


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Complete generic T3 y1/y2 one-loop group-factor coefficients."
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
