from __future__ import annotations

"""Analytic one-loop group factors for beta_y2 in the T3 models.

Base/vector-like branch:
    16*pi^2 beta_y2 =
        G_S2 Tr(y2 y2^\dagger) y2
      + 1/2 (G_L2 + G_F2) y2 y2^\dagger y2
      + 1/2 G_L1 y1 y1^\dagger y2
      + 1/2 Ye Ye^\dagger y2
      - 3 g2^2 [C2(L)+C2(F)] y2
      - 3 gY^2 [Y_L^2+Y_F^2] y2.

For the physical self-conjugate branch (Y_F=0 and odd dF), RGBeta exposes
two additional nongauge structures:
        + G_S2 Tr(y2 y1^\dagger) y1
        + 1/2 G_F2 y2 y1^\dagger y1.

The current Matchete-side physical branch criterion is therefore
    alpha == -1 and dF odd.

Even-dF alpha=-1 RGBeta points use a different neutral-field convention and
are not exact validations of the current physical T3 model.
"""

from dataclasses import asdict, dataclass
from fractions import Fraction
import argparse
import json


def _txt(x: Fraction) -> str:
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def _c2(d: int) -> Fraction:
    return Fraction(d * d - 1, 4)


def _leg_factors(dS: int, dF: int):
    dgt = max(dS, dF)
    return {
        "GS": Fraction(dgt, dS),
        "GL": Fraction(dgt, 2),
        "GF": Fraction(dgt, dF),
    }


def physical_self_conjugate_f(dF: int, alpha: int) -> bool:
    return int(alpha) == -1 and int(dF) % 2 == 1


def _validate(d1: int, d2: int, dF: int) -> None:
    if any(d not in (1, 2, 3) for d in (d1, d2, dF)):
        raise ValueError("Current T3 implementation supports d in {1,2,3}.")
    if abs(d1 - dF) != 1 or abs(d2 - dF) != 1:
        raise ValueError("T3 Yukawa invariance requires dSi=dF+/-1.")


@dataclass(frozen=True)
class Yukawa2BetaGroupFactors:
    dS1: int
    dS2: int
    dF: int
    alpha: int
    self_conjugate_F: bool
    base_nongauge: dict[str, str]
    self_conjugate_extra: dict[str, str]
    gauge: dict[str, str]


def yukawa2_beta_group_factors(
    dS1: int,
    dS2: int,
    dF: int,
    alpha: int,
) -> Yukawa2BetaGroupFactors:
    d1, d2, dF, alpha = map(int, (dS1, dS2, dF, alpha))
    _validate(d1, d2, dF)

    f1 = _leg_factors(d1, dF)
    f2 = _leg_factors(d2, dF)

    YL = Fraction(-1, 2)
    YF = Fraction(alpha + 1, 2)
    C2L = Fraction(3, 4)
    C2F = _c2(dF)

    base = {
        "Tr_y2*y2": _txt(f2["GS"]),
        "y2_y2dag_y2": _txt(Fraction(1, 2) * (f2["GL"] + f2["GF"])),
        "y1_y1dag_y2": _txt(Fraction(1, 2) * f1["GL"]),
        "Ye_Yedag_y2": "1/2",
    }

    self_conj = {}
    if physical_self_conjugate_f(dF, alpha):
        self_conj = {
            "Tr_y2_y1dag*y1": _txt(f2["GS"]),
            "y2_y1dag_y1": _txt(Fraction(1, 2) * f2["GF"]),
        }

    gauge = {
        "g2_sq*y2": _txt(-3 * (C2L + C2F)),
        "gY_sq*y2": _txt(-3 * (YL * YL + YF * YF)),
    }

    return Yukawa2BetaGroupFactors(
        dS1=d1,
        dS2=d2,
        dF=dF,
        alpha=alpha,
        self_conjugate_F=physical_self_conjugate_f(dF, alpha),
        base_nongauge=base,
        self_conjugate_extra=self_conj,
        gauge=gauge,
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("dS1", type=int)
    parser.add_argument("dS2", type=int)
    parser.add_argument("dF", type=int)
    parser.add_argument("alpha", type=int)
    args = parser.parse_args()

    result = yukawa2_beta_group_factors(
        args.dS1, args.dS2, args.dF, args.alpha
    )
    print(json.dumps(asdict(result), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
