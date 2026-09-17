from __future__ import annotations

"""Analytic SU(2) recoupling factor for direct C12 -> Weinberg mixing.

For the T3 operator chain

    (L L)(S1 S2)  --lambdaT3-->  (L L)(H H),

the one-loop direct EFT1 mixing can be written

    16*pi^2 beta_kappa = R_W * lambdaT3 * C12,

with

    R_W =
      (4/3) * eta_F
      * sqrt(3 dS1 dS2 dF)
      * Wigner6j(1/2, 1/2, 1; jS2, jS1, jF).

Here j=(d-1)/2.

The Wigner 6j contains the invariant SU(2) recoupling geometry.  eta_F is a
basis-phase convention inherited from the project's Yukawa/CG ordering.  In
the currently supported d<=3 T3 convention,

    eta_F = +1 for dF=1,2
            -1 for dF=3.

Only the overall sign is convention-dependent; the magnitude is fixed by the
recoupling invariant.  The implementation does not use A--E model labels.
"""

from dataclasses import asdict, dataclass
import argparse
import json

import sympy as sp
from sympy.physics.wigner import wigner_6j


def _validate(d1: int, d2: int, dF: int) -> None:
    if any(d not in (1, 2, 3) for d in (d1, d2, dF)):
        raise ValueError("Current T3 implementation supports d in {1,2,3}.")
    if abs(d1 - dF) != 1 or abs(d2 - dF) != 1:
        raise ValueError("T3 Yukawa invariance requires dSi=dF+/-1.")


def _j(d: int) -> sp.Rational:
    return sp.Rational(d - 1, 2)


def project_cg_phase(dF: int) -> int:
    """Project CG convention phase for the currently supported dF<=3 basis."""
    dF = int(dF)
    if dF in (1, 2):
        return 1
    if dF == 3:
        return -1
    raise ValueError("Unsupported dF for current project CG phase convention.")


@dataclass(frozen=True)
class DirectWeinbergGroupFactor:
    dS1: int
    dS2: int
    dF: int
    jS1: str
    jS2: str
    jF: str
    wigner6j: str
    normalization: str
    cg_phase: int
    reduced_factor: str


def direct_weinberg_group_factor(
    dS1: int,
    dS2: int,
    dF: int,
) -> DirectWeinbergGroupFactor:
    d1, d2, dF = map(int, (dS1, dS2, dF))
    _validate(d1, d2, dF)

    j1, j2, jf = _j(d1), _j(d2), _j(dF)
    sixj = sp.simplify(
        wigner_6j(
            sp.Rational(1, 2),
            sp.Rational(1, 2),
            1,
            j2,
            j1,
            jf,
        )
    )

    norm = sp.sqrt(3 * d1 * d2 * dF)
    phase = project_cg_phase(dF)
    reduced = sp.factor(sp.Rational(4, 3) * phase * norm * sixj)

    return DirectWeinbergGroupFactor(
        dS1=d1,
        dS2=d2,
        dF=dF,
        jS1=sp.sstr(j1),
        jS2=sp.sstr(j2),
        jF=sp.sstr(jf),
        wigner6j=sp.sstr(sixj),
        normalization=sp.sstr(norm),
        cg_phase=phase,
        reduced_factor=sp.sstr(reduced),
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("dS1", type=int)
    parser.add_argument("dS2", type=int)
    parser.add_argument("dF", type=int)
    args = parser.parse_args()

    result = direct_weinberg_group_factor(args.dS1, args.dS2, args.dF)
    print(json.dumps(asdict(result), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
