from __future__ import annotations

import sys
from pathlib import Path

import sympy as sp

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from RGE.running.RGBetaT3Running import run_rgbeta_t3


MODELS = [
    ("A", 1, 3, 2, 0),
    ("B", 2, 2, 1, -1),
    ("C", 2, 2, 3, -1),
    ("D", 3, 1, 2, -2),
    ("E", 3, 3, 2, 0),
]


def main() -> None:
    required = {
        "gY",
        "g2",
        "g3",
        "yu",
        "yd",
        "ye",
        "y1",
        "y2",
        "MF",
        "mS1Sq",
        "mS2Sq",
        "lambdaH",
        "lambdaS1",
        "lambdaS2",
        "lambdaH1",
        "lambdaH2",
        "lambda12",
        "lambdaT3",
    }

    for label, d_s1, d_s2, d_f, alpha in MODELS:
        result = run_rgbeta_t3(d_s1, d_s2, d_f, alpha)

        missing = required.difference(result.betas)

        print(
            f"T3-{label}: status={result.status}, "
            f"betas={len(result.betas)}, missing={sorted(missing)}"
        )

        if missing:
            raise AssertionError(f"T3-{label} is missing beta functions: {missing}")

    shared = run_rgbeta_t3(2, 2, 1, -1, shared_scalar=True)
    shared_required = {
        "gY", "g2", "g3", "yu", "yd", "ye", "h", "MF",
        "mSSq", "lambdaH", "lambdaS", "lambda3", "lambda4", "lambda5",
    }
    shared_missing = shared_required.difference(shared.betas)
    print(
        f"shared dS=2,dF=1: status={shared.status}, "
        f"betas={len(shared.betas)}, missing={sorted(shared_missing)}"
    )
    if shared_missing:
        raise AssertionError(
            f"Shared-scalar scotogenic model is missing beta functions: {shared_missing}"
        )

    def parse_report_beta(text: str):
        return sp.sympify(
            str(text).replace("^", "**"),
            locals={"gY": sp.Symbol("gY"), "g2": sp.Symbol("g2"), "g3": sp.Symbol("g3")},
        )

    gY, g2, g3 = sp.symbols("gY g2 g3")
    assert sp.simplify(parse_report_beta(shared.report_betas["gY"]) - 7 * gY**3) == 0
    assert sp.simplify(parse_report_beta(shared.report_betas["g2"]) + 3 * g2**3) == 0
    assert sp.simplify(parse_report_beta(shared.report_betas["g3"]) + 7 * g3**3) == 0

    print("RGBeta Python integration smoke test: PASS")


if __name__ == "__main__":
    main()
