from __future__ import annotations

import sys
from pathlib import Path

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

    print("RGBeta Python integration smoke test: PASS")


if __name__ == "__main__":
    main()
