from __future__ import annotations

import argparse
from pathlib import Path

import sympy as sp

from RGE.matching.FlavorMatchedC5 import (
    build_flavor_matched_c5,
    extract_t3_loop_kernel,
)
from RGE.matching.MatchedEFTRGE import parse_matchete_c5


def default_c5_path() -> Path:
    return (
        Path(__file__).resolve().parent
        / "wolfram"
        / "output"
        / "T3_B_alpha_m1"
        / "c5_coefficient.txt"
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--c5",
        type=Path,
        default=default_c5_path(),
    )
    args = parser.parse_args()

    kappa_1g = parse_matchete_c5(
        args.c5.read_text(encoding="utf-8")
    )
    kernel = extract_t3_loop_kernel(kappa_1g)

    y1 = sp.Symbol("y1")
    y2 = sp.Symbol("y2")

    Y1 = sp.Matrix([[y1]])
    Y2 = sp.Matrix([[y2]])

    K1 = build_flavor_matched_c5(
        kernel,
        Y1,
        Y2,
        heavy_masses=None,
    )

    one_generation_difference = sp.simplify(
        K1[0, 0] - kappa_1g
    )

    a11, a12, a21, a22 = sp.symbols(
        "a11 a12 a21 a22"
    )
    b11, b12, b21, b22 = sp.symbols(
        "b11 b12 b21 b22"
    )

    Y1_2 = sp.Matrix([
        [a11, a12],
        [a21, a22],
    ])
    Y2_2 = sp.Matrix([
        [b11, b12],
        [b21, b22],
    ])

    K2 = build_flavor_matched_c5(
        kernel,
        Y1_2,
        Y2_2,
        heavy_masses=[
            sp.Symbol("MF1"),
            sp.Symbol("MF2"),
        ],
    )

    symmetry_difference = (
        K2 - K2.T
    ).applyfunc(sp.simplify)

    print("=" * 72)
    print("T3 FLAVOR-LIFT MATCHING REGRESSION")
    print("=" * 72)
    print()
    print("Original one-generation Matchete C5:")
    print(kappa_1g)
    print()
    print("Extracted loop kernel:")
    print(kernel)
    print()
    print("One-generation reconstructed C5:")
    print(K1[0, 0])
    print()
    print("One-generation difference:")
    print(one_generation_difference)
    print()
    print("Two-flavor symmetry check K-K^T:")
    print(symmetry_difference)
    print()

    if one_generation_difference != 0:
        print("FAIL: one-generation reconstruction")
        return 1

    if any(x != 0 for x in symmetry_difference):
        print("FAIL: flavor C5 is not symmetric")
        return 1

    print("PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
