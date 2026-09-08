from __future__ import annotations

"""Consolidated Wolfram->Python T3 tensor/flavor bridge regression."""

import argparse
import json
from pathlib import Path

import sympy as sp

from RGE.matching.FlavorMatchedC5 import (
    build_flavor_matched_c5,
    extract_t3_loop_kernel,
)
from RGE.matching.MatchedEFTRGE import parse_matchete_c5
from tests.python.Weinberg.T3YukawaAdapter import (
    T3WeylConvention,
    real_yukawa_tensor_from_exchange,
)


TEST_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_EXCHANGE = (
    TEST_ROOT / "wolfram" / "output" / "T3_B_alpha_m1" / "rge_tensor_exchange.json"
)
DEFAULT_C5 = (
    TEST_ROOT / "wolfram" / "output" / "T3_B_alpha_m1" / "c5_coefficient.txt"
)


def check_yukawa_exchange(exchange_path: Path) -> None:
    data = json.loads(exchange_path.read_text(encoding="utf-8"))

    scalar_model, fermion_basis, tensor = real_yukawa_tensor_from_exchange(
        data,
        T3WeylConvention(
            lepton_multiplicity=1,
            heavy_fermion_multiplicity=1,
            symmetrize_fermion_indices=True,
        ),
    )

    assert scalar_model.total_real_scalar_dimension > 0
    assert fermion_basis.dimension > 0
    assert tensor

    # Adapter contract: y_ija is symmetric in the two Weyl slots.
    for (i, j, a), value in tensor.items():
        assert sp.simplify(tensor.get((j, i, a), 0) - value) == 0


def check_flavor_lift(c5_path: Path) -> None:
    kappa_1g = parse_matchete_c5(c5_path.read_text(encoding="utf-8"))
    kernel = extract_t3_loop_kernel(kappa_1g)

    y1, y2 = sp.symbols("y1 y2")
    K1 = build_flavor_matched_c5(
        kernel,
        sp.Matrix([[y1]]),
        sp.Matrix([[y2]]),
        heavy_masses=None,
    )
    assert sp.simplify(K1[0, 0] - kappa_1g) == 0

    a11, a12, a21, a22 = sp.symbols("a11 a12 a21 a22")
    b11, b12, b21, b22 = sp.symbols("b11 b12 b21 b22")

    K2 = build_flavor_matched_c5(
        kernel,
        sp.Matrix([[a11, a12], [a21, a22]]),
        sp.Matrix([[b11, b12], [b21, b22]]),
        heavy_masses=[sp.Symbol("MF1"), sp.Symbol("MF2")],
    )

    symmetry = (K2 - K2.T).applyfunc(sp.simplify)
    assert all(x == 0 for x in symmetry)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--exchange", type=Path, default=DEFAULT_EXCHANGE)
    parser.add_argument("--c5", type=Path, default=DEFAULT_C5)
    args = parser.parse_args()

    if not args.exchange.exists():
        raise FileNotFoundError(
            f"Missing {args.exchange}. Run the consolidated Wolfram RGE export test first."
        )
    if not args.c5.exists():
        raise FileNotFoundError(
            f"Missing {args.c5}. Run the normal matching/C5 pipeline first."
        )

    check_yukawa_exchange(args.exchange)
    check_flavor_lift(args.c5)

    print("PASS: Wolfram Yukawa exchange -> Python y_ija")
    print("PASS: one-generation C5 reconstruction")
    print("PASS: multi-flavor C5 symmetry")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
