from __future__ import annotations

"""Exact scalar-quartic recoupling factors in beta_lambdaT3.

This module uses the existing Matchete -> real-scalar quartic tensor adapter.
No A--E coefficient table is encoded.

Convention
----------
The adapter constructs the fully symmetric real tensor lambda_abcd from

    V4 = (1/4!) lambda_abcd phi_a phi_b phi_c phi_d.

For this convention the pure-scalar one-loop contribution is

    beta_abcd|_{lambda^2}
      = sum_{e,f} [
          lambda_abef lambda_efcd
        + lambda_acef lambda_efbd
        + lambda_adef lambda_efbc
        ].

This is equivalent to (1/8) times the full 24-permutation sum.

For a target real quartic X, write

    lambda_abcd = lambdaT3 * T_abcd + X * X_abcd + ...

The coefficient linear in lambdaT3 * X is then built exactly from the
cross-contractions T*X + X*T in the three channels and projected back onto T.

The output reports
- the exact projection coefficient R_X,
- the projection residual,
- whether the full generated tensor is proportional to the lambdaT3 tensor.

The lambdaT3 coupling is kept holomorphic: conjugate(lambdaT3) is treated as
an independent direction.  The adjoint/cross quartics tested here are real, so
their conjugates are identified with themselves before extracting their basis
tensors.
"""

import argparse
from collections import Counter
from dataclasses import dataclass, asdict
from fractions import Fraction
from itertools import combinations_with_replacement
import json
from math import factorial
from pathlib import Path
import sys
from typing import Iterable

import sympy as sp


PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from RGE.running.EFT1QuarticAdapter import (
    SparseQuarticTensor,
    load_and_build_eft1_quartic_tensor,
)


REAL_TARGETS = (
    "lambdaH1Adj",
    "lambdaH2Adj",
    "lambda12Adj",
    "lambda12Cross",
)


@dataclass(frozen=True)
class ProjectionResult:
    coupling: str
    present_in_seed: bool
    t3_nonzero_components: int
    generated_nonzero_components: int
    projection_coefficient: str | None
    norm_t3: str
    max_residual: str
    residual_nonzero_components: int
    proportional_to_t3: bool


def _text(expr: sp.Expr) -> str:
    return sp.sstr(sp.factor(sp.simplify(expr)))


def _multiplicity_weight(key: tuple[int, int, int, int]) -> int:
    """Number of ordered index tuples represented by one sorted key."""
    counts = Counter(key)
    weight = factorial(4)
    for count in counts.values():
        weight //= factorial(count)
    return weight


def _basis_tensor(
    full: SparseQuarticTensor,
    coupling_name: str,
    *,
    identify_conjugate: bool,
) -> SparseQuarticTensor:
    """Coefficient tensor of one coupling in the full quartic tensor."""
    symbol = sp.Symbol(coupling_name)
    conjugate = sp.conjugate(symbol)

    entries: dict[tuple[int, int, int, int], sp.Expr] = {}
    for key, raw in full.nonzero_items():
        expr = sp.expand(raw)
        if identify_conjugate:
            expr = sp.expand(expr.xreplace({conjugate: symbol}))
        coefficient = sp.simplify(expr.coeff(symbol))
        if coefficient != 0:
            entries[key] = coefficient

    return SparseQuarticTensor(entries)


def _scalar_dimension(full: SparseQuarticTensor) -> int:
    return max((max(key) for key in full.entries), default=0)


def _cross_component(
    t3: SparseQuarticTensor,
    other: SparseQuarticTensor,
    a: int,
    b: int,
    c: int,
    d: int,
    scalar_dimension: int,
) -> sp.Expr:
    """Coefficient of lambdaT3*X in beta_abcd from pure scalar quartics."""
    total = sp.S.Zero

    pairings = (
        ((a, b), (c, d)),
        ((a, c), (b, d)),
        ((a, d), (b, c)),
    )

    for (p, q), (r, s) in pairings:
        for e in range(1, scalar_dimension + 1):
            for f in range(1, scalar_dimension + 1):
                total += (
                    t3[p, q, e, f] * other[e, f, r, s]
                    + other[p, q, e, f] * t3[e, f, r, s]
                )

    return sp.simplify(total)


def _build_cross_tensor(
    t3: SparseQuarticTensor,
    other: SparseQuarticTensor,
    scalar_dimension: int,
) -> SparseQuarticTensor:
    entries: dict[tuple[int, int, int, int], sp.Expr] = {}

    # The largest T3 model here has 16 real scalars, hence C(19,4)=3876 sorted
    # output components.  Scanning the complete symmetric output space gives a
    # genuine closure/residual check rather than checking only T3 support.
    for key in combinations_with_replacement(
        range(1, scalar_dimension + 1),
        4,
    ):
        value = _cross_component(
            t3,
            other,
            *key,
            scalar_dimension,
        )
        if value != 0:
            entries[key] = value

    return SparseQuarticTensor(entries)


def _inner_product(
    left: SparseQuarticTensor,
    right: SparseQuarticTensor,
) -> sp.Expr:
    """Hermitian tensor inner product summed over all ordered a,b,c,d."""
    keys = set(left.entries) | set(right.entries)
    total = sp.S.Zero
    for key in keys:
        total += (
            _multiplicity_weight(key)
            * sp.conjugate(left[key])
            * right[key]
        )
    return sp.simplify(total)


def _project(
    generated: SparseQuarticTensor,
    t3: SparseQuarticTensor,
) -> tuple[sp.Expr, SparseQuarticTensor, sp.Expr]:
    norm = sp.simplify(_inner_product(t3, t3))
    if norm == 0:
        raise ValueError("lambdaT3 tensor has zero norm.")

    coefficient = sp.simplify(_inner_product(t3, generated) / norm)

    residual_entries: dict[tuple[int, int, int, int], sp.Expr] = {}
    for key in set(generated.entries) | set(t3.entries):
        residual = sp.simplify(generated[key] - coefficient * t3[key])
        if residual != 0:
            residual_entries[key] = residual

    return coefficient, SparseQuarticTensor(residual_entries), norm


def _max_residual_text(residual: SparseQuarticTensor) -> str:
    if not residual.entries:
        return "0"

    # Exact expressions do not have a useful universal ordering.  Return a
    # representative simplified nonzero residual rather than a floating norm.
    key = sorted(residual.entries)[0]
    return _text(residual.entries[key])


def recoupling_for_seed(
    quartic_seed_path: Path,
    rgbeta_path: Path,
    *,
    targets: Iterable[str] = REAL_TARGETS,
) -> dict:
    full = load_and_build_eft1_quartic_tensor(
        quartic_seed_path,
        rgbeta_path,
    )
    scalar_dimension = _scalar_dimension(full)

    # lambdaT3 is complex.  Keep only its holomorphic tensor direction.
    t3 = _basis_tensor(
        full,
        "lambdaT3",
        identify_conjugate=False,
    )
    if not t3.entries:
        raise ValueError(
            f"No holomorphic lambdaT3 tensor found in {quartic_seed_path}."
        )

    results: list[ProjectionResult] = []

    for coupling_name in targets:
        other = _basis_tensor(
            full,
            coupling_name,
            identify_conjugate=True,
        )

        if not other.entries:
            results.append(
                ProjectionResult(
                    coupling=coupling_name,
                    present_in_seed=False,
                    t3_nonzero_components=len(t3.entries),
                    generated_nonzero_components=0,
                    projection_coefficient=None,
                    norm_t3=_text(_inner_product(t3, t3)),
                    max_residual="0",
                    residual_nonzero_components=0,
                    proportional_to_t3=True,
                )
            )
            continue

        generated = _build_cross_tensor(
            t3,
            other,
            scalar_dimension,
        )
        coefficient, residual, norm = _project(generated, t3)

        results.append(
            ProjectionResult(
                coupling=coupling_name,
                present_in_seed=True,
                t3_nonzero_components=len(t3.entries),
                generated_nonzero_components=len(generated.entries),
                projection_coefficient=_text(coefficient),
                norm_t3=_text(norm),
                max_residual=_max_residual_text(residual),
                residual_nonzero_components=len(residual.entries),
                proportional_to_t3=not residual.entries,
            )
        )

    return {
        "status": (
            "Success"
            if all(item.proportional_to_t3 for item in results)
            else "Failed"
        ),
        "quartic_seed": str(quartic_seed_path),
        "rgbeta_file": str(rgbeta_path),
        "scalar_dimension": scalar_dimension,
        "convention": (
            "V4=lambda_abcd phi_a phi_b phi_c phi_d/4!; "
            "beta_lambda|scalar^2 has the three pair-contraction channels"
        ),
        "lambdaT3_holomorphic_nonzero_components": len(t3.entries),
        "results": [asdict(item) for item in results],
    }


def _default_rgbeta(seed: Path) -> Path:
    candidate = seed.with_name("eft1_rgbeta_rge.json")
    if not candidate.is_file():
        raise FileNotFoundError(
            f"Could not infer EFT1 RGBeta metadata next to seed: {candidate}"
        )
    return candidate


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Derive lambdaT3 non-singlet scalar-quartic recoupling factors "
            "from exact Matchete real-scalar tensors."
        )
    )
    parser.add_argument(
        "quartic_seed_json",
        type=Path,
        help="eft1_after_F_scalar_quartic_seed.json",
    )
    parser.add_argument(
        "--rgbeta",
        type=Path,
        default=None,
        help="eft1_rgbeta_rge.json; inferred from the seed directory by default.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Optional JSON output path.",
    )
    args = parser.parse_args()

    rgbeta = args.rgbeta or _default_rgbeta(args.quartic_seed_json)
    payload = recoupling_for_seed(
        args.quartic_seed_json,
        rgbeta,
    )

    print(f"seed: {args.quartic_seed_json}")
    print(f"real scalar dimension: {payload['scalar_dimension']}")
    print(
        "lambdaT3 holomorphic tensor components: "
        f"{payload['lambdaT3_holomorphic_nonzero_components']}"
    )
    print()
    print(
        f"{'coupling':<20} {'present':<8} {'R_X':>12} "
        f"{'residual nz':>12} {'closure':>9}"
    )
    print("-" * 68)

    for item in payload["results"]:
        coefficient = item["projection_coefficient"]
        print(
            f"{item['coupling']:<20} "
            f"{str(item['present_in_seed']):<8} "
            f"{str(coefficient):>12} "
            f"{item['residual_nonzero_components']:>12} "
            f"{('PASS' if item['proportional_to_t3'] else 'FAIL'):>9}"
        )

    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(
            json.dumps(payload, indent=2),
            encoding="utf-8",
        )
        print()
        print(f"JSON summary: {args.output}")

    print(f"OVERALL: {payload['status']}")
    return 0 if payload["status"] == "Success" else 1


if __name__ == "__main__":
    raise SystemExit(main())
