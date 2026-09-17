from __future__ import annotations

"""Inspect the actual Matchete scalar-quartic basis relevant to beta_lambdaT3.

This is a diagnostic companion to MixingQuarticRecoupling.py.  It does not
assume RGBeta coupling names such as lambdaH1Adj.  Instead it:

1. reads the exact EFT1 Matchete scalar-quartic seed;
2. lists every coupling actually present in ScalarQuarticTerms;
3. builds the exact real-scalar lambda_abcd tensor;
4. extracts each coupling-direction tensor;
5. computes its one-loop pure-scalar cross-contraction with the holomorphic
   lambdaT3 tensor;
6. projects back onto the lambdaT3 direction and reports the exact coefficient
   and closure residual.

This establishes which Matchete invariant-basis directions correspond to the
non-singlet recouplings before any basis conversion to RGBeta's Adj/Cross
couplings is attempted.
"""

import argparse
from collections import Counter
from itertools import combinations_with_replacement
import json
from math import factorial
from pathlib import Path
import re
import sys

import sympy as sp


PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from RGE.running.EFT1QuarticAdapter import (
    SparseQuarticTensor,
    load_and_build_eft1_quartic_tensor,
)


def _text(x: sp.Expr) -> str:
    return sp.sstr(sp.factor(sp.simplify(x)))


def _matching_bracket(text: str, open_index: int) -> int:
    depth = 0
    for i in range(open_index, len(text)):
        if text[i] == "[":
            depth += 1
        elif text[i] == "]":
            depth -= 1
            if depth == 0:
                return i
    raise ValueError("unbalanced brackets")


def _coupling_name(term: str) -> str:
    start = term.find("Coupling[")
    if start < 0:
        raise ValueError(f"No Coupling[...] in term: {term}")
    open_index = start + len("Coupling")
    close = _matching_bracket(term, open_index)
    body = term[open_index + 1:close]
    # Coupling[name, {}, 0]; the first comma is top level here.
    return body.split(",", 1)[0].strip()


def _basis_tensor(
    full: SparseQuarticTensor,
    coupling_name: str,
    *,
    identify_conjugate: bool = True,
) -> SparseQuarticTensor:
    symbol = sp.Symbol(coupling_name)
    conjugate = sp.conjugate(symbol)
    entries = {}

    for key, raw in full.nonzero_items():
        expr = sp.expand(raw)
        if identify_conjugate:
            expr = sp.expand(expr.xreplace({conjugate: symbol}))
        coefficient = sp.simplify(expr.coeff(symbol))
        if coefficient != 0:
            entries[key] = coefficient

    return SparseQuarticTensor(entries)


def _dimension(full: SparseQuarticTensor) -> int:
    return max((max(key) for key in full.entries), default=0)


def _weight(key: tuple[int, int, int, int]) -> int:
    counts = Counter(key)
    value = factorial(4)
    for n in counts.values():
        value //= factorial(n)
    return value


def _inner(a: SparseQuarticTensor, b: SparseQuarticTensor) -> sp.Expr:
    total = sp.S.Zero
    for key in set(a.entries) | set(b.entries):
        total += _weight(key) * sp.conjugate(a[key]) * b[key]
    return sp.simplify(total)


def _cross_component(
    t: SparseQuarticTensor,
    x: SparseQuarticTensor,
    a: int, b: int, c: int, d: int,
    n: int,
) -> sp.Expr:
    total = sp.S.Zero
    for (p, q), (r, s) in (
        ((a, b), (c, d)),
        ((a, c), (b, d)),
        ((a, d), (b, c)),
    ):
        for e in range(1, n + 1):
            for f in range(1, n + 1):
                total += (
                    t[p, q, e, f] * x[e, f, r, s]
                    + x[p, q, e, f] * t[e, f, r, s]
                )
    return sp.simplify(total)


def _cross_tensor(
    t: SparseQuarticTensor,
    x: SparseQuarticTensor,
    n: int,
) -> SparseQuarticTensor:
    entries = {}
    for key in combinations_with_replacement(range(1, n + 1), 4):
        value = _cross_component(t, x, *key, n)
        if value != 0:
            entries[key] = value
    return SparseQuarticTensor(entries)


def _project(
    generated: SparseQuarticTensor,
    target: SparseQuarticTensor,
) -> tuple[sp.Expr, int]:
    norm = sp.simplify(_inner(target, target))
    if norm == 0:
        raise ValueError("lambdaT3 tensor has zero norm")

    coeff = sp.simplify(_inner(target, generated) / norm)
    residual_count = 0
    for key in set(generated.entries) | set(target.entries):
        if sp.simplify(generated[key] - coeff * target[key]) != 0:
            residual_count += 1
    return coeff, residual_count


def inspect(seed_path: Path, rgbeta_path: Path) -> dict:
    seed = json.loads(seed_path.read_text(encoding="utf-8"))
    names = sorted({
        _coupling_name(record["TermInputForm"])
        for record in seed.get("ScalarQuarticTerms", [])
    })

    full = load_and_build_eft1_quartic_tensor(seed_path, rgbeta_path)
    n = _dimension(full)

    t3 = _basis_tensor(full, "lambdaT3", identify_conjugate=False)
    if not t3.entries:
        raise ValueError("No holomorphic lambdaT3 tensor found")

    rows = []
    for name in names:
        if name == "lambdaT3":
            continue

        direction = _basis_tensor(full, name, identify_conjugate=True)
        if not direction.entries:
            continue

        generated = _cross_tensor(t3, direction, n)
        coeff, residual_count = _project(generated, t3)

        rows.append({
            "coupling": name,
            "direction_components": len(direction.entries),
            "projection_coefficient": _text(coeff),
            "residual_nonzero_components": residual_count,
            "closes_on_lambdaT3": residual_count == 0,
        })

    return {
        "seed": str(seed_path),
        "rgbeta": str(rgbeta_path),
        "couplings_in_seed": names,
        "real_scalar_dimension": n,
        "lambdaT3_components": len(t3.entries),
        "recoupling_rows": rows,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("quartic_seed_json", type=Path)
    parser.add_argument("--rgbeta", type=Path, default=None)
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()

    rgbeta = args.rgbeta or args.quartic_seed_json.with_name("eft1_rgbeta_rge.json")
    result = inspect(args.quartic_seed_json, rgbeta)

    print("Couplings present in exact Matchete scalar seed:")
    for name in result["couplings_in_seed"]:
        print(f"  {name}")

    print()
    print("Pure-scalar recoupling onto holomorphic lambdaT3:")
    print(f"{'coupling':<22} {'R':>12} {'residual nz':>12} {'closure':>9}")
    print("-" * 61)
    for row in result["recoupling_rows"]:
        print(
            f"{row['coupling']:<22} "
            f"{row['projection_coefficient']:>12} "
            f"{row['residual_nonzero_components']:>12} "
            f"{('PASS' if row['closes_on_lambdaT3'] else 'FAIL'):>9}"
        )

    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2), encoding="utf-8")
        print()
        print(f"JSON summary: {args.output}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
