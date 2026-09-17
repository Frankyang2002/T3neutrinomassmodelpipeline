from __future__ import annotations

"""Validate completeness of mixed scalar-quartic invariant sectors in a Matchete seed.

For two SU(2) multiplets A and B with dimensions d_A and d_B, the quartics

    A^\dagger A B^\dagger B

span one invariant for every common total isospin J appearing in

    A^\dagger \otimes A  and  B^\dagger \otimes B.

For SU(2), j \otimes j contains J = 0,1,...,2j once, hence

    N_invariants(A,B) = min(d_A, d_B).

This gives the expected dimensions:
    H(2)-singlet(1): 1
    H(2)-doublet(2): 2
    H(2)-triplet(3): 2
    doublet-doublet: 2
    triplet-triplet: 3
    singlet-triplet: 1

The script compares this representation-theory count against the distinct
couplings actually exported in ScalarQuarticTerms.  It deliberately does not
assume coupling names such as Inv1/Inv2 beyond extracting the unique coupling
symbol from each quartic term.
"""

import argparse
from collections import Counter
import json
from pathlib import Path


FIELD_DIMS_KEY = {
    "H": "H",
    "NewScalar1": "S1",
    "NewScalar2": "S2",
}


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


def _split_top_level(text: str) -> list[str]:
    out = []
    start = 0
    sq = cu = pa = 0
    for i, ch in enumerate(text):
        if ch == "[":
            sq += 1
        elif ch == "]":
            sq -= 1
        elif ch == "{":
            cu += 1
        elif ch == "}":
            cu -= 1
        elif ch == "(":
            pa += 1
        elif ch == ")":
            pa -= 1
        elif ch == "," and sq == 0 and cu == 0 and pa == 0:
            out.append(text[start:i].strip())
            start = i + 1
    out.append(text[start:].strip())
    return out


def _coupling_name(term: str) -> str:
    start = term.find("Coupling[")
    if start < 0:
        raise ValueError("quartic term contains no Coupling[...]")
    open_index = start + len("Coupling")
    close = _matching_bracket(term, open_index)
    return _split_top_level(term[open_index + 1:close])[0]


def _fields(term: str) -> tuple[str, ...]:
    names = []
    pos = 0
    while True:
        start = term.find("Field[", pos)
        if start < 0:
            break
        open_index = start + len("Field")
        close = _matching_bracket(term, open_index)
        args = _split_top_level(term[open_index + 1:close])
        if len(args) >= 2 and args[1] == "Scalar" and args[0] in FIELD_DIMS_KEY:
            names.append(args[0])
        pos = close + 1
    return tuple(names)


def _sector(fields: tuple[str, ...]) -> str | None:
    c = Counter(fields)
    if c == Counter({"H": 2, "NewScalar1": 2}):
        return "H1"
    if c == Counter({"H": 2, "NewScalar2": 2}):
        return "H2"
    if c == Counter({"NewScalar1": 2, "NewScalar2": 2}):
        return "12"
    return None


def expected_invariant_count(d_a: int, d_b: int) -> int:
    return min(int(d_a), int(d_b))


def validate(seed_path: Path, d_s1: int, d_s2: int) -> dict:
    seed = json.loads(seed_path.read_text(encoding="utf-8"))
    dims = {"H": 2, "S1": int(d_s1), "S2": int(d_s2)}
    sector_fields = {
        "H1": ("H", "S1"),
        "H2": ("H", "S2"),
        "12": ("S1", "S2"),
    }

    couplings = {name: set() for name in sector_fields}

    for record in seed.get("ScalarQuarticTerms", []):
        term = str(record.get("TermInputForm", ""))
        sec = _sector(_fields(term))
        if sec is not None:
            couplings[sec].add(_coupling_name(term))

    rows = []
    for sec, (a, b) in sector_fields.items():
        expected = expected_invariant_count(dims[a], dims[b])
        observed_names = sorted(couplings[sec])
        observed = len(observed_names)
        rows.append({
            "sector": sec,
            "fields": [a, b],
            "dimensions": [dims[a], dims[b]],
            "expected_invariant_count": expected,
            "observed_distinct_coupling_count": observed,
            "observed_couplings": observed_names,
            "missing_count": max(0, expected - observed),
            "status": "PASS" if expected == observed else "FAIL",
        })

    return {
        "status": "Success" if all(r["status"] == "PASS" for r in rows) else "Incomplete",
        "seed": str(seed_path),
        "dimensions": {"dH": 2, "dS1": int(d_s1), "dS2": int(d_s2)},
        "representation_theory_rule": "N(A^dagger A B^dagger B)=min(d_A,d_B) for SU(2)",
        "sectors": rows,
    }


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("quartic_seed_json", type=Path)
    p.add_argument("dS1", type=int)
    p.add_argument("dS2", type=int)
    p.add_argument("--output", type=Path, default=None)
    a = p.parse_args()

    payload = validate(a.quartic_seed_json, a.dS1, a.dS2)

    print(f"seed: {a.quartic_seed_json}")
    print("rule: N_invariants = min(d_A, d_B)")
    print()
    print(f"{'sector':<7} {'dims':<8} {'expected':>8} {'observed':>8} {'missing':>8}  couplings")
    print("-" * 80)
    for row in payload["sectors"]:
        dims = "x".join(map(str, row["dimensions"]))
        print(
            f"{row['sector']:<7} {dims:<8} "
            f"{row['expected_invariant_count']:>8} "
            f"{row['observed_distinct_coupling_count']:>8} "
            f"{row['missing_count']:>8}  "
            f"{row['observed_couplings']}  {row['status']}"
        )

    if a.output:
        a.output.parent.mkdir(parents=True, exist_ok=True)
        a.output.write_text(json.dumps(payload, indent=2), encoding="utf-8")
        print()
        print(f"JSON summary: {a.output}")

    print(f"OVERALL: {payload['status']}")
    return 0 if payload["status"] == "Success" else 1


if __name__ == "__main__":
    raise SystemExit(main())
