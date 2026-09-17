from __future__ import annotations

"""Inspect exact Matchete scalar-quartic seed sectors.

Purpose
-------
ScalarQuarticBasisMap.py showed that some triplet portal sectors apparently
contain only one coupling direction, while representation theory requires two.
Before changing any basis map, inspect the actual exported ScalarQuarticTerms
without assuming coupling names.

For every quartic term this script extracts
- coupling name,
- field multiset including conjugation,
- CG tensor names,
- exact TermInputForm.

It then groups terms into H1, H2 and 12 sectors by field content.
"""

import argparse
from collections import Counter
import json
from pathlib import Path
import re


FIELD_NAMES = ("H", "NewScalar1", "NewScalar2")


def _matching_bracket(text: str, open_index: int) -> int:
    if text[open_index] != "[":
        raise ValueError("expected '['")
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
        return "<none>"
    open_index = start + len("Coupling")
    close = _matching_bracket(term, open_index)
    args = _split_top_level(term[open_index + 1:close])
    return args[0].strip()


def _fields(term: str) -> tuple[str, ...]:
    found = []
    pos = 0
    while True:
        start = term.find("Field[", pos)
        if start < 0:
            break
        open_index = start + len("Field")
        close = _matching_bracket(term, open_index)
        args = _split_top_level(term[open_index + 1:close])
        name = args[0].strip()

        conjugated = False
        if start >= 4 and term[start - 4:start] == "Bar[":
            try:
                conjugated = _matching_bracket(term, start - 1) == close + 1
            except Exception:
                pass

        if name in FIELD_NAMES:
            found.append(("Bar[" + name + "]") if conjugated else name)
        pos = close + 1

    return tuple(found)


def _cg_names(term: str) -> tuple[str, ...]:
    names = []
    pos = 0
    while True:
        start = term.find("CG[", pos)
        if start < 0:
            break
        open_index = start + len("CG")
        close = _matching_bracket(term, open_index)
        args = _split_top_level(term[open_index + 1:close])
        names.append(args[0].strip())
        pos = close + 1
    return tuple(names)


def _base_field(label: str) -> str:
    return label[4:-1] if label.startswith("Bar[") else label


def _sector(fields: tuple[str, ...]) -> str | None:
    counts = Counter(_base_field(x) for x in fields)
    if counts == Counter({"H": 2, "NewScalar1": 2}):
        return "H1"
    if counts == Counter({"H": 2, "NewScalar2": 2}):
        return "H2"
    if counts == Counter({"NewScalar1": 2, "NewScalar2": 2}):
        return "12"
    return None


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("quartic_seed_json", type=Path)
    parser.add_argument("--json-output", type=Path, default=None)
    args = parser.parse_args()

    payload = json.loads(args.quartic_seed_json.read_text(encoding="utf-8"))
    records = []

    for index, record in enumerate(payload.get("ScalarQuarticTerms", []), 1):
        term = str(record.get("TermInputForm", ""))
        fields = _fields(term)
        sec = _sector(fields)
        if sec is None:
            continue

        records.append({
            "index": index,
            "sector": sec,
            "coupling": _coupling_name(term),
            "fields": list(fields),
            "cg_names": list(_cg_names(term)),
            "term_input_form": term,
        })

    print(f"seed: {args.quartic_seed_json}")
    print()

    for sec in ("H1", "H2", "12"):
        subset = [r for r in records if r["sector"] == sec]
        print(f"[{sec}] {len(subset)} exported term(s)")
        names = sorted(set(r["coupling"] for r in subset))
        print("  couplings:", names)
        print("  independent coupling count:", len(names))
        for r in subset:
            print(
                f"  term {r['index']}: coupling={r['coupling']} "
                f"fields={r['fields']} CG={r['cg_names']}"
            )
            print("    " + r["term_input_form"])
        print()

    if args.json_output:
        args.json_output.parent.mkdir(parents=True, exist_ok=True)
        args.json_output.write_text(
            json.dumps({"seed": str(args.quartic_seed_json), "records": records}, indent=2),
            encoding="utf-8",
        )
        print(f"JSON summary: {args.json_output}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
