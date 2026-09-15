from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path
import sympy as sp

from RGE.running.EFT1WilsonAdapter import build_eft1_wilson_tensor


def canonical_items(tensor):
    items = {}
    for key, raw in tensor.nonzero_items():
        value = sp.simplify(raw)
        if value != 0:
            items[tuple(key)] = value
    return items


def check_scaled_term(seed_payload, term, d_s1, d_s2, scale=sp.Integer(7)):
    one = copy.deepcopy(seed_payload)
    one["TreeWilsonTerms"] = [copy.deepcopy(term)]
    one["TreeWilsonTermCount"] = 1

    baseline = build_eft1_wilson_tensor(
        one,
        d_s1=d_s1,
        d_s2=d_s2,
        chirality="PL",
        project_operator_symmetry=True,
    )
    baseline_items = canonical_items(baseline)

    scaled = copy.deepcopy(one)
    scaled_term = scaled["TreeWilsonTerms"][0]
    scaled_term["TermInputForm"] = (
        f"({sp.sstr(scale)})*({scaled_term['TermInputForm']})"
    )

    scaled_tensor = build_eft1_wilson_tensor(
        scaled,
        d_s1=d_s1,
        d_s2=d_s2,
        chirality="PL",
        project_operator_symmetry=True,
    )
    scaled_items = canonical_items(scaled_tensor)

    keys = sorted(set(baseline_items) | set(scaled_items))
    failures = []
    observed_ratios = set()

    for key in keys:
        a = sp.simplify(baseline_items.get(key, 0))
        b = sp.simplify(scaled_items.get(key, 0))

        if a == 0:
            if b != 0:
                failures.append((key, a, b, "new component appeared"))
            continue

        ratio = sp.simplify(b / a)
        observed_ratios.add(sp.sstr(ratio))
        if sp.simplify(b - scale * a) != 0:
            failures.append((key, a, b, f"ratio={ratio}"))

    return {
        "term_index": int(term["Index"]),
        "cg_names": term.get("CGNames", []),
        "baseline_nonzero": len(baseline_items),
        "scaled_nonzero": len(scaled_items),
        "expected_scale": sp.sstr(scale),
        "observed_ratios": sorted(observed_ratios),
        "passes_exact_scaling": not failures,
        "failure_count": len(failures),
        "first_failures": [
            {
                "component": list(key),
                "baseline": sp.sstr(a),
                "scaled": sp.sstr(b),
                "detail": detail,
            }
            for key, a, b, detail in failures[:10]
        ],
    }


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Check whether EFT1WilsonAdapter preserves arbitrary outer "
            "numerical prefactors from Matchete TermInputForm."
        )
    )
    parser.add_argument("wilson_seed_json", type=Path)
    parser.add_argument("--d-s1", type=int, required=True)
    parser.add_argument("--d-s2", type=int, required=True)
    parser.add_argument("--scale", type=int, default=7)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    seed = json.loads(args.wilson_seed_json.read_text(encoding="utf-8"))

    pl_terms = [
        term for term in seed.get("TreeWilsonTerms", [])
        if term.get("Chirality") == "PL"
    ]
    if not pl_terms:
        raise RuntimeError("No PL terms found in the Wilson seed.")

    results = [
        check_scaled_term(
            seed,
            term,
            d_s1=args.d_s1,
            d_s2=args.d_s2,
            scale=sp.Integer(args.scale),
        )
        for term in pl_terms
    ]

    payload = {
        "status": "Success",
        "seed": str(args.wilson_seed_json),
        "all_terms_preserve_outer_prefactor": all(
            r["passes_exact_scaling"] for r in results
        ),
        "terms": results,
    }

    text = json.dumps(payload, indent=2)
    print(text)

    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")

    return 0 if payload["all_terms_preserve_outer_prefactor"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
