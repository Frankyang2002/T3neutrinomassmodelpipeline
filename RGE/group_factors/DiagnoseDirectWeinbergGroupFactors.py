from __future__ import annotations

"""Extract the direct EFT1 C12 -> Weinberg mixing coefficient across T3 models.

This is the first C5 group-factor diagnostic.  It does not fit a formula yet.

Authoritative fixed-order counting:
    tree LLSS boundary               O(hbar^0)
    direct LLSS -> O5 running        O(hbar^1)
    LLSS self-running                O(hbar^1)
    self-running + scalar matching   O(hbar^2), excluded here

For each successful run we read
    data/eft1_wilson_flavor_at_S_threshold.json

and report the stored flavor-blind prefactor in

    16*pi^2 beta_kappa = prefactor * C12[p,q].

We also divide out lambdaT3 when possible, because the remaining coefficient is
the SU(2) operator-geometry/group factor that we want to reformulate.
"""

import argparse
import json
from pathlib import Path
import re

import sympy as sp


MODEL_DIMS = {
    "A": (1, 3, 2),
    "B": (2, 2, 1),
    "C": (2, 2, 3),
    "D": (3, 1, 2),
    "E": (3, 3, 2),
}


def model_alpha(path: Path):
    name = path.parents[1].name
    match = re.fullmatch(r"T3_([A-E])_alpha_([mp])(\d+)", name)
    if not match:
        raise ValueError(name)
    model, sign, mag = match.groups()
    alpha = int(mag)
    return model, (-alpha if sign == "m" else alpha)


def expr(raw):
    return sp.sympify(
        str(raw).replace("^", "**"),
        locals={
            "lambdaT3": sp.Symbol("lambdaT3"),
            "sqrt": sp.sqrt,
            "log": sp.log,
            "conjugate": sp.conjugate,
        },
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--root",
        type=Path,
        default=Path("output/full/hypercharge"),
    )
    parser.add_argument(
        "--json-output",
        type=Path,
        default=Path(
            "output/group_factors/direct_weinberg_group_factor_diagnostic.json"
        ),
    )
    args = parser.parse_args()

    paths = sorted(
        args.root.glob(
            "T3_*_alpha_*/data/eft1_wilson_flavor_at_S_threshold.json"
        )
    )
    if not paths:
        raise SystemExit(
            "No eft1_wilson_flavor_at_S_threshold.json files found under "
            f"{args.root}"
        )

    lam = sp.Symbol("lambdaT3")
    rows = []

    print("Direct EFT1 C12 -> Weinberg group-factor diagnostic")
    print("model alpha dS1 dS2 dF  prefactor                     prefactor/lambdaT3")
    print("-" * 92)

    for path in paths:
        payload = json.loads(path.read_text(encoding="utf-8-sig"))
        if payload.get("status") != "Success":
            continue

        model, alpha = model_alpha(path)
        d1, d2, dF = MODEL_DIMS[model]

        direct = payload.get("direct_weinberg_beta") or {}
        raw_prefactor = direct.get("flavor_blind_prefactor")
        if raw_prefactor is None:
            continue

        prefactor = sp.factor(expr(raw_prefactor))
        reduced = sp.factor(prefactor / lam) if prefactor.has(lam) else None

        print(
            f"T3-{model} {alpha:>5} {d1:>3} {d2:>3} {dF:>2}  "
            f"{sp.sstr(prefactor):<28} "
            f"{sp.sstr(reduced) if reduced is not None else '<no lambdaT3>'}"
        )

        rows.append({
            "model": model,
            "alpha": alpha,
            "dS1": d1,
            "dS2": d2,
            "dF": dF,
            "prefactor": sp.sstr(prefactor),
            "prefactor_over_lambdaT3": (
                sp.sstr(reduced) if reduced is not None else None
            ),
            "source": str(path),
            "one_generation_beta_kappa_16pi2": direct.get(
                "beta_kappa_16pi2_one_generation"
            ),
            "canonical_c12_one_generation": direct.get(
                "canonical_c12_one_generation"
            ),
        })

    # Check alpha-independence model by model.  Direct scalar-quartic mixing
    # should depend on SU(2) geometry, not hypercharge.
    grouped = {}
    for row in rows:
        grouped.setdefault(row["model"], set()).add(
            row["prefactor_over_lambdaT3"]
        )

    alpha_independent = all(len(values) == 1 for values in grouped.values())

    print()
    print("Distinct reduced factors by model:")
    for model in sorted(grouped):
        print(f"  T3-{model}: {sorted(grouped[model])}")
    print(f"alpha-independent: {alpha_independent}")

    result = {
        "status": "Success" if rows else "NoData",
        "alpha_independent": alpha_independent,
        "distinct_reduced_factors": {
            model: sorted(values)
            for model, values in sorted(grouped.items())
        },
        "rows": rows,
    }

    args.json_output.parent.mkdir(parents=True, exist_ok=True)
    args.json_output.write_text(json.dumps(result, indent=2), encoding="utf-8")

    print(f"JSON summary: {args.json_output}")
    return 0 if rows else 1


if __name__ == "__main__":
    raise SystemExit(main())
