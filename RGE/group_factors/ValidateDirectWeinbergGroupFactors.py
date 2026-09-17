from __future__ import annotations

"""Validate the analytic direct-Weinberg SU(2) recoupling formula.

Input:
    output/group_factors/direct_weinberg_group_factor_diagnostic.json

The diagnostic was extracted from the authoritative EFT1 full-flavor
C12 -> Weinberg running object.  This validator compares its reduced
prefactor/lambdaT3 against the Wigner-6j formula.
"""

import argparse
import json
from pathlib import Path
import sys

import sympy as sp

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from RGE.group_factors.DirectWeinbergGroupFactors import (
    direct_weinberg_group_factor,
)


def E(text):
    return sp.sympify(
        str(text),
        locals={"sqrt": sp.sqrt},
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        type=Path,
        default=Path(
            "output/group_factors/direct_weinberg_group_factor_diagnostic.json"
        ),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "output/group_factors/direct_weinberg_group_factor_validation.json"
        ),
    )
    args = parser.parse_args()

    payload = json.loads(args.input.read_text(encoding="utf-8"))
    rows = payload.get("rows") or []
    if not rows:
        raise SystemExit(f"No diagnostic rows in {args.input}")

    out_rows = []
    overall = True

    print("Direct C12 -> Weinberg analytic group-factor validation")
    print("model alpha  predicted             extracted             status")
    print("-" * 74)

    for row in rows:
        pred = direct_weinberg_group_factor(
            int(row["dS1"]),
            int(row["dS2"]),
            int(row["dF"]),
        )
        predicted = sp.factor(E(pred.reduced_factor))
        extracted = sp.factor(E(row["prefactor_over_lambdaT3"]))
        residual = sp.simplify(predicted - extracted)
        ok = residual == 0
        overall &= ok

        print(
            f"T3-{row['model']} {int(row['alpha']):>5}  "
            f"{sp.sstr(predicted):<21} "
            f"{sp.sstr(extracted):<21} "
            f"{'PASS' if ok else 'FAIL'}"
        )

        out_rows.append({
            "model": row["model"],
            "alpha": row["alpha"],
            "dS1": row["dS1"],
            "dS2": row["dS2"],
            "dF": row["dF"],
            "wigner6j": pred.wigner6j,
            "normalization": pred.normalization,
            "cg_phase": pred.cg_phase,
            "predicted": sp.sstr(predicted),
            "extracted": sp.sstr(extracted),
            "residual": sp.sstr(residual),
            "match": ok,
        })

    result = {
        "status": "Success" if overall else "Failed",
        "formula": (
            "(4/3)*etaF*sqrt(3*dS1*dS2*dF)"
            "*Wigner6j(1/2,1/2,1;jS2,jS1,jF)"
        ),
        "phase_convention": {
            "dF=1": 1,
            "dF=2": 1,
            "dF=3": -1,
        },
        "rows": out_rows,
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2), encoding="utf-8")

    print()
    print(f"JSON summary: {args.output}")
    print(f"OVERALL:      {result['status']}")
    return 0 if overall else 1


if __name__ == "__main__":
    raise SystemExit(main())
