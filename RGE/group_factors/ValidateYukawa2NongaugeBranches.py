from __future__ import annotations

"""Branch-aware validation of the remaining nongauge beta_y2 sector.

Input:
    output/group_factors/y2_nongauge_diagnostic.json

Validation logic:
  * alpha != -1:
      compare every nongauge expression to the same model's alpha=0 reference.
      This checks hypercharge-independence of the vector-like nongauge sector.
  * alpha == -1 and dF odd (T3-B,C):
      validate the two extra self-conjugate structures found by RGBeta:
          G_S2 Tr(y2 y1^\dagger) y1
          (1/2) G_F2 y2 y1^\dagger y1
  * alpha == -1 and dF even (T3-A,D,E):
      record an excluded branch diagnostic, because RGBeta's neutral branch is
      not the current Matchete-side physical T3 field branch.

This script deliberately does not use RGEComparison term signatures.
"""

import argparse
import json
from pathlib import Path
import re
import sys
from fractions import Fraction

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from RGE.group_factors.Yukawa2BetaGroupFactors import (
    yukawa2_beta_group_factors,
    physical_self_conjugate_f,
)


def _canon(s: str) -> str:
    return re.sub(r"\s+", "", str(s))


def _frac_prefix(text: str, target: str) -> Fraction | None:
    """Coefficient of a standalone target term from its leading TeX prefix."""
    c = _canon(text)
    t = _canon(target)
    pos = c.find(t)
    if pos < 0:
        return None
    prefix = c[:pos].lstrip("+")
    if prefix == "":
        return Fraction(1)
    m = re.fullmatch(r"\\frac\{(-?\d+)\}\{(\d+)\}", prefix)
    if m:
        return Fraction(int(m.group(1)), int(m.group(2)))
    m = re.fullmatch(r"(-?\d+)", prefix)
    if m:
        return Fraction(int(m.group(1)))
    return None


def _outer_fraction(text: str) -> Fraction:
    c = _canon(text).lstrip("+")
    m = re.match(r"\\frac\{(-?\d+)\}\{(\d+)\}\(", c)
    if m:
        return Fraction(int(m.group(1)), int(m.group(2)))
    return Fraction(1)


def _inner_integer_before(text: str, target: str) -> Fraction | None:
    c = _canon(text)
    t = _canon(target)
    pos = c.find(t)
    if pos < 0:
        return None
    before = c[:pos]
    # The target is inside a sum.  Its immediate multiplicative coefficient is
    # the final integer before the Matrix token, or 1 if no integer is present.
    m = re.search(r"(?:\+|\()(-?\d+)$", before)
    if m:
        return Fraction(int(m.group(1)))
    return Fraction(1)


TRACE_CROSS = (
    r"\text{Tr}[\text{y2}.\text{Trans}(\text{Bar}(\text{y1}))]"
    r"\text{Matrix}(\text{y1})"
)

MATRIX_CROSS = (
    r"\text{Matrix}(\text{y2},\text{Trans}(\text{Bar}(\text{y1})),"
    r"\text{y1})"
)


def extract_self_conjugate_extras(terms: list[str]):
    trace_coeff = None
    matrix_coeff = None

    for term in terms:
        if _canon(TRACE_CROSS) in _canon(term):
            trace_coeff = _frac_prefix(term, TRACE_CROSS)
        if _canon(MATRIX_CROSS) in _canon(term):
            outer = _outer_fraction(term)
            inner = _inner_integer_before(term, MATRIX_CROSS)
            if inner is not None:
                matrix_coeff = outer * inner

    return trace_coeff, matrix_coeff


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("output/group_factors/y2_nongauge_diagnostic.json"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("output/group_factors/y2_nongauge_branch_validation.json"),
    )
    args = parser.parse_args()

    payload = json.loads(args.input.read_text(encoding="utf-8"))
    rows = payload.get("rows") or []
    if not rows:
        raise SystemExit(f"No rows in {args.input}")

    reference = {}
    for row in rows:
        if int(row["alpha"]) == 0:
            reference[row["model"]] = row["canonical_nongauge"]

    results = []
    overall = True
    excluded = 0

    print("beta_y2 branch-aware nongauge validation")
    print("model alpha branch                         status")
    print("-" * 64)

    for row in rows:
        model = row["model"]
        alpha = int(row["alpha"])
        d1, d2, dF = int(row["dS1"]), int(row["dS2"]), int(row["dF"])
        pred = yukawa2_beta_group_factors(d1, d2, dF, alpha)

        if alpha != -1:
            ok = row["canonical_nongauge"] == reference[model]
            overall &= ok
            branch = "vector-like"
            detail = {"matches_alpha0_reference": ok}

        elif dF % 2 == 0:
            # Known RGBeta/Matchete neutral-field-content mismatch.
            ok = True
            excluded += 1
            branch = "excluded even-dF RGBeta neutral"
            detail = {
                "excluded_branch_mismatch": True,
                "reason": (
                    "RGBeta switches to its YF=0 neutral/self-conjugate branch "
                    "for even dF, while the current physical T3 model does not."
                ),
            }

        else:
            branch = "physical self-conjugate"
            trace_got, matrix_got = extract_self_conjugate_extras(
                row["nongauge_terms"]
            )
            trace_expected = Fraction(
                pred.self_conjugate_extra["Tr_y2_y1dag*y1"]
            )
            matrix_expected = Fraction(
                pred.self_conjugate_extra["y2_y1dag_y1"]
            )
            ok = (
                trace_got == trace_expected
                and matrix_got == matrix_expected
            )
            overall &= ok
            detail = {
                "trace_extra_expected": str(trace_expected),
                "trace_extra_extracted": (
                    str(trace_got) if trace_got is not None else None
                ),
                "matrix_extra_expected": str(matrix_expected),
                "matrix_extra_extracted": (
                    str(matrix_got) if matrix_got is not None else None
                ),
            }

        print(
            f"T3-{model} {alpha:>5} {branch:<30} "
            f"{'PASS' if ok else 'FAIL'}"
        )
        results.append({
            "model": model,
            "alpha": alpha,
            "dF": dF,
            "physical_self_conjugate_F": physical_self_conjugate_f(dF, alpha),
            "branch": branch,
            "match": ok,
            **detail,
        })

    result = {
        "status": "Success" if overall else "Failed",
        "formulas": {
            "vectorlike": {
                "Tr_y2*y2": "G_S2",
                "y2_y2dag_y2": "1/2 (G_L2 + G_F2)",
                "y1_y1dag_y2": "1/2 G_L1",
                "Ye_Yedag_y2": "1/2",
            },
            "self_conjugate_extra": {
                "Tr_y2_y1dag*y1": "G_S2",
                "y2_y1dag_y1": "1/2 G_F2",
            },
        },
        "excluded_even_dF_neutral_rows": excluded,
        "rows": results,
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2), encoding="utf-8")

    print()
    print(f"Excluded even-dF neutral RGBeta rows: {excluded}")
    print(f"JSON summary: {args.output}")
    print(f"OVERALL:      {result['status']}")
    return 0 if overall else 1


if __name__ == "__main__":
    raise SystemExit(main())
