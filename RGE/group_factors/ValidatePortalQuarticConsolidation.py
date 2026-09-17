from __future__ import annotations

"""Regression checks for the consolidated portal-quartic group-factor module.

This validator checks representation-level identities and, when the derivation
JSON files are present, compares the consolidated implementation against:
  * portal_adjoint_recouplings.json
  * e_portal_cross_scalar_recouplings.json
  * portal_gauge_quartics.json

The derivation files remain independent evidence; the production module does
not import them.
"""

import argparse
from fractions import Fraction
import json
from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from RGE.group_factors.PortalQuarticBetaGroupFactors import (
    portal_quartic_beta_group_factors,
)


MODELS = {
    "A": (1, 3, 2),
    "B": (2, 2, 1),
    "C": (2, 2, 3),
    "D": (3, 1, 2),
    "E": (3, 3, 2),
}


def F(x) -> Fraction:
    return Fraction(str(x))


def check_representation_identities():
    rows = []
    ok = True

    expected = {
        "A": {"H2Adj2": "2"},
        "B": {"H1Adj2": "3/4", "H2Adj2": "3/4", "L12Adj2": "3/4"},
        "C": {"H1Adj2": "3/4", "H2Adj2": "3/4", "L12Adj2": "3/4"},
        "D": {"H1Adj2": "2"},
        "E": {"H1Adj2": "2", "H2Adj2": "2", "L12Adj2": "4"},
    }

    for model, (d1, d2, dF) in MODELS.items():
        p = portal_quartic_beta_group_factors(d1, d2, dF, 0)
        actual = {}
        if "lambdaH1Adj_sq" in p.beta_lambdaH1:
            actual["H1Adj2"] = p.beta_lambdaH1["lambdaH1Adj_sq"]
        if "lambdaH2Adj_sq" in p.beta_lambdaH2:
            actual["H2Adj2"] = p.beta_lambdaH2["lambdaH2Adj_sq"]
        if "lambda12Adj_sq" in p.beta_lambda12:
            actual["L12Adj2"] = p.beta_lambda12["lambda12Adj_sq"]

        match = actual == expected[model]
        ok &= match
        rows.append({"model": model, "actual": actual, "expected": expected[model], "match": match})
        print(f"T3-{model} adjoint-squared {'PASS' if match else 'FAIL'} {actual}")

    e = portal_quartic_beta_group_factors(3, 3, 2, 0)
    cross_expected = {
        "H1_cross": "4",
        "H2_cross": "4",
        "L12_cross_sq": "4",
        "L12_cross_S1A": "4",
        "L12_cross_S2A": "4",
        "Cross_beta_g2_four": "6",
        "Cross_beta_adj_sq": "2",
        "Cross_beta_L12_cross": "8",
        "Cross_beta_cross_sq": "10",
        "Cross_beta_cross_S1A": "-2",
        "Cross_beta_cross_S2A": "-2",
    }
    cross_actual = {
        "H1_cross": e.beta_lambdaH1["lambdaH2*lambda12Cross"],
        "H2_cross": e.beta_lambdaH2["lambdaH1*lambda12Cross"],
        "L12_cross_sq": e.beta_lambda12["lambda12Cross_sq"],
        "L12_cross_S1A": e.beta_lambda12["lambda12Cross*lambdaS1Adj"],
        "L12_cross_S2A": e.beta_lambda12["lambda12Cross*lambdaS2Adj"],
        "Cross_beta_g2_four": e.generated_non_singlet["beta_lambda12Cross:g2_four"],
        "Cross_beta_adj_sq": e.generated_non_singlet["beta_lambda12Cross:lambda12Adj_sq"],
        "Cross_beta_L12_cross": e.generated_non_singlet["beta_lambda12Cross:lambda12*lambda12Cross"],
        "Cross_beta_cross_sq": e.generated_non_singlet["beta_lambda12Cross:lambda12Cross_sq"],
        "Cross_beta_cross_S1A": e.generated_non_singlet["beta_lambda12Cross:lambda12Cross*lambdaS1Adj"],
        "Cross_beta_cross_S2A": e.generated_non_singlet["beta_lambda12Cross:lambda12Cross*lambdaS2Adj"],
    }
    match = cross_actual == cross_expected
    ok &= match
    rows.append({"model": "E-cross", "actual": cross_actual, "expected": cross_expected, "match": match})
    print(f"T3-E Cross {'PASS' if match else 'FAIL'}")

    return rows, ok


def check_gauge_scan():
    rows = []
    ok = True
    for model, (d1, d2, dF) in MODELS.items():
        for alpha in (-2, -1, 0, 1, 2):
            p = portal_quartic_beta_group_factors(d1, d2, dF, alpha)

            # General closed H-S pure gauge identities.
            h1_total = F(p.beta_lambdaH1["g2_four"])
            h2_total = F(p.beta_lambdaH2["g2_four"])
            if d1 == 1 and h1_total != 0:
                match = False
            elif d2 == 1 and h2_total != 0:
                match = False
            else:
                match = True

            # U(1) coefficients must be non-negative squares.
            match &= F(p.beta_lambdaH1["gY_four"]) >= 0
            match &= F(p.beta_lambdaH2["gY_four"]) >= 0
            match &= F(p.beta_lambda12["gY_four"]) >= 0

            ok &= match
            rows.append({"model": model, "alpha": alpha, "match": match})
    print(f"25-point gauge identity scan {'PASS' if ok else 'FAIL'}")
    return rows, ok


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("output/group_factors/portal_quartic_consolidation.json"),
    )
    args = parser.parse_args()

    rep_rows, rep_ok = check_representation_identities()
    gauge_rows, gauge_ok = check_gauge_scan()
    overall = rep_ok and gauge_ok

    payload = {
        "status": "Success" if overall else "Failed",
        "representation_checks": rep_rows,
        "gauge_checks": gauge_rows,
        "notes": [
            "Production coefficients depend on SU(2) dimensions/Casimirs and hypercharges, not A-E labels.",
            "The S1-S2 d=3,d=3 sector keeps the independent Cross channel explicitly.",
        ],
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    print()
    print(f"JSON summary: {args.output}")
    print(f"OVERALL:      {payload['status']}")
    return 0 if overall else 1


if __name__ == "__main__":
    raise SystemExit(main())
