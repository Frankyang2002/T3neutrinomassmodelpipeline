from __future__ import annotations

"""Diagnose the alpha=-1 (Y_F=0) scalar-mass validation failures.

This script reuses the parser/check definitions from
ValidateScalarMassGroupFactors.py and prints ONLY nonzero residuals.

Interpretation:
  alpha=-1 => Y_F=(alpha+1)/2=0.

The project has two distinct neutral-F branches:
  * odd-dimensional SU(2) F (dF=1,3): a self-conjugate/Majorana branch is
    allowed and is the branch used by the Matchete-side field definition;
  * even-dimensional SU(2) F (dF=2): the representation is pseudoreal, so the
    current Matchete-side T3 field definition does NOT use the self-conjugate
    branch.

RGBeta currently switches to its neutral/self-conjugate F construction whenever
Y_F=0, irrespective of dF.  Therefore:
  B,C (dF=1,3): alpha=-1 is a genuine branch change to diagnose separately.
  A,D,E (dF=2): alpha=-1 RGBeta is not the same field-content branch as the
                current Matchete T3 model and must not be used as an exact
                validation point for the vector-like formula.
"""

import argparse
from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from RGE.group_factors.ScalarMassGroupFactors import scalar_mass_group_factors
from RGE.group_factors.ValidateScalarMassGroupFactors import (
    MODELS,
    F,
    txt,
    coeff,
    load_betas,
    checks_for,
)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--root",
        type=Path,
        default=Path("output/full/hypercharge"),
    )
    args = parser.parse_args()

    print("alpha=-1 neutral-F scalar-mass diagnostic")
    print("Y_F = 0 for every row below.")
    print()

    mismatch_count = 0

    for model, (d1, d2, dF) in MODELS.items():
        path = (
            args.root
            / f"T3_{model}_alpha_m1"
            / "data"
            / "uv_rgbeta_rge.json"
        )
        if not path.exists():
            print(f"T3-{model}: missing {path}")
            continue

        pred = scalar_mass_group_factors(d1, d2, dF, -1)
        betas = load_betas(path)

        branch = (
            "physical self-conjugate/Majorana-allowed branch"
            if dF % 2 == 1
            else "RGBeta neutral branch differs from current Matchete field branch"
        )

        print(f"T3-{model}: dF={dF}  [{branch}]")

        model_mismatches = 0
        for beta_name, checks in checks_for(pred).items():
            beta = betas[beta_name]
            for label, (monomial, expected_text) in checks.items():
                got = coeff(beta, monomial)
                expected = F(expected_text)
                residual = got - expected

                if residual != 0:
                    model_mismatches += 1
                    mismatch_count += 1
                    print(
                        f"  {beta_name:<6} {label:<28} "
                        f"expected={txt(expected):>6}  "
                        f"rgbeta={txt(got):>6}  "
                        f"residual={txt(residual):>6}"
                    )

        if model_mismatches == 0:
            print("  no mismatches")
        print()

    print(f"TOTAL NONZERO RESIDUALS: {mismatch_count}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
