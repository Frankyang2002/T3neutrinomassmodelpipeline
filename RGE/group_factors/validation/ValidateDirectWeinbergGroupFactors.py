from __future__ import annotations

"""Extract and validate the direct EFT1 C12 -> Weinberg SU(2) group factor.

This file combines the former extraction diagnostic and analytic validator.

Authoritative fixed-order counting:
    tree LLSS boundary               O(hbar^0)
    direct LLSS -> O5 running        O(hbar^1)
    LLSS self-running                O(hbar^1)
    self-running + scalar matching   O(hbar^2), excluded here

For every successful T3 run under ``--root`` we read

    data/eft1_wilson_flavor_at_S_threshold.json

and extract the stored flavor-blind prefactor in

    16*pi^2 beta_kappa = prefactor * C12[p,q].

After dividing out lambdaT3, the extracted SU(2) recoupling factor is compared
with the analytic Wigner-6j expression defined directly in this validator.

Important
---------
The direct C12 -> C5 Wigner-6j phase/normalization remains a known physics
issue in this project. This validator preserves the existing comparison; it
does not claim to resolve that issue or change the formula.
"""

import argparse
from dataclasses import dataclass
import json
from pathlib import Path
import re
from typing import Any

import sympy as sp
from sympy.physics.wigner import wigner_6j


MODEL_DIMS = {
    "A": (1, 3, 2),
    "B": (2, 2, 1),
    "C": (2, 2, 3),
    "D": (3, 1, 2),
    "E": (3, 3, 2),
}


def _validate_t3_dimensions(d1: int, d2: int, d_f: int) -> None:
    if any(d not in (1, 2, 3) for d in (d1, d2, d_f)):
        raise ValueError("Current T3 implementation supports d in {1,2,3}.")
    if abs(d1 - d_f) != 1 or abs(d2 - d_f) != 1:
        raise ValueError("T3 Yukawa invariance requires dSi=dF+/-1.")


def _j(dimension: int) -> sp.Rational:
    return sp.Rational(dimension - 1, 2)


def project_cg_phase(d_f: int) -> int:
    """Project CG convention phase for the currently supported dF<=3 basis."""
    d_f = int(d_f)
    if d_f in (1, 2):
        return 1
    if d_f == 3:
        return -1
    raise ValueError(
        "Unsupported dF for current project CG phase convention."
    )


@dataclass(frozen=True)
class DirectWeinbergGroupFactor:
    dS1: int
    dS2: int
    dF: int
    jS1: str
    jS2: str
    jF: str
    wigner6j: str
    normalization: str
    cg_phase: int
    reduced_factor: str


def direct_weinberg_group_factor(
    dS1: int,
    dS2: int,
    dF: int,
) -> DirectWeinbergGroupFactor:
    """Return the preserved analytic C12 -> C5 Wigner-6j prediction."""

    d1, d2, d_f = map(int, (dS1, dS2, dF))
    _validate_t3_dimensions(d1, d2, d_f)

    j1, j2, jf = _j(d1), _j(d2), _j(d_f)
    sixj = sp.simplify(
        wigner_6j(
            sp.Rational(1, 2),
            sp.Rational(1, 2),
            1,
            j2,
            j1,
            jf,
        )
    )

    normalization = sp.sqrt(3 * d1 * d2 * d_f)
    phase = project_cg_phase(d_f)
    reduced = sp.factor(
        sp.Rational(4, 3)
        * phase
        * normalization
        * sixj
    )

    return DirectWeinbergGroupFactor(
        dS1=d1,
        dS2=d2,
        dF=d_f,
        jS1=sp.sstr(j1),
        jS2=sp.sstr(j2),
        jF=sp.sstr(jf),
        wigner6j=sp.sstr(sixj),
        normalization=sp.sstr(normalization),
        cg_phase=phase,
        reduced_factor=sp.sstr(reduced),
    )


def model_alpha(path: Path) -> tuple[str, int]:
    name = path.parents[1].name
    match = re.fullmatch(r"T3_([A-E])_alpha_([mp])(\d+)", name)
    if not match:
        raise ValueError(name)

    model, sign, magnitude = match.groups()
    alpha = int(magnitude)
    return model, (-alpha if sign == "m" else alpha)


def parse_expression(raw: Any) -> sp.Expr:
    return sp.sympify(
        str(raw).replace("^", "**"),
        locals={
            "lambdaT3": sp.Symbol("lambdaT3"),
            "sqrt": sp.sqrt,
            "log": sp.log,
            "conjugate": sp.conjugate,
        },
    )


def parse_reduced_factor(raw: Any) -> sp.Expr:
    return sp.sympify(
        str(raw),
        locals={"sqrt": sp.sqrt},
    )


def extract_rows(root: Path) -> tuple[list[dict[str, Any]], bool, dict[str, list[str]]]:
    paths = sorted(
        root.glob(
            "T3_*_alpha_*/data/eft1_wilson_flavor_at_S_threshold.json"
        )
    )
    if not paths:
        raise FileNotFoundError(
            "No eft1_wilson_flavor_at_S_threshold.json files found under "
            f"{root}"
        )

    lambda_t3 = sp.Symbol("lambdaT3")
    rows: list[dict[str, Any]] = []

    for path in paths:
        payload = json.loads(path.read_text(encoding="utf-8-sig"))
        if payload.get("status") != "Success":
            continue

        model, alpha = model_alpha(path)
        d_s1, d_s2, d_f = MODEL_DIMS[model]

        direct = payload.get("direct_weinberg_beta") or {}
        raw_prefactor = direct.get("flavor_blind_prefactor")
        if raw_prefactor is None:
            continue

        prefactor = sp.factor(parse_expression(raw_prefactor))
        reduced = (
            sp.factor(prefactor / lambda_t3)
            if prefactor.has(lambda_t3)
            else None
        )

        rows.append(
            {
                "model": model,
                "alpha": alpha,
                "dS1": d_s1,
                "dS2": d_s2,
                "dF": d_f,
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
            }
        )

    grouped: dict[str, set[str | None]] = {}
    for row in rows:
        grouped.setdefault(str(row["model"]), set()).add(
            row["prefactor_over_lambdaT3"]
        )

    alpha_independent = all(len(values) == 1 for values in grouped.values())
    distinct = {
        model: sorted(
            str(value) if value is not None else "<no lambdaT3>"
            for value in values
        )
        for model, values in sorted(grouped.items())
    }

    return rows, alpha_independent, distinct


def load_rows_from_diagnostic(
    path: Path,
) -> tuple[list[dict[str, Any]], bool | None, dict[str, list[str]] | None]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    rows = payload.get("rows") or []
    if not rows:
        raise ValueError(f"No diagnostic rows in {path}")

    return (
        rows,
        payload.get("alpha_independent"),
        payload.get("distinct_reduced_factors"),
    )


def validate_rows(rows: list[dict[str, Any]]) -> tuple[list[dict[str, Any]], bool]:
    out_rows: list[dict[str, Any]] = []
    overall = True

    for row in rows:
        pred = direct_weinberg_group_factor(
            int(row["dS1"]),
            int(row["dS2"]),
            int(row["dF"]),
        )

        reduced_text = row.get("prefactor_over_lambdaT3")
        if reduced_text is None:
            predicted = None
            extracted = None
            residual = None
            ok = False
        else:
            predicted = sp.factor(parse_reduced_factor(pred.reduced_factor))
            extracted = sp.factor(parse_reduced_factor(reduced_text))
            residual = sp.simplify(predicted - extracted)
            ok = residual == 0

        overall &= ok

        out_rows.append(
            {
                "model": row["model"],
                "alpha": row["alpha"],
                "dS1": row["dS1"],
                "dS2": row["dS2"],
                "dF": row["dF"],
                "source": row.get("source"),
                "prefactor": row.get("prefactor"),
                "prefactor_over_lambdaT3": reduced_text,
                "one_generation_beta_kappa_16pi2": row.get(
                    "one_generation_beta_kappa_16pi2"
                ),
                "canonical_c12_one_generation": row.get(
                    "canonical_c12_one_generation"
                ),
                "wigner6j": pred.wigner6j,
                "normalization": pred.normalization,
                "cg_phase": pred.cg_phase,
                "predicted": (
                    sp.sstr(predicted) if predicted is not None else None
                ),
                "extracted": (
                    sp.sstr(extracted) if extracted is not None else None
                ),
                "residual": (
                    sp.sstr(residual) if residual is not None else None
                ),
                "match": ok,
            }
        )

    return out_rows, overall


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--root",
        type=Path,
        default=Path("output/full/hypercharge"),
        help=(
            "Root containing T3_*_alpha_*/data/"
            "eft1_wilson_flavor_at_S_threshold.json."
        ),
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=None,
        help=(
            "Optional legacy diagnostic JSON. If supplied, validate its rows "
            "instead of extracting directly from --root."
        ),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "output/group_factors/direct_weinberg_group_factor_validation.json"
        ),
    )
    parser.add_argument(
        "--diagnostic-output",
        type=Path,
        default=None,
        help=(
            "Optional path for writing the former extraction-only diagnostic "
            "payload as a compatibility/debug artifact."
        ),
    )
    args = parser.parse_args()

    if args.input is not None:
        rows, alpha_independent, distinct = load_rows_from_diagnostic(args.input)
        source_mode = "legacy_diagnostic_json"
    else:
        try:
            rows, alpha_independent, distinct = extract_rows(args.root)
        except FileNotFoundError as exc:
            raise SystemExit(str(exc)) from exc
        source_mode = "direct_extraction"

    if not rows:
        raise SystemExit("No successful direct-Weinberg rows were extracted.")

    print("Direct EFT1 C12 -> Weinberg group-factor validation")
    print()
    print("Extracted prefactors:")
    print(
        "model alpha dS1 dS2 dF  "
        "prefactor                     prefactor/lambdaT3"
    )
    print("-" * 92)

    for row in rows:
        print(
            f"T3-{row['model']} {int(row['alpha']):>5} "
            f"{int(row['dS1']):>3} {int(row['dS2']):>3} "
            f"{int(row['dF']):>2}  "
            f"{str(row.get('prefactor')):<28} "
            f"{str(row.get('prefactor_over_lambdaT3'))}"
        )

    print()
    print("Distinct reduced factors by model:")
    if distinct is not None:
        for model in sorted(distinct):
            print(f"  T3-{model}: {distinct[model]}")
    print(f"alpha-independent: {alpha_independent}")

    if args.diagnostic_output is not None:
        diagnostic = {
            "status": "Success",
            "alpha_independent": alpha_independent,
            "distinct_reduced_factors": distinct,
            "rows": rows,
        }
        args.diagnostic_output.parent.mkdir(parents=True, exist_ok=True)
        args.diagnostic_output.write_text(
            json.dumps(diagnostic, indent=2),
            encoding="utf-8",
        )
        print(f"Diagnostic JSON: {args.diagnostic_output}")

    out_rows, overall = validate_rows(rows)

    print()
    print("Analytic Wigner-6j comparison:")
    print("model alpha  predicted             extracted             status")
    print("-" * 74)

    for row in out_rows:
        print(
            f"T3-{row['model']} {int(row['alpha']):>5}  "
            f"{str(row['predicted']):<21} "
            f"{str(row['extracted']):<21} "
            f"{'PASS' if row['match'] else 'FAIL'}"
        )

    result = {
        "status": "Success" if overall else "Failed",
        "source_mode": source_mode,
        "root": str(args.root) if args.input is None else None,
        "input": str(args.input) if args.input is not None else None,
        "alpha_independent": alpha_independent,
        "distinct_reduced_factors": distinct,
        "formula": (
            "(4/3)*etaF*sqrt(3*dS1*dS2*dF)"
            "*Wigner6j(1/2,1/2,1;jS2,jS1,jF)"
        ),
        "phase_convention": {
            "dF=1": 1,
            "dF=2": 1,
            "dF=3": -1,
        },
        "known_issue": (
            "Direct C12->C5 Wigner-6j phase/normalization remains an "
            "unresolved project physics issue; this file only preserves the "
            "existing extraction and comparison."
        ),
        "rows": out_rows,
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2), encoding="utf-8")

    print()
    print(f"Validation JSON: {args.output}")
    print(f"OVERALL:         {result['status']}")
    return 0 if overall else 1


if __name__ == "__main__":
    raise SystemExit(main())
