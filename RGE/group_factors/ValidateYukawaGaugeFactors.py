from __future__ import annotations

"""Validate analytic T3 Yukawa gauge group factors against saved RGBeta JSON.

The validator reads the existing hypercharge-comparison outputs under

    output/full/hypercharge/T3_<A-E>_alpha_<mN|pN>/data/uv_rgbeta_rge.json

and compares the gY^2 y_i and g2^2 y_i coefficients in beta_y1 and beta_y2
against

    A_Y  = -3 [Y_L^2 + Y_F^2]
         = -3/4 [1 + (alpha+1)^2],

    A_2  = -3 [C2(L) + C2(F)]
         = -3/4 (d_F^2 + 2).

It intentionally validates the saved RGBeta output, rather than comparing the
analytic formula to itself.
"""

import argparse
from dataclasses import dataclass, asdict
from fractions import Fraction
import json
from pathlib import Path
import re
import sys
from typing import Any

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from Reports.RGEComparison import _normalise_rgbeta_latex
from Reports.ReportGeneration import split_latex_terms


MODEL_DIMS = {
    "A": (1, 3, 2),
    "B": (2, 2, 1),
    "C": (2, 2, 3),
    "D": (3, 1, 2),
    "E": (3, 3, 2),
}


@dataclass(frozen=True)
class GaugeCheck:
    model: str
    alpha: int
    yukawa: str
    expected_u1: str
    rgbeta_u1: str | None
    u1_residual: str | None
    u1_match: bool
    expected_su2: str
    rgbeta_su2: str | None
    su2_residual: str | None
    su2_match: bool
    source: str


def _fraction_text(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def _alpha_from_dirname(name: str) -> int:
    match = re.fullmatch(r"T3_[A-E]_alpha_([mp])(\d+)", name)
    if not match:
        raise ValueError(f"Unrecognised hypercharge directory name: {name}")
    sign, magnitude = match.groups()
    value = int(magnitude)
    return -value if sign == "m" else value


def _model_from_dirname(name: str) -> str:
    match = re.fullmatch(r"T3_([A-E])_alpha_[mp]\d+", name)
    if not match:
        raise ValueError(f"Unrecognised hypercharge directory name: {name}")
    return match.group(1)


def _expected_u1(alpha: int) -> Fraction:
    # -3/4 [1 + (alpha+1)^2]
    return -Fraction(3, 4) * (1 + (alpha + 1) ** 2)


def _expected_su2(dF: int) -> Fraction:
    # C2(L)=3/4, C2(F)=(dF^2-1)/4
    return -Fraction(3, 4) * (dF * dF + 2)


def _compact_latex(text: str) -> str:
    text = text.replace(r"\,", "")
    text = text.replace(r"\!", "")
    text = text.replace(r"\left", "")
    text = text.replace(r"\right", "")
    text = text.replace(" ", "")
    return text


def _leading_rational(term: str) -> Fraction:
    """Extract the leading numerical coefficient of a normalised LaTeX term."""
    s = _compact_latex(term)

    sign = 1
    if s.startswith("+"):
        s = s[1:]
    elif s.startswith("-"):
        sign = -1
        s = s[1:]

    # \frac{n}{d}
    match = re.match(r"\\frac\{(\d+)\}\{(\d+)\}", s)
    if match:
        return sign * Fraction(int(match.group(1)), int(match.group(2)))

    # plain integer/rational integer prefix
    match = re.match(r"(\d+)", s)
    if match:
        return sign * Fraction(int(match.group(1)), 1)

    # no explicit magnitude -> +/-1
    return Fraction(sign, 1)


def _contains_gauge_yukawa(term: str, gauge: str, yukawa: str) -> bool:
    s = _compact_latex(term)

    if gauge == "gY":
        gauge_present = (
            r"g_Y" in s
            or r"g_{Y}" in s
            or "gY" in s
        )
    elif gauge == "g2":
        gauge_present = (
            r"g_2" in s
            or r"g_{2}" in s
            or "g2" in s
        )
    else:
        raise ValueError(gauge)

    y_index = "1" if yukawa == "y1" else "2"
    yukawa_present = (
        rf"y_{y_index}" in s
        or rf"y_{{{y_index}}}" in s
        or yukawa in s
    )

    # Require a squared gauge coupling. RGBeta report normalisation may render
    # this as ^2 or ^{2}.
    squared = "^2" in s or "^{2}" in s

    return gauge_present and yukawa_present and squared


def _extract_from_latex(payload: dict[str, Any], yukawa: str, gauge: str) -> Fraction:
    latex_map = payload.get("report_beta_latex", {}) or {}
    raw = latex_map.get(yukawa)
    if not raw:
        raise KeyError(f"Missing report_beta_latex[{yukawa!r}]")

    cleaned = _normalise_rgbeta_latex(str(raw))
    terms = split_latex_terms(cleaned) or [cleaned]

    candidates = [
        term for term in terms
        if _contains_gauge_yukawa(term, gauge, yukawa)
    ]

    if len(candidates) != 1:
        diagnostic = "\n".join(f"  {item}" for item in candidates) or "  (none)"
        raise ValueError(
            f"Expected exactly one {gauge}^2 {yukawa} term, found "
            f"{len(candidates)}:\n{diagnostic}\n"
            f"Full normalised beta:\n{cleaned}"
        )

    return _leading_rational(candidates[0])


def _load_payload(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8-sig"))
    if payload.get("status") != "Success":
        raise ValueError(f"RGBeta status is {payload.get('status')!r}: {path}")
    return payload


def validate_one(path: Path) -> list[GaugeCheck]:
    run_dir = path.parents[1]
    model = _model_from_dirname(run_dir.name)
    alpha = _alpha_from_dirname(run_dir.name)
    _dS1, _dS2, dF = MODEL_DIMS[model]

    payload = _load_payload(path)

    expected_u1 = _expected_u1(alpha)
    expected_su2 = _expected_su2(dF)

    rows: list[GaugeCheck] = []
    for yukawa in ("y1", "y2"):
        try:
            got_u1 = _extract_from_latex(payload, yukawa, "gY")
            u1_error = None
        except Exception as exc:
            got_u1 = None
            u1_error = str(exc)

        try:
            got_su2 = _extract_from_latex(payload, yukawa, "g2")
            su2_error = None
        except Exception as exc:
            got_su2 = None
            su2_error = str(exc)

        u1_residual = None if got_u1 is None else got_u1 - expected_u1
        su2_residual = None if got_su2 is None else got_su2 - expected_su2

        row = GaugeCheck(
            model=model,
            alpha=alpha,
            yukawa=yukawa,
            expected_u1=_fraction_text(expected_u1),
            rgbeta_u1=None if got_u1 is None else _fraction_text(got_u1),
            u1_residual=None if u1_residual is None else _fraction_text(u1_residual),
            u1_match=(u1_residual == 0),
            expected_su2=_fraction_text(expected_su2),
            rgbeta_su2=None if got_su2 is None else _fraction_text(got_su2),
            su2_residual=None if su2_residual is None else _fraction_text(su2_residual),
            su2_match=(su2_residual == 0),
            source=str(path),
        )
        rows.append(row)

        if u1_error:
            print(f"[parse warning] {run_dir.name} {yukawa} U(1): {u1_error}")
        if su2_error:
            print(f"[parse warning] {run_dir.name} {yukawa} SU(2): {su2_error}")

    return rows


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Validate T3 y1/y2 gauge group factors against saved RGBeta results."
    )
    parser.add_argument(
        "--root",
        type=Path,
        default=Path("output/full/hypercharge"),
        help="Hypercharge-comparison output root.",
    )
    parser.add_argument(
        "--json-output",
        type=Path,
        default=Path("output/group_factors/yukawa_gauge_validation.json"),
    )
    args = parser.parse_args()

    paths = sorted(args.root.glob("T3_*_alpha_*/data/uv_rgbeta_rge.json"))
    if not paths:
        raise FileNotFoundError(
            f"No uv_rgbeta_rge.json files found below {args.root}"
        )

    rows: list[GaugeCheck] = []
    for path in paths:
        rows.extend(validate_one(path))

    rows.sort(key=lambda item: (item.model, item.alpha, item.yukawa))

    print(
        f"{'model':<6} {'alpha':>5} {'y':>3} "
        f"{'U1 exp':>8} {'U1 RGB':>8} {'dU1':>6} "
        f"{'SU2 exp':>8} {'SU2 RGB':>8} {'dSU2':>6} status"
    )
    print("-" * 86)

    for row in rows:
        status = "PASS" if row.u1_match and row.su2_match else "FAIL"
        print(
            f"T3-{row.model:<3} {row.alpha:>5} {row.yukawa:>3} "
            f"{row.expected_u1:>8} {str(row.rgbeta_u1):>8} "
            f"{str(row.u1_residual):>6} "
            f"{row.expected_su2:>8} {str(row.rgbeta_su2):>8} "
            f"{str(row.su2_residual):>6} {status}"
        )

    payload = {
        "status": (
            "Success"
            if all(row.u1_match and row.su2_match for row in rows)
            else "Failed"
        ),
        "file_count": len(paths),
        "check_count": len(rows),
        "u1_match_count": sum(row.u1_match for row in rows),
        "su2_match_count": sum(row.su2_match for row in rows),
        "all_u1_match": all(row.u1_match for row in rows),
        "all_su2_match": all(row.su2_match for row in rows),
        "formulae": {
            "u1": "-3*(YL^2+YF^2) = -3/4*(1+(alpha+1)^2)",
            "su2": "-3*(C2(L)+C2(F)) = -3/4*(dF^2+2)",
        },
        "checks": [asdict(row) for row in rows],
    }

    args.json_output.parent.mkdir(parents=True, exist_ok=True)
    args.json_output.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    print()
    print(f"Validated RGBeta files: {len(paths)}")
    print(f"Yukawa gauge checks:     {len(rows)}")
    print(f"U(1) exact matches:      {payload['u1_match_count']}/{len(rows)}")
    print(f"SU(2) exact matches:     {payload['su2_match_count']}/{len(rows)}")
    print(f"JSON summary:            {args.json_output}")
    print(f"OVERALL:                 {payload['status']}")

    return 0 if payload["status"] == "Success" else 1


if __name__ == "__main__":
    raise SystemExit(main())
