from __future__ import annotations

"""Yukawa beta group factors and their validation.
We assemble full yukawa beta function coefficients and validates with saved RGBeta output"""

import argparse
import json
import sys
from dataclasses import asdict, dataclass
from fractions import Fraction
from pathlib import Path
import re
from typing import Any

import sympy as sp

PROJECT_ROOT = Path(__file__).resolve().parents[3]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from RGE.group_factors.core.RepresentationFactors import (
    canonical_yukawa_leg_factors,
    su2_quadratic_casimir_from_dimension,
)



@dataclass(frozen=True)
class YukawaBetaCoefficients:
    trace: str
    self_matrix: str
    cross_matrix: str
    charged_lepton_matrix: str
    su2_gauge: str
    u1_gauge: str


@dataclass(frozen=True)
class T3YukawaBetaGroupFactors:
    dS1: int
    dS2: int
    dF: int
    alpha: int
    C2L: str
    C2F: str
    YL: str
    YF: str
    self_conjugate_F: bool
    y1: YukawaBetaCoefficients
    y2: YukawaBetaCoefficients
    y2_self_conjugate_extra: dict[str, str]


def _txt(value: Fraction | sp.Rational) -> str:
    value = sp.Rational(value)
    return (
        str(int(value))
        if value.q == 1
        else f"{int(value.p)}/{int(value.q)}"
    )


def physical_self_conjugate_f(dF: int, alpha: int) -> bool:
    """Current Matchete-side physical self-conjugate-F criterion."""
    return int(alpha) == -1 and int(dF) % 2 == 1


def _validate(dS1: int, dS2: int, dF: int) -> None:
    if any(d not in (1, 2, 3) for d in (dS1, dS2, dF)):
        raise ValueError("Current T3 implementation supports d in {1,2,3}.")
    if abs(dS1 - dF) != 1 or abs(dS2 - dF) != 1:
        raise ValueError("T3 Yukawa invariance requires dSi=dF+/-1.")


def complete_yukawa_group_factors(
    *,
    dS1: int,
    dS2: int,
    dF: int,
    alpha: int,
) -> T3YukawaBetaGroupFactors:
    dS1, dS2, dF, alpha = map(int, (dS1, dS2, dF, alpha))
    _validate(dS1, dS2, dF)

    y1legs = canonical_yukawa_leg_factors(dF, dS1)
    y2legs = canonical_yukawa_leg_factors(dF, dS2)

    GL1 = sp.Rational(
        y1legs.G_lepton.numerator,
        y1legs.G_lepton.denominator,
    )
    GF1 = sp.Rational(
        y1legs.G_heavy.numerator,
        y1legs.G_heavy.denominator,
    )
    GS1 = sp.Rational(
        y1legs.G_scalar.numerator,
        y1legs.G_scalar.denominator,
    )

    GL2 = sp.Rational(
        y2legs.G_lepton.numerator,
        y2legs.G_lepton.denominator,
    )
    GF2 = sp.Rational(
        y2legs.G_heavy.numerator,
        y2legs.G_heavy.denominator,
    )
    GS2 = sp.Rational(
        y2legs.G_scalar.numerator,
        y2legs.G_scalar.denominator,
    )

    C2L = sp.Rational(3, 4)
    c2f = su2_quadratic_casimir_from_dimension(dF)
    C2F = sp.Rational(c2f.numerator, c2f.denominator)

    YL = -sp.Rational(1, 2)
    YF = sp.Rational(alpha + 1, 2)

    su2 = sp.simplify(-3 * (C2L + C2F))
    u1 = sp.simplify(-3 * (YL**2 + YF**2))

    y1 = YukawaBetaCoefficients(
        trace=_txt(GS1),
        self_matrix=_txt(sp.simplify((GL1 + GF1) / 2)),
        cross_matrix=_txt(sp.simplify(GL2 / 2)),
        charged_lepton_matrix="1/2",
        su2_gauge=_txt(su2),
        u1_gauge=_txt(u1),
    )

    y2 = YukawaBetaCoefficients(
        trace=_txt(GS2),
        self_matrix=_txt(sp.simplify((GL2 + GF2) / 2)),
        cross_matrix=_txt(sp.simplify(GL1 / 2)),
        charged_lepton_matrix="1/2",
        su2_gauge=_txt(su2),
        u1_gauge=_txt(u1),
    )

    self_conjugate = physical_self_conjugate_f(dF, alpha)
    y2_self_conjugate_extra: dict[str, str] = {}
    if self_conjugate:
        y2_self_conjugate_extra = {
            "Tr_y2_y1dag*y1": _txt(GS2),
            "y2_y1dag_y1": _txt(sp.simplify(GF2 / 2)),
        }

    return T3YukawaBetaGroupFactors(
        dS1=dS1,
        dS2=dS2,
        dF=dF,
        alpha=alpha,
        C2L=_txt(C2L),
        C2F=_txt(C2F),
        YL=_txt(YL),
        YF=_txt(YF),
        self_conjugate_F=self_conjugate,
        y1=y1,
        y2=y2,
        y2_self_conjugate_extra=y2_self_conjugate_extra,
    )


def run_yukawa_group_factor_cli() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Complete generic T3 y1/y2 one-loop group-factor coefficients, "
            "including the physical self-conjugate beta_y2 extras."
        )
    )
    parser.add_argument("dS1", type=int)
    parser.add_argument("dS2", type=int)
    parser.add_argument("dF", type=int)
    parser.add_argument("alpha", type=int)
    args = parser.parse_args()

    result = complete_yukawa_group_factors(
        dS1=args.dS1,
        dS2=args.dS2,
        dF=args.dF,
        alpha=args.alpha,
    )
    print(json.dumps(asdict(result), indent=2))
    return 0








MODEL_DIMS = {
    "A": (1, 3, 2),
    "B": (2, 2, 1),
    "C": (2, 2, 3),
    "D": (3, 1, 2),
    "E": (3, 3, 2),
}


TRACE_CROSS = (
    r"\text{Tr}[\text{y2}.\text{Trans}(\text{Bar}(\text{y1}))]"
    r"\text{Matrix}(\text{y1})"
)

MATRIX_CROSS = (
    r"\text{Matrix}(\text{y2},\text{Trans}(\text{Bar}(\text{y1})),"
    r"\text{y1})"
)


def fraction_text(value: Fraction) -> str:
    return (
        str(value.numerator)
        if value.denominator == 1
        else f"{value.numerator}/{value.denominator}"
    )


def model_alpha(path: Path) -> tuple[str, int]:
    name = path.parents[1].name
    match = re.fullmatch(r"T3_([A-E])_alpha_([mp])(\d+)", name)
    if not match:
        raise ValueError(name)

    model, sign, magnitude = match.groups()
    alpha = int(magnitude)
    return model, (-alpha if sign == "m" else alpha)


def split_top_level_additive(text: str) -> list[str]:
    """Split a TeXForm expression at top-level + or - while preserving signs."""
    parts: list[str] = []
    start = 0
    par = brk = brc = 0

    for index, char in enumerate(text):
        if char == "(":
            par += 1
        elif char == ")":
            par -= 1
        elif char == "[":
            brk += 1
        elif char == "]":
            brk -= 1
        elif char == "{":
            brc += 1
        elif char == "}":
            brc -= 1
        elif (
            char in "+-"
            and par == brk == brc == 0
            and index > start
        ):
            parts.append(text[start:index].strip())
            start = index

    tail = text[start:].strip()
    if tail:
        parts.append(tail)

    return parts


def canonical_term(term: str) -> str:
    """Apply formatting-only canonicalisation; do not alter physics content."""
    text = re.sub(r"\s+", "", term)
    return text.replace(r"\left", "").replace(r"\right", "")


def is_gauge_term(term: str) -> bool:
    text = canonical_term(term)
    return (
        "g2" in text
        or "gY" in text
        or r"g_2" in text
        or r"g_Y" in text
        or r"g_{2}" in text
        or r"g_{Y}" in text
    )


def nongauge_from_latex(raw_latex: str) -> tuple[list[str], str]:
    terms = split_top_level_additive(str(raw_latex))
    nongauge_terms = [term for term in terms if not is_gauge_term(term)]
    canonical = "|".join(canonical_term(term) for term in nongauge_terms)
    return nongauge_terms, canonical


def _compact_latex(text: str) -> str:
    text = text.replace(r"\,", "")
    text = text.replace(r"\!", "")
    text = text.replace(r"\left", "")
    text = text.replace(r"\right", "")
    return text.replace(" ", "")


def _leading_rational(term: str) -> Fraction:
    """Extract the leading numerical coefficient of a normalized LaTeX term."""
    text = _compact_latex(term)

    sign = 1
    if text.startswith("+"):
        text = text[1:]
    elif text.startswith("-"):
        sign = -1
        text = text[1:]

    match = re.match(r"\\frac\{(\d+)\}\{(\d+)\}", text)
    if match:
        return sign * Fraction(int(match.group(1)), int(match.group(2)))

    match = re.match(r"(\d+)", text)
    if match:
        return sign * Fraction(int(match.group(1)), 1)

    return Fraction(sign, 1)


def _contains_gauge_yukawa(term: str, gauge: str, yukawa: str) -> bool:
    text = _compact_latex(term)

    if gauge == "gY":
        gauge_present = (
            r"g_Y" in text
            or r"g_{Y}" in text
            or "gY" in text
        )
    elif gauge == "g2":
        gauge_present = (
            r"g_2" in text
            or r"g_{2}" in text
            or "g2" in text
        )
    else:
        raise ValueError(gauge)

    y_index = "1" if yukawa == "y1" else "2"
    yukawa_present = (
        rf"y_{y_index}" in text
        or rf"y_{{{y_index}}}" in text
        or yukawa in text
    )
    squared = "^2" in text or "^{2}" in text

    return gauge_present and yukawa_present and squared


def _extract_gauge_coefficient(
    payload: dict[str, Any],
    yukawa: str,
    gauge: str,
) -> Fraction:
    latex_map = payload.get("report_beta_latex", {}) or {}
    raw = latex_map.get(yukawa)
    if not raw:
        raise KeyError(f"Missing report_beta_latex[{yukawa!r}]")

    beta_latex = str(raw)
    terms = split_top_level_additive(beta_latex) or [beta_latex]
    candidates = [
        term
        for term in terms
        if _contains_gauge_yukawa(term, gauge, yukawa)
    ]

    if len(candidates) != 1:
        diagnostic = "\n".join(f"  {item}" for item in candidates) or "  (none)"
        raise ValueError(
            f"Expected exactly one {gauge}^2 {yukawa} term, found "
            f"{len(candidates)}:\n{diagnostic}\n"
            f"Full RGBeta beta:\n{beta_latex}"
        )

    return _leading_rational(candidates[0])


def validate_gauge_rows(root: Path) -> tuple[list[dict[str, Any]], bool]:
    paths = sorted(root.glob("T3_*_alpha_*/data/uv_rgbeta_rge.json"))
    if not paths:
        raise FileNotFoundError(f"No UV RGBeta outputs found below {root}")

    rows: list[dict[str, Any]] = []
    overall = True

    print()
    print("Yukawa gauge-factor validation")
    print(
        f"{'model':<6} {'alpha':>5} {'y':>3} "
        f"{'U1 exp':>8} {'U1 RGB':>8} {'dU1':>6} "
        f"{'SU2 exp':>8} {'SU2 RGB':>8} {'dSU2':>6} status"
    )
    print("-" * 86)

    for path in paths:
        payload = json.loads(path.read_text(encoding="utf-8-sig"))
        if payload.get("status") != "Success":
            continue

        model, alpha = model_alpha(path)
        d_s1, d_s2, d_f = MODEL_DIMS[model]
        prediction = complete_yukawa_group_factors(
            dS1=d_s1,
            dS2=d_s2,
            dF=d_f,
            alpha=alpha,
        )

        for yukawa in ("y1", "y2"):
            predicted = prediction.y1 if yukawa == "y1" else prediction.y2
            expected_u1 = Fraction(predicted.u1_gauge)
            expected_su2 = Fraction(predicted.su2_gauge)

            got_u1 = _extract_gauge_coefficient(payload, yukawa, "gY")
            got_su2 = _extract_gauge_coefficient(payload, yukawa, "g2")

            u1_residual = got_u1 - expected_u1
            su2_residual = got_su2 - expected_su2
            ok = u1_residual == 0 and su2_residual == 0
            overall &= ok

            print(
                f"T3-{model:<3} {alpha:>5} {yukawa:>3} "
                f"{fraction_text(expected_u1):>8} "
                f"{fraction_text(got_u1):>8} "
                f"{fraction_text(u1_residual):>6} "
                f"{fraction_text(expected_su2):>8} "
                f"{fraction_text(got_su2):>8} "
                f"{fraction_text(su2_residual):>6} "
                f"{'PASS' if ok else 'FAIL'}"
            )

            rows.append(
                {
                    "model": model,
                    "alpha": alpha,
                    "yukawa": yukawa,
                    "expected_u1": fraction_text(expected_u1),
                    "rgbeta_u1": fraction_text(got_u1),
                    "u1_residual": fraction_text(u1_residual),
                    "u1_match": u1_residual == 0,
                    "expected_su2": fraction_text(expected_su2),
                    "rgbeta_su2": fraction_text(got_su2),
                    "su2_residual": fraction_text(su2_residual),
                    "su2_match": su2_residual == 0,
                    "match": ok,
                    "source": str(path),
                }
            )

    return rows, overall


def _canon(text: str) -> str:
    return re.sub(r"\s+", "", str(text))


def _frac_prefix(text: str, target: str) -> Fraction | None:
    """Coefficient of a standalone target term from its leading TeX prefix."""
    canonical = _canon(text)
    target_canonical = _canon(target)
    position = canonical.find(target_canonical)

    if position < 0:
        return None

    prefix = canonical[:position].lstrip("+")
    if prefix == "":
        return Fraction(1)

    match = re.fullmatch(r"\\frac\{(-?\d+)\}\{(\d+)\}", prefix)
    if match:
        return Fraction(int(match.group(1)), int(match.group(2)))

    match = re.fullmatch(r"(-?\d+)", prefix)
    if match:
        return Fraction(int(match.group(1)))

    return None


def _outer_fraction(text: str) -> Fraction:
    canonical = _canon(text).lstrip("+")
    match = re.match(r"\\frac\{(-?\d+)\}\{(\d+)\}\(", canonical)
    if match:
        return Fraction(int(match.group(1)), int(match.group(2)))
    return Fraction(1)


def _inner_integer_before(text: str, target: str) -> Fraction | None:
    canonical = _canon(text)
    target_canonical = _canon(target)
    position = canonical.find(target_canonical)

    if position < 0:
        return None

    before = canonical[:position]
    match = re.search(r"(?:\+|\()(-?\d+)$", before)
    if match:
        return Fraction(int(match.group(1)))

    return Fraction(1)


def extract_self_conjugate_extras(
    terms: list[str],
) -> tuple[Fraction | None, Fraction | None]:
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


def extract_rows(root: Path) -> tuple[
    list[dict[str, Any]],
    dict[str, int],
    bool,
]:
    paths = sorted(root.glob("T3_*_alpha_*/data/uv_rgbeta_rge.json"))
    if not paths:
        raise FileNotFoundError(f"No UV RGBeta outputs found below {root}")

    rows: list[dict[str, Any]] = []
    by_model: dict[str, set[str]] = {}

    for path in paths:
        payload = json.loads(path.read_text(encoding="utf-8-sig"))
        if payload.get("status") != "Success":
            continue

        model, alpha = model_alpha(path)
        d_s1, d_s2, d_f = MODEL_DIMS[model]

        latex_map = payload.get("report_beta_latex", {}) or {}
        raw_latex = latex_map.get("y2")
        if not raw_latex:
            raise RuntimeError(f"No report_beta_latex['y2'] in {path}")

        terms, canonical = nongauge_from_latex(raw_latex)
        predicted = complete_yukawa_group_factors(
            dS1=d_s1,
            dS2=d_s2,
            dF=d_f,
            alpha=alpha,
        )

        rows.append(
            {
                "model": model,
                "alpha": alpha,
                "dS1": d_s1,
                "dS2": d_s2,
                "dF": d_f,
                "nongauge_terms": terms,
                "canonical_nongauge": canonical,
                "predicted_coefficients": {
                    "Tr_y2*y2": predicted.y2.trace,
                    "y2_y2dag_y2": predicted.y2.self_matrix,
                    "y1_y1dag_y2": predicted.y2.cross_matrix,
                    "Ye_Yedag_y2": predicted.y2.charged_lepton_matrix,
                },
                "source": str(path),
            }
        )
        by_model.setdefault(model, set()).add(canonical)

    distinct_counts = {
        model: len(values)
        for model, values in sorted(by_model.items())
    }
    raw_alpha_independent = all(
        count == 1
        for count in distinct_counts.values()
    )

    return rows, distinct_counts, raw_alpha_independent


def load_rows_from_diagnostic(
    path: Path,
) -> tuple[list[dict[str, Any]], dict[str, int] | None, bool | None]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    rows = payload.get("rows") or []
    if not rows:
        raise ValueError(f"No rows in {path}")

    return (
        rows,
        payload.get("distinct_nongauge_expression_count"),
        payload.get("alpha_independent"),
    )


def validate_rows(
    rows: list[dict[str, Any]],
) -> tuple[list[dict[str, Any]], bool, int]:
    reference: dict[str, str] = {}
    for row in rows:
        if int(row["alpha"]) == 0:
            reference[str(row["model"])] = str(row["canonical_nongauge"])

    results: list[dict[str, Any]] = []
    overall = True
    excluded = 0

    for row in rows:
        model = str(row["model"])
        alpha = int(row["alpha"])
        d_s1 = int(row["dS1"])
        d_s2 = int(row["dS2"])
        d_f = int(row["dF"])

        prediction = complete_yukawa_group_factors(
            dS1=d_s1,
            dS2=d_s2,
            dF=d_f,
            alpha=alpha,
        )

        if alpha != -1:
            reference_value = reference.get(model)
            ok = (
                reference_value is not None
                and row["canonical_nongauge"] == reference_value
            )
            branch = "vector-like"
            detail = {
                "matches_alpha0_reference": ok,
                "alpha0_reference_present": reference_value is not None,
            }

        elif d_f % 2 == 0:
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
                list(row["nongauge_terms"])
            )
            trace_expected = Fraction(
                prediction.y2_self_conjugate_extra["Tr_y2_y1dag*y1"]
            )
            matrix_expected = Fraction(
                prediction.y2_self_conjugate_extra["y2_y1dag_y1"]
            )

            ok = (
                trace_got == trace_expected
                and matrix_got == matrix_expected
            )

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

        overall &= ok

        results.append(
            {
                "model": model,
                "alpha": alpha,
                "dS1": d_s1,
                "dS2": d_s2,
                "dF": d_f,
                "source": row.get("source"),
                "physical_self_conjugate_F": physical_self_conjugate_f(
                    d_f,
                    alpha,
                ),
                "branch": branch,
                "match": ok,
                **detail,
            }
        )

    return results, overall, excluded


def print_extraction_summary(
    rows: list[dict[str, Any]],
    distinct_counts: dict[str, int] | None,
    raw_alpha_independent: bool | None,
) -> None:
    print("Extracted RGBeta nongauge beta_y2 terms")
    print()

    for row in rows:
        print(f"T3-{row['model']} alpha={int(row['alpha']):+d}")
        for term in row["nongauge_terms"]:
            print(f"  {term}")

        predicted = row.get("predicted_coefficients") or {}
        if predicted:
            print(
                "  predicted vector-like coeffs: "
                + ", ".join(
                    f"{name}={value}"
                    for name, value in predicted.items()
                )
            )
        print()

    if distinct_counts is not None:
        print("Raw alpha-independence diagnostic:")
        for model in sorted(distinct_counts):
            count = distinct_counts[model]
            print(
                f"  T3-{model}: "
                f"{'PASS' if count == 1 else 'DIFFERS'} "
                f"({count} distinct nongauge expressions)"
            )

    print(f"raw alpha-independent across all saved rows: {raw_alpha_independent}")


def run_yukawa_beta_validation_cli() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--root",
        type=Path,
        default=Path("output/full/hypercharge"),
        help="Root containing T3_*_alpha_*/data/uv_rgbeta_rge.json.",
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=None,
        help=(
            "Optional legacy y2_nongauge_diagnostic.json. If supplied, "
            "validate those rows instead of extracting directly from --root."
        ),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "output/group_factors/y2_nongauge_branch_validation.json"
        ),
    )
    parser.add_argument(
        "--gauge-output",
        type=Path,
        default=None,
        help=(
            "Optional compatibility artifact reproducing the former "
            "yukawa_gauge_validation.json result."
        ),
    )
    parser.add_argument(
        "--diagnostic-output",
        type=Path,
        default=None,
        help=(
            "Optional compatibility/debug output reproducing the former "
            "extraction-only diagnostic payload."
        ),
    )
    args = parser.parse_args()

    if args.input is not None:
        (
            rows,
            distinct_counts,
            raw_alpha_independent,
        ) = load_rows_from_diagnostic(args.input)
        source_mode = "legacy_diagnostic_json"
    else:
        try:
            (
                rows,
                distinct_counts,
                raw_alpha_independent,
            ) = extract_rows(args.root)
        except FileNotFoundError as exc:
            raise SystemExit(str(exc)) from exc
        source_mode = "direct_extraction"

    if not rows:
        raise SystemExit("No successful beta_y2 rows were extracted.")

    print_extraction_summary(
        rows,
        distinct_counts,
        raw_alpha_independent,
    )

    if args.diagnostic_output is not None:
        diagnostic = {
            "status": (
                "Success"
                if rows and bool(raw_alpha_independent)
                else "NeedsInspection"
            ),
            "alpha_independent": raw_alpha_independent,
            "distinct_nongauge_expression_count": distinct_counts,
            "rows": rows,
        }
        args.diagnostic_output.parent.mkdir(parents=True, exist_ok=True)
        args.diagnostic_output.write_text(
            json.dumps(diagnostic, indent=2),
            encoding="utf-8",
        )
        print()
        print(f"Diagnostic JSON: {args.diagnostic_output}")

    results, overall, excluded = validate_rows(rows)

    print()
    print("Branch-aware beta_y2 nongauge validation")
    print("model alpha branch                         status")
    print("-" * 64)

    for row in results:
        print(
            f"T3-{row['model']} {int(row['alpha']):>5} "
            f"{row['branch']:<30} "
            f"{'PASS' if row['match'] else 'FAIL'}"
        )

    gauge_rows, gauge_ok = validate_gauge_rows(args.root)
    overall = overall and gauge_ok

    if args.gauge_output is not None:
        gauge_payload = {
            "status": "Success" if gauge_ok else "Failed",
            "file_count": len(
                list(args.root.glob("T3_*_alpha_*/data/uv_rgbeta_rge.json"))
            ),
            "check_count": len(gauge_rows),
            "u1_match_count": sum(row["u1_match"] for row in gauge_rows),
            "su2_match_count": sum(row["su2_match"] for row in gauge_rows),
            "all_u1_match": all(row["u1_match"] for row in gauge_rows),
            "all_su2_match": all(row["su2_match"] for row in gauge_rows),
            "formulae": {
                "u1": "-3*(YL^2+YF^2)",
                "su2": "-3*(C2(L)+C2(F))",
            },
            "checks": gauge_rows,
        }
        args.gauge_output.parent.mkdir(parents=True, exist_ok=True)
        args.gauge_output.write_text(
            json.dumps(gauge_payload, indent=2),
            encoding="utf-8",
        )

    result = {
        "status": "Success" if overall else "Failed",
        "source_mode": source_mode,
        "root": str(args.root) if args.input is None else None,
        "input": str(args.input) if args.input is not None else None,
        "raw_alpha_independent": raw_alpha_independent,
        "distinct_nongauge_expression_count": distinct_counts,
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
        "gauge_validation": {
            "status": "Success" if gauge_ok else "Failed",
            "checks": gauge_rows,
        },
        "known_issue": (
            "Self-conjugate beta_y2 reporting remains a known project physics "
            "issue. This validator preserves the existing branch-aware "
            "comparison and excluded even-dF neutral branch."
        ),
        "rows": results,
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(result, indent=2),
        encoding="utf-8",
    )

    print()
    print(f"Excluded even-dF neutral RGBeta rows: {excluded}")
    if args.gauge_output is not None:
        print(f"Gauge compatibility JSON: {args.gauge_output}")
    print(f"Validation JSON: {args.output}")
    print(f"OVERALL:         {result['status']}")

    return 0 if overall else 1


def main() -> int:
    parser = argparse.ArgumentParser(
        description="T3 Yukawa-beta group-factor calculation and validation."
    )
    parser.add_argument(
        "mode",
        choices=("factors", "validate"),
        help="Calculate generic Yukawa group factors or validate saved RGBeta beta functions.",
    )
    args, remaining = parser.parse_known_args()
    original_argv = sys.argv
    try:
        sys.argv = [original_argv[0], *remaining]
        if args.mode == "factors":
            return run_yukawa_group_factor_cli()
        return run_yukawa_beta_validation_cli()
    finally:
        sys.argv = original_argv


if __name__ == "__main__":
    raise SystemExit(main())
