"""Exact, independent one-loop gauge-factor projections for the general T3 basis.

This module is intentionally independent of Matchete, RGBeta, and repository
state. It predicts universal one-loop gauge/group factors from quantum numbers
and extracts the corresponding rational coefficient from *native* RGBeta beta
expressions. It does NOT establish an inter-package physical coupling map.

Convention: beta_X = 16*pi^2*dX/dln(mu), Q=T3+Y, T(2)=1/2, three F flavours.
"""
from __future__ import annotations

from fractions import Fraction
import re
from typing import Iterable

import sympy as sp
from sympy.parsing.mathematica import parse_mathematica


class ProjectionError(ValueError):
    """Expression could not be projected without assumptions."""


def C2(dim: int) -> Fraction:
    if type(dim) is not int or dim not in (1, 2, 3):
        raise ProjectionError(f"Unsupported SU(2) representation dimension: {dim}")
    return Fraction(dim * dim - 1, 4)


def T2(dim: int) -> Fraction:
    C2(dim)
    return Fraction(dim * (dim * dim - 1), 12)


def field_data(dims: tuple[int, int, int], alpha: int) -> dict[str, tuple[int, Fraction]]:
    d1, d2, dF = dims
    return {
        "H": (2, Fraction(1, 2)),
        "S1": (d1, Fraction(alpha, 2)),
        "S2": (d2, Fraction(alpha + 2, 2)),
        "F": (dF, Fraction(alpha + 1, 2)),
    }


def expected_gauge_betas(dims: tuple[int, int, int], alpha: int,
                         stage: str) -> dict[str, Fraction]:
    """b_i in beta(g_i)=b_i*g_i^3; SM with two complex scalars and 3 heavy flavours.

    For UV, F is 3 Majorana Weyl flavours when YF=0 (B/C alpha=-1),
    otherwise three Dirac/vector-like flavours. EFT1 contains no F.
    The new scalar fields are *complex* even when their SU(2) irreps are real.
    """
    if stage not in ("UV", "EFT1"):
        raise ProjectionError(f"Unknown RGE stage: {stage}")
    f = field_data(dims, alpha)
    b1 = Fraction(41, 6)
    b2 = Fraction(-19, 6)
    for species in ("S1", "S2"):
        d, y = f[species]
        b1 += Fraction(d, 3) * y * y
        b2 += Fraction(1, 3) * T2(d)
    if stage == "UV":
        dF, yF = f["F"]
        # Two Weyl flavours per vector-like fermion, one per real Majorana.
        mult = 1 if yF == 0 else 2
        b1 += Fraction(2 * 3 * mult * dF, 3) * yF * yF
        b2 += Fraction(2 * 3 * mult, 3) * T2(dF)
    return {"gY": b1, "g2": b2, "g3": Fraction(-7)}


def expected_scalar_gauge_factors(
    counts: Iterable[int], dims: tuple[int, int, int], alpha: int,
) -> dict[str, Fraction]:
    """-3 sum_i C2(R_i) and -3 sum_i Y_i^2 for a scalar interaction.

    `counts` orders H,Hdag,S1,S1dag,S2,S2dag and applies equally to scalar
    quartics (4 fields) and bilinears (2 fields). These diagonal gauge-linear
    coefficients follow from scalar wavefunctions in the chosen convention.
    They are not the full beta functions, which also mix scalar invariants.
    """
    c = tuple(counts)
    if len(c) != 6 or any(type(x) is not int or x < 0 for x in c) or sum(c) not in (2, 4):
        raise ProjectionError(f"Invalid scalar operator counts: {counts}")
    f = field_data(dims, alpha)
    dH, yH = f["H"]
    d1, y1 = f["S1"]
    d2, y2 = f["S2"]
    su2 = -3 * (sum(c[:2]) * C2(dH) + sum(c[2:4]) * C2(d1) + sum(c[4:]) * C2(d2))
    u1 = -3 * (sum(c[:2]) * yH**2 + sum(c[2:4]) * y1**2 + sum(c[4:]) * y2**2)
    return {"g2": su2, "gY": u1}


def expected_yukawa_gauge_factors(dims: tuple[int, int, int], alpha: int) -> dict[str, Fraction]:
    dF, yF = field_data(dims, alpha)["F"]
    return {"g2": -3 * (C2(2) + C2(dF)),
            "gY": -3 * (Fraction(1, 4) + yF**2)}


def _matching_bracket(text: str, opened: int) -> int:
    if text[opened] != "[":
        raise ProjectionError("Expected [")
    depth = 0
    for i in range(opened, len(text)):
        if text[i] == "[":
            depth += 1
        elif text[i] == "]":
            depth -= 1
            if depth == 0:
                return i
    raise ProjectionError("Unbalanced Mathematica brackets")


def _matrix_replacement(text: str, matrix_target: str) -> str:
    """Replace matrix-valued factors atomically, preserving only the target.

    RGBeta expresses Yukawa beta matrices as Matrix[y1][...] products. Only
    a linear target matrix may contribute to the diagonal gauge projection.
    Every other Matrix[...][...] is an independent dummy source set to zero.
    """
    target = "Matrix[" + matrix_target + "]"
    out: list[str] = []
    pos = 0
    while True:
        match = re.search(r"(?<![A-Za-z0-9])Matrix\[", text[pos:])
        if match is None:
            out.append(text[pos:])
            break
        start = pos + match.start()
        out.append(text[pos:start])
        opener = start + len("Matrix")
        close = _matching_bracket(text, opener)
        head = text[start:close + 1].replace(" ", "")
        end = close + 1
        # RGBeta optionally appends output flavour indices after the matrix.
        while end < len(text) and text[end] == "[":
            end = _matching_bracket(text, end) + 1
        out.append("TARGETM" if head == target else "OTHERMATRIX")
        pos = end
    return "".join(out)


def _top_level_terms(text: str) -> list[str]:
    """Split native InputForm at additive operators outside all delimiters."""
    terms: list[str] = []
    start = 0
    p = b = c = 0
    for i, char in enumerate(text):
        if char == "(": p += 1
        elif char == ")": p -= 1
        elif char == "[": b += 1
        elif char == "]": b -= 1
        elif char == "{": c += 1
        elif char == "}": c -= 1
        elif char in "+-" and i > start and p == b == c == 0:
            prev = text[:i].rstrip()
            if prev.endswith(("*", "/", "^", "E", "e")):
                continue
            terms.append(text[start:i].strip())
            start = i
        if min(p, b, c) < 0:
            raise ProjectionError("Unbalanced native beta expression")
    if p or b or c:
        raise ProjectionError("Unbalanced native beta expression")
    if text[start:].strip():
        terms.append(text[start:].strip())
    return terms


def _symbol_match(text: str, token: str) -> bool:
    return re.search(r"(?<![A-Za-z0-9_])" + re.escape(token) + r"(?![A-Za-z0-9_])", text) is not None


def _rational_coefficient(expr: sp.Expr, monomial: sp.Expr) -> Fraction:
    """Return the exact coefficient of `monomial` at zero other couplings.

    Never evaluate non-polynomial, symbolic or matrix-valued remainders as 0
    by implicit approximation; return an explicit `ProjectionError` instead.
    """
    if not isinstance(expr, sp.Expr) or not isinstance(monomial, sp.Expr):
        raise ProjectionError("Expression parser did not return SymPy algebra")
    # All non-target couplings and unidentified traces are independent variables.
    for symbol in expr.free_symbols - monomial.free_symbols:
        expr = expr.subs(symbol, 0)
    # If function heads survive, they have not been shown to vanish. Fail closed.
    if expr.atoms(sp.Function):
        raise ProjectionError("Unresolved function in projected gauge coefficient")
    expr = sp.expand(expr)
    try:
        coefficient = sp.Poly(expr, *sorted(monomial.free_symbols, key=str)).coeff_monomial(monomial)
    except (sp.PolynomialError, ValueError, TypeError) as exc:
        raise ProjectionError(f"Cannot isolate polynomial gauge coefficient: {exc}") from exc
    coefficient = sp.cancel(coefficient)
    if coefficient.is_Rational is not True:
        raise ProjectionError(f"Gauge coefficient not a rational number: {coefficient}")
    return Fraction(int(coefficient.p), int(coefficient.q))


def project_native_beta(
    raw: str, *, gauge: str, coupling: str | None = None,
    matrix_coupling: str | None = None, zero_is_valid: bool = False,
) -> Fraction:
    """Project b_g or a diagonal -g^2*coupling coefficient from native RGBeta.

    When coupling is None, the expected monomial is g^3.
    When coupling is given, the expected monomial is g^2*coupling.
    When matrix_coupling is given, it is g^2*Matrix[coupling] (Yukawa).
    Selects only syntactically relevant summands, reducing parse cost and
    preventing unrelated complicated scalar invariant functions from hiding a
    successful coefficient projection. Full expression remains unchanged.
    """
    if gauge not in ("gY", "g2", "g3"):
        raise ProjectionError(f"Unsupported gauge symbol: {gauge}")
    if (coupling is not None) and (matrix_coupling is not None):
        raise ProjectionError("Choose scalar or matrix coupling, not both")
    if not isinstance(raw, str) or not raw.strip():
        raise ProjectionError("Empty RGBeta beta expression")
    g = sp.Symbol(gauge)
    if matrix_coupling is not None:
        target = sp.Symbol("TARGETM")
    elif coupling is not None:
        target = sp.Symbol(coupling)
    else:
        target = None
    monomial = g**3 if target is None else g**2 * target
    selected: list[str] = []
    for term in _top_level_terms(raw):
        if not _symbol_match(term, gauge):
            continue
        if coupling is not None and not _symbol_match(term, coupling):
            continue
        if matrix_coupling is not None and not _symbol_match(term, "Matrix"):
            continue
        if target is None and ("Matrix[" in term or "Tr[" in term):
            # Gauge beta must be a pure gauge polynomial at one loop.
            raise ProjectionError("Unexpected trace or matrix in gauge beta")
        selected.append(_matrix_replacement(term, matrix_coupling) if matrix_coupling else term)
    if not selected:
        if zero_is_valid:
            return Fraction(0)
        raise ProjectionError(f"Missing gauge monomial {gauge}, target={coupling or matrix_coupling}")
    try:
        parsed = sp.Add(*(parse_mathematica(term) for term in selected))
    except (ValueError, TypeError, SyntaxError, IndexError, RecursionError) as exc:
        raise ProjectionError(f"Mathematica parser rejected relevant gauge terms: {exc}") from exc
    result = _rational_coefficient(parsed, monomial)
    if result == 0 and not zero_is_valid:
        raise ProjectionError("Target linear gauge monomial absent after exact polynomial projection")
    return result
