"""Parsing helpers for Matchete-exported Weinberg coefficients.

This module owns only the textual Matchete -> SymPy compatibility boundary.
No RGE equation, numerical evolution, or stage/output orchestration belongs
here.
"""

from __future__ import annotations

import sympy as sp


def _matching_square_bracket(text: str, open_index: int) -> int:
    if text[open_index] != "[":
        raise ValueError("Expected '[' at open_index.")

    depth = 0
    for index in range(open_index, len(text)):
        char = text[index]
        if char == "[":
            depth += 1
        elif char == "]":
            depth -= 1
            if depth == 0:
                return index

    raise ValueError("Unbalanced Mathematica square brackets.")


def _first_top_level_argument(arguments: str) -> str:
    depth = 0

    for pos, char in enumerate(arguments):
        if char in "[{(":
            depth += 1
        elif char in "]})":
            depth -= 1
        elif char == "," and depth == 0:
            return arguments[:pos].strip()

    return arguments.strip()


def _replace_matchete_barred_couplings(text: str) -> str:
    marker = "Bar[Coupling["
    result = text

    while marker in result:
        start = result.index(marker)
        bar_open = start + len("Bar")
        bar_close = _matching_square_bracket(result, bar_open)
        inner = result[bar_open + 1 : bar_close]

        if not inner.startswith("Coupling["):
            raise ValueError(
                f"Unexpected Bar contents while parsing Matchete C5: {inner}"
            )

        coupling_open = inner.index("[")
        coupling_close = _matching_square_bracket(inner, coupling_open)
        arguments = inner[coupling_open + 1 : coupling_close]
        name = _first_top_level_argument(arguments)

        result = (
            result[:start]
            + f"conjugate({name})"
            + result[bar_close + 1 :]
        )

    return result


def _replace_matchete_couplings(text: str) -> str:
    marker = "Coupling["
    result = text

    while marker in result:
        start = result.index(marker)
        open_index = start + len("Coupling")
        close_index = _matching_square_bracket(result, open_index)
        arguments = result[open_index + 1 : close_index]
        name = _first_top_level_argument(arguments)
        result = result[:start] + name + result[close_index + 1 :]

    return result


def parse_matchete_c5(text: str) -> sp.Expr:
    """Convert one Matchete C5 expression into the project's SymPy form."""

    stripped = text.strip()
    if not stripped:
        raise ValueError("The matched C5 coefficient file is empty.")

    cleaned = _replace_matchete_barred_couplings(stripped)
    cleaned = _replace_matchete_couplings(cleaned)
    cleaned = cleaned.replace("Conjugate[", "conjugate(")
    cleaned = cleaned.replace("Log[", "log(")
    cleaned = cleaned.replace("Sqrt[", "sqrt(")
    cleaned = cleaned.replace("]", ")")
    cleaned = cleaned.replace("^", "**")

    return sp.simplify(
        sp.sympify(
            cleaned,
            locals={
                "log": sp.log,
                "sqrt": sp.sqrt,
                "conjugate": sp.conjugate,
            },
        )
    )
