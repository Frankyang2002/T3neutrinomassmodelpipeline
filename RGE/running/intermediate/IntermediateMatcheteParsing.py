from __future__ import annotations

"""Small Matchete ``InputForm`` parsing helpers for intermediate-EFT adapters.

Converts the subset of Matchete Mathematica expressions used by the scalar-only
matching/RGE bridge into exact SymPy objects.
"""

import re
from dataclasses import dataclass

import sympy as sp


def matching_bracket(
    text: str,
    open_index: int,
    left: str = "[",
    right: str = "]",
) -> int:
    """Return the matching closing delimiter for one nested expression."""
    if text[open_index] != left:
        raise ValueError(f"Expected {left!r} at index {open_index}.")

    depth = 0

    for pos in range(open_index, len(text)):
        if text[pos] == left:
            depth += 1
        elif text[pos] == right:
            depth -= 1
            if depth == 0:
                return pos

    raise ValueError(f"Unbalanced {left}{right} brackets.")


def matching_brace(text: str, open_index: int) -> int:
    """Return the matching closing Mathematica curly brace."""
    return matching_bracket(text, open_index, "{", "}")


def split_top_level(text: str, separator: str = ",") -> list[str]:
    """Split on separators that are outside [], {}, and ()."""
    result: list[str] = []
    start = 0
    square = curly = paren = 0

    for pos, char in enumerate(text):
        if char == "[":
            square += 1
        elif char == "]":
            square -= 1
        elif char == "{":
            curly += 1
        elif char == "}":
            curly -= 1
        elif char == "(":
            paren += 1
        elif char == ")":
            paren -= 1
        elif (
            char == separator
            and square == 0
            and curly == 0
            and paren == 0
        ):
            result.append(text[start:pos].strip())
            start = pos + 1

    result.append(text[start:].strip())
    return result


def strip_bar(text: str) -> tuple[str, bool]:
    """Strip one outer Matchete Bar[...] wrapper, if present."""
    stripped = text.strip()

    if stripped.startswith("Bar[") and stripped.endswith("]"):
        open_index = stripped.index("[")
        if matching_bracket(stripped, open_index) == len(stripped) - 1:
            return stripped[open_index + 1 : -1].strip(), True

    return stripped, False


def parse_index(text: str) -> tuple[str, str]:
    """Parse Index[dummy, representation] into its two labels."""
    inner, _ = strip_bar(text)

    if not inner.startswith("Index["):
        raise ValueError(f"Expected Index[...] but got {text!r}.")

    open_index = inner.index("[")
    close_index = matching_bracket(inner, open_index)
    args = split_top_level(inner[open_index + 1 : close_index])

    if len(args) != 2:
        raise ValueError(f"Unexpected index expression {text!r}.")

    dummy = args[0].strip()
    rep, _ = strip_bar(args[1].strip())

    return dummy, rep


def parse_reps(text: str) -> tuple[str, ...]:
    """Parse a Mathematica representation list, removing outer Bar wrappers."""
    stripped = text.strip()

    if not (stripped.startswith("{") and stripped.endswith("}")):
        raise ValueError(f"Unexpected representation list {text!r}.")

    reps: list[str] = []

    for item in split_top_level(stripped[1:-1]):
        rep, _ = strip_bar(item)
        reps.append(rep)

    return tuple(reps)


def parse_mathematica_scalar(text: str) -> sp.Expr:
    """Parse the scalar-expression subset used in exported Matchete data."""
    cleaned = text.strip()
    sqrt_pattern = re.compile(r"Sqrt\[([^\[\]]+)\]")

    while sqrt_pattern.search(cleaned):
        cleaned = sqrt_pattern.sub(r"sqrt(\1)", cleaned)

    cleaned = cleaned.replace("^", "**")

    return sp.simplify(
        sp.sympify(
            cleaned,
            locals={
                "sqrt": sp.sqrt,
                "I": sp.I,
            },
        )
    )


def parse_sparse_array(text: str) -> sp.MutableDenseNDimArray:
    """Parse Matchete's arbitrary-rank CSR-style SparseArray InputForm."""
    s = text.strip()
    prefix = "SparseArray[Automatic,"

    if not s.startswith(prefix):
        raise ValueError(
            f"Unsupported CG tensor representation: {s[:100]}"
        )

    dim_open = s.find("{", len(prefix))

    if dim_open < 0:
        raise ValueError("SparseArray dimensions were not found.")

    dim_close = matching_brace(s, dim_open)
    dims = tuple(
        int(piece.strip())
        for piece in split_top_level(s[dim_open + 1 : dim_close])
    )

    if len(dims) < 1:
        raise ValueError("SparseArray must have rank at least one.")

    payload_marker = "{1, {{"
    payload_start = s.find(payload_marker, dim_close)

    if payload_start < 0:
        raise ValueError("SparseArray CSR payload was not found.")

    row_start = payload_start + len("{1, {")

    if s[row_start] != "{":
        raise ValueError("Malformed SparseArray row-pointer list.")

    row_end = matching_brace(s, row_start)
    row_ptr = [
        int(piece.strip())
        for piece in split_top_level(s[row_start + 1 : row_end])
        if piece.strip()
    ]

    coord_start = s.find("{", row_end + 1)

    if coord_start < 0:
        raise ValueError("SparseArray coordinate list was not found.")

    coord_end = matching_brace(s, coord_start)
    coord_text = s[coord_start + 1 : coord_end].strip()

    coordinates: list[tuple[int, ...]] = []

    if coord_text:
        pos = 0

        while pos < len(coord_text):
            while pos < len(coord_text) and coord_text[pos] in " ,\t\r\n":
                pos += 1

            if pos >= len(coord_text):
                break

            if coord_text[pos] != "{":
                raise ValueError("Malformed SparseArray coordinate entry.")

            close = matching_brace(coord_text, pos)
            coords = tuple(
                int(piece.strip())
                for piece in split_top_level(coord_text[pos + 1 : close])
                if piece.strip()
            )
            coordinates.append(coords)
            pos = close + 1

    inner_pair_end = matching_brace(s, row_start - 1)
    values_start = s.find("{", inner_pair_end + 1)

    if values_start < 0:
        raise ValueError("SparseArray value list was not found.")

    values_end = matching_brace(s, values_start)
    value_text = s[values_start + 1 : values_end].strip()
    values = (
        []
        if not value_text
        else [
            parse_mathematica_scalar(piece)
            for piece in split_top_level(value_text)
        ]
    )

    if len(row_ptr) != dims[0] + 1:
        raise ValueError(
            "SparseArray row-pointer length is inconsistent with "
            f"first dimension {dims[0]}: {row_ptr}"
        )

    if len(coordinates) != len(values):
        raise ValueError(
            "SparseArray coordinate/value lengths do not agree: "
            f"{len(coordinates)} != {len(values)}."
        )

    expected_tail_rank = len(dims) - 1

    if any(len(coords) != expected_tail_rank for coords in coordinates):
        raise ValueError(
            "SparseArray coordinate rank does not match tensor rank."
        )

    data = sp.MutableDenseNDimArray.zeros(*dims)

    for first_index in range(dims[0]):
        lo = row_ptr[first_index]
        hi = row_ptr[first_index + 1]

        if not (0 <= lo <= hi <= len(values)):
            raise ValueError(
                "SparseArray row pointers are outside the nonzero-entry range."
            )

        for pointer in range(lo, hi):
            tail = coordinates[pointer]

            for axis, component in enumerate(tail, start=1):
                if not 1 <= component <= dims[axis]:
                    raise IndexError(
                        f"SparseArray component {component} is outside "
                        f"dimension {dims[axis]} on axis {axis + 1}."
                    )

            index = (first_index,) + tuple(
                component - 1
                for component in tail
            )
            data[index] = values[pointer]

    return data


@dataclass(frozen=True)
class CGTensor:
    name: str
    reps: tuple[str, ...]
    tensor: sp.MutableDenseNDimArray


def load_cg_registry(seed: dict) -> dict[str, CGTensor]:
    """Parse a Matchete ``CGRegistry`` payload into exact tensors."""
    result: dict[str, CGTensor] = {}
    for item in seed.get("CGRegistry", []):
        result[item["Name"]] = CGTensor(
            name=item["Name"],
            reps=parse_reps(item["RepsInputForm"]),
            tensor=parse_sparse_array(item["TensorInputForm"]),
        )
    return result
