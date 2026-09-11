from __future__ import annotations

"""Convert the exported EFT1 scalar potential to real-basis lambda_abcd.

Input:
    eft1_after_F_scalar_quartic_seed.json
    eft1_rgbeta_rge.json

Output convention:
    V4 = (1/4!) lambda_abcd phi_a phi_b phi_c phi_d

The Matchete export contains Lagrangian terms, so this adapter first forms
V4 = -L4.  Complex multiplets are then converted as

    Phi_m = (R_m + i I_m)/sqrt(2),
    Phi_m^* = (R_m - i I_m)/sqrt(2),

with the same interleaved real-scalar ordering used by RGEModel:

    H : 1..4
    S1: 5..(4 + 2 dS1)
    S2: next 2 dS2 entries.

The resulting lookup is fully symmetric in all four scalar indices.
"""

from dataclasses import dataclass
from itertools import product
import json
from math import factorial
from pathlib import Path
import re

import sympy as sp


@dataclass(frozen=True)
class ScalarLeg:
    name: str
    dummy: str | None
    rep: str | None
    conjugated: bool


@dataclass(frozen=True)
class CGTensor:
    name: str
    reps: tuple[str, ...]
    tensor: sp.MutableDenseNDimArray


@dataclass(frozen=True)
class CGCall:
    name: str
    conjugated: bool
    indices: tuple[tuple[str, str], ...]


@dataclass(frozen=True)
class ScalarLayout:
    d_s1: int
    d_s2: int

    @property
    def dimensions(self) -> dict[str, int]:
        return {
            "H": 2,
            "NewScalar1": self.d_s1,
            "NewScalar2": self.d_s2,
        }

    def global_real_index(
        self,
        field_name: str,
        component: int,
        *,
        imaginary: bool,
    ) -> int:
        if field_name == "H":
            first = 1
            dimension = 2
        elif field_name == "NewScalar1":
            first = 5
            dimension = self.d_s1
        elif field_name == "NewScalar2":
            first = 5 + 2 * self.d_s1
            dimension = self.d_s2
        else:
            raise KeyError(f"Unknown scalar field {field_name!r}.")

        if not 1 <= component <= dimension:
            raise IndexError(
                f"{field_name} component {component} outside 1,...,{dimension}."
            )

        return first + 2 * (component - 1) + (1 if imaginary else 0)


class SparseQuarticTensor:
    """Fully symmetric real-scalar lambda_abcd lookup."""

    def __init__(self, entries: dict[tuple[int, int, int, int], sp.Expr]):
        self.entries = {
            tuple(sorted(map(int, key))): sp.simplify(value)
            for key, value in entries.items()
            if sp.simplify(value) != 0
        }

    def __getitem__(self, key: tuple[int, int, int, int]) -> sp.Expr:
        if len(key) != 4:
            raise KeyError("Quartic tensor requires four scalar indices.")
        return self.entries.get(tuple(sorted(map(int, key))), sp.S.Zero)

    def __call__(self, a: int, b: int, c: int, d: int) -> sp.Expr:
        return self[a, b, c, d]

    def nonzero_items(self):
        return self.entries.items()


def _matching_bracket(
    text: str,
    open_index: int,
    left: str = "[",
    right: str = "]",
) -> int:
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


def _matching_brace(text: str, open_index: int) -> int:
    return _matching_bracket(text, open_index, "{", "}")


def _split_top_level(text: str, separator: str = ",") -> list[str]:
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


def _strip_bar(text: str) -> tuple[str, bool]:
    stripped = text.strip()
    if stripped.startswith("Bar[") and stripped.endswith("]"):
        if _matching_bracket(stripped, 3) == len(stripped) - 1:
            return stripped[4:-1].strip(), True
    return stripped, False


def _parse_index(text: str) -> tuple[str, str]:
    inner, _ = _strip_bar(text)
    if not inner.startswith("Index["):
        raise ValueError(f"Expected Index[...] but got {text!r}.")
    close_index = _matching_bracket(inner, 5)
    args = _split_top_level(inner[6:close_index])
    if len(args) != 2:
        raise ValueError(f"Unexpected index expression {text!r}.")
    rep, _ = _strip_bar(args[1])
    return args[0].strip(), rep.strip()


def _parse_reps(text: str) -> tuple[str, ...]:
    stripped = text.strip()
    if not (stripped.startswith("{") and stripped.endswith("}")):
        raise ValueError(f"Unexpected representation list {text!r}.")
    reps: list[str] = []
    for item in _split_top_level(stripped[1:-1]):
        rep, _ = _strip_bar(item)
        reps.append(rep)
    return tuple(reps)


def _parse_mathematica_scalar(text: str) -> sp.Expr:
    cleaned = text.strip().replace("^", "**")
    sqrt_pattern = re.compile(r"Sqrt\[([^\[\]]+)\]")
    while sqrt_pattern.search(cleaned):
        cleaned = sqrt_pattern.sub(r"sqrt(\1)", cleaned)

    return sp.simplify(
        sp.sympify(
            cleaned,
            locals={"sqrt": sp.sqrt, "I": sp.I},
        )
    )


def _parse_sparse_array(text: str) -> sp.MutableDenseNDimArray:
    """Parse Matchete's arbitrary-rank CSR-style SparseArray InputForm."""

    s = text.strip()
    prefix = "SparseArray[Automatic,"
    if not s.startswith(prefix):
        raise ValueError(f"Unsupported CG tensor representation: {s[:100]}")

    dim_open = s.find("{", len(prefix))
    dim_close = _matching_brace(s, dim_open)
    dims = tuple(
        int(piece.strip())
        for piece in _split_top_level(s[dim_open + 1 : dim_close])
    )

    payload_start = s.find("{1, {{", dim_close)
    if payload_start < 0:
        raise ValueError("SparseArray CSR payload was not found.")

    row_start = payload_start + len("{1, {")
    row_end = _matching_brace(s, row_start)
    row_ptr = [
        int(piece.strip())
        for piece in _split_top_level(s[row_start + 1 : row_end])
        if piece.strip()
    ]

    coord_start = s.find("{", row_end + 1)
    coord_end = _matching_brace(s, coord_start)
    coord_text = s[coord_start + 1 : coord_end].strip()

    coordinates: list[tuple[int, ...]] = []
    pos = 0
    while pos < len(coord_text):
        while pos < len(coord_text) and coord_text[pos] in " ,\t\r\n":
            pos += 1
        if pos >= len(coord_text):
            break
        close = _matching_brace(coord_text, pos)
        coordinates.append(
            tuple(
                int(piece.strip())
                for piece in _split_top_level(coord_text[pos + 1 : close])
            )
        )
        pos = close + 1

    inner_pair_start = row_start - 1
    inner_pair_end = _matching_brace(s, inner_pair_start)
    values_start = s.find("{", inner_pair_end + 1)
    values_end = _matching_brace(s, values_start)
    values_text = s[values_start + 1 : values_end].strip()
    values = (
        []
        if not values_text
        else [
            _parse_mathematica_scalar(piece)
            for piece in _split_top_level(values_text)
        ]
    )

    if len(row_ptr) != dims[0] + 1:
        raise ValueError("SparseArray row-pointer length is inconsistent.")
    if len(coordinates) != len(values):
        raise ValueError("SparseArray coordinate/value lengths do not agree.")

    tensor = sp.MutableDenseNDimArray.zeros(*dims)

    for first in range(dims[0]):
        for pointer in range(row_ptr[first], row_ptr[first + 1]):
            tail = coordinates[pointer]
            if len(tail) != len(dims) - 1:
                raise ValueError("SparseArray coordinate rank mismatch.")
            index = (first,) + tuple(component - 1 for component in tail)
            tensor[index] = values[pointer]

    return tensor


def _load_cg_registry(seed: dict) -> dict[str, CGTensor]:
    result: dict[str, CGTensor] = {}
    for item in seed.get("CGRegistry", []):
        result[item["Name"]] = CGTensor(
            name=item["Name"],
            reps=_parse_reps(item["RepsInputForm"]),
            tensor=_parse_sparse_array(item["TensorInputForm"]),
        )
    return result


def _find_scalar_legs(term: str) -> tuple[ScalarLeg, ...]:
    """Read actual scalar Field[...] factors, including outer Bar[...] wrappers."""

    legs: list[ScalarLeg] = []
    pos = 0

    while True:
        start = term.find("Field[", pos)
        if start < 0:
            break

        open_index = start + len("Field")
        close_index = _matching_bracket(term, open_index)
        args = _split_top_level(term[open_index + 1 : close_index])

        if len(args) >= 3 and args[1].strip() == "Scalar":
            field_name = args[0].strip()

            indices = args[2].strip()
            if not (indices.startswith("{") and indices.endswith("}")):
                raise ValueError("Scalar Field index argument is not a list.")

            index_body = indices[1:-1].strip()

            # SU(2) singlets are printed by Matchete with no representation
            # index at all: Field[NewScalar1, Scalar, {}, {}].  They still
            # have one complex component, but there is no dummy index to
            # parse or contract through a CG tensor.
            if not index_body:
                dummy = None
                rep = None
            else:
                index_items = _split_top_level(index_body)
                if len(index_items) != 1:
                    raise ValueError(
                        "EFT1 scalar adapter currently expects at most one "
                        "SU(2) index per complex scalar multiplet."
                    )
                dummy, rep = _parse_index(index_items[0])

            conjugated = False
            if start >= 4 and term[start - 4 : start] == "Bar[":
                bar_open = start - 1
                conjugated = (
                    _matching_bracket(term, bar_open) == close_index + 1
                )

            legs.append(
                ScalarLeg(
                    name=field_name,
                    dummy=dummy,
                    rep=rep,
                    conjugated=conjugated,
                )
            )

        pos = close_index + 1

    if len(legs) != 4:
        raise ValueError(f"Expected four scalar legs, found {len(legs)}.")

    return tuple(legs)


def _find_cg_calls(term: str) -> tuple[CGCall, ...]:
    calls: list[CGCall] = []
    pos = 0

    while True:
        start = term.find("CG[", pos)
        if start < 0:
            break

        open_index = start + 2
        close_index = _matching_bracket(term, open_index)
        args = _split_top_level(term[open_index + 1 : close_index])

        if len(args) != 2:
            raise ValueError("Unexpected CG call structure.")

        name, conjugated = _strip_bar(args[0])
        index_list = args[1].strip()

        if not (index_list.startswith("{") and index_list.endswith("}")):
            raise ValueError("CG index argument is not a list.")

        calls.append(
            CGCall(
                name=name,
                conjugated=conjugated,
                indices=tuple(
                    _parse_index(piece)
                    for piece in _split_top_level(index_list[1:-1])
                ),
            )
        )

        pos = close_index + 1

    return tuple(calls)


def _find_coupling(term: str) -> tuple[str, bool]:
    starts: list[int] = []
    pos = 0
    while True:
        start = term.find("Coupling[", pos)
        if start < 0:
            break
        starts.append(start)
        pos = start + 1

    if len(starts) != 1:
        raise ValueError(
            f"Expected exactly one quartic coupling, found {len(starts)}."
        )

    start = starts[0]
    open_index = start + len("Coupling")
    close_index = _matching_bracket(term, open_index)
    args = _split_top_level(term[open_index + 1 : close_index])
    name = args[0].strip()

    conjugated = False
    if start >= 4 and term[start - 4 : start] == "Bar[":
        conjugated = (
            _matching_bracket(term, start - 1) == close_index + 1
        )

    return name, conjugated


def _object_spans(term: str) -> list[tuple[int, int]]:
    """Spans of Field, CG and Coupling objects for prefactor extraction."""

    spans: list[tuple[int, int]] = []

    for marker in ("Field[", "CG[", "Coupling["):
        pos = 0
        while True:
            start = term.find(marker, pos)
            if start < 0:
                break

            open_index = start + marker.index("[")
            close_index = _matching_bracket(term, open_index)
            span_start = start
            span_end = close_index + 1

            if start >= 4 and term[start - 4 : start] == "Bar[":
                bar_open = start - 1
                if _matching_bracket(term, bar_open) == close_index + 1:
                    span_start = start - 4
                    span_end = close_index + 2

            spans.append((span_start, span_end))
            pos = close_index + 1

    # Nested objects such as Coupling inside Bar[...] can produce overlapping
    # spans only with their own wrapper; merge conservatively.
    spans.sort()
    merged: list[list[int]] = []
    for start, end in spans:
        if not merged or start >= merged[-1][1]:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)

    return [(start, end) for start, end in merged]


def _numeric_prefactor(term: str) -> sp.Expr:
    stripped = term
    for start, end in reversed(_object_spans(term)):
        stripped = stripped[:start] + "1" + stripped[end:]

    return _parse_mathematica_scalar(stripped)


def _cg_value(
    call: CGCall,
    registry: dict[str, CGTensor],
    assignment: dict[tuple[str, str], int],
) -> sp.Expr:
    if call.name not in registry:
        raise KeyError(f"CG {call.name!r} missing from quartic seed.")

    tensor = registry[call.name]

    if len(call.indices) != len(tensor.reps):
        raise ValueError(
            f"{call.name}: call rank and registered tensor rank disagree."
        )

    components: list[int] = []
    for dummy, rep in call.indices:
        key = (dummy, rep)
        if key not in assignment:
            candidates = [
                value
                for (label, _), value in assignment.items()
                if label == dummy
            ]
            if len(candidates) != 1:
                raise KeyError(
                    f"No unambiguous component assignment for {dummy}, {rep}."
                )
            components.append(candidates[0])
        else:
            components.append(assignment[key])

    value = tensor.tensor[tuple(component - 1 for component in components)]
    return sp.conjugate(value) if call.conjugated else value


def _real_expansion(
    leg: ScalarLeg,
    component: int,
    layout: ScalarLayout,
) -> tuple[tuple[int, sp.Expr], tuple[int, sp.Expr]]:
    root2 = sp.sqrt(2)

    real_index = layout.global_real_index(
        leg.name,
        component,
        imaginary=False,
    )
    imag_index = layout.global_real_index(
        leg.name,
        component,
        imaginary=True,
    )

    imag_factor = (
        -sp.I / root2 if leg.conjugated else sp.I / root2
    )

    return (
        (real_index, 1 / root2),
        (imag_index, imag_factor),
    )


def build_eft1_quartic_tensor(
    seed: dict,
    *,
    d_s1: int,
    d_s2: int,
) -> SparseQuarticTensor:
    """Build lambda_abcd from the exact Matchete scalar-quartic seed."""

    registry = _load_cg_registry(seed)
    layout = ScalarLayout(d_s1=int(d_s1), d_s2=int(d_s2))
    dimensions = layout.dimensions

    # Polynomial coefficient of each sorted real-field monomial in V4.
    polynomial: dict[tuple[int, int, int, int], sp.Expr] = {}

    for record in seed.get("ScalarQuarticTerms", []):
        term = record["TermInputForm"]
        legs = _find_scalar_legs(term)
        cg_calls = _find_cg_calls(term)
        coupling_name, coupling_conjugated = _find_coupling(term)

        coupling = sp.Symbol(coupling_name)
        if coupling_conjugated:
            coupling = sp.conjugate(coupling)

        # Exported expressions are Lagrangian terms.  The master-RGE quartic
        # convention is defined through V4 = -L4.
        prefactor = -_numeric_prefactor(term) * coupling

        # Repeated Mathematica dummy indices denote the same summed component.
        # Enumerate unique (dummy, representation) labels rather than four
        # field legs independently.
        unique_keys: list[tuple[str, str]] = []
        key_dimension: dict[tuple[str, str], int] = {}

        for leg in legs:
            # SU(2) singlets have no dummy representation index. Their only
            # complex component is fixed to component 1 and therefore does
            # not participate in the summed-index assignment.
            if leg.dummy is None:
                if dimensions[leg.name] != 1:
                    raise ValueError(
                        f"Missing SU(2) index for non-singlet {leg.name}."
                    )
                continue

            key = (leg.dummy, leg.rep)
            dimension = dimensions[leg.name]

            if key not in unique_keys:
                unique_keys.append(key)
                key_dimension[key] = dimension
            elif key_dimension[key] != dimension:
                raise ValueError(
                    f"Inconsistent dimensions for repeated index {key}."
                )

        ranges = [
            range(1, key_dimension[key] + 1)
            for key in unique_keys
        ]

        for values in product(*ranges):
            assignment = dict(zip(unique_keys, values))

            cg_factor = sp.S.One
            for call in cg_calls:
                cg_factor *= _cg_value(call, registry, assignment)
            cg_factor = sp.simplify(cg_factor)

            if cg_factor == 0:
                continue

            expansions = [
                _real_expansion(
                    leg,
                    (
                        1
                        if leg.dummy is None
                        else assignment[(leg.dummy, leg.rep)]
                    ),
                    layout,
                )
                for leg in legs
            ]

            for choices in product(*expansions):
                real_indices = tuple(index for index, _ in choices)
                basis_factor = sp.prod(factor for _, factor in choices)

                key = tuple(sorted(real_indices))
                polynomial[key] = sp.simplify(
                    polynomial.get(key, sp.S.Zero)
                    + prefactor * cg_factor * basis_factor
                )

    # If
    #
    #   V4 = c * phi_1^m1 phi_2^m2 ...
    #
    # and lambda_abcd is fully symmetric with V4 = lambda_abcd phi^4 / 4!,
    # then lambda for that multiset is c * m1! m2! ...
    entries: dict[tuple[int, int, int, int], sp.Expr] = {}

    for key, coefficient in polynomial.items():
        multiplicity_factor = sp.S.One
        for index in set(key):
            multiplicity_factor *= factorial(key.count(index))

        value = sp.simplify(coefficient * multiplicity_factor)
        if value != 0:
            entries[key] = value

    return SparseQuarticTensor(entries)


def load_and_build_eft1_quartic_tensor(
    quartic_seed_path: Path,
    rgbeta_path: Path,
) -> SparseQuarticTensor:
    seed = json.loads(Path(quartic_seed_path).read_text(encoding="utf-8"))
    rgbeta = json.loads(Path(rgbeta_path).read_text(encoding="utf-8"))
    metadata = rgbeta["metadata"]

    return build_eft1_quartic_tensor(
        seed,
        d_s1=int(metadata["dS1"]),
        d_s2=int(metadata["dS2"]),
    )


def higgs_block_diagnostics(tensor: SparseQuarticTensor) -> dict:
    """Validate the H-only block against the project's Higgs convention.

    The exported Matchete term is

        L4,H = -(lambda/2) (H† H)^2,

    hence

        V4,H = +(lambda/2) (H† H)^2
             = lambda/8 (phi_a phi_a)^2

    for the four real Higgs fields.  Therefore

        lambda_abcd =
            lambda (
                delta_ab delta_cd
              + delta_ac delta_bd
              + delta_ad delta_bc
            ).
    """

    lam = sp.Symbol("λ")
    failures: dict[tuple[int, int, int, int], sp.Expr] = {}

    for a, b, c, d in product(range(1, 5), repeat=4):
        expected = lam * (
            int(a == b) * int(c == d)
            + int(a == c) * int(b == d)
            + int(a == d) * int(b == c)
        )
        residual = sp.simplify(tensor[a, b, c, d] - expected)
        if residual != 0:
            failures[(a, b, c, d)] = residual

    return {
        "higgs_block_matches_known_convention": not failures,
        "higgs_block_failure_count": len(failures),
        "lambda_1111": str(tensor[1, 1, 1, 1]),
        "lambda_1122": str(tensor[1, 1, 2, 2]),
        "lambda_1133": str(tensor[1, 1, 3, 3]),
        "lambda_1234": str(tensor[1, 2, 3, 4]),
    }


def tensor_diagnostics(tensor: SparseQuarticTensor) -> dict:
    max_index = max(
        (max(key) for key in tensor.entries),
        default=0,
    )

    permutation_failures = 0
    for key, value in tensor.entries.items():
        a, b, c, d = key
        tests = (
            (b, a, c, d),
            (a, c, b, d),
            (d, c, b, a),
        )
        for test_key in tests:
            if sp.simplify(tensor[test_key] - value) != 0:
                permutation_failures += 1

    return {
        "nonzero_symmetric_components": len(tensor.entries),
        "largest_real_scalar_index": max_index,
        "permutation_failure_count": permutation_failures,
        "fully_symmetric_lookup": permutation_failures == 0,
        **higgs_block_diagnostics(tensor),
    }


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description=(
            "Convert an EFT1 Matchete scalar-quartic seed to real-basis "
            "lambda_abcd."
        )
    )
    parser.add_argument("quartic_seed_json", type=Path)
    parser.add_argument("rgbeta_json", type=Path)
    args = parser.parse_args()

    quartic = load_and_build_eft1_quartic_tensor(
        args.quartic_seed_json,
        args.rgbeta_json,
    )

    print(json.dumps(tensor_diagnostics(quartic), indent=2))

    for key, value in sorted(quartic.nonzero_items()):
        print(f"{key}: {value}")
