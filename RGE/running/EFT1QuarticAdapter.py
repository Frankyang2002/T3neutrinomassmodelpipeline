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

import sympy as sp

from RGE.running.MatcheteParsing import (
    CGTensor,
    load_cg_registry,
    matching_brace as _matching_brace,
    matching_bracket as _matching_bracket,
    parse_index as _parse_index,
    parse_mathematica_scalar as _parse_mathematica_scalar,
    parse_reps as _parse_reps,
    parse_sparse_array as _parse_sparse_array,
    split_top_level as _split_top_level,
    strip_bar as _strip_bar,
)


@dataclass(frozen=True)
class ScalarLeg:
    name: str
    dummy: str | None
    rep: str | None
    conjugated: bool


@dataclass(frozen=True)
class CGCall:
    name: str
    conjugated: bool
    indices: tuple[tuple[str, str], ...]


@dataclass(frozen=True)
class ScalarLayout:
    d_s1: int
    d_s2: int
    shared_scalar: bool = False

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
            first = 5 if self.shared_scalar else 5 + 2 * self.d_s1
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


def _shared_scalar_component_map(
    leg: ScalarLeg,
    component: int,
    dimension: int,
) -> tuple[int, bool, sp.Expr]:
    if leg.name == "NewScalar2":
        return component, leg.conjugated, sp.S.One
    if leg.name != "NewScalar1":
        return component, leg.conjugated, sp.S.One
    physical_component = dimension + 1 - component
    phase = sp.Integer(-1) ** (component - 1)
    return physical_component, (not leg.conjugated), phase


def _real_expansion(
    leg: ScalarLeg,
    component: int,
    layout: ScalarLayout,
) -> tuple[tuple[int, sp.Expr], tuple[int, sp.Expr]]:
    root2 = sp.sqrt(2)
    physical_component = component
    conjugated = leg.conjugated
    phase = sp.S.One
    if layout.shared_scalar and leg.name in {"NewScalar1", "NewScalar2"}:
        physical_component, conjugated, phase = _shared_scalar_component_map(
            leg, component, layout.d_s1
        )
    real_index = layout.global_real_index(leg.name, physical_component, imaginary=False)
    imag_index = layout.global_real_index(leg.name, physical_component, imaginary=True)
    imag_factor = -sp.I / root2 if conjugated else sp.I / root2
    return (
        (real_index, sp.simplify(phase / root2)),
        (imag_index, sp.simplify(phase * imag_factor)),
    )


def build_eft1_quartic_tensor(
    seed: dict,
    *,
    d_s1: int,
    d_s2: int,
    shared_scalar: bool = False,
) -> SparseQuarticTensor:
    """Build lambda_abcd from the exact Matchete scalar-quartic seed."""

    registry = load_cg_registry(seed)
    layout = ScalarLayout(
        d_s1=int(d_s1), d_s2=int(d_s2), shared_scalar=bool(shared_scalar)
    )
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
        shared_scalar=bool(metadata.get("SharedScalar", False)),
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
