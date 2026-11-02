"""Ordered, numerical Matchete LLSS -> real-scalar component expansion.

This module is intentionally NOT the symmetric canonical Wilson C_ijab:
  - the Matchete lepton/Yukawa order is retained per source term;
  - 3x3 flavour coefficients are kept on each ordered SU(2)/scalar component;
  - no automatic 1/2, 1/4, charge-conjugation or EFT1 coupling-map factor
    is introduced.

It handles original S1/S2, the six self-bilinear real fields, and the B/C
X/Y complex-doublet basis. All CG numbers are read from a Matchete CGRegistry.
Unknown fields/couplings/CGs fail closed. The scalar EFT1 anomalous dimension
is not supplied by this module.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
import re
from typing import Any, Mapping

import numpy as np
import sympy as sp


REPRESENTATIONS = {
    "A": (1, 3, 2), "B": (2, 2, 1), "C": (2, 2, 3),
    "D": (3, 1, 2), "E": (3, 3, 2),
}
MIXED = frozenset(("B_-1", "C_-1"))
SELF_BILINEAR = frozenset(("A_-2", "A_+0", "D_-2", "D_+0", "E_-2", "E_+0"))
_FIELD_START = "Field["


@dataclass(frozen=True)
class ScalarLeg:
    name: str
    dummy: str | None
    rep: str | None
    barred: bool


@dataclass(frozen=True)
class RealScalarBasis:
    """Global real-scalar indices 1-4 = four real Higgs components."""
    model: str
    dimensions: Mapping[str, int]
    real_fields: frozenset[str]
    offsets: Mapping[str, int]
    total_dimension: int

    def expansion(self, leg: ScalarLeg, component: int) -> tuple[tuple[int, complex], ...]:
        name = leg.name
        if name not in self.dimensions:
            raise ValueError(f"Unrecognised scalar leg {name}")
        dim = self.dimensions[name]
        if not isinstance(component, int) or not 1 <= component <= dim:
            raise ValueError(f"Invalid scalar component {name}[{component}]")
        base = self.offsets[name]
        if name in self.real_fields:
            return ((base + component - 1, 1.0 + 0j),)
        k = base + 2 * (component - 1)
        return ((k, complex(1 / math.sqrt(2))),
                (k + 1, complex(-1j if leg.barred else 1j) / math.sqrt(2)))


def scalar_basis(model: str) -> RealScalarBasis:
    family = model.split("_")[0]
    if family not in REPRESENTATIONS:
        raise ValueError(f"Unknown model {model}")
    d1, d2, _ = REPRESENTATIONS[family]
    if model in MIXED:
        fields = (("NewScalarX", 2, False), ("NewScalarY", 2, False))
    elif model in SELF_BILINEAR:
        if model.endswith("_-2"):
            fields = (("NewScalar1", d1, False),
                      ("T3Scalar2Minus", d2, True), ("T3Scalar2Plus", d2, True))
        else:
            fields = (("T3Scalar1Minus", d1, True), ("T3Scalar1Plus", d1, True),
                      ("NewScalar2", d2, False))
    else:
        fields = (("NewScalar1", d1, False), ("NewScalar2", d2, False))
    next_index = 5
    offsets = {}
    for name, dimension, real in fields:
        offsets[name] = next_index
        next_index += dimension if real else 2 * dimension
    return RealScalarBasis(model, {x: d for x, d, _ in fields},
                           frozenset(x for x, _, real in fields if real),
                           offsets, next_index - 1)


def _closing(text: str, start: int, left: str = "[", right: str = "]") -> int:
    if start < 0 or text[start] != left:
        raise ValueError("Malformed Matchete bracket")
    depth = 0
    for pos in range(start, len(text)):
        if text[pos] == left:
            depth += 1
        elif text[pos] == right:
            depth -= 1
            if depth == 0:
                return pos
    raise ValueError("Unclosed Matchete bracket")


def _split(text: str) -> list[str]:
    parts, last, sq, curly, paren = [], 0, 0, 0, 0
    for i, ch in enumerate(text):
        if ch == "[": sq += 1
        elif ch == "]": sq -= 1
        elif ch == "{": curly += 1
        elif ch == "}": curly -= 1
        elif ch == "(": paren += 1
        elif ch == ")": paren -= 1
        elif ch == "," and sq == curly == paren == 0:
            parts.append(text[last:i].strip()); last = i + 1
        if min(sq, curly, paren) < 0:
            raise ValueError("Unbalanced Matchete expression")
    if (sq, curly, paren) != (0, 0, 0):
        raise ValueError("Unbalanced Matchete expression")
    parts.append(text[last:].strip())
    return parts


def scalar_legs(term: str, basis: RealScalarBasis) -> tuple[ScalarLeg, ScalarLeg]:
    """Read exact light-scalar leg names and SU(2) index labels."""
    result = []
    index = 0
    while True:
        pos = term.find(_FIELD_START, index)
        if pos < 0:
            break
        end = _closing(term, pos + len("Field"))
        args = _split(term[pos + len(_FIELD_START):end])
        index = end + 1
        if len(args) < 3 or args[1] != "Scalar":
            continue
        name = args[0]
        if name not in basis.dimensions:
            raise ValueError(f"Unexpected scalar field {name} (model {basis.model})")
        is_barred = False
        if term[max(pos - 4, 0):pos] == "Bar[":
            is_barred = _closing(term, pos - 1) == end + 1
        ix = args[2]
        if not (ix.startswith("{") and ix.endswith("}")):
            raise ValueError("Malformed scalar index list")
        indices = _split(ix[1:-1]) if ix[1:-1].strip() else []
        if len(indices) > 1:
            raise ValueError("Expected at most one SU(2) scalar index")
        if indices:
            raw = indices[0].strip()
            if not raw.startswith("Index[") or not raw.endswith("]"):
                raise ValueError(f"Invalid scalar index {raw}")
            pieces = _split(raw[len("Index["):-1])
            if len(pieces) != 2:
                raise ValueError(f"Invalid scalar index {raw}")
            dummy, rep = pieces
            if rep.startswith("Bar[") and rep.endswith("]"):
                rep = rep[4:-1]
        else:
            dummy = rep = None
        result.append(ScalarLeg(name, dummy, rep, is_barred))
    if len(result) != 2:
        raise ValueError(f"Expected 2 light-scalar legs, found {len(result)}")
    for leg in result:
        if basis.dimensions[leg.name] == 1 and leg.dummy is not None:
            raise ValueError("SU(2) singlet unexpectedly indexed")
        if basis.dimensions[leg.name] > 1 and leg.dummy is None:
            raise ValueError("Non-singlet scalar missing SU(2) index")
    return result[0], result[1]


def _number(z: Any) -> complex:
    if isinstance(z, (float, int, complex, np.number)):
        out = complex(z)
    else:
        out = complex(sp.N(z, 25))
    if not np.isfinite(out.real) or not np.isfinite(out.imag):
        raise ValueError("Nonfinite CG/prefactor value")
    return out


def _matrix(raw: Mapping[str, Any]) -> np.ndarray:
    r, i = raw["real"], raw["imag"]
    result = np.asarray(r, float) + 1j * np.asarray(i, float)
    if result.shape != (3, 3) or not np.isfinite(result).all():
        raise ValueError("Expected finite ordered 3x3 flavour matrix")
    return result


def _encode_matrix(mat: np.ndarray) -> dict:
    z = np.asarray(mat, dtype=complex)
    if z.shape != (3, 3) or not np.isfinite(z).all():
        raise ValueError("Invalid tensor coefficient matrix")
    return {"real": z.real.tolist(), "imag": z.imag.tolist()}


def project_ordered_term(term: Mapping[str, Any], registry: list[dict],
                         model: str, *, max_entries: int = 3000) -> dict:
    """Exact scalar/CG expansion of an already source-factorised PL term.

    The external flavour matrix is kept ordered. A tensor entry is a 3x3
    matrix rather than a one-generation symbolic Yukawa monomial.
    """
    from RGE.running.intermediate.ScalarOnlyTensorAdapters import (
        load_wilson_cg_registry, _find_wilson_cg_calls,
        _find_lepton_dummies, _wilson_cg_value, _term_scalar_expression,
    )
    from RGE.running.intermediate.IntermediateMatcheteParsing import parse_mathematica_scalar

    basis = scalar_basis(model)
    source = term["operator_template_unit_flavor_inputform"]
    if "Coupling[" in source:
        raise ValueError("Unfactorised coupling remains in template")
    if "GammaCC" not in source:
        raise ValueError("Not a Majorana-type Matchete source")
    legs = scalar_legs(source, basis)
    order = _find_lepton_dummies(source)
    calls = _find_wilson_cg_calls(source)
    seed = {"CGRegistry": registry}
    cg_registry = load_wilson_cg_registry(seed)
    if any(c.name not in cg_registry for c in calls):
        raise ValueError("Missing exact source CG tensor")
    template_for_prefactor = {"Index": term.get("source_index", 0),
                              "TermInputForm": source}
    prefactor_raw = _term_scalar_expression(template_for_prefactor)
    prefactor = _number(parse_mathematica_scalar(prefactor_raw))
    w = _matrix(term["ordered_flavour_matrix_gev_inv"])
    # A repeated Matchete dummy denotes ONE summed SU(2) component.
    # Its occurrences must agree; they must not be treated as independent
    # scalar indices (nor counted more than once in the contraction).
    dimensions_by_dummy = {(order[0], "SU2L[fund]"): 2,
                           (order[1], "SU2L[fund]"): 2}
    for leg in legs:
        if leg.dummy is None:
            continue
        key = (leg.dummy, leg.rep)
        dimension = basis.dimensions[leg.name]
        if key in dimensions_by_dummy and dimensions_by_dummy[key] != dimension:
            raise ValueError(f"Inconsistent dimensions for repeated SU(2) index {key}")
        dimensions_by_dummy[key] = dimension

    sparse = {}
    for i in (1, 2):
        for j in (1, 2):
            for aa in range(1, basis.dimensions[legs[0].name] + 1):
                for bb in range(1, basis.dimensions[legs[1].name] + 1):
                    assignment = {(order[0], "SU2L[fund]"): i,
                                  (order[1], "SU2L[fund]"): j}
                    consistent = True
                    for leg, comp in ((legs[0], aa), (legs[1], bb)):
                        if leg.dummy is not None:
                            key = (leg.dummy, leg.rep)
                            if key in assignment and assignment[key] != comp:
                                consistent = False
                                break
                            assignment[key] = comp
                    if not consistent:
                        # This trial tuple violates the repeated-index equality;
                        # it contributes zero rather than invalidating the term.
                        continue
                    cg = sp.S.One
                    for call in calls:
                        cg *= _wilson_cg_value(call, cg_registry, assignment)
                    zcg = _number(cg)
                    if abs(zcg) == 0:
                        continue
                    for a, ca in basis.expansion(legs[0], aa):
                        for b, cb in basis.expansion(legs[1], bb):
                            key = (i, j, a, b)
                            sparse[key] = sparse.get(key, 0j) + prefactor * zcg * ca * cb
    if len(sparse) > max_entries:
        raise ValueError("Tensor operator unexpectedly exceeds supported size")
    entries = [{"lepton_su2": [i, j], "scalar_real": [a, b],
                "ordered_flavour_coefficient_gev_inv": _encode_matrix(w * z)}
               for (i, j, a, b), z in sorted(sparse.items()) if abs(z) > 1e-15]
    if not entries:
        raise ValueError("No nonzero real-component Wilson terms")
    return {
        "source_index": term["source_index"],
        "basis_scalar_fields": dict(basis.dimensions),
        "real_scalar_dimension_including_higgs": basis.total_dimension,
        "flavour_matrix_axes_external_indices": term.get("flavour_matrix_axes_external_indices"),
        "ordered_not_canonical_symmetrized": True,
        "cg_registry_used": sorted({c.name for c in calls}),
        "component_count": len(entries), "components": entries,
    }
