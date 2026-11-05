"""Exact real-field adapters for the D(-2) fermion-first intermediate EFT.

Read the matched Matchete InputForm and same-kernel CG tensor registries.
Scalar ordering: (H_R,H_I) interleaved (1..4), S1 triplet real
coordinates (5..10), the PHYSICAL neutral singlets T3S2R=11, T3S2I=12.

This is not the old complex-NewScalar2 adapter.  The two neutral singlets
are treated as independent REAL field components, with no 1/sqrt(2) factor
added a second time. Phase spurions are kept symbolic and unnormalised.
"""
from __future__ import annotations

from dataclasses import dataclass
from itertools import product
from math import factorial
import re
from typing import Any

import sympy as sp

FIELD_DIMS = {"H": 2, "NewScalar1": 3, "T3S2R": 1, "T3S2I": 1}
REAL_SINGLET_INDICES = {"T3S2R": 11, "T3S2I": 12}
PHASE = "T3S2Phase"
YUKAWAS = {"y1", "y2", "y2Conj"}


@dataclass(frozen=True)
class FieldLeg:
    name: str
    dummy: str | None
    rep: str | None
    barred: bool


def bracket_end(source: str, open_at: int) -> int:
    if source[open_at] != "[":
        raise ValueError("Expected opening bracket")
    count = 0
    for at in range(open_at, len(source)):
        if source[at] == "[":
            count += 1
        elif source[at] == "]":
            count -= 1
            if count == 0:
                return at
    raise ValueError("Unbalanced Mathematica brackets")


def split_args(text: str) -> list[str]:
    result: list[str] = []
    start = 0
    sq = brace = par = 0
    for i, ch in enumerate(text):
        if ch == "[": sq += 1
        elif ch == "]": sq -= 1
        elif ch == "{": brace += 1
        elif ch == "}": brace -= 1
        elif ch == "(": par += 1
        elif ch == ")": par -= 1
        elif ch == "," and not (sq or brace or par):
            result.append(text[start:i].strip())
            start = i + 1
    result.append(text[start:].strip())
    return result


def unwrap_bar(text: str) -> tuple[str, bool]:
    text = text.strip()
    if text.startswith("Bar[") and bracket_end(text, 3) == len(text) - 1:
        return text[4:-1], True
    return text, False


def read_index(item: str) -> tuple[str, str]:
    item, _ = unwrap_bar(item)
    if not item.startswith("Index[") or bracket_end(item, 5) != len(item)-1:
        raise ValueError(f"Not an Index[...]: {item[:100]}")
    parts = split_args(item[6:-1])
    if len(parts) != 2:
        raise ValueError("Expected two Index arguments")
    rep, _ = unwrap_bar(parts[1])
    return parts[0], rep


def scalar_legs(text: str, count: int) -> tuple[FieldLeg, ...]:
    """Find full physical Field[...,Scalar,...] factors, including Bar[Field]."""
    legs: list[FieldLeg] = []
    pos = 0
    while (at := text.find("Field[", pos)) != -1:
        stop = bracket_end(text, at + 5)
        args = split_args(text[at + 6:stop])
        if len(args) >= 3 and args[1] == "Scalar":
            name = args[0]
            if name not in FIELD_DIMS:
                raise ValueError(f"Unknown scalar in matched EFT1: {name}")
            inds = args[2]
            if not inds.startswith("{") or not inds.endswith("}"):
                raise ValueError("Non-list scalar indices")
            contents = inds[1:-1].strip()
            if contents:
                index_list = split_args(contents)
                if len(index_list) != 1:
                    raise ValueError("Expected at most one scalar SU(2) index")
                dummy, rep = read_index(index_list[0])
            else:
                dummy = rep = None
            barred = False
            if at >= 4 and text[at-4:at] == "Bar[":
                barred = bracket_end(text, at-1) == stop+1
            if name in REAL_SINGLET_INDICES and dummy is not None:
                raise ValueError("Real D(-2) singlet incorrectly has SU(2) index")
            if name in ("H", "NewScalar1") and dummy is None:
                raise ValueError("Missing gauge index for complex multiplet")
            # Matchete prints repeated identical fields as powers, e.g.
            # Field[T3S2R, Scalar, {}, {}]^4.  One textual Field is then
            # FOUR tensor legs.  Likewise for Bar[Field[...]]^n.
            outer_stop = stop + (1 if barred else 0)
            power_match = re.match(r"\^([0-9]+)", text[outer_stop + 1:])
            multiplicity = int(power_match.group(1)) if power_match else 1
            if not 1 <= multiplicity <= 4:
                raise ValueError(f"Invalid physical scalar monomial power {multiplicity}")
            legs.extend([FieldLeg(name, dummy, rep, barred)] * multiplicity)
        pos = stop + 1
    if len(legs) != count:
        raise ValueError(f"Expected {count} scalar legs, found {len(legs)}")
    return tuple(legs)


def scalar_coordinates(leg: FieldLeg, component: int) -> tuple[tuple[int, sp.Expr], ...]:
    """Expand one physical scalar field in the RGE's real-coordinate basis."""
    if leg.name in REAL_SINGLET_INDICES:
        if component != 1:
            raise ValueError("Real singlet has only component 1")
        return ((REAL_SINGLET_INDICES[leg.name], sp.S.One),)
    d = FIELD_DIMS[leg.name]
    if not (1 <= component <= d):
        raise ValueError("Scalar component outside SU(2) irrep")
    first = 1 if leg.name == "H" else 5
    r = first + 2*(component-1)
    root2 = sp.sqrt(2)
    return ((r, 1/root2), (r+1, (-sp.I if leg.barred else sp.I)/root2))


def _coupling_name(text: str) -> str:
    if not text.startswith("Coupling[") or bracket_end(text, 8) != len(text)-1:
        raise ValueError(f"Unrecognized coupling term: {text[:100]}")
    args = split_args(text[9:-1])
    if len(args) < 1 or not re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*|[\w\\\[\]]+", args[0]):
        raise ValueError(f"Bad Matchete coupling name: {args[0] if args else text}")
    name = args[0]
    # SymPy cannot read Mathematica's explicit escaped symbol spelling.
    return {r"\[Lambda]": "lambdaSM", "λ": "lambdaSM"}.get(name, name)


def scalar_prefactor(text: str) -> tuple[sp.Expr, tuple[str, ...]]:
    """Remove Field/CG/NCM, retaining every Coupling, Bar and phase exactly.

    Reject any other unparsed Mathematica function instead of silently
    assigning a numerical normalization to an unknown term.
    """
    out: list[str] = []
    seen: list[str] = []
    at = 0
    while at < len(text):
        consumed = False
        for prefix, placeholder, is_bar in (
            ("Bar[Coupling[", "coupling", True),
            ("Coupling[", "coupling", False),
            ("Bar[CG[", "one", True),
            ("CG[", "one", False),
            ("Bar[Field[", "one", True),
            ("Field[", "one", False),
            ("NCM[", "one", False),
        ):
            if not text.startswith(prefix, at):
                continue
            outer_open = at + (3 if is_bar else prefix.index("["))
            stop = bracket_end(text, outer_open)
            if is_bar:
                inside = text[outer_open+1:stop]
                # Bar wraps the complete object, not just its first argument.
                if inside.startswith("Coupling["):
                    name = _coupling_name(inside)
                    out.append(f"conjugate({name})")
                    seen.append(name)
                else:
                    out.append("1")
            else:
                if placeholder == "coupling":
                    name = _coupling_name(text[at:stop+1])
                    seen.append(name)
                    out.append(name)
                else:
                    out.append("1")
            at = stop + 1
            consumed = True
            break
        if not consumed:
            out.append(text[at])
            at += 1
    cleaned = "".join(out)
    if "Field[" in cleaned or "CG[" in cleaned or "Coupling[" in cleaned:
        raise ValueError("Unremoved Matchete object in prefactor")
    if "T3S2Phase[" in cleaned:
        raise ValueError("Phase not converted to a registered Matchete Coupling")
    while (match := re.search(r"Sqrt\[([^\[\]]+)\]", cleaned)):
        cleaned = cleaned[:match.start()] + f"sqrt({match.group(1)})" + cleaned[match.end():]
    cleaned = cleaned.replace("^", "**")
    try:
        value = sp.sympify(cleaned, locals={"conjugate": sp.conjugate, "sqrt": sp.sqrt, "I": sp.I})
    except (ValueError, TypeError, SyntaxError, sp.SympifyError) as exc:
        raise ValueError(f"Cannot parse physical EFT1 coefficient: {cleaned[:160]}") from exc
    return sp.simplify(value), tuple(seen)


def _expanded_cg_calls(text: str, calls: tuple | list) -> tuple:
    """Retain multiplicities of Matchete CG[...]^n factors.

    The project CG parsers return each textual CG occurrence once. Powers of
    contractions therefore need expansion just like powers of scalar fields.
    An outer Bar[CG[...]] conjugates the entire CG (as opposed to CG[Bar[name],...]).
    """
    from dataclasses import replace

    positions: list[tuple[int, int]] = []
    at = 0
    while (start := text.find("CG[", at)) != -1:
        stop = bracket_end(text, start + 2)
        outer_bar = (
            start >= 4 and text[start - 4:start] == "Bar["
            and stop + 1 < len(text) and text[stop + 1] == "]"
        )
        last = stop + (1 if outer_bar else 0)
        match = re.match(r"\^([0-9]+)", text[last + 1:])
        multiplicity = int(match.group(1)) if match else 1
        if not 1 <= multiplicity <= 4:
            raise ValueError(f"Invalid CG tensor power {multiplicity}")
        positions.append((multiplicity, outer_bar))
        at = stop + 1
    if len(positions) != len(calls):
        raise ValueError("CG parser and physical-basis CG multiplicity disagree")
    out = []
    for call, (multiplicity, outer_bar) in zip(calls, positions):
        adjusted = replace(call, conjugated=not call.conjugated) if outer_bar else call
        out.extend([adjusted] * multiplicity)
    return tuple(out)


def _l_indices(term: str) -> tuple[str, str]:
    hits = re.findall(
        r"(?:Bar\[)?Field\[l,\s*Fermion,\s*\{Index\[([^,\]]+),\s*SU2L\[fund\]\]",
        term,
    )
    dummies = list(dict.fromkeys(hits))
    if len(dummies) != 2:
        raise ValueError(f"Two distinct external lepton SU(2) dummies required; got {dummies}")
    return dummies[0], dummies[1]


def _assignment_ranges(legs: tuple[FieldLeg, ...]) -> tuple[list[tuple[str,str]], list[range]]:
    keys: list[tuple[str, str]] = []
    ranges: list[range] = []
    for leg in legs:
        if leg.dummy is None:
            continue
        key = (leg.dummy, leg.rep)
        if key not in keys:
            keys.append(key)
            ranges.append(range(1, FIELD_DIMS[leg.name] + 1))
        elif FIELD_DIMS[leg.name] != len(ranges[keys.index(key)]):
            raise ValueError(f"Inconsistent SU(2) component dimensions for {key}")
    return keys, ranges


def build_real_quartic(seed: dict[str, Any]):
    """Return project SparseQuarticTensor with V4=lambda_abcd phi^4/4!."""
    from RGE.running.intermediate.ScalarOnlyTensorAdapters import (
        SparseQuarticTensor, _find_cg_calls, _cg_value,
    )
    from RGE.running.intermediate.IntermediateMatcheteParsing import load_cg_registry
    registry = load_cg_registry(seed)
    poly: dict[tuple[int, int, int, int], sp.Expr] = {}
    terms = seed.get("ScalarQuarticTerms", [])
    if not terms or seed.get("ScalarQuarticTermCount") != len(terms):
        raise ValueError("Missing or inconsistent physical scalar quartic seed")
    for record in terms:
        term = record["TermInputForm"]
        legs = scalar_legs(term, 4)
        cg_calls = _expanded_cg_calls(term, _find_cg_calls(term))
        prefactor, names = scalar_prefactor(term)
        nonphase = [name for name in names if name != PHASE]
        if len(nonphase) != 1:
            raise ValueError(f"Quartic term {record['Index']} must contain exactly one potential coupling; got {names}")
        keys, ranges = _assignment_ranges(legs)
        for values in product(*ranges):
            assign = dict(zip(keys, values))
            cg = sp.S.One
            for call in cg_calls:
                cg *= _cg_value(call, registry, assign)
            if cg == 0:
                continue
            leg_expansions = [
                scalar_coordinates(leg, 1 if leg.dummy is None else assign[(leg.dummy, leg.rep)])
                for leg in legs
            ]
            for opts in product(*leg_expansions):
                key = tuple(sorted(index for index, _ in opts))
                # Quartic InputForm is from L_4, while beta uses V_4=-L_4.
                contribution = -prefactor * cg * sp.prod(factor for _, factor in opts)
                poly[key] = poly.get(key, sp.S.Zero) + contribution
    entries = {}
    for key, val in poly.items():
        v = sp.simplify(val * sp.prod(factorial(key.count(x)) for x in set(key)))
        if v != 0:
            entries[key] = v
    return SparseQuarticTensor(entries)


def build_real_wilson(seed: dict[str, Any], lepton_indices: tuple[int, int]):
    """Return projected C_ijab for the matched PL two-lepton operators."""
    from RGE.running.intermediate.ScalarOnlyTensorAdapters import (
        SparseWilsonTensor, _find_wilson_cg_calls, _wilson_cg_value,
        load_wilson_cg_registry,
    )
    registry = load_wilson_cg_registry(seed)
    terms = seed.get("TreeWilsonTerms", [])
    if not terms or seed.get("TreeWilsonTermCount") != len(terms):
        raise ValueError("Missing or inconsistent real-basis Wilson seed")
    if any(rec.get("Chirality") not in {"PL", "PR"} for rec in terms):
        raise ValueError("Unclassified chirality in D(-2) Wilson seed")
    selected = [rec for rec in terms if rec["Chirality"] == "PL"]
    if not selected:
        raise ValueError("No left-chiral matched D(-2) Wilson seed")
    entries: dict[tuple[int,int,int,int], sp.Expr] = {}
    used_yukawas: set[str] = set()
    for rec in selected:
        term = rec["TermInputForm"]
        legs = scalar_legs(term, 2)
        labels = _l_indices(term)
        cg_calls = _expanded_cg_calls(term, _find_wilson_cg_calls(term))
        pref, couplings = scalar_prefactor(term)
        yukawas = [name for name in couplings if name in YUKAWAS]
        if len(yukawas) != 2 or any(x not in YUKAWAS|{PHASE,"MF"} for x in couplings):
            raise ValueError(f"Unrecognised PL Wilson prefactor in term {rec['Index']}: {couplings}")
        mf = sp.Symbol("MF")
        if not pref.has(mf) or sp.simplify(pref*mf).has(mf):
            raise ValueError("Matched Wilson term does not scale as precisely 1/MF")
        used_yukawas.update(yukawas)
        keys, ranges = _assignment_ranges(legs)
        for l1 in (1, 2):
            for l2 in (1, 2):
                for vals in product(*ranges):
                    assign = { (labels[0], "SU2L[fund]"): l1,
                               (labels[1], "SU2L[fund]"): l2 }
                    for key, value in zip(keys, vals):
                        if key in assign and assign[key] != value:
                            raise ValueError("A lepton SU(2) index reused inconsistently")
                        assign[key] = value
                    cg = sp.S.One
                    for call in cg_calls:
                        cg *= _wilson_cg_value(call, registry, assign)
                    if cg == 0:
                        continue
                    expansions = [scalar_coordinates(leg, 1 if leg.dummy is None else assign[(leg.dummy, leg.rep)]) for leg in legs]
                    for opts in product(*expansions):
                        i, j = lepton_indices[l1-1], lepton_indices[l2-1]
                        a, b = opts[0][0], opts[1][0]
                        key = (i,j,a,b)
                        entries[key] = entries.get(key,sp.S.Zero) + pref*cg*opts[0][1]*opts[1][1]
    return SparseWilsonTensor(entries).symmetrized(), used_yukawas
