from __future__ import annotations

"""
We get our python objects of matchete EFT operators
and we convert them into RGE tensors

Specifically
Wilson C_ijab for LLSS
Scalar potential lambda_abcd
"""

from dataclasses import dataclass
from itertools import product
import json
from math import factorial
from pathlib import Path
import re
import sys

import sympy as sp

PROJECT_ROOT = Path(__file__).resolve().parents[3]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from RGE.running.eft1.MatcheteParsing import (
    CGTensor,
    load_cg_registry,
    matching_bracket as _matching_bracket,
    parse_index as _parse_index,
    parse_mathematica_scalar as _parse_mathematica_scalar,
    parse_reps as _parse_reps,
    parse_sparse_array as _parse_sparse_array,
    split_top_level as _split_top_level,
    strip_bar as _strip_bar,
)


# ---------------------------------------------------------------------------
# Wilson-seed -> C_ijab adapter
# ---------------------------------------------------------------------------

class WilsonCGTensor:
    name: str
    reps: tuple[str, ...]
    matrix: sp.MutableDenseNDimArray


@dataclass(frozen=True)
class WilsonScalarLeg:
    name: str
    dummy: str | None
    rep: str | None
    conjugated: bool


@dataclass(frozen=True)
class WilsonCGCall:
    name: str
    conjugated: bool
    indices: tuple[tuple[str, str], ...]


@dataclass(frozen=True)
class WilsonScalarLayout:
    """We find where EFT1 scalars exist in our global real basis
    like where 1234 -> H"""

    d_s1: int
    d_s2: int
    include_higgs: bool = True
    shared_scalar: bool = False

    # Higgs takes up 1-4, so if we have higgs the s1 is at 5
    @property
    def s1_first(self) -> int:
        return 5 if self.include_higgs else 1

    @property
    def s2_first(self) -> int:
        if self.shared_scalar:
            return self.s1_first
        return self.s1_first + 2 * self.d_s1

    def global_real_index(
        self,
        scalar_name: str,
        complex_component: int,
        *,
        imaginary: bool,
    ) -> int:
        if scalar_name == "S1":
            first = self.s1_first
            dimension = self.d_s1
        elif scalar_name == "S2":
            first = self.s2_first
            dimension = self.d_s2
        else:
            raise KeyError(f"Unknown scalar {scalar_name!r}.")

        if not 1 <= complex_component <= dimension:
            raise IndexError(
                f"{scalar_name} component {complex_component} outside 1,...,{dimension}."
            )

        local = 2 * complex_component - (0 if imaginary else 1)
        return first + local - 1


class SparseWilsonTensor:
    """Sparse C_ijab lookup compatible with WilsonTensorRGE call sites."""

    def __init__(self, entries: dict[tuple[int, int, int, int], sp.Expr]):
        self.entries = {
            tuple(map(int, key)): sp.simplify(value)
            for key, value in entries.items()
            if sp.simplify(value) != 0
        }

    def __getitem__(self, key: tuple[int, int, int, int]) -> sp.Expr:
        return self.entries.get(tuple(key), sp.S.Zero)

    def __call__(self, i: int, j: int, a: int, b: int) -> sp.Expr:
        return self[i, j, a, b]

    def nonzero_items(self):
        return self.entries.items()


    def symmetrized(self) -> "SparseWilsonTensor":
        """We have ij and ab symmetry, we have it symmetrised as 
            Csym_ijab = 1/4 (
                C_ijab + C_jiab + C_ijba + C_jiba
            ).
            So antisymmetrised components cancel
        """

        keys = set(self.entries)
        for i, j, a, b in tuple(keys):
            keys.update(
                {
                    (j, i, a, b),
                    (i, j, b, a),
                    (j, i, b, a),
                }
            )

        result: dict[tuple[int, int, int, int], sp.Expr] = {}
        for i, j, a, b in keys:
            value = sp.simplify(
                (
                    self[i, j, a, b]
                    + self[j, i, a, b]
                    + self[i, j, b, a]
                    + self[j, i, b, a]
                )
                / 4
            )
            if value != 0:
                result[(i, j, a, b)] = value

        return SparseWilsonTensor(result)

    def symmetry_residuals(self) -> dict[str, dict[tuple[int, int, int, int], sp.Expr]]:
        """Return exact residuals for the two operator-index symmetries.
        We can check violations of our symmetries here"""

        keys = set(self.entries)
        for i, j, a, b in tuple(keys):
            keys.update({(j, i, a, b), (i, j, b, a)})

        fermion: dict[tuple[int, int, int, int], sp.Expr] = {}
        scalar: dict[tuple[int, int, int, int], sp.Expr] = {}

        for i, j, a, b in keys:
            rf = sp.simplify(self[i, j, a, b] - self[j, i, a, b])
            rs = sp.simplify(self[i, j, a, b] - self[i, j, b, a])
            if rf != 0:
                fermion[(i, j, a, b)] = rf
            if rs != 0:
                scalar[(i, j, a, b)] = rs

        return {
            "fermion_swap": fermion,
            "scalar_swap": scalar,
        }


def load_wilson_cg_registry(seed: dict) -> dict[str, WilsonCGTensor]:
    registry: dict[str, WilsonCGTensor] = {}
    for item in seed.get("CGRegistry", []):
        registry[item["Name"]] = WilsonCGTensor(
            name=item["Name"],
            reps=_parse_reps(item["RepsInputForm"]),
            matrix=_parse_sparse_array(item["TensorInputForm"]),
        )
    return registry


def _find_wilson_cg_calls(term: str) -> tuple[WilsonCGCall, ...]:
    calls: list[WilsonCGCall] = []
    pos = 0
    while True:
        start = term.find("CG[", pos)
        if start < 0:
            break
        open_index = start + 2
        close_index = _matching_bracket(term, open_index)
        inside = term[open_index + 1 : close_index]
        args = _split_top_level(inside)
        if len(args) != 2:
            raise ValueError(f"Unexpected CG call: {term[start:close_index + 1]}")

        name, conjugated = _strip_bar(args[0])
        index_list = args[1].strip()
        if not (index_list.startswith("{") and index_list.endswith("}")):
            raise ValueError("CG index argument is not a list.")
        indices = tuple(
            _parse_index(piece)
            for piece in _split_top_level(index_list[1:-1])
        )
        calls.append(
            WilsonCGCall(
                name=name,
                conjugated=conjugated,
                indices=indices,
            )
        )
        pos = close_index + 1
    return tuple(calls)


def _find_wilson_scalar_legs(term: str) -> tuple[WilsonScalarLeg, ...]:
    """
    We get the 2 scalar fields from the Matchete output into our fields
    """

    legs: list[WilsonScalarLeg] = []
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

            if field_name not in {"NewScalar1", "NewScalar2"}:
                pos = close_index + 1
                continue

            scalar_name = "S1" if field_name == "NewScalar1" else "S2"

            indices = args[2].strip()
            if not (indices.startswith("{") and indices.endswith("}")):
                raise ValueError(
                    f"Scalar Field index argument is not a list: {indices!r}."
                )

            index_body = indices[1:-1].strip()

            if not index_body:
                dummy = None
                rep = None
            else:
                index_items = _split_top_level(index_body)
                if len(index_items) != 1:
                    raise ValueError(
                        "EFT1 Wilson adapter currently expects at most one "
                        "SU(2) index per scalar multiplet."
                    )
                dummy, rep = _parse_index(index_items[0])
                rep, _ = _strip_bar(rep)

            conjugated = False
            if start >= 4 and term[start - 4 : start] == "Bar[":
                bar_open = start - 1
                conjugated = (
                    _matching_bracket(term, bar_open) == close_index + 1
                )

            legs.append(
                WilsonScalarLeg(
                    name=scalar_name,
                    dummy=dummy,
                    rep=rep,
                    conjugated=conjugated,
                )
            )

        pos = close_index + 1

    if len(legs) != 2:
        raise ValueError(
            f"Expected exactly two scalar legs, found {len(legs)} in term."
        )

    return tuple(legs)


_LEPTON_RE = re.compile(
    r"(?:Bar\[)?Field\[l,\s*Fermion,\s*"
    r"\{Index\[(?P<dummy>[^,\]]+),\s*SU2L\[fund\]\]"
)


def _find_lepton_dummies(term: str) -> tuple[str, str]:
    '''Find our Lepton Doublet indices and return it'''
    labels: list[str] = []
    for match in _LEPTON_RE.finditer(term):
        label = match.group("dummy").strip()
        if label not in labels:
            labels.append(label)
    if len(labels) != 2:
        raise ValueError(
            f"Expected two external lepton SU(2) labels, found {labels}."
        )
    return labels[0], labels[1]


def _term_scalar_expression(term: dict) -> str:
    """Reduce a Matchete Wilson term to its exact scalar prefactor.

    Fields, CG calls and spinor chains are replaced by 1.  Yukawa couplings,
    MF, signs and all numerical factors remain.  Therefore factors such as
    1/Sqrt[3] and Sqrt[2/3] are preserved automatically.
    """

    text = term["TermInputForm"]
    out: list[str] = []
    pos = 0

    while pos < len(text):
        if text.startswith("Bar[Coupling[", pos):
            bar_open = pos + len("Bar")
            bar_close = _matching_bracket(text, bar_open)
            inner = text[bar_open + 1 : bar_close]

            if inner.startswith("Coupling[y1,"):
                out.append("yb1")
            elif inner.startswith("Coupling[y2,"):
                out.append("yb2")
            else:
                raise ValueError(
                    f"Unsupported barred coupling in Wilson term {term['Index']}: "
                    f"{inner[:80]}"
                )

            pos = bar_close + 1
            continue

        if text.startswith("Coupling[", pos):
            open_index = pos + len("Coupling")
            close_index = _matching_bracket(text, open_index)
            inner = text[open_index + 1 : close_index]
            args = _split_top_level(inner)
            name = args[0].strip() if args else ""

            if name in {"y1", "y2", "MF"}:
                out.append(name)
            else:
                raise ValueError(
                    f"Unsupported coupling {name!r} in Wilson term {term['Index']}."
                )

            pos = close_index + 1
            continue

        if text.startswith("CG[", pos):
            open_index = pos + len("CG")
            close_index = _matching_bracket(text, open_index)
            out.append("1")
            pos = close_index + 1
            continue

        if text.startswith("Bar[Field[", pos):
            bar_open = pos + len("Bar")
            bar_close = _matching_bracket(text, bar_open)
            out.append("1")
            pos = bar_close + 1
            continue

        if text.startswith("Field[", pos):
            open_index = pos + len("Field")
            close_index = _matching_bracket(text, open_index)
            out.append("1")
            pos = close_index + 1
            continue

        if text.startswith("NCM[", pos):
            open_index = pos + len("NCM")
            close_index = _matching_bracket(text, open_index)
            out.append("1")
            pos = close_index + 1
            continue

        out.append(text[pos])
        pos += 1

    return "".join(out)


def _term_prefactor(term: dict) -> sp.Expr:
    '''We reconstruct scalar coefficient by
    multiplying fields and CG contractions.
    So no more fields, CG objects and spinors butretaining
    yukawa and fermion mass'''
    scalar_text = _term_scalar_expression(term)

    cleaned = scalar_text.strip()
    sqrt_pattern = re.compile(r"Sqrt\[([^\[\]]+)\]")
    while sqrt_pattern.search(cleaned):
        cleaned = sqrt_pattern.sub(r"sqrt(\1)", cleaned)
    cleaned = cleaned.replace("^", "**")

    y1, y2 = sp.symbols("y1 y2")
    yb1, yb2 = sp.symbols("yb1 yb2")
    mf = sp.Symbol("MF", nonzero=True)

    raw = sp.simplify(
        sp.sympify(
            cleaned,
            locals={
                "sqrt": sp.sqrt,
                "I": sp.I,
                "y1": y1,
                "y2": y2,
                "yb1": yb1,
                "yb2": yb2,
                "MF": mf,
            },
        )
    )

    # Cheap structural regressions before replacing barred Yukawa placeholders.
    powers = raw.as_powers_dict()
    yukawa_power = sum(
        int(powers.get(symbol, 0))
        for symbol in (y1, y2, yb1, yb2)
    )
    if yukawa_power != 2:
        raise ValueError(
            f"Expected two Yukawa factors in Wilson term {term['Index']}, "
            f"reconstructed scalar expression={raw}."
        )

    if powers.get(mf, 0) != -1:
        raise ValueError(
            f"Expected exactly one inverse MF in Wilson term {term['Index']}, "
            f"reconstructed scalar expression={raw}."
        )

    result = raw.xreplace(
        {
            yb1: sp.conjugate(y1),
            yb2: sp.conjugate(y2),
        }
    )
    return sp.simplify(result)


def _wilson_cg_value(
    call: WilsonCGCall,
    registry: dict[str, WilsonCGTensor],
    assignment: dict[tuple[str, str], int],
) -> sp.Expr:
    '''For assigned SU2 Components we evaluate the CG'''
    if call.name not in registry:
        raise KeyError(f"CG {call.name!r} is absent from the exported registry.")
    tensor = registry[call.name]

    if len(call.indices) != len(tensor.reps):
        raise ValueError(
            f"{call.name}: call rank {len(call.indices)} does not match "
            f"registry rank {len(tensor.reps)}."
        )

    components: list[int] = []
    for (dummy, rep_from_call), rep_from_registry in zip(call.indices, tensor.reps):
        # Barred representations were stripped when parsing.  The component
        # numbering itself is unchanged by conjugation.
        key = (dummy, rep_from_call)
        if key not in assignment:
            # Matchete can print the same representation with a Bar only in
            # one of the two places.  Fall back to matching dummy label when
            # that label is unique.
            candidates = [
                value for (label, _), value in assignment.items()
                if label == dummy
            ]
            if len(candidates) != 1:
                raise KeyError(
                    f"No unambiguous component assignment for {dummy}, {rep_from_call}."
                )
            components.append(candidates[0])
        else:
            components.append(assignment[key])

    value = tensor.matrix[tuple(c - 1 for c in components)]
    return sp.conjugate(value) if call.conjugated else value


def _wilson_shared_scalar_component_map(
    leg: WilsonScalarLeg,
    component: int,
    dimension: int,
) -> tuple[int, bool, sp.Expr]:
    """Map formal S1/S2 legs to one physical scalar S.

    S1_a = C_ab S_b^*, 
    S2_a = S_a with the standard SU(2)
    C_ab is the charge-conjugation metric in descending-m ordering.
    
    This is for the shared scalar mode where they are counted the same
    """
    if leg.name == "S2":
        return component, leg.conjugated, sp.S.One
    if leg.name != "S1":
        raise KeyError(f"Unknown shared-scalar formal leg {leg.name!r}.")
    physical_component = dimension + 1 - component
    phase = sp.Integer(-1) ** (component - 1)
    return physical_component, (not leg.conjugated), phase


def _scalar_real_expansion(
    leg: WilsonScalarLeg,
    component: int,
    layout: WilsonScalarLayout,
) -> tuple[tuple[int, sp.Expr], tuple[int, sp.Expr]]:
    '''Get complex leg -> real coordinates'''
    root2 = sp.sqrt(2)
    physical_component = component
    conjugated = leg.conjugated
    phase = sp.S.One
    if layout.shared_scalar:
        physical_component, conjugated, phase = _wilson_shared_scalar_component_map(
            leg, component, layout.d_s1
        )
    r = layout.global_real_index(leg.name, physical_component, imaginary=False)
    im = layout.global_real_index(leg.name, physical_component, imaginary=True)
    imag_factor = -sp.I / root2 if conjugated else sp.I / root2
    return (
        (r, sp.simplify(phase / root2)),
        (im, sp.simplify(phase * imag_factor)),
    )


def build_eft1_wilson_tensor(
    seed: dict,
    *,
    d_s1: int,
    d_s2: int,
    chirality: str = "PL",
    lepton_indices: tuple[int, int] = (1, 2),
    include_higgs: bool = True,
    project_operator_symmetry: bool = False,
    shared_scalar: bool = False,
) -> SparseWilsonTensor:
    """Tree Wilson Terms becomes C_ijab
    We get our scalar and lepton legs and get our CG contraction
    We get the coefficients in terms of yukawa, 1/MF etc 
    Then we evaluate our CG and make everything in terms of real
    This gives us C_ijab with coefficient prefactor
    """

    if chirality not in {"PL", "PR"}:
        raise ValueError("chirality must be 'PL' or 'PR'.")

    registry = load_wilson_cg_registry(seed)
    layout = WilsonScalarLayout(
        d_s1=int(d_s1),
        d_s2=int(d_s2),
        include_higgs=include_higgs,
        shared_scalar=bool(shared_scalar),
    )
    scalar_dims = {"S1": int(d_s1), "S2": int(d_s2)}

    entries: dict[tuple[int, int, int, int], sp.Expr] = {}

    selected = [
        term for term in seed.get("TreeWilsonTerms", [])
        if term.get("Chirality") == chirality
    ]
    if not selected:
        raise ValueError(f"No {chirality} Wilson terms were found.")

    for term in selected:
        text = term["TermInputForm"]
        scalar_legs = _find_wilson_scalar_legs(text)
        lepton_dummies = _find_lepton_dummies(text)
        cg_calls = _find_wilson_cg_calls(text)
        prefactor = _term_prefactor(term)

        d1 = scalar_dims[scalar_legs[0].name]
        d2 = scalar_dims[scalar_legs[1].name]

        for lcomp1 in range(1, 3):
            for lcomp2 in range(1, 3):
                for scomp1 in range(1, d1 + 1):
                    for scomp2 in range(1, d2 + 1):
                        assignment: dict[tuple[str, str], int] = {}

                        assignment[(lepton_dummies[0], "SU2L[fund]")] = lcomp1
                        assignment[(lepton_dummies[1], "SU2L[fund]")] = lcomp2

                        for leg, component, dimension in (
                            (scalar_legs[0], scomp1, d1),
                            (scalar_legs[1], scomp2, d2),
                        ):
                            if leg.dummy is None:
                                if dimension != 1:
                                    raise ValueError(
                                        f"Missing SU(2) index for non-singlet {leg.name}."
                                    )
                                if component != 1:
                                    raise ValueError(
                                        f"Invalid singlet component {component} for {leg.name}."
                                    )
                                continue

                            assignment[(leg.dummy, leg.rep)] = component

                        cg_factor = sp.S.One
                        for call in cg_calls:
                            cg_factor *= _wilson_cg_value(call, registry, assignment)
                        cg_factor = sp.simplify(cg_factor)
                        if cg_factor == 0:
                            continue

                        i = lepton_indices[lcomp1 - 1]
                        j = lepton_indices[lcomp2 - 1]

                        for a, fa in _scalar_real_expansion(
                            scalar_legs[0], scomp1, layout
                        ):
                            for b, fb in _scalar_real_expansion(
                                scalar_legs[1], scomp2, layout
                            ):
                                key = (i, j, a, b)
                                entries[key] = sp.simplify(
                                    entries.get(key, sp.S.Zero)
                                    + prefactor * cg_factor * fa * fb
                                )

    tensor = SparseWilsonTensor(entries)

    if project_operator_symmetry:
        tensor = tensor.symmetrized()

    return tensor


def load_and_build_eft1_wilson_tensor(
    seed_path: Path,
    rgbeta_path: Path,
    *,
    chirality: str = "PL",
    lepton_indices: tuple[int, int] = (1, 2),
    project_operator_symmetry: bool = False,
) -> SparseWilsonTensor:
    """Build the tensor using the two files emitted by the current pipeline
    The RGbeta and the matchete EFT1 seed."""

    seed = json.loads(Path(seed_path).read_text(encoding="utf-8"))
    rgbeta = json.loads(Path(rgbeta_path).read_text(encoding="utf-8"))
    meta = rgbeta["metadata"]

    return build_eft1_wilson_tensor(
        seed,
        d_s1=int(meta["dS1"]),
        d_s2=int(meta["dS2"]),
        chirality=chirality,
        lepton_indices=lepton_indices,
        project_operator_symmetry=project_operator_symmetry,
        shared_scalar=bool(meta.get("SharedScalar", False)),
    )


def wilson_tensor_diagnostics(tensor: SparseWilsonTensor) -> dict:
    """Compares Matchete tensor to our new projected one."""

    residuals = tensor.symmetry_residuals()
    projected = tensor.symmetrized()
    projected_residuals = projected.symmetry_residuals()

    return {
        "raw_nonzero_components": len(tensor.entries),
        "raw_fermion_swap_failure_count": len(
            residuals["fermion_swap"]
        ),
        "raw_scalar_swap_failure_count": len(
            residuals["scalar_swap"]
        ),
        "projected_nonzero_components": len(projected.entries),
        "projected_fermion_swap_failure_count": len(
            projected_residuals["fermion_swap"]
        ),
        "projected_scalar_swap_failure_count": len(
            projected_residuals["scalar_swap"]
        ),
        "projected_operator_symmetry": (
            not projected_residuals["fermion_swap"]
            and not projected_residuals["scalar_swap"]
        ),
    }


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Convert a Matchete EFT1 Wilson seed to the real-scalar tensor."
    )
    parser.add_argument("seed_json", type=Path)
    parser.add_argument("rgbeta_json", type=Path)
    parser.add_argument("--chirality", choices=("PL", "PR"), default="PL")
    parser.add_argument(
        "--project-operator-symmetry",
        action="store_true",
        help="Return the separately (ij)- and (ab)-symmetric Wilson tensor.",
    )
    args = parser.parse_args()

    tensor = load_and_build_eft1_wilson_tensor(
        args.seed_json,
        args.rgbeta_json,
        chirality=args.chirality,
        project_operator_symmetry=args.project_operator_symmetry,
    )
    print(json.dumps(wilson_tensor_diagnostics(tensor), indent=2))
    for key, value in sorted(tensor.nonzero_items()):
        print(f"{key}: {value}")


# ---------------------------------------------------------------------------
# Scalar-potential -> lambda_abcd adapter
# ---------------------------------------------------------------------------

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
    """This is fully symmetric."""

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
    """Read actual scalar Field[...] factors, including outer Bar[...] wrappers.
    from our matchete"""

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
    '''get CG from matchete output'''
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
    '''Get coupling from matchete output'''
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
    """
    We find all fields, CG and couplings and note their positions
    Then we can remove them to just get the coefficient
    ."""

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
    '''Replace all fields and couplings and cg with 1'''
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
