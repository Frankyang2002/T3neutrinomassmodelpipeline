from __future__ import annotations

"""Adapter from Matchete EFT1 Wilson-seed JSON to the real-scalar C_ijab basis.

The Matchete matching stage exports exact tree-level psi^2 phi^2 terms and the
Clebsch-Gordan tensors they reference.  This module reconstructs those gauge
contractions and converts complex scalar components

    phi_m = (R_m + i I_m)/sqrt(2),
    phi_m^* = (R_m - i I_m)/sqrt(2)

into the interleaved real-scalar basis used by RGE.general.RGEModel and
GaugeGenerators:

    (R_1, I_1, R_2, I_2, ...).

Only the PL Matchete terms are used by default because they contain unbarred
left-handed lepton fields and therefore correspond directly to the project's
left-handed Weyl-fermion convention.  The PR terms are their Hermitian
conjugates.

This file deliberately performs only the Wilson-tensor conversion.  It does
not yet run the master RGE or claim that the overall operator-normalisation
factor has been calibrated against the existing Weinberg adapter.
"""

from dataclasses import dataclass
import json
from pathlib import Path
import re
from typing import Iterable

import sympy as sp

from RGE.running.MatcheteParsing import (
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
class CGTensor:
    name: str
    reps: tuple[str, ...]
    matrix: sp.MutableDenseNDimArray


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
    """Global real-scalar offsets for ordinary or shared-scalar T3 EFT1."""

    d_s1: int
    d_s2: int
    include_higgs: bool = True
    shared_scalar: bool = False

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
    """Sparse C_ijab lookup compatible with MasterWeinbergRGE call sites."""

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
        """Project onto C_(ij)(ab) without changing the represented operator.

        For left-handed two-component fermions the Lorentz scalar psi_i psi_j
        is symmetric under i<->j, while ordinary scalar fields commute under
        a<->b.  Therefore only the separately symmetric part of C_ijab is
        physically observable in the Lagrangian.

        The projector is

            Csym_ijab = 1/4 (
                C_ijab + C_jiab + C_ijba + C_jiba
            ).

        Averaging, rather than summing, is essential: summing would multiply
        the Lagrangian by combinatorial factors.
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
        """Return exact residuals for the two operator-index symmetries."""

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


def load_cg_registry(seed: dict) -> dict[str, CGTensor]:
    registry: dict[str, CGTensor] = {}
    for item in seed.get("CGRegistry", []):
        registry[item["Name"]] = CGTensor(
            name=item["Name"],
            reps=_parse_reps(item["RepsInputForm"]),
            matrix=_parse_sparse_array(item["TensorInputForm"]),
        )
    return registry


def _find_cg_calls(term: str) -> tuple[CGCall, ...]:
    calls: list[CGCall] = []
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
            CGCall(
                name=name,
                conjugated=conjugated,
                indices=indices,
            )
        )
        pos = close_index + 1
    return tuple(calls)


def _find_scalar_legs(term: str) -> tuple[ScalarLeg, ...]:
    """Read the two actual NewScalar1/NewScalar2 Field[...] factors.

    Matchete omits the SU(2) representation index entirely for singlets, e.g.

        Field[NewScalar1, Scalar, {}, {}]

    so a regex that requires Index[...] misses that scalar. Parse Field[...]
    structurally instead and preserve multiplicity of repeated scalar fields.
    """

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
                ScalarLeg(
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


def prefactor_scaling_regression(
    term: dict,
    scale: sp.Expr = sp.Integer(7),
) -> bool:
    """Cheap check that an overall term rescaling survives parsing exactly."""

    scaled = dict(term)
    scaled["TermInputForm"] = f"({sp.sstr(scale)})*({term['TermInputForm']})"
    lhs = sp.simplify(_term_prefactor(scaled))
    rhs = sp.simplify(scale * _term_prefactor(term))
    return sp.simplify(lhs - rhs) == 0


def _cg_value(
    call: CGCall,
    registry: dict[str, CGTensor],
    assignment: dict[tuple[str, str], int],
) -> sp.Expr:
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


def _shared_scalar_component_map(
    leg: ScalarLeg,
    component: int,
    dimension: int,
) -> tuple[int, bool, sp.Expr]:
    """Map formal S1/S2 legs to one physical scalar S.

    S1_a = C_ab S_b^*, S2_a = S_a with the standard SU(2)
    charge-conjugation metric in descending-m ordering.
    """
    if leg.name == "S2":
        return component, leg.conjugated, sp.S.One
    if leg.name != "S1":
        raise KeyError(f"Unknown shared-scalar formal leg {leg.name!r}.")
    physical_component = dimension + 1 - component
    phase = sp.Integer(-1) ** (component - 1)
    return physical_component, (not leg.conjugated), phase


def _scalar_real_expansion(
    leg: ScalarLeg,
    component: int,
    layout: ScalarLayout,
) -> tuple[tuple[int, sp.Expr], tuple[int, sp.Expr]]:
    root2 = sp.sqrt(2)
    physical_component = component
    conjugated = leg.conjugated
    phase = sp.S.One
    if layout.shared_scalar:
        physical_component, conjugated, phase = _shared_scalar_component_map(
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
    """Convert exported tree Wilson terms to the real-scalar C_ijab basis.

    ``lepton_indices`` gives the global fermion-basis positions of the two
    SU(2) components of L.  In the current FermionBasis construction L is the
    first multiplet, so the default is (1, 2).

    The returned tensor is in the direct Lagrangian coefficient normalisation
    obtained from the Matchete term.  A later calibration step should compare
    this convention with ``build_weinberg_wilson_tensor`` before using it for
    a final numerical RGE.
    """

    if chirality not in {"PL", "PR"}:
        raise ValueError("chirality must be 'PL' or 'PR'.")

    registry = load_cg_registry(seed)
    layout = ScalarLayout(
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
        scalar_legs = _find_scalar_legs(text)
        lepton_dummies = _find_lepton_dummies(text)
        cg_calls = _find_cg_calls(text)
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
                            cg_factor *= _cg_value(call, registry, assignment)
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
    """Build the tensor using the two files emitted by the current pipeline."""

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


def tensor_diagnostics(tensor: SparseWilsonTensor) -> dict:
    """Serialisable raw-vs-projected symmetry diagnostics."""

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
    print(json.dumps(tensor_diagnostics(tensor), indent=2))
    for key, value in sorted(tensor.nonzero_items()):
        print(f"{key}: {value}")
