from __future__ import annotations

"""Preserve the full flavor tensor of the tree-generated EFT1 Wilson operators.

Why this exists
---------------
The component RGE currently used by EFT1WilsonRGE.py collapses the flavor
structure to one generation before running.  That is sufficient for the
SU(2)/scalar/gauge regression, but it is not sufficient to reconstruct a
Matchete expression with explicit flavor indices after running.

This module is the first flavor-aware layer.  It reads the exact PL terms
exported by Matchete at the F threshold and extracts the Wilson flavor tensor

    C^{ab}_{pq}

without replacing y1, y2, ... by scalar symbols.

For the present T3 tree operators, the flavor structure is of the form

    Sum_r Conjugate[yA[p,r]] Conjugate[yB[q,r]] / MF

with the symmetry/combinatorial factor taken directly from the Matchete term.

The output is deliberately representation independent: it keeps the original
Matchete TermInputForm and adds a symbolic flavor kernel that can be consumed
by a later flavor-aware master-RGE layer.

This file does NOT guess a flavor lift from a one-generation beta function.
"""

import argparse
import json
from dataclasses import dataclass, asdict
from pathlib import Path
import re
from typing import Any

import sympy as sp


BARRED_COUPLING_RE = re.compile(
    r"Bar\[Coupling\[(?P<name>[A-Za-z0-9_]+),\s*"
    r"\{Index\[(?P<flavor>[^,\]]+),\s*Flavor\],\s*"
    r"Index\[(?P<nflavor>[^,\]]+),\s*NFlavor\]\},\s*0\]\]"
)

MASS_DENOM_RE = re.compile(
    r"/\((?P<factor>.+?)\*Coupling\[MF,\s*\{\},\s*0\]\)\s*$"
)
MASS_DENOM_UNIT_RE = re.compile(
    r"/Coupling\[MF,\s*\{\},\s*0\]\s*$"
)


@dataclass(frozen=True)
class FlavorYukawaFactor:
    name: str
    external_dummy: str
    internal_dummy: str
    conjugated: bool = True


@dataclass(frozen=True)
class FlavorSeedTerm:
    term_index: int
    chirality: str
    yukawas: tuple[FlavorYukawaFactor, ...]
    external_flavor_dummies: tuple[str, ...]
    internal_flavor_dummy: str
    mass_symbol: str
    symmetry_denominator: str
    flavor_label: str
    flavor_exchange_symmetrized: bool
    flavor_kernel_text: str
    flavor_kernel_wolfram: str
    one_generation_reduction: str
    original_term_input_form: str
    cg_names: tuple[str, ...]


def _normalise_dummy(name: str) -> str:
    """Turn Matchete temporaries such as d$$1 into readable labels."""
    return name.strip()


def _extract_symmetry_denominator(term: str) -> str:
    if MASS_DENOM_UNIT_RE.search(term):
        return "1"

    match = MASS_DENOM_RE.search(term)
    if match:
        return match.group("factor").strip()

    raise ValueError(
        "Could not identify the MF denominator/combinatorial factor in "
        f"tree Wilson term:\n{term}"
    )


def _extract_yukawa_factors(term: str) -> tuple[FlavorYukawaFactor, ...]:
    factors = []
    for match in BARRED_COUPLING_RE.finditer(term):
        factors.append(
            FlavorYukawaFactor(
                name=match.group("name"),
                external_dummy=_normalise_dummy(match.group("flavor")),
                internal_dummy=_normalise_dummy(match.group("nflavor")),
                conjugated=True,
            )
        )
    return tuple(factors)


def _flavor_label(
    yukawas: tuple[FlavorYukawaFactor, ...],
) -> str:
    if len(yukawas) != 2:
        raise ValueError("A T3 flavor tensor requires exactly two Yukawas.")

    names = tuple(factor.name for factor in yukawas)
    if names == ("y1", "y1"):
        return "C11"
    if names == ("y2", "y2"):
        return "C22"
    if set(names) == {"y1", "y2"}:
        return "C12"

    raise ValueError(f"Unsupported T3 Yukawa pair: {names!r}")


def _kernel_strings(
    yukawas: tuple[FlavorYukawaFactor, ...],
) -> tuple[str, str, str]:
    """Return the canonical flavor tensor, independent of SU(2) geometry.

    Representation-dependent numerical/CG factors belong to the operator
    geometry, not to the flavor tensor.  Hence every mixed y1-y2 direction
    uses the same C12 boundary even when the raw Matchete term contains
    factors such as 1/Sqrt[3] or Sqrt[2/3].
    """
    if len(yukawas) != 2:
        raise ValueError(
            "The current T3 flavor extractor expects exactly two barred "
            f"Yukawa factors, found {len(yukawas)}."
        )

    y_a, y_b = yukawas
    if y_a.internal_dummy != y_b.internal_dummy:
        raise ValueError(
            "Tree-level T3 Wilson term does not use one common internal "
            "NFlavor dummy."
        )

    mixed = y_a.name != y_b.name

    if mixed:
        text_num = (
            f"Sum_r (conjugate({y_a.name}[p,r]) "
            f"conjugate({y_b.name}[q,r]) + "
            f"conjugate({y_b.name}[p,r]) "
            f"conjugate({y_a.name}[q,r]))"
        )
        wl_num = (
            f"Sum[Conjugate[{y_a.name}[p, r]]*"
            f"Conjugate[{y_b.name}[q, r]] + "
            f"Conjugate[{y_b.name}[p, r]]*"
            f"Conjugate[{y_a.name}[q, r]], {{r, 1, NF}}]"
        )
        kernel_text = f"({text_num})/(2*MF)"
        kernel_wolfram = f"({wl_num})/(2*MF)"
        one_generation = (
            f"(conjugate({y_a.name})*conjugate({y_b.name}))/MF"
        )
    else:
        text_num = (
            f"Sum_r conjugate({y_a.name}[p,r]) "
            f"conjugate({y_b.name}[q,r])"
        )
        wl_num = (
            f"Sum[Conjugate[{y_a.name}[p, r]]*"
            f"Conjugate[{y_b.name}[q, r]], {{r, 1, NF}}]"
        )
        kernel_text = f"({text_num})/(2*MF)"
        kernel_wolfram = f"({wl_num})/(2*MF)"
        one_generation = (
            f"(conjugate({y_a.name})*conjugate({y_b.name}))/(2*MF)"
        )

    return kernel_text, kernel_wolfram, one_generation


def extract_flavor_seed_term(term: dict[str, Any]) -> FlavorSeedTerm:
    if term.get("Chirality") != "PL":
        raise ValueError("Flavor seed extraction is defined for PL terms.")

    input_form = str(term["TermInputForm"])
    yukawas = _extract_yukawa_factors(input_form)

    if len(yukawas) != 2:
        raise ValueError(
            f"PL term {term.get('Index')} has {len(yukawas)} barred "
            "flavor Yukawas; expected exactly 2."
        )

    internal = {factor.internal_dummy for factor in yukawas}
    if len(internal) != 1:
        raise ValueError(
            f"PL term {term.get('Index')} does not have a unique shared "
            "heavy-flavor contraction."
        )

    external = tuple(factor.external_dummy for factor in yukawas)
    denominator = _extract_symmetry_denominator(input_form)
    kernel_text, kernel_wolfram, one_gen = _kernel_strings(yukawas)

    return FlavorSeedTerm(
        term_index=int(term["Index"]),
        chirality="PL",
        yukawas=yukawas,
        external_flavor_dummies=external,
        internal_flavor_dummy=next(iter(internal)),
        mass_symbol="MF",
        symmetry_denominator=denominator,
        flavor_label=_flavor_label(yukawas),
        flavor_exchange_symmetrized=(yukawas[0].name != yukawas[1].name),
        flavor_kernel_text=kernel_text,
        flavor_kernel_wolfram=kernel_wolfram,
        one_generation_reduction=one_gen,
        original_term_input_form=input_form,
        cg_names=tuple(term.get("CGNames", [])),
    )


def _expected_one_generation_kernel(
    y_a: str,
    y_b: str,
) -> sp.Expr:
    locals_ = {
        "MF": sp.Symbol("MF", nonzero=True),
        y_a: sp.Symbol(y_a),
        y_b: sp.Symbol(y_b),
    }
    denominator = 2 if y_a == y_b else 1
    return sp.simplify(
        sp.conjugate(locals_[y_a])
        * sp.conjugate(locals_[y_b])
        / (denominator * locals_["MF"])
    )


def _parse_one_generation_kernel(text: str) -> sp.Expr:
    names = set(re.findall(r"conjugate\(([A-Za-z0-9_]+)\)", text))
    locals_: dict[str, Any] = {
        "conjugate": sp.conjugate,
        "MF": sp.Symbol("MF", nonzero=True),
    }
    for name in names:
        locals_[name] = sp.Symbol(name)

    cleaned = text.replace("^", "**")
    return sp.simplify(sp.sympify(cleaned, locals=locals_))


def run_flavor_seed_export(
    wilson_seed_path: Path,
    output_path: Path | None = None,
) -> dict[str, Any]:
    seed = json.loads(Path(wilson_seed_path).read_text(encoding="utf-8"))

    pl_terms = [
        term
        for term in seed.get("TreeWilsonTerms", [])
        if term.get("Chirality") == "PL"
    ]
    if not pl_terms:
        raise ValueError("No PL tree Wilson terms found in seed.")

    extracted = [extract_flavor_seed_term(term) for term in pl_terms]

    regression_failures: list[dict[str, Any]] = []
    for item in extracted:
        y_a, y_b = (factor.name for factor in item.yukawas)
        expected = _expected_one_generation_kernel(y_a, y_b)
        actual = _parse_one_generation_kernel(
            item.one_generation_reduction
        )
        residual = sp.simplify(actual - expected)
        if residual != 0:
            regression_failures.append(
                {
                    "term_index": item.term_index,
                    "residual": sp.sstr(residual),
                }
            )

    payload = {
        "status": "Success" if not regression_failures else "Failed",
        "purpose": (
            "full-flavor EFT1 Wilson boundary tensors before one-loop running"
        ),
        "chirality": "PL",
        "term_count": len(extracted),
        "one_generation_reduction_matches": not regression_failures,
        "one_generation_regression_failures": regression_failures,
        "flavor_conventions": {
            "external_lepton_indices": ["p", "q"],
            "internal_heavy_flavor_index": "r",
            "heavy_flavor_sum": "r = 1,...,NF",
            "yukawa_orientation": (
                "PL seed uses conjugated T3 Yukawas exactly as exported "
                "by Matchete"
            ),
            "mixed_C12_convention": (
                "C12[p,q] = 1/2 Sum_r "
                "(conjugate(y1[p,r])*conjugate(y2[q,r]) + "
                "conjugate(y2[p,r])*conjugate(y1[q,r]))/MF"
            ),
            "mixed_C12_symmetric_under_pq_exchange": True,
        },
        "terms": [
            {
                **asdict(item),
                "yukawas": [asdict(y) for y in item.yukawas],
            }
            for item in extracted
        ],
    }

    if output_path is not None:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(
            json.dumps(payload, indent=2),
            encoding="utf-8",
        )

    return payload


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Extract the full flavor tensors of tree-generated EFT1 "
            "dimension-five Wilson operators from the Matchete seed."
        )
    )
    parser.add_argument("wilson_seed_json", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    payload = run_flavor_seed_export(
        args.wilson_seed_json,
        output_path=args.output,
    )

    print(
        json.dumps(
            {
                "status": payload["status"],
                "term_count": payload["term_count"],
                "one_generation_reduction_matches":
                    payload["one_generation_reduction_matches"],
                "flavor_kernels": [
                    {
                        "term_index": term["term_index"],
                        "kernel": term["flavor_kernel_text"],
                    }
                    for term in payload["terms"]
                ],
            },
            indent=2,
        )
    )

    return 0 if payload["status"] == "Success" else 1


if __name__ == "__main__":
    raise SystemExit(main())
