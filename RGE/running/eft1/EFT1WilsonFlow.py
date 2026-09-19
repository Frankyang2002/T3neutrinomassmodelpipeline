from __future__ import annotations

"""EFT1 Wilson-coefficient transport and full-flavor seed extraction.

This module combines the two adjacent EFT1 flow steps while keeping their
calculation bodies separate: transport of the dimension-five Wilson tensor
between thresholds and extraction of the full-flavor seed from Matchete.
"""

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path
import re
import sys
from typing import Any

import sympy as sp

# Allow direct execution:
#     python RGE/running/eft1/EFT1WilsonFlow.py ...
PROJECT_ROOT = Path(__file__).resolve().parents[3]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from RGE.running.eft1.EFT1TensorAdapters import load_and_build_eft1_wilson_tensor

# ---------------------------------------------------------------------------
# Wilson-coefficient transport
# ---------------------------------------------------------------------------


HBAR = sp.Symbol("hbar")


def _expr(raw: str | int | float | sp.Expr) -> sp.Expr:
    """Parse expressions stored by the symbolic RGE pipeline."""
    if isinstance(raw, sp.Expr):
        return raw
    return sp.sympify(
        str(raw),
        locals={
            "conjugate": sp.conjugate,
            "sqrt": sp.sqrt,
            "log": sp.log,
            "I": sp.I,
            "pi": sp.pi,
        },
    )


def _component_key(component: tuple[int, int, int, int]) -> str:
    return ",".join(str(index) for index in component)


def _parse_component_key(key: str) -> tuple[int, int, int, int]:
    pieces = tuple(int(piece) for piece in key.split(","))
    if len(pieces) != 4:
        raise ValueError(f"Invalid Wilson component key: {key!r}")
    return pieces


def _canonical_component(
    component: tuple[int, int, int, int],
) -> tuple[int, int, int, int]:
    i, j, a, b = component
    if i > j:
        i, j = j, i
    if a > b:
        a, b = b, a
    return i, j, a, b


def _text(expression: sp.Expr) -> str:
    return sp.sstr(sp.factor(sp.simplify(expression)))


def _load_rge_payload(path: Path) -> dict[str, Any]:
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    if payload.get("status") != "Success":
        raise ValueError(
            f"EFT1 Wilson-RGE payload is not successful: "
            f"{payload.get('status')!r}"
        )
    if "betas" not in payload:
        raise ValueError("EFT1 Wilson-RGE payload contains no beta functions.")
    return payload


def _tree_boundary_components(
    wilson_seed_path: Path,
    rgbeta_path: Path,
) -> dict[tuple[int, int, int, int], sp.Expr]:
    """Reconstruct the same projected C^(0) tensor used by EFT1WilsonRGE."""
    tensor = load_and_build_eft1_wilson_tensor(
        wilson_seed_path,
        rgbeta_path,
        chirality="PL",
        project_operator_symmetry=True,
    )

    result: dict[tuple[int, int, int, int], sp.Expr] = {}
    for component, value in tensor.nonzero_items():
        canonical = _canonical_component(component)
        result[canonical] = sp.simplify(
            result.get(canonical, sp.S.Zero) + value
        )

    return {
        component: sp.simplify(value)
        for component, value in result.items()
        if sp.simplify(value) != 0
    }


def _beta_components(
    rge_payload: dict[str, Any],
) -> dict[tuple[int, int, int, int], sp.Expr]:
    result: dict[tuple[int, int, int, int], sp.Expr] = {}

    for key, entry in rge_payload["betas"].items():
        if not isinstance(entry, dict) or "total" not in entry:
            continue

        component = _canonical_component(_parse_component_key(key))
        value = sp.simplify(_expr(entry["total"]))
        if value == 0:
            continue

        result[component] = sp.simplify(
            result.get(component, sp.S.Zero) + value
        )

    return result


def run_eft1_wilson_transport(
    wilson_seed_path: Path,
    wilson_rge_path: Path,
    rgbeta_path: Path,
    mu_high: str | sp.Expr,
    mu_low: str | sp.Expr,
    output_path: Path | None = None,
) -> dict[str, Any]:
    """Transport EFT1 Wilson coefficients between two thresholds at O(hbar).

    Convention
    ----------
    The master-RGE payload stores beta^(1) defined by

        16*pi^2 dC/dln(mu) = beta^(1)[C^(0)].

    The EFT expansion is kept in the form

        C(mu_low) = C^(0)
                    + hbar * C_run^(1)
                    + O(hbar^2),

        hbar = 1/(16*pi^2),

    with

        C_run^(1)
          = log(mu_low/mu_high) * beta^(1)[C^(0)].

    Crucially, C^(0) itself is not overwritten.  This separation is what the
    next threshold needs for fixed-order matching:

        M_2^(0)[C^(0)]                       tree piece
        M_2^(1)[C^(0)]                       threshold one-loop piece
        M_2^(0)[hbar * C_run^(1)]            inherited running piece

    The stage-1 one-loop matching boundary term is deliberately not evolved
    here; doing so would first contribute at O(hbar^2).
    """

    mu_high_expr = sp.simplify(_expr(mu_high))
    mu_low_expr = sp.simplify(_expr(mu_low))

    if mu_high_expr == 0 or mu_low_expr == 0:
        raise ValueError("Threshold scales must be nonzero.")

    rge_payload = _load_rge_payload(wilson_rge_path)
    tree = _tree_boundary_components(wilson_seed_path, rgbeta_path)
    beta = _beta_components(rge_payload)

    log_ratio = sp.log(mu_low_expr / mu_high_expr)

    components = sorted(set(tree) | set(beta))

    running_one_loop: dict[tuple[int, int, int, int], sp.Expr] = {}
    full_fixed_order: dict[tuple[int, int, int, int], sp.Expr] = {}

    for component in components:
        c0 = sp.simplify(tree.get(component, sp.S.Zero))
        b1 = sp.simplify(beta.get(component, sp.S.Zero))
        c1_run = sp.simplify(log_ratio * b1)
        full = sp.simplify(c0 + HBAR * c1_run)

        if c1_run != 0:
            running_one_loop[component] = c1_run
        if full != 0:
            full_fixed_order[component] = full

    generated_by_running = sorted(
        component
        for component, value in running_one_loop.items()
        if component not in tree and sp.simplify(value) != 0
    )

    # This is useful for later consistency checks: if the two thresholds
    # coincide, the entire running correction must vanish identically.
    zero_log_check = all(
        sp.simplify(value.subs(mu_low_expr, mu_high_expr)) == 0
        if mu_low_expr != mu_high_expr
        else sp.simplify(value) == 0
        for value in running_one_loop.values()
    )

    payload: dict[str, Any] = {
        "status": "Success",
        "theory": rge_payload.get(
            "theory",
            "SM + S1 + S2 + tree-generated dimension-5 psi2phi2",
        ),
        "order": "fixed overall one loop",
        "convention": {
            "hbar": "1/(16*pi^2)",
            "rge": "16*pi^2*dC/dln(mu)=beta^(1)",
            "transport": (
                "C(mu_low)=C^(0)+hbar*log(mu_low/mu_high)"
                "*beta^(1)[C^(0)]+O(hbar^2)"
            ),
        },
        "scales": {
            "mu_high": _text(mu_high_expr),
            "mu_low": _text(mu_low_expr),
            "log_ratio": _text(log_ratio),
        },
        "metadata": rge_payload.get("metadata", {}),
        "tree_boundary_component_count": len(tree),
        "beta_component_count": len(beta),
        "running_correction_component_count": len(running_one_loop),
        "generated_by_running_component_count": len(generated_by_running),
        "generated_by_running_components": [
            list(component) for component in generated_by_running
        ],
        "equal_scale_running_vanishes": zero_log_check,
        # C^(0): the tree boundary tensor entering the interval.
        "tree_boundary_components": {
            _component_key(component): _text(value)
            for component, value in sorted(tree.items())
        },
        # C_run^(1): coefficient multiplying hbar in the inherited O(hbar)
        # running correction.  This is the form intended for stage-2
        # tree-level propagation.
        "one_loop_running_components": {
            _component_key(component): _text(value)
            for component, value in sorted(running_one_loop.items())
        },
        # Convenience view only.  The perturbatively separated fields above
        # are the authoritative representation for later matching.
        "full_fixed_order_components": {
            _component_key(component): _text(value)
            for component, value in sorted(full_fixed_order.items())
        },
        "weinberg_subspace_validation": rge_payload.get(
            "weinberg_subspace_validation"
        ),
    }

    if output_path is not None:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(
            json.dumps(payload, indent=2),
            encoding="utf-8",
        )

    return payload


def transport_main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Transport tree-generated EFT1 dimension-five Wilson "
            "coefficients between two thresholds at fixed one-loop order."
        )
    )
    parser.add_argument("wilson_seed_json", type=Path)
    parser.add_argument("wilson_rge_json", type=Path)
    parser.add_argument("rgbeta_json", type=Path)
    parser.add_argument(
        "--mu-high",
        required=True,
        help="Upper threshold scale, e.g. MF or 1e10.",
    )
    parser.add_argument(
        "--mu-low",
        required=True,
        help="Lower threshold scale, e.g. MS or 1e8.",
    )
    parser.add_argument("--output", type=Path)

    args = parser.parse_args()

    payload = run_eft1_wilson_transport(
        wilson_seed_path=args.wilson_seed_json,
        wilson_rge_path=args.wilson_rge_json,
        rgbeta_path=args.rgbeta_json,
        mu_high=args.mu_high,
        mu_low=args.mu_low,
        output_path=args.output,
    )

    summary = {
        "status": payload["status"],
        "mu_high": payload["scales"]["mu_high"],
        "mu_low": payload["scales"]["mu_low"],
        "log_ratio": payload["scales"]["log_ratio"],
        "tree_boundary_component_count": payload[
            "tree_boundary_component_count"
        ],
        "beta_component_count": payload["beta_component_count"],
        "running_correction_component_count": payload[
            "running_correction_component_count"
        ],
        "generated_by_running_component_count": payload[
            "generated_by_running_component_count"
        ],
        "equal_scale_running_vanishes": payload[
            "equal_scale_running_vanishes"
        ],
    }
    print(json.dumps(summary, indent=2))

    return 0


# ---------------------------------------------------------------------------
# Flavor-seed extraction
# ---------------------------------------------------------------------------

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


def flavor_seed_main() -> int:
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

def main() -> int:
    parser = argparse.ArgumentParser(
        description="EFT1 Wilson flow utilities."
    )
    parser.add_argument(
        "mode",
        choices=("transport", "flavor-seed"),
        help="Run Wilson transport or flavor-seed extraction.",
    )
    args, remaining = parser.parse_known_args()
    original_argv = sys.argv
    try:
        sys.argv = [original_argv[0], *remaining]
        if args.mode == "transport":
            return transport_main()
        return flavor_seed_main()
    finally:
        sys.argv = original_argv


if __name__ == "__main__":
    raise SystemExit(main())
