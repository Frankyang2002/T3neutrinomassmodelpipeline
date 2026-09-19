from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
from typing import Any

import sympy as sp

# Allow direct execution from the repository root:
#     python RGE/running/EFT1WilsonRunning.py ...
PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from RGE.running.EFT1TensorAdapters import load_and_build_eft1_wilson_tensor


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


def main() -> int:
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


if __name__ == "__main__":
    raise SystemExit(main())
