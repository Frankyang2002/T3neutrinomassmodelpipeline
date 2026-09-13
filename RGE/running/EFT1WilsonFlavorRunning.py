from __future__ import annotations

"""Leading-log transport for the flavor-aware EFT1 Wilson tensors.

At fixed overall one-loop order,

    C_X(mu_low) =
        C_X(mu_high)
        + hbar * log(mu_low/mu_high) * beta_X[C(mu_high)],

with hbar = 1/(16*pi^2).

The flavor-aware RGE file stores

    16*pi^2 beta_X[p,q]
      = sum_Y A_XY C_Y[p,q]
      + sum_Y B_XY (He.C_Y + C_Y.He^T)[p,q],

where He = Ye Ye^dagger.

This module inserts the full-flavor boundary tensors exported by
EFT1WilsonFlavorSeed.py and writes the transported O(hbar) corrections in a
representation that can be translated to Matchete without collapsing flavor.

It intentionally does not run the stage-1 one-loop matching boundary again.
"""

import argparse
import json
from pathlib import Path
from typing import Any

import sympy as sp


HBAR = sp.Symbol("hbar")


def _expr(raw: Any) -> sp.Expr:
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


def _text(expr: sp.Expr) -> str:
    return sp.sstr(sp.factor(sp.simplify(expr)))


def _load_json(path: Path, *, require_success: bool = False) -> dict[str, Any]:
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    if require_success and payload.get("status") != "Success":
        raise ValueError(
            f"{path} has status {payload.get('status')!r}, not 'Success'."
        )
    return payload


def _canonical_boundary_map(
    flavor_seed: dict[str, Any],
) -> dict[str, dict[str, Any]]:
    """Return one canonical C11/C12/C22 flavor definition per Yukawa pair."""
    result: dict[str, dict[str, Any]] = {}
    for term in flavor_seed.get("terms") or []:
        label = str(term["flavor_label"])
        if label not in result:
            result[label] = term
        else:
            # Every operator contraction carrying one flavor label must use the
            # exact same canonical full-flavor boundary.
            for key in (
                "flavor_kernel_text",
                "flavor_kernel_wolfram",
                "one_generation_reduction",
            ):
                if str(result[label][key]) != str(term[key]):
                    raise ValueError(
                        f"Inconsistent canonical {label} flavor boundary."
                    )
    if not result:
        raise ValueError("Flavor seed contains no canonical boundaries.")
    return result


def _basis_metadata(flavor_rge: dict[str, Any]) -> dict[str, dict[str, Any]]:
    basis = flavor_rge.get("basis") or []
    metadata = flavor_rge.get("basis_metadata") or {}
    if not basis or any(name not in metadata for name in basis):
        raise ValueError("Flavor RGE operator-basis metadata is incomplete.")
    return {name: metadata[name] for name in basis}


def _he_action_text(source: str) -> str:
    return (
        f"Sum_s(He[p,s]*{source}[s,q] + "
        f"{source}[p,s]*He[q,s])"
    )


def _he_action_wolfram(source: str) -> str:
    return (
        f"Sum[He[p, s]*{source}[s, q] + "
        f"{source}[p, s]*He[q, s], {{s, 1, 3}}]"
    )


def _boundary_definition_text(name: str, term: dict[str, Any]) -> str:
    return str(term["flavor_kernel_text"]).replace("Sum_r ", "Sum_r ")


def _boundary_definition_wolfram(name: str, term: dict[str, Any]) -> str:
    # The flavor-seed exporter used symbolic yA[p,r] syntax. Keep that
    # convention here; the later Matchete exporter maps these to Coupling[].
    return str(term["flavor_kernel_wolfram"])


def _one_generation_boundary(term: dict[str, Any]) -> sp.Expr:
    return _expr(term["one_generation_reduction"])


def run_flavor_wilson_transport(
    flavor_seed_path: Path,
    flavor_rge_path: Path,
    mu_high: str | sp.Expr,
    mu_low: str | sp.Expr,
    output_path: Path | None = None,
) -> dict[str, Any]:
    flavor_seed = _load_json(flavor_seed_path, require_success=True)
    flavor_rge = _load_json(flavor_rge_path, require_success=True)

    if not flavor_rge.get("one_generation_reduction_matches"):
        raise ValueError(
            "Flavor RGE one-generation regression is not successful."
        )

    boundary = _canonical_boundary_map(flavor_seed)
    basis_meta = _basis_metadata(flavor_rge)
    operator_rge = flavor_rge.get("operator_rge") or {}

    mu_high_expr = _expr(mu_high)
    mu_low_expr = _expr(mu_low)
    log_ratio = sp.simplify(sp.log(mu_low_expr / mu_high_expr))

    corrections: dict[str, Any] = {}
    one_gen_regression: dict[str, Any] = {}
    all_one_gen_match = True

    ye = sp.Symbol("ye")
    ye2 = ye * sp.conjugate(ye)

    projected_one_gen = {
        name: _expr(expr)
        for name, expr in (
            flavor_rge.get("projected_one_generation_beta") or {}
        ).items()
    }

    canonical_one_gen = {
        label: _one_generation_boundary(term)
        for label, term in boundary.items()
    }

    operator_boundary_substitutions = {
        sp.Symbol(name): (
            _expr(meta["boundary_weight"])
            * canonical_one_gen[str(meta["flavor_label"])]
        )
        for name, meta in basis_meta.items()
    }

    for target, sources in operator_rge.items():
        # At fixed overall one-loop order, beta is evaluated on the tree-level
        # UV boundary.  Therefore independent operator coefficients are replaced
        # by boundary_weight * canonical flavor tensor before transport.
        aggregated: dict[tuple[str, str], sp.Expr] = {}

        for source, pieces in sources.items():
            meta = basis_meta[source]
            label = str(meta["flavor_label"])
            weight = sp.simplify(_expr(meta["boundary_weight"]))

            a = sp.simplify(_expr(pieces.get("flavor_blind", "0")) * weight)
            b = sp.simplify(_expr(pieces.get("He_two_leg_action", "0")) * weight)

            if a != 0:
                key = (label, "flavor_blind")
                aggregated[key] = sp.simplify(
                    aggregated.get(key, sp.S.Zero) + a
                )
            if b != 0:
                key = (label, "He_two_leg_action")
                aggregated[key] = sp.simplify(
                    aggregated.get(key, sp.S.Zero) + b
                )

        term_records: list[dict[str, Any]] = []
        one_gen_beta = sp.S.Zero

        for (label, kind), coefficient in sorted(aggregated.items()):
            coefficient = sp.simplify(coefficient)
            if coefficient == 0:
                continue

            if kind == "flavor_blind":
                tensor_text = f"{label}[p,q]"
                tensor_wolfram = f"{label}[p, q]"
                one_gen_beta += coefficient * canonical_one_gen[label]
            else:
                tensor_text = _he_action_text(label)
                tensor_wolfram = _he_action_wolfram(label)
                one_gen_beta += (
                    2 * coefficient * ye2 * canonical_one_gen[label]
                )

            term_records.append(
                {
                    "source": label,
                    "kind": kind,
                    "beta_coefficient": _text(coefficient),
                    "running_coefficient": _text(log_ratio * coefficient),
                    "tensor_text": tensor_text,
                    "tensor_wolfram": tensor_wolfram,
                }
            )

        expected_operator = projected_one_gen.get(target, sp.S.Zero)
        expected = sp.simplify(
            expected_operator.subs(operator_boundary_substitutions)
        )
        residual = sp.simplify(one_gen_beta - expected)
        matches = residual == 0
        all_one_gen_match = all_one_gen_match and matches
        one_gen_regression[target] = {
            "matches": matches,
            "residual": _text(residual),
        }

        corrections[target] = {
            "terms": term_records,
            "beta_tensor_text": (
                " + ".join(
                    f"({record['beta_coefficient']})*{record['tensor_text']}"
                    for record in term_records
                )
                or "0"
            ),
            "running_tensor_text": (
                " + ".join(
                    f"({record['running_coefficient']})*"
                    f"{record['tensor_text']}"
                    for record in term_records
                )
                or "0"
            ),
            "running_tensor_wolfram": (
                " + ".join(
                    f"({record['running_coefficient']})*"
                    f"{record['tensor_wolfram']}"
                    for record in term_records
                )
                or "0"
            ),
        }

    equal_scale_vanishes = all(
        sp.simplify(
            _expr(record["running_coefficient"]).subs(mu_low_expr, mu_high_expr)
        ) == 0
        for target in corrections.values()
        for record in target["terms"]
    )

    payload = {
        "status": (
            "Success"
            if all_one_gen_match and equal_scale_vanishes
            else "RegressionFailed"
        ),
        "order": "fixed overall one loop",
        "scales": {
            "mu_high": _text(mu_high_expr),
            "mu_low": _text(mu_low_expr),
            "log_ratio": _text(log_ratio),
        },
        "boundary_tensors": {
            name: {
                "flavor_label": name,
                "definition_text": _boundary_definition_text(name, term),
                "definition_wolfram": _boundary_definition_wolfram(name, term),
                "one_generation_reduction":
                    term["one_generation_reduction"],
                "manifestly_symmetric_under_pq": True,
            }
            for name, term in boundary.items()
        },
        "operator_directions": {
            name: {
                **meta,
                "boundary_weight": _text(_expr(meta["boundary_weight"])),
            }
            for name, meta in basis_meta.items()
        },
        "definitions": {
            "He": "Ye.Ye^dagger",
            "running_rule": (
                "C(mu_low)=C(mu_high)+hbar*log(mu_low/mu_high)*beta16pi2"
            ),
            "hbar": "1/(16*pi^2)",
        },
        "running_corrections": corrections,
        "one_generation_regression": one_gen_regression,
        "one_generation_reduction_matches": all_one_gen_match,
        "equal_scale_running_vanishes": equal_scale_vanishes,
        "bookkeeping": {
            "stage1_one_loop_boundary_evolved": False,
            "stage2_heavy_running_piece": list(basis_meta),
            "stage2_direct_low_energy_piece": ["Weinberg"],
        },
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
            "Transport the full-flavor EFT1 Wilson tensors between two "
            "thresholds at fixed overall one-loop order."
        )
    )
    parser.add_argument("flavor_seed_json", type=Path)
    parser.add_argument("flavor_rge_json", type=Path)
    parser.add_argument("mu_high")
    parser.add_argument("mu_low")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    result = run_flavor_wilson_transport(
        flavor_seed_path=args.flavor_seed_json,
        flavor_rge_path=args.flavor_rge_json,
        mu_high=args.mu_high,
        mu_low=args.mu_low,
        output_path=args.output,
    )

    print(
        json.dumps(
            {
                "status": result["status"],
                "mu_high": result["scales"]["mu_high"],
                "mu_low": result["scales"]["mu_low"],
                "log_ratio": result["scales"]["log_ratio"],
                "transported_operator_count": len(
                    result["running_corrections"]
                ),
                "one_generation_reduction_matches":
                    result["one_generation_reduction_matches"],
                "equal_scale_running_vanishes":
                    result["equal_scale_running_vanishes"],
            },
            indent=2,
        )
    )

    return 0 if result["status"] == "Success" else 1


if __name__ == "__main__":
    raise SystemExit(main())
