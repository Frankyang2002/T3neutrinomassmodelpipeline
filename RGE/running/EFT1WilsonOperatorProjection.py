from __future__ import annotations

"""Project transported EFT1 Wilson components back onto operator directions.

The component-level master RGE is useful for calculation, but the next
Matchete threshold needs gauge-invariant operators rather than arbitrary
real-basis entries C_{ijab}.  This module performs the first inverse map:

  1. reconstruct each tree-generated Matchete operator as its own real-basis
     tensor direction using EFT1WilsonAdapter;
  2. reconstruct the generated Weinberg direction from the already validated
     LLHH beta subspace in eft1_wilson_rge.json;
  3. solve exactly for the coefficients multiplying those operator
     directions in the transported O(hbar) correction;
  4. report any residual component that is not in this operator span.

The residual check is essential.  We do not silently discard component-level
mixing that has not yet been identified with a gauge-invariant operator.
"""

import argparse
import json
from pathlib import Path
import sys
from typing import Any

import sympy as sp


PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from RGE.running.EFT1WilsonAdapter import build_eft1_wilson_tensor


Component = tuple[int, int, int, int]


def _expr(raw: Any) -> sp.Expr:
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


def _text(expression: sp.Expr) -> str:
    return sp.sstr(sp.factor(sp.simplify(expression)))


def _parse_component(key: str) -> Component:
    component = tuple(int(piece) for piece in key.split(","))
    if len(component) != 4:
        raise ValueError(f"Invalid Wilson component key {key!r}.")
    return component  # type: ignore[return-value]


def _component_key(component: Component) -> str:
    return ",".join(str(value) for value in component)


def _canonical(component: Component) -> Component:
    i, j, a, b = component
    return min(i, j), max(i, j), min(a, b), max(a, b)


def _canonical_tensor(entries) -> dict[Component, sp.Expr]:
    """Collapse a symmetric tensor to independent (i<=j,a<=b) entries.

    Symmetry-related entries must agree.  We assign rather than sum, since
    SparseWilsonTensor.symmetrized() explicitly stores those equivalent
    permutations.
    """

    result: dict[Component, sp.Expr] = {}

    for component, raw_value in entries:
        component = _canonical(component)
        value = sp.simplify(raw_value)
        if value == 0:
            continue

        if component in result:
            if sp.simplify(result[component] - value) != 0:
                raise ValueError(
                    "Symmetry-related Wilson components disagree for "
                    f"{component}: {result[component]} vs {value}."
                )
        else:
            result[component] = value

    return result


def _load_json(path: Path, *, require_success: bool = False) -> dict[str, Any]:
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    if require_success and payload.get("status") != "Success":
        raise ValueError(
            f"{path} has status {payload.get('status')!r}, not 'Success'."
        )
    return payload


def _tree_operator_directions(
    seed: dict[str, Any],
    rgbeta: dict[str, Any],
    *,
    chirality: str = "PL",
) -> list[dict[str, Any]]:
    """Build one projected tensor direction for every tree Matchete term."""

    metadata = rgbeta["metadata"]
    d_s1 = int(metadata["dS1"])
    d_s2 = int(metadata["dS2"])

    directions: list[dict[str, Any]] = []

    for term in seed.get("TreeWilsonTerms", []):
        if term.get("Chirality") != chirality:
            continue

        one_term_seed = dict(seed)
        one_term_seed["TreeWilsonTerms"] = [term]

        tensor = build_eft1_wilson_tensor(
            one_term_seed,
            d_s1=d_s1,
            d_s2=d_s2,
            chirality=chirality,
            project_operator_symmetry=True,
        )

        entries = _canonical_tensor(tensor.nonzero_items())
        if not entries:
            raise ValueError(
                f"Tree Wilson term {term.get('Index')} maps to an empty tensor."
            )

        directions.append(
            {
                "kind": "tree",
                "name": f"TreeWilsonTerm{term['Index']}",
                "term_index": int(term["Index"]),
                "term_input_form": term["TermInputForm"],
                "cg_names": list(term.get("CGNames", [])),
                "entries": entries,
            }
        )

    if not directions:
        raise ValueError(f"No {chirality} tree Wilson directions were found.")

    return directions


def _rge_beta_components(
    rge: dict[str, Any],
) -> dict[Component, sp.Expr]:
    result: dict[Component, sp.Expr] = {}

    for key, entry in (rge.get("betas") or {}).items():
        if not isinstance(entry, dict) or "total" not in entry:
            continue
        value = sp.simplify(_expr(entry["total"]))
        if value == 0:
            continue
        component = _canonical(_parse_component(key))

        if component in result:
            if sp.simplify(result[component] - value) != 0:
                raise ValueError(
                    f"RGE symmetry-related components disagree at {component}."
                )
        else:
            result[component] = value

    return result


def _weinberg_direction(
    rge: dict[str, Any],
) -> dict[str, Any] | None:
    """Recover the normalized Weinberg tensor from validated LLHH running.

    EFT1WilsonRGE already validates that the LLHH beta is exactly proportional
    to the benchmark Weinberg tensor and stores the proportionality coefficient
    beta_kappa.  Therefore

        W_ijab = beta_ijab / beta_kappa

    on the Higgs-Higgs subspace gives the normalized operator direction with
    kappa=1, without duplicating a second normalization convention here.
    """

    validation = rge.get("weinberg_subspace_validation")
    if not isinstance(validation, dict):
        return None
    if not validation.get("matches_weinberg_subspace"):
        return None

    beta_kappa_raw = validation.get("beta_kappa_16pi2")
    if beta_kappa_raw in (None, "", "0"):
        return None

    beta_kappa = sp.simplify(_expr(beta_kappa_raw))
    if beta_kappa == 0:
        return None

    beta = _rge_beta_components(rge)
    entries: dict[Component, sp.Expr] = {}

    # ScalarLayout always places the Higgs real components first:
    # H = complex doublet -> four real components 1..4.
    # L is the first fermion multiplet -> components 1,2.
    for component, value in beta.items():
        i, j, a, b = component
        if i in (1, 2) and j in (1, 2) and a <= 4 and b <= 4:
            direction_value = sp.simplify(value / beta_kappa)
            if direction_value != 0:
                entries[component] = direction_value

    if not entries:
        return None

    return {
        "kind": "weinberg",
        "name": "Weinberg",
        "beta_kappa_16pi2": beta_kappa,
        "reference_component": validation.get("reference_component"),
        "entries": entries,
    }


def _transport_target(
    transport: dict[str, Any],
) -> dict[Component, sp.Expr]:
    result: dict[Component, sp.Expr] = {}

    for key, raw_value in (
        transport.get("one_loop_running_components") or {}
    ).items():
        value = sp.simplify(_expr(raw_value))
        if value == 0:
            continue
        component = _canonical(_parse_component(key))

        if component in result:
            if sp.simplify(result[component] - value) != 0:
                raise ValueError(
                    "Transport symmetry-related components disagree at "
                    f"{component}."
                )
        else:
            result[component] = value

    return result


def _solve_projection(
    target: dict[Component, sp.Expr],
    directions: list[dict[str, Any]],
) -> tuple[list[sp.Expr], dict[Component, sp.Expr], int]:
    """Solve target = sum_r x_r direction_r exactly.

    The solve is performed on the union of all target/basis support.  We use
    SymPy's exact linear solver, so square roots, I, symbolic couplings and
    logarithms are retained exactly.
    """

    components = sorted(
        set(target).union(
            *(
                set(direction["entries"])
                for direction in directions
            )
        )
    )

    if not components:
        return [sp.S.Zero] * len(directions), {}, 0

    matrix = sp.Matrix(
        [
            [
                sp.simplify(direction["entries"].get(component, sp.S.Zero))
                for direction in directions
            ]
            for component in components
        ]
    )
    vector = sp.Matrix(
        [sp.simplify(target.get(component, sp.S.Zero)) for component in components]
    )

    rank = int(matrix.rank())

    solution_set = sp.linsolve((matrix, vector))
    if solution_set is sp.EmptySet or not solution_set:
        # No exact solution.  Compute coefficients from an independent row
        # subset where possible, then expose the full residual rather than
        # pretending the basis is complete.
        independent_rows = list(sp.Matrix(matrix.T).rref()[1])
        # rref of M^T gives independent rows of M via pivot columns.
        row_indices = independent_rows[:rank]

        if rank == 0:
            coefficients = [sp.S.Zero] * len(directions)
        else:
            reduced_matrix = matrix[row_indices, :]
            reduced_vector = vector[row_indices, :]
            reduced_solution = sp.linsolve((reduced_matrix, reduced_vector))
            if not reduced_solution:
                coefficients = [sp.S.Zero] * len(directions)
            else:
                raw = next(iter(reduced_solution))
                free_symbols = sorted(
                    set().union(*(entry.free_symbols for entry in raw)),
                    key=str,
                )
                # Set solver-introduced free parameters to zero for one
                # deterministic representative. Physical symbols already
                # occur in the matrix/vector and are not solver parameters.
                physical_symbols = set()
                for entry in matrix:
                    physical_symbols.update(entry.free_symbols)
                for entry in vector:
                    physical_symbols.update(entry.free_symbols)
                substitutions = {
                    symbol: sp.S.Zero
                    for symbol in free_symbols
                    if symbol not in physical_symbols
                }
                coefficients = [
                    sp.simplify(entry.subs(substitutions))
                    for entry in raw
                ]
    else:
        raw = next(iter(solution_set))

        # If the basis is linearly dependent linsolve can introduce free
        # parameters.  Choose the deterministic representative with those
        # solver-only parameters set to zero.
        physical_symbols = set()
        for entry in matrix:
            physical_symbols.update(entry.free_symbols)
        for entry in vector:
            physical_symbols.update(entry.free_symbols)

        free_symbols = set().union(*(entry.free_symbols for entry in raw))
        substitutions = {
            symbol: sp.S.Zero
            for symbol in free_symbols
            if symbol not in physical_symbols
        }
        coefficients = [
            sp.simplify(entry.subs(substitutions))
            for entry in raw
        ]

    residual: dict[Component, sp.Expr] = {}
    for row, component in enumerate(components):
        reconstructed = sp.simplify(
            sum(
                coefficients[column] * matrix[row, column]
                for column in range(len(directions))
            )
        )
        difference = sp.simplify(vector[row] - reconstructed)
        if difference != 0:
            residual[component] = difference

    return coefficients, residual, rank


def project_eft1_wilson_running(
    wilson_seed_path: Path,
    rgbeta_path: Path,
    wilson_rge_path: Path,
    transport_path: Path,
    output_path: Path | None = None,
) -> dict[str, Any]:
    """Project the transported O(hbar) Wilson correction onto operators."""

    seed = _load_json(wilson_seed_path)
    rgbeta = _load_json(rgbeta_path, require_success=True)
    rge = _load_json(wilson_rge_path, require_success=True)
    transport = _load_json(transport_path, require_success=True)

    directions = _tree_operator_directions(seed, rgbeta, chirality="PL")

    weinberg = _weinberg_direction(rge)
    if weinberg is not None:
        directions.append(weinberg)

    target = _transport_target(transport)
    coefficients, residual, rank = _solve_projection(target, directions)

    log_ratio = sp.simplify(
        _expr((transport.get("scales") or {}).get("log_ratio", "1"))
    )

    operator_coefficients = []
    for direction, coefficient in zip(directions, coefficients):
        entry: dict[str, Any] = {
            "kind": direction["kind"],
            "name": direction["name"],
            "running_multiplier": _text(coefficient),
        }

        if log_ratio != 0:
            entry["beta_multiplier_16pi2"] = _text(
                sp.simplify(coefficient / log_ratio)
            )

        if direction["kind"] == "tree":
            entry.update(
                {
                    "term_index": direction["term_index"],
                    "cg_names": direction["cg_names"],
                    "term_input_form": direction["term_input_form"],
                    "interpretation": (
                        "Multiply this original Matchete tree operator by "
                        "hbar * running_multiplier before stage-2 tree matching."
                    ),
                }
            )
        elif direction["kind"] == "weinberg":
            entry.update(
                {
                    "beta_kappa_16pi2": _text(
                        direction["beta_kappa_16pi2"]
                    ),
                    "reference_component": direction["reference_component"],
                    "interpretation": (
                        "Normalized Weinberg operator with kappa=1; its "
                        "O(hbar) coefficient is running_multiplier."
                    ),
                }
            )

        operator_coefficients.append(entry)

    payload: dict[str, Any] = {
        "status": "Success" if not residual else "IncompleteBasis",
        "order": "fixed overall one loop",
        "chirality": "PL",
        "scales": transport.get("scales", {}),
        "target_running_component_count": len(target),
        "operator_direction_count": len(directions),
        "operator_basis_rank": rank,
        "exact_operator_projection": not residual,
        "residual_component_count": len(residual),
        "operator_coefficients": operator_coefficients,
        "residual_components": {
            _component_key(component): _text(value)
            for component, value in sorted(residual.items())
        },
        "bookkeeping": {
            "insert_into_stage2_as": (
                "hbar * [sum tree running_multiplier * original tree operator "
                "+ Weinberg running_multiplier * O5]"
            ),
            "do_not_run_stage1_one_loop_boundary_again": True,
        },
    }

    if output_path is not None:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    return payload


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Project transported EFT1 Wilson components back onto the "
            "tree-generated gauge-invariant operators plus Weinberg."
        )
    )
    parser.add_argument("wilson_seed_json", type=Path)
    parser.add_argument("rgbeta_json", type=Path)
    parser.add_argument("wilson_rge_json", type=Path)
    parser.add_argument("transport_json", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    payload = project_eft1_wilson_running(
        wilson_seed_path=args.wilson_seed_json,
        rgbeta_path=args.rgbeta_json,
        wilson_rge_path=args.wilson_rge_json,
        transport_path=args.transport_json,
        output_path=args.output,
    )

    summary = {
        "status": payload["status"],
        "target_running_component_count": payload[
            "target_running_component_count"
        ],
        "operator_direction_count": payload["operator_direction_count"],
        "operator_basis_rank": payload["operator_basis_rank"],
        "exact_operator_projection": payload["exact_operator_projection"],
        "residual_component_count": payload["residual_component_count"],
    }
    print(json.dumps(summary, indent=2))

    # IncompleteBasis is a scientifically meaningful diagnostic, not a parser
    # crash, but return nonzero so automated pipeline tests cannot overlook it.
    return 0 if payload["exact_operator_projection"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
