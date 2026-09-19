from __future__ import annotations

"""Flavor-aware operator-space RGE for the EFT after integrating out F.

This module separates full-flavor tensors from gauge/operator geometry.

The canonical flavor tensors are C11[p,q], C12[p,q], C22[p,q].  A given
SU(2) representation can contain several raw Matchete contractions carrying
the same flavor tensor.  Those contractions are rank-reduced before the
anomalous-dimension projection.

Strategy
--------
1. Reconstruct the three PL real-scalar tensor directions from the Matchete
   Wilson seed.
2. Divide out the one-generation flavor kernel exported by
   EFT1WilsonFlavorSeed.py.  This leaves pure gauge/scalar operator geometry.
3. Seed the master RGE with independent symbols c11,c12,c22 multiplying those
   geometry tensors.
4. Project the one-loop beta tensor back onto the same three operator
   directions plus the Weinberg direction.
5. Split each projected coefficient into

       A * C[p,q]
       + B * (He.C + C.He^T)[p,q]

   where He = Ye Ye^dagger.

The coefficient B is fixed by the exact one-generation reduction: in one
generation

    (He.C + C.He^T) -> 2 |ye|^2 C.

This removes y1 and y2 from the anomalous-dimension matrix entirely; they now
live only in the boundary tensors C11,C12,C22, where their full flavor
contractions belong.

The output includes an exact one-generation regression against the current
component-level EFT1 Wilson RGE.  We do not proceed if that reduction fails.
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

from RGE.general.AnomalousDimensions import calculate_complete_master_rge
from RGE.running.EFT1WilsonRGE import (
    build_eft1_rge_context,
    closure_output_components,
)
from RGE.matching.SMEFTBasis import build_sm_eft
from RGE.matching.WeinbergWilsonAdapter import (
    build_weinberg_wilson_tensor,
)
from RGE.running.EFT1WilsonAdapter import SparseWilsonTensor
from RGE.running.EFT1WilsonFlavorSeed import run_flavor_seed_export


YE = sp.Symbol("ye")
YE2 = sp.Symbol("Ye2")


def _parse_expr(raw: Any) -> sp.Expr:
    return sp.sympify(
        str(raw),
        locals={
            "conjugate": sp.conjugate,
            "sqrt": sp.sqrt,
            "log": sp.log,
            "I": sp.I,
            "pi": sp.pi,
            "MF": sp.Symbol("MF", nonzero=True),
            "y1": sp.Symbol("y1"),
            "y2": sp.Symbol("y2"),
            "ye": YE,
        },
    )


def _canonical(component: tuple[int, int, int, int]):
    i, j, a, b = component
    return min(i, j), max(i, j), min(a, b), max(a, b)


def _canonical_entries(tensor) -> dict[tuple[int, int, int, int], sp.Expr]:
    result: dict[tuple[int, int, int, int], sp.Expr] = {}
    for component, raw in tensor.nonzero_items():
        key = _canonical(component)
        value = sp.simplify(raw)
        if value == 0:
            continue
        if key in result:
            if sp.simplify(result[key] - value) != 0:
                raise ValueError(
                    f"Symmetry-related tensor entries disagree at {key}."
                )
        else:
            result[key] = value
    return result


def _load_flavor_seed(
    flavor_seed_path: Path | None,
    wilson_seed_path: Path,
) -> dict[str, Any]:
    if flavor_seed_path is not None:
        payload = json.loads(
            Path(flavor_seed_path).read_text(encoding="utf-8")
        )
    else:
        payload = run_flavor_seed_export(wilson_seed_path)

    if payload.get("status") != "Success":
        raise ValueError("Flavor seed is not successful.")
    if not payload.get("one_generation_reduction_matches"):
        raise ValueError("Flavor seed one-generation regression failed.")
    if int(payload.get("term_count", 0)) < 1:
        raise ValueError("Flavor seed contains no PL seed terms.")
    return payload


def _flavor_kernel_expr(text: str) -> sp.Expr:
    """Parse the already-recorded one-generation reduction."""
    return _parse_expr(text)


def _entries_vector(entries, components):
    return sp.Matrix([entries.get(component, sp.S.Zero) for component in components])


def _reduce_geometry_group(flavor_label, raw_directions):
    """Rank-reduce raw contractions while preserving the complete UV boundary."""
    components = sorted(set().union(*(set(item["entries"]) for item in raw_directions)))
    matrix = sp.Matrix.hstack(*[
        _entries_vector(item["entries"], components)
        for item in raw_directions
    ])
    pivots = list(matrix.rref()[1])
    if not pivots:
        raise ValueError(f"No nonzero operator geometry for {flavor_label}.")

    independent = matrix[:, pivots]
    total_boundary = matrix * sp.ones(len(raw_directions), 1)
    solution = sp.linsolve((independent, total_boundary))
    if not solution:
        raise ValueError(
            f"Could not express the {flavor_label} tree boundary in its independent basis."
        )
    weights = next(iter(solution))
    if any(sym.name.startswith("tau") for expr in weights for sym in expr.free_symbols):
        raise ValueError(f"Rank-reduced {flavor_label} boundary is underdetermined.")

    rank = len(pivots)
    all_terms = [int(item["term_index"]) for item in raw_directions]
    result=[]
    for ordinal,(pivot,weight) in enumerate(zip(pivots,weights,strict=True),start=1):
        raw=raw_directions[pivot]
        name = flavor_label if rank == 1 else f"{flavor_label}_{ordinal}"
        result.append({
            **raw,
            "name": name,
            "flavor_label": flavor_label,
            "boundary_weight": sp.simplify(weight),
            "source_term_indices": all_terms,
            "raw_group_count": len(raw_directions),
            "group_rank": rank,
        })
    return result


def _geometry_directions(context, flavor_seed):
    """Build a rank-reduced operator basis grouped by canonical flavor tensor."""
    terms_by_index={int(term["term_index"]):term for term in flavor_seed["terms"]}
    seed_payload=json.loads(Path(context.metadata["_wilson_seed_path"]).read_text(encoding="utf-8"))
    from RGE.running.EFT1WilsonAdapter import build_eft1_wilson_tensor

    d_s1=int(context.metadata["dS1"])
    d_s2=int(context.metadata["dS2"])
    grouped={}
    for term_index in sorted(terms_by_index):
        term_meta=terms_by_index[term_index]
        label=str(term_meta.get("flavor_label") or _infer_flavor_label(term_meta))
        original_term=next(term for term in seed_payload["TreeWilsonTerms"] if int(term["Index"])==term_index)
        one_term_seed=dict(seed_payload)
        one_term_seed["TreeWilsonTerms"]=[original_term]
        tensor=build_eft1_wilson_tensor(
            one_term_seed,d_s1=d_s1,d_s2=d_s2,chirality="PL",project_operator_symmetry=True
        )
        kernel=_flavor_kernel_expr(term_meta["one_generation_reduction"])
        if sp.simplify(kernel)==0:
            raise ValueError(f"Zero flavor kernel for term {term_index}.")
        entries={key:sp.simplify(value/kernel) for key,value in _canonical_entries(tensor).items()}
        grouped.setdefault(label,[]).append({
            "term_index":term_index,
            "entries":entries,
            "original_term_input_form":term_meta["original_term_input_form"],
            "raw_mass_denominator_factor":term_meta.get("symmetry_denominator","1"),
            "cg_names":term_meta.get("cg_names",[]),
        })

    directions=[]
    for label in ("C11","C12","C22"):
        if label in grouped:
            directions.extend(_reduce_geometry_group(label,grouped[label]))
    if not directions:
        raise ValueError("No independent PL operator directions were found.")
    return directions


def _infer_flavor_label(term_meta):
    names=tuple(item["name"] for item in term_meta.get("yukawas",[]))
    if names==("y1","y1"):
        return "C11"
    if names==("y2","y2"):
        return "C22"
    if set(names)=={"y1","y2"}:
        return "C12"
    raise ValueError(f"Unsupported T3 Yukawa pair {names!r}.")


def _combined_coefficient(directions):
    entries={}
    symbols=tuple(sp.Symbol(direction["name"]) for direction in directions)
    for symbol,direction in zip(symbols,directions,strict=True):
        for key,value in direction["entries"].items():
            entries[key]=sp.simplify(entries.get(key,sp.S.Zero)+symbol*value)
            i,j,a,b=key
            for permuted in {(j,i,a,b),(i,j,b,a),(j,i,b,a)}:
                entries[permuted]=entries[key]
    return SparseWilsonTensor(entries), symbols


def _weinberg_geometry(context) -> dict[tuple[int, int, int, int], sp.Expr]:
    # Use exactly the same benchmark construction as EFT1WilsonRGE.py.
    # build_weinberg_wilson_tensor requires the scalar model and fermion
    # basis explicitly; kappa=1 fixes the normalization of the direction.
    sm_scalar_model, sm_fermion_basis = build_sm_eft()

    if sm_fermion_basis.dimension != context.fermion_basis.dimension:
        raise RuntimeError(
            "SM benchmark fermion basis does not match EFT1 basis."
        )

    tensor_dict = build_weinberg_wilson_tensor(
        sm_scalar_model,
        sm_fermion_basis,
        sp.S.One,
    )

    return {
        _canonical(component): sp.simplify(value)
        for component, value in tensor_dict.items()
        if sp.simplify(value) != 0
    }


def _project(
    beta: dict[tuple[int, int, int, int], sp.Expr],
    directions: list[dict[str, Any]],
    weinberg: dict[tuple[int, int, int, int], sp.Expr],
) -> dict[str, sp.Expr]:
    all_dirs = directions + [{"name": "Weinberg", "entries": weinberg}]
    components = sorted(
        set(beta).union(
            *(set(direction["entries"]) for direction in all_dirs)
        )
    )

    matrix = sp.Matrix(
        [
            [
                direction["entries"].get(component, sp.S.Zero)
                for direction in all_dirs
            ]
            for component in components
        ]
    )
    vector = sp.Matrix(
        [beta.get(component, sp.S.Zero) for component in components]
    )

    solution = sp.linsolve((matrix, vector))
    if not solution:
        basis = ", ".join(direction["name"] for direction in directions)
        raise ValueError(
            "Flavor-aware beta does not close on the rank-reduced operator "
            f"basis {{{basis}, Weinberg}}."
        )

    raw = next(iter(solution))
    if any(
        symbol.name.startswith("tau")
        for expr in raw
        for symbol in expr.free_symbols
    ):
        raise ValueError("Operator projection is underdetermined.")

    reconstructed = matrix * sp.Matrix(raw)
    residual = [
        sp.simplify(vector[i] - reconstructed[i])
        for i in range(len(components))
    ]
    if any(value != 0 for value in residual):
        raise ValueError("Nonzero residual after flavor operator projection.")

    return {
        direction["name"]: sp.simplify(value)
        for direction, value in zip(all_dirs, raw, strict=True)
    }


def _calculate_projected_beta(
    context,
    coefficient: SparseWilsonTensor,
    directions,
) -> dict[str, sp.Expr]:
    beta_components: dict[tuple[int, int, int, int], sp.Expr] = {}

    for component in closure_output_components(context):
        contributions = calculate_complete_master_rge(
            model=context.model,
            inputs=context.inputs,
            output_component=component,
            coefficient=coefficient,
            simplify_each=False,
        )
        value = sp.simplify(contributions["total"])
        if value != 0:
            beta_components[_canonical(component)] = value

    return _project(
        beta_components,
        directions,
        _weinberg_geometry(context),
    )


def _replace_abs_ye_squared(expr: sp.Expr) -> sp.Expr:
    """Replace ye*conjugate(ye) by a single marker Ye2."""
    expanded = sp.expand(expr)
    return sp.simplify(
        expanded.xreplace(
            {
                YE * sp.conjugate(YE): YE2,
                sp.conjugate(YE) * YE: YE2,
            }
        )
    )


def _split_flavor_action(expr: sp.Expr) -> tuple[sp.Expr, sp.Expr]:
    """Return flavor-blind and two-leg He-action coefficients.

    The input is linear in the rank-reduced operator-basis coefficients.  For each source coefficient Cx:

        coeff = a + b*|ye|^2

    is lifted to

        a Cx[p,q]
        + (b/2) (He.Cx + Cx.He^T)[p,q].

    The factor 1/2 is fixed because the latter reduces to
    2 |ye|^2 Cx in one generation.
    """

    expr = _replace_abs_ye_squared(expr)
    blind = sp.S.Zero
    action = sp.S.Zero

    # This function is called after selecting one source basis symbol,
    # so the remaining expression should be affine in Ye2.
    poly = sp.Poly(sp.expand(expr), YE2, domain="EX")
    if poly.degree() > 1:
        raise ValueError(
            "Unexpected higher-than-linear charged-lepton Yukawa dependence."
        )

    blind = sp.simplify(poly.coeff_monomial(1))
    ye2_coeff = sp.simplify(poly.coeff_monomial(YE2))
    action = sp.simplify(ye2_coeff / 2)

    return blind, action


def _operator_matrix(projected, basis_names, basis_symbols):
    result={}
    for target_name,expression in projected.items():
        target={}
        residual=sp.simplify(expression-sum(sp.diff(expression,sym)*sym for sym in basis_symbols))
        if residual!=0:
            raise ValueError(f"{target_name} beta is not linear in boundary C tensors: {residual}")
        for source_name,source_symbol in zip(basis_names,basis_symbols,strict=True):
            source_coeff=sp.simplify(sp.diff(expression,source_symbol))
            blind,action=_split_flavor_action(source_coeff)
            if blind!=0 or action!=0:
                target[source_name]={"flavor_blind":blind,"He_two_leg_action":action}
        result[target_name]=target
    return result


def _one_generation_reconstruct(matrix, basis_names, basis_symbols):
    result={}
    ye2=YE*sp.conjugate(YE)
    symbol_by_name=dict(zip(basis_names,basis_symbols,strict=True))
    for target,sources in matrix.items():
        expr=sp.S.Zero
        for source,pieces in sources.items():
            c=symbol_by_name[source]
            expr += pieces["flavor_blind"]*c
            expr += 2*pieces["He_two_leg_action"]*ye2*c
        result[target]=sp.simplify(expr)
    return result


def _serialise_matrix(
    matrix: dict[str, dict[str, dict[str, sp.Expr]]]
) -> dict[str, Any]:
    return {
        target: {
            source: {
                name: sp.sstr(sp.factor(value))
                for name, value in pieces.items()
            }
            for source, pieces in sources.items()
        }
        for target, sources in matrix.items()
    }


def run_flavor_rge(
    wilson_seed_path: Path,
    quartic_seed_path: Path,
    rgbeta_path: Path,
    flavor_seed_path: Path | None = None,
    output_path: Path | None = None,
) -> dict[str, Any]:
    context = build_eft1_rge_context(
        wilson_seed_path,
        quartic_seed_path,
        rgbeta_path,
    )

    # Keep these paths only as local reconstruction metadata.  They are not
    # part of the physics model.
    context.metadata["_wilson_seed_path"] = str(wilson_seed_path)
    context.metadata["_rgbeta_path"] = str(rgbeta_path)

    flavor_seed = _load_flavor_seed(
        flavor_seed_path,
        wilson_seed_path,
    )
    directions = _geometry_directions(context, flavor_seed)
    coefficient, basis_symbols = _combined_coefficient(directions)
    basis_names = tuple(direction["name"] for direction in directions)

    projected = _calculate_projected_beta(
        context,
        coefficient,
        directions,
    )
    matrix = _operator_matrix(projected, basis_names, basis_symbols)
    reconstructed = _one_generation_reconstruct(
        matrix, basis_names, basis_symbols
    )

    regression = {}
    all_match = True
    for target, exact in projected.items():
        residual = sp.simplify(
            _replace_abs_ye_squared(exact)
            - _replace_abs_ye_squared(reconstructed[target])
        )
        ok = residual == 0
        regression[target] = {
            "matches": ok,
            "residual": sp.sstr(residual),
        }
        all_match = all_match and ok

    payload = {
        "status": "Success" if all_match else "RegressionFailed",
        "basis": list(basis_names),
        "basis_terms": {
            direction["name"]: int(direction["term_index"])
            for direction in directions
        },
        "basis_metadata": {
            direction["name"]: {
                "flavor_label": direction["flavor_label"],
                "boundary_weight": sp.sstr(
                    sp.factor(direction["boundary_weight"])
                ),
                "representative_term_index": int(direction["term_index"]),
                "source_term_indices": [
                    int(i) for i in direction["source_term_indices"]
                ],
                "raw_group_count": int(direction["raw_group_count"]),
                "group_rank": int(direction["group_rank"]),
                "raw_mass_denominator_factor":
                    direction["raw_mass_denominator_factor"],
                "original_term_input_form":
                    direction["original_term_input_form"],
                "cg_names": direction.get("cg_names", []),
            }
            for direction in directions
        },
        "generated_basis": ["Weinberg"],
        "flavor_equation": (
            "16*pi^2 beta_X[p,q] = "
            "sum_Y A_XY C_Y[p,q] + "
            "sum_Y B_XY (He.C_Y + C_Y.He^T)[p,q], "
            "He = Ye Ye^dagger"
        ),
        "operator_rge": _serialise_matrix(matrix),
        "one_generation_regression": regression,
        "one_generation_reduction_matches": all_match,
        "projected_one_generation_beta": {
            name: sp.sstr(sp.factor(expr))
            for name, expr in projected.items()
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
            "Build the flavor-aware EFT1 T3 Wilson operator-space RGE and "
            "verify its one-generation reduction."
        )
    )
    parser.add_argument("wilson_seed_json", type=Path)
    parser.add_argument("quartic_seed_json", type=Path)
    parser.add_argument("rgbeta_json", type=Path)
    parser.add_argument("--flavor-seed", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    result = run_flavor_rge(
        args.wilson_seed_json,
        args.quartic_seed_json,
        args.rgbeta_json,
        flavor_seed_path=args.flavor_seed,
        output_path=args.output,
    )

    print(
        json.dumps(
            {
                "status": result["status"],
                "basis": result["basis"],
                "generated_basis": result["generated_basis"],
                "one_generation_reduction_matches":
                    result["one_generation_reduction_matches"],
                "nonzero_channels": sum(
                    len(sources)
                    for sources in result["operator_rge"].values()
                ),
            },
            indent=2,
        )
    )

    return 0 if result["status"] == "Success" else 1


if __name__ == "__main__":
    raise SystemExit(main())
