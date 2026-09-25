from __future__ import annotations

"""Complete one-loop ``psi^2 phi^2`` RGE in the scalar-only T3 EFT.

Seed support means nonzero matched tensor components at the fermion threshold.
Closure scans the allowed output components so operator mixing can generate
components that vanish at the boundary.
"""

from dataclasses import dataclass
from itertools import combinations_with_replacement
import json
from pathlib import Path
import sys
import time


PROJECT_ROOT = Path(__file__).resolve().parents[3]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import sympy as sp

from RGE.general.WilsonTensorRGE import (
    WilsonRGEInputs,
    calculate_complete_wilson_tensor_rge,
)
from RGE.general.ScalarBasis import ScalarBasis
from RGE.general.FermionBasis import build_gauge_sectors
from RGE.matching.WeinbergTensorAdapter import (
    SparseWilsonLookup,
    build_sm_eft,
    build_sm_yukawa,
    build_weinberg_wilson_tensor,
)
from RGE.running.intermediate.ScalarOnlyTensorAdapters import (
    load_and_build_scalar_only_quartic_tensor,
    load_and_build_scalar_only_wilson_tensor,
)


@dataclass(frozen=True)
class ScalarOnlyRGEContext:
    model: ScalarBasis
    fermion_basis: object
    coefficient: object
    quartic: object
    inputs: WilsonRGEInputs
    metadata: dict


def _parse_rational(value) -> sp.Rational:
    """Parse Matchete/RGBeta metadata values such as ``-1/2`` exactly."""
    if isinstance(value, int):
        return sp.Rational(value)
    return sp.Rational(str(value))


def _canonical_component(
    component: tuple[int, int, int, int],
) -> tuple[int, int, int, int]:
    i, j, a, b = component
    return (min(i, j), max(i, j), min(a, b), max(a, b))


def _symmetry_residual_count(coefficient) -> int:
    failures = 0
    keys = set(coefficient.entries)
    for i, j, a, b in tuple(keys):
        keys.update(
            {
                (j, i, a, b),
                (i, j, b, a),
                (j, i, b, a),
            }
        )
    for i, j, a, b in keys:
        reference = coefficient[i, j, a, b]
        if sp.simplify(reference - coefficient[j, i, a, b]) != 0:
            failures += 1
        if sp.simplify(reference - coefficient[i, j, b, a]) != 0:
            failures += 1
    return failures


def build_scalar_only_rge_context(
    wilson_seed_path: Path,
    quartic_seed_path: Path,
    rgbeta_path: Path,
) -> ScalarOnlyRGEContext:
    """Build the scalar-only intermediate EFT tensor-RGE context."""
    rgbeta = json.loads(Path(rgbeta_path).read_text(encoding="utf-8"))
    metadata = rgbeta["metadata"]

    d_s1 = int(metadata["dS1"])
    d_s2 = int(metadata["dS2"])
    y_s1 = _parse_rational(metadata["YS1"])
    y_s2 = _parse_rational(metadata["YS2"])
    shared_scalar = bool(metadata.get("SharedScalar", False))

    if shared_scalar:
        scalar_model = ScalarBasis.t3_shared(
            d_s=d_s2, y_s=y_s2, include_higgs=True
        )
    else:
        scalar_model = ScalarBasis.t3(
            d_s1=d_s1,
            y_s1=y_s1,
            d_s2=d_s2,
            y_s2=y_s2,
            include_higgs=True,
        )

    # F has already been removed, so the active fermion basis is the validated
    # one-generation SM Weyl basis.
    _, fermion_basis = build_sm_eft()

    l1 = fermion_basis.global_index("L", 1)
    l2 = fermion_basis.global_index("L", 2)

    coefficient = load_and_build_scalar_only_wilson_tensor(
        wilson_seed_path,
        rgbeta_path,
        chirality="PL",
        lepton_indices=(l1, l2),
        project_operator_symmetry=True,
    )

    residual_count = _symmetry_residual_count(coefficient)
    if residual_count:
        raise RuntimeError(
            "Projected scalar-only Wilson tensor is not separately symmetric: "
            f"{residual_count} residual(s)."
        )

    quartic = load_and_build_scalar_only_quartic_tensor(
        quartic_seed_path,
        rgbeta_path,
    )

    yukawa = build_sm_yukawa(
        scalar_model,
        fermion_basis,
    )

    inputs = WilsonRGEInputs(
        fermion_dimension=fermion_basis.dimension,
        yukawa=yukawa,
        quartic=quartic,
        gauge_sectors=build_gauge_sectors(
            scalar_model,
            fermion_basis,
        ),
    )

    return ScalarOnlyRGEContext(
        model=scalar_model,
        fermion_basis=fermion_basis,
        coefficient=coefficient,
        quartic=quartic,
        inputs=inputs,
        metadata=metadata,
    )


def seed_support_components(context: ScalarOnlyRGEContext):
    """Return independent components that are nonzero at the boundary."""
    return sorted(
        {
            _canonical_component(component)
            for component, value in context.coefficient.nonzero_items()
            if sp.simplify(value) != 0
        }
    )


def closure_output_components(context: ScalarOnlyRGEContext):
    """Return the L/eC one-loop closure component set."""
    fb = context.fermion_basis
    fermions = sorted(
        {
            fb.global_index("L", 1),
            fb.global_index("L", 2),
            fb.global_index("eC", 1),
        }
    )
    scalars = range(1, context.model.total_real_scalar_dimension + 1)

    return [
        (i, j, a, b)
        for i, j in combinations_with_replacement(fermions, 2)
        for a, b in combinations_with_replacement(scalars, 2)
    ]


def calculate_component_betas(
    context: ScalarOnlyRGEContext,
    *,
    seed_support_only: bool = False,
    simplify_each: bool = False,
) -> dict[tuple[int, int, int, int], dict[str, sp.Expr]]:
    """Evaluate the complete one-loop beta tensor on an independent component set."""

    components = (
        seed_support_components(context)
        if seed_support_only
        else closure_output_components(context)
    )

    result: dict[
        tuple[int, int, int, int],
        dict[str, sp.Expr],
    ] = {}

    for component in components:
        contributions = calculate_complete_wilson_tensor_rge(
            model=context.model,
            inputs=context.inputs,
            output_component=component,
            coefficient=context.coefficient,
            simplify_each=simplify_each,
        )

        total = sp.simplify(contributions["total"])
        if total != 0:
            result[component] = {
                **contributions,
                "total": total,
            }

    return result


def _expr_text(expression: sp.Expr) -> str:
    return sp.sstr(sp.factor(sp.simplify(expression)))


def _serialise_component(component) -> str:
    return ",".join(str(index) for index in component)


def weinberg_subspace_diagnostics(
    context: ScalarOnlyRGEContext,
    betas: dict[tuple[int, int, int, int], dict[str, sp.Expr]],
) -> dict:
    """Check whether the generated LL-HH beta is exactly a Weinberg operator.

    With template kappa=1, a proportionality factor ``r`` means

        16*pi^2 d(kappa)/dln(mu) = r

    for the part generated by scalar-only intermediate running.
    """

    sm_scalar_model, sm_fermion_basis = build_sm_eft()

    if sm_fermion_basis.dimension != context.fermion_basis.dimension:
        raise RuntimeError(
            "SM benchmark fermion basis does not match the intermediate basis."
        )

    template_dict = build_weinberg_wilson_tensor(
        sm_scalar_model,
        sm_fermion_basis,
        sp.S.One,
    )
    template = SparseWilsonLookup(template_dict)

    l1 = context.fermion_basis.global_index("L", 1)
    l2 = context.fermion_basis.global_index("L", 2)
    h_indices = tuple(sm_scalar_model.block("H").indices)

    components = [
        (i, j, a, b)
        for i, j in combinations_with_replacement((l1, l2), 2)
        for a, b in combinations_with_replacement(h_indices, 2)
    ]

    def beta_value(component):
        return sp.simplify(
            betas.get(_canonical_component(component), {}).get(
                "total",
                sp.S.Zero,
            )
        )

    template_nonzero = [
        component
        for component in components
        if sp.simplify(template[component]) != 0
    ]

    if not template_nonzero:
        raise RuntimeError(
            "Weinberg normalization template has no nonzero components."
        )

    reference_component = template_nonzero[0]
    beta_kappa = sp.simplify(
        beta_value(reference_component) / template[reference_component]
    )

    failures = {}
    for component in components:
        expected = sp.simplify(beta_kappa * template[component])
        residual = sp.simplify(beta_value(component) - expected)
        if residual != 0:
            failures[_serialise_component(component)] = _expr_text(residual)

    return {
        "matches_weinberg_subspace": not failures,
        "failure_count": len(failures),
        "template_nonzero_component_count": len(template_nonzero),
        "reference_component": list(reference_component),
        "beta_kappa_16pi2": _expr_text(beta_kappa),
        "residuals": failures,
    }


def run_scalar_only_wilson_rge(
    wilson_seed_path: Path,
    quartic_seed_path: Path,
    rgbeta_path: Path,
    output_path: Path | None = None,
    *,
    seed_support_only: bool = False,
) -> dict:
    """Evaluate and serialize the scalar-only dimension-five Wilson RGE."""
    started = time.perf_counter()

    context = build_scalar_only_rge_context(
        wilson_seed_path=wilson_seed_path,
        quartic_seed_path=quartic_seed_path,
        rgbeta_path=rgbeta_path,
    )

    candidates = (
        seed_support_components(context)
        if seed_support_only
        else closure_output_components(context)
    )

    betas = calculate_component_betas(
        context,
        seed_support_only=seed_support_only,
        simplify_each=False,
    )

    initial_support = {
        _canonical_component(component)
        for component, value in context.coefficient.nonzero_items()
        if sp.simplify(value) != 0
    }

    generated = sorted(
        component
        for component in betas
        if component not in initial_support
    )

    weinberg_diagnostics = (
        None
        if seed_support_only
        else weinberg_subspace_diagnostics(context, betas)
    )

    elapsed = time.perf_counter() - started

    payload = {
        "status": "Success",
        "theory": (
            "SM + S + tree-generated dimension-5 psi2phi2"
            if bool(context.metadata.get("SharedScalar", False))
            else "SM + S1 + S2 + tree-generated dimension-5 psi2phi2"
        ),
        "loop_order": 1,
        # Historical serialized wording retained for output compatibility.
        "fixed_order_input": "tree-level EFT1 Wilson coefficients only",
        "scan_mode": (
            "seed_support_only"
            if seed_support_only
            else "L/eC one-loop closure"
        ),
        "metadata": {
            "SharedScalar": bool(context.metadata.get("SharedScalar", False)),
            "dS": context.metadata.get("dS"),
            "dS1": int(context.metadata["dS1"]),
            "dS2": int(context.metadata["dS2"]),
            "YS1": context.metadata["YS1"],
            "YS2": context.metadata["YS2"],
            "YS": context.metadata.get("YS"),
            "real_scalar_dimension": context.model.total_real_scalar_dimension,
            "fermion_dimension": context.fermion_basis.dimension,
        },
        "candidate_component_count": len(candidates),
        "initial_independent_component_count": len(initial_support),
        "nonzero_beta_component_count": len(betas),
        "generated_component_count": len(generated),
        "generated_components": [list(item) for item in generated],
        "weinberg_subspace_validation": weinberg_diagnostics,
        "elapsed_seconds": elapsed,
        "betas": {
            _serialise_component(component): {
                name: _expr_text(value)
                for name, value in contributions.items()
            }
            for component, contributions in sorted(betas.items())
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


# Compatibility alias for callers that still use the historical function name.
run_eft1_wilson_rge = run_scalar_only_wilson_rge
build_eft1_rge_context = build_scalar_only_rge_context
EFT1RGEContext = ScalarOnlyRGEContext


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description=(
            "Evaluate the complete one-loop scalar-only psi^2 phi^2 master RGE."
        )
    )
    parser.add_argument("wilson_seed_json", type=Path)
    parser.add_argument("quartic_seed_json", type=Path)
    parser.add_argument("rgbeta_json", type=Path)
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Optional JSON output path.",
    )
    parser.add_argument(
        "--seed-support-only",
        action="store_true",
        help=(
            "Fast plumbing test: evaluate only components nonzero at the "
            "matching boundary instead of the full L/eC closure."
        ),
    )
    args = parser.parse_args()

    payload = run_scalar_only_wilson_rge(
        wilson_seed_path=args.wilson_seed_json,
        quartic_seed_path=args.quartic_seed_json,
        rgbeta_path=args.rgbeta_json,
        output_path=args.output,
        seed_support_only=args.seed_support_only,
    )

    summary = {
        key: payload[key]
        for key in (
            "status",
            "scan_mode",
            "candidate_component_count",
            "initial_independent_component_count",
            "nonzero_beta_component_count",
            "generated_component_count",
            "elapsed_seconds",
        )
    }
    print(json.dumps(summary, indent=2))

    if payload["weinberg_subspace_validation"] is not None:
        print(
            json.dumps(
                {
                    "weinberg_subspace_validation":
                        payload["weinberg_subspace_validation"]
                },
                indent=2,
            )
        )

    for component, contributions in list(payload["betas"].items())[:8]:
        print(f"{component}: {contributions['total']}")
