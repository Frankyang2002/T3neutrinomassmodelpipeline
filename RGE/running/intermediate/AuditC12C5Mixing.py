"""Independent two-seed audit of dimension-five operator mixing in scalar T3 EFT.

Diagnostic only: never changes production matching, RGEs, or numerical results.

Usage:
 python -m RGE.running.intermediate.AuditC12C5Mixing --data-dir PATH

The context is the existing *one-generation* symbolic representation. The C12
seed is the mixed-y1-y2 part of the Matchete-derived LLSS tensor, while the C5
seed is an independent Weinberg operator embedded in the same real-scalar basis.
A nonzero cross block means the corresponding source generates a component of
the other operator class; it does NOT establish a complete closed 2x2 EFT basis.
"""
from __future__ import annotations

import argparse
from dataclasses import replace
import json
from pathlib import Path

import sympy as sp

from RGE.matching.WeinbergTensorAdapter import (
    SparseWilsonLookup,
    build_weinberg_wilson_tensor,
)
from RGE.running.intermediate.ScalarOnlyWilsonTensorRGE import (
    build_scalar_only_rge_context,
    calculate_component_betas,
)


def canonical(key: tuple[int, int, int, int]) -> tuple[int, int, int, int]:
    i, j, a, b = key
    return min(i, j), max(i, j), min(a, b), max(a, b)


def complete_symmetric_entries(mapping):
    """Canonicalise a separately symmetric fermion/scalar tensor, checking clashes."""
    canonical_values = {}
    for key, raw in mapping.items():
        key = canonical(tuple(key))
        value = sp.simplify(raw)
        if key in canonical_values:
            if sp.simplify(canonical_values[key] - value) != 0:
                raise ValueError(f"Tensor symmetry conflict at {key}")
        else:
            canonical_values[key] = value
    full = {}
    for (i, j, a, b), v in canonical_values.items():
        if v == 0:
            continue
        for key in ((i, j, a, b), (j, i, a, b),
                    (i, j, b, a), (j, i, b, a)):
            full[key] = v
    return full


def mixed_yukawa_part(expr):
    """Remove terms surviving when either y1 or y2 is set to zero."""
    y1, y2 = sp.symbols("y1 y2")
    e = sp.expand(expr)
    zero1 = {y1: 0, sp.conjugate(y1): 0}
    zero2 = {y2: 0, sp.conjugate(y2): 0}
    return sp.simplify(
        e - e.subs(zero1) - e.subs(zero2) + e.subs({**zero1, **zero2})
    )


def make_seeds(context):
    raw = getattr(context.coefficient, "entries", None)
    if raw is None:
        raise TypeError("Expected Matchete-derived sparse Wilson entries")
    mixed = {
        k: v
        for k, original in raw.items()
        if (v := mixed_yukawa_part(original)) != 0
    }
    c12 = complete_symmetric_entries(mixed)
    kappa = sp.Symbol("kappa")
    c5 = complete_symmetric_entries(
        build_weinberg_wilson_tensor(
            context.model, context.fermion_basis, kappa
        )
    )
    if not c12:
        raise ValueError("No mixed y1-y2 LLSS tensor components found")
    if not c5:
        raise ValueError("No embedded Weinberg tensor components found")
    return c12, c5


def operator_components(context):
    fb = context.fermion_basis
    l = (fb.global_index("L", 1), fb.global_index("L", 2))
    h = set(context.model.block("H").indices)
    scalars = range(1, context.model.total_real_scalar_dimension + 1)
    c5 = []
    llss = []
    other = []
    for i in l:
        for j in l:
            if j < i:
                continue
            for a in scalars:
                for b in range(a, context.model.total_real_scalar_dimension + 1):
                    key = (i, j, a, b)
                    if a in h and b in h:
                        c5.append(key)
                    elif a not in h and b not in h:
                        llss.append(key)
                    else:
                        other.append(key)
    return {"C5_HH": c5, "LLSS": llss, "mixed_HS": other}


def seed_projection_report(source, target, beta_map, target_keys):
    """Compute component-level projection, plus failures of exact proportionality.

    For the LLSS class, a single C12 basis tensor need not span every allowed
    invariant. We report that distinction rather than discarding residuals.
    """
    generated = {
        ",".join(map(str, k)): sp.sstr(sp.factor(sp.simplify(
            beta_map.get(k, {}).get("total", 0)
        )))
        for k in target_keys
        if sp.simplify(beta_map.get(k, {}).get("total", 0)) != 0
    }
    reference = [k for k in target_keys
                 if sp.simplify(target.get(k, 0)) != 0]
    ratio = None
    residuals = {}
    if reference:
        pivot = reference[0]
        b = sp.simplify(beta_map.get(pivot, {}).get("total", 0))
        ratio = sp.simplify(b / target[pivot])
        for k in target_keys:
            value = sp.simplify(beta_map.get(k, {}).get("total", 0))
            expected = ratio * target.get(k, 0)
            residual = sp.simplify(value - expected)
            if residual != 0:
                residuals[",".join(map(str,k))] = sp.sstr(residual)
    return {
        "source": source,
        "nonzero_target_components": len(generated),
        "generated_components": generated,
        "reference_tensor_component_count": len(reference),
        "proportionality_factor": sp.sstr(ratio) if ratio is not None else None,
        "proportionality_residual_count": len(residuals),
        "proportionality_residuals": residuals,
    }


def audit(context):
    c12, c5 = make_seeds(context)
    groups = operator_components(context)
    results = {}
    for label, seed in (("C12_only", c12), ("C5_only", c5)):
        specialised = replace(context, coefficient=SparseWilsonLookup(seed))
        # Restrict output to external LL channels. No automatic restriction to
        # initially nonzero support: this catches cross-operator generation.
        from RGE.general.WilsonTensorRGE import calculate_complete_wilson_tensor_rge
        betas = {}
        for keys in groups.values():
            for key in keys:
                contributions = calculate_complete_wilson_tensor_rge(
                    model=specialised.model, inputs=specialised.inputs,
                    output_component=key, coefficient=specialised.coefficient,
                    simplify_each=False,
                )
                if sp.simplify(contributions["total"]) != 0:
                    betas[key] = contributions
        results[label] = {}
        for name, keys in groups.items():
            template = c5 if name == "C5_HH" else c12 if name == "LLSS" else {}
            results[label][name] = seed_projection_report(
                label, template, betas, keys
            )
    return {
        "status": "AuditComplete",
        "scope": "one-generation component-level audit, not a full-flavor RGE",
        "c12_seed_components": len(c12),
        "c5_seed_components": len(c5),
        "mixing": results,
        "warning": ("An LLSS component need not lie in the one-dimensional C12 "
                    "subspace. Check proportionality residuals before interpreting "
                    "the entries as a closed 2x2 anomalous-dimension matrix."),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, required=True,
                        help="Directory containing eft1_after_F_*_seed.json and eft1_rgbeta_rge.json")
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()
    d = args.data_dir
    context = build_scalar_only_rge_context(
        d / "eft1_after_F_wilson_seed.json",
        d / "eft1_after_F_scalar_quartic_seed.json",
        d / "eft1_rgbeta_rge.json",
    )
    result = audit(context)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps({
        "status": result["status"],
        "c12_seed_components": result["c12_seed_components"],
        "c5_seed_components": result["c5_seed_components"],
        "counts": {s: {t: v["nonzero_target_components"]
                      for t, v in targets.items()}
                   for s, targets in result["mixing"].items()},
        "proportionality_residual_counts": {
            s: {t: v["proportionality_residual_count"]
                for t, v in targets.items()}
            for s, targets in result["mixing"].items()},
        "full_output": str(args.output) if args.output else "not saved",
    }, indent=2))


if __name__ == "__main__":
    main()
