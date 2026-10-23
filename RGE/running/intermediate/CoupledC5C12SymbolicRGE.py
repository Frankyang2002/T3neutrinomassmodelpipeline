"""Exact one-loop EFT1 C5 <-> LLSS symbolic mixing.

This is an independent symbolic block calculation, not the fixed-order
production threshold flow. It reuses the existing WilsonTensorRGE kernel.
"""
from __future__ import annotations

import argparse
from dataclasses import replace
import json
from pathlib import Path
import sympy as sp

from RGE.general.WilsonTensorRGE import calculate_complete_wilson_tensor_rge
from RGE.matching.WeinbergTensorAdapter import build_weinberg_wilson_tensor, SparseWilsonLookup
from RGE.running.intermediate.ScalarOnlyWilsonTensorRGE import (
    build_scalar_only_rge_context,
)


def canon(key):
    i, j, a, b = key
    return (min(i, j), max(i, j), min(a, b), max(a, b))


def simp(value):
    return sp.factor(sp.simplify(value))


def expr(value):
    return sp.sstr(simp(value))


def canonical_entries(tensor):
    output = {}
    for key, value in tensor.items():
        key = canon(key)
        if key in output and simp(output[key] - value) != 0:
            raise ValueError(f"Inconsistent symmetric components for {key}")
        output[key] = simp(value)
    return {k: v for k, v in output.items() if v != 0}


def matched_components(ctx):
    return canonical_entries(ctx.coefficient.entries)


def evaluate_block(ctx, source, components):
    inputs = replace(ctx.inputs)
    # The Wilson coefficient is supplied independently of the renormalisable
    # inputs, so beta is linear in this source (plus its complex conjugate).
    terms = {}
    for key in components:
        res = calculate_complete_wilson_tensor_rge(
            model=ctx.model, inputs=inputs, output_component=key,
            coefficient=source, simplify_each=False,
        )
        terms[key] = {name: simp(val) for name, val in res.items()}
    return terms


def projection(template, beta, components):
    """Exact reference-component projection with an explicit closure residual.

    The result uses the template's original Matchete normalisation.
    Nonzero residuals mean the one-dimensional subspace is not closed.
    """
    support = [k for k in components if template.get(k, 0) != 0]
    if not support:
        raise ValueError("Projection target has no nonzero components")
    ref = support[0]
    ratio = simp(beta.get(ref, sp.S.Zero) / template[ref])
    residuals = {
        ",".join(map(str, k)): expr(
            beta.get(k, sp.S.Zero) - ratio * template.get(k, sp.S.Zero)
        )
        for k in components
        if simp(beta.get(k, sp.S.Zero) -
                ratio * template.get(k, sp.S.Zero)) != 0
    }
    return {
        "reference_component": list(ref),
        "coefficient": expr(ratio),
        "closed": not residuals,
        "residual_component_count": len(residuals),
        "residuals": residuals,
    }


def compute(wilson_seed, quartic_seed, rgbeta, *, only_mixed=False):
    ctx = build_scalar_only_rge_context(wilson_seed, quartic_seed, rgbeta)
    c12 = matched_components(ctx)
    c5 = canonical_entries(build_weinberg_wilson_tensor(
        ctx.model, ctx.fermion_basis, sp.Integer(1)))
    if set(c12) & set(c5):
        raise ValueError("Weinberg and matched LLSS templates overlap")

    if only_mixed:
        outputs = sorted(set(c12) | set(c5))
    else:
        from RGE.running.intermediate.ScalarOnlyWilsonTensorRGE import (
            closure_output_components,
        )
        outputs = closure_output_components(ctx)

    result = {
        "status": "Success",
        "kind": "symbolic_operator_mixing_blocks",
        "convention": "16*pi^2 dC/dln(mu) = beta[C]",
        "normalisation": (
            "C5 template kappa=1; C12 template is the Matchete PL seed "
            "including its symbolic Yukawa/MF prefactor"
        ),
        "limitations": (
            "The C12 input is a one-generation matched direction, not "
            "an independently normalised complete LLSS operator basis. "
            "Residuals identify non-closure."
        ),
        "metadata": ctx.metadata,
        "input_templates": {
            "C5": {",".join(map(str,k)): expr(v) for k,v in sorted(c5.items())},
            "C12": {",".join(map(str,k)): expr(v) for k,v in sorted(c12.items())},
        },
        "blocks": {},
    }

    targets = {"C5": c5, "C12": c12}
    for source_name, entries in targets.items():
        expanded = {}
        for (i,j,a,b), value in entries.items():
            for key in ((i,j,a,b),(j,i,a,b),(i,j,b,a),(j,i,b,a)):
                expanded[key] = value
        source = SparseWilsonLookup(expanded)
        betas = evaluate_block(ctx, source, outputs)
        serialized = {
            ",".join(map(str,k)): {
                name: expr(value) for name,value in d.items() if value != 0
            }
            for k,d in betas.items() if d["total"] != 0
        }
        result["blocks"][source_name] = {
            "beta_components": serialized,
            "projections": {}
        }
        for target_name, template in targets.items():
            # Project within the target scalar sector, but retain all canonical
            # L/eC components for the complete closure diagnostic.
            region = set(k for k in outputs if k[2] in
                         ({1,2,3,4} if target_name == "C5" else
                          set(range(5,ctx.model.total_real_scalar_dimension+1)))
                         and k[3] in
                         ({1,2,3,4} if target_name == "C5" else
                          set(range(5,ctx.model.total_real_scalar_dimension+1))))
            region |= set(template)
            region = sorted(region)
            totals = {k: betas.get(k,{}).get("total",sp.S.Zero) for k in region}
            result["blocks"][source_name]["projections"][target_name] = (
                projection(template, totals, region))
    return result


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("wilson_seed", type=Path)
    ap.add_argument("quartic_seed", type=Path)
    ap.add_argument("eft1_rgbeta", type=Path)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--only-mixed", action="store_true",
                    help="Fast check of C5 and matched C12 output supports only")
    args = ap.parse_args()
    data = compute(args.wilson_seed, args.quartic_seed, args.eft1_rgbeta,
                   only_mixed=args.only_mixed)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(data, indent=2)+"\\n", encoding="utf-8")
    print("Saved:",args.output)
    for src, block in data["blocks"].items():
        for tgt, proj in block["projections"].items():
            print(f"{src} -> {tgt}: {proj['coefficient']} "
                  f"closed={proj['closed']} residuals={proj['residual_component_count']}")


if __name__ == "__main__":
    main()
