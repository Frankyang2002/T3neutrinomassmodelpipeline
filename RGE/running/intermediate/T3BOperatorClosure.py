"""T3-B one-generation operator-subspace closure via the production master RGE.

This is an exact symbolic test **within** the psi^2 phi^2 real-tensor sector
implemented by RGE.general.WilsonTensorRGE, not a proof of all EFT operators
or of 3-generation closure.

Start from independently matched LLSS and embedded LLHH Weinberg tensors,
evolve their one-loop beta tensors, and test whether either result belongs to
the two-template span. Optional --feedback computes a further beta application
to the *unretained* residual and tests return to LLSS/LLHH support.

No production code or numerical benchmarks are changed.
"""
from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import replace
import json
from pathlib import Path
import sys
import time
from typing import Callable, Mapping

import sympy as sp

Key = tuple[int, int, int, int]

def canonical(key: Key) -> Key:
    i, j, a, b = key
    return (min(i,j), max(i,j), min(a,b), max(a,b))


def reduced(expr: sp.Expr) -> sp.Expr:
    """Zero testing uses the same exact SymPy simplification as master RGE."""
    return sp.simplify(expr)


def canonical_values(tensor: Mapping[Key, sp.Expr]) -> dict[Key, sp.Expr]:
    """Require tensor consistency under independent fermion/scalar swaps."""
    unique: dict[Key, sp.Expr] = {}
    for raw_key, value in tensor.items():
        key = canonical(tuple(map(int, raw_key)))
        value = reduced(sp.sympify(value))
        if key in unique and reduced(unique[key] - value) != 0:
            raise ValueError(f"Inconsistent tensor symmetry at {key}")
        unique[key] = value
    return {key: value for key, value in unique.items() if value != 0}


def project_onto_two_templates(
    beta: Mapping[Key, sp.Expr],
    templates: Mapping[str, Mapping[Key, sp.Expr]],
    components: list[Key],
) -> tuple[dict[str, sp.Expr], dict[Key, sp.Expr]]:
    """Project by disjoint LLSS/LLHH pivot entries; retain *all* residuals.

    This is an algebraic-coordinate projection: nonzero residual => not closed.
    A nonzero coefficient alone is not a physical projection when residuals
    remain; the full operator basis must then be determined.
    """
    keys = [canonical(k) for k in components]
    if len(set(keys)) != len(keys):
        raise ValueError("components must have unique canonical entries")
    clean = {name: canonical_values(t) for name, t in templates.items()}
    if set(clean) != {"C12", "C5"}:
        raise ValueError("Expected exactly C12 and C5 template tensors")
    if set(clean["C12"]) & set(clean["C5"]):
        raise ValueError("LLSS and LLHH templates must have disjoint supports")
    for name, t in clean.items():
        if not t:
            raise ValueError(f"Empty {name} template")
        if not set(t).issubset(keys):
            raise ValueError(f"{name} support is not in the scanned component set")
    normalized_beta = canonical_values(beta)
    if not set(normalized_beta).issubset(keys):
        raise ValueError("Beta contains entries outside scanned components")

    factors = {}
    for name in ("C12", "C5"):
        pivot = next(key for key in keys if key in clean[name])
        factors[name] = reduced(normalized_beta.get(pivot, 0) / clean[name][pivot])

    residual = {}
    for key in keys:
        value = reduced(
            normalized_beta.get(key, 0)
            - factors["C12"] * clean["C12"].get(key, 0)
            - factors["C5"] * clean["C5"].get(key, 0)
        )
        if value != 0:
            residual[key] = value
    return factors, residual


def classify(key: Key, context) -> str:
    """Classify the external fields; no invariant-operator claim is implied."""
    i, j, a, b = key
    fb = context.fermion_basis
    lepton = {fb.global_index("L", 1), fb.global_index("L", 2)}
    if i not in lepton or j not in lepton:
        return "has_eC_external"
    h = set(context.model.block("H").indices)
    if a in h and b in h:
        return "LL_HH"
    if a not in h and b not in h:
        return "LL_SS"
    return "LL_HS"


def evaluated_betas(context, coefficient: Mapping[Key, sp.Expr],
                    components: list[Key], *, report_every: int = 30,
                    label: str = "") -> dict[Key, sp.Expr]:
    """Evaluate the full existing production tensor beta (not a surrogate)."""
    from RGE.running.intermediate.AuditC12C5Mixing import complete_symmetric_entries
    from RGE.matching.WeinbergTensorAdapter import SparseWilsonLookup
    from RGE.general.WilsonTensorRGE import calculate_complete_wilson_tensor_rge
    source = SparseWilsonLookup(complete_symmetric_entries(coefficient))
    scoped = replace(context, coefficient=source)
    betas: dict[Key, sp.Expr] = {}
    t0 = time.monotonic()
    for index, output in enumerate(components, start=1):
        result = calculate_complete_wilson_tensor_rge(
            model=scoped.model, inputs=scoped.inputs,
            output_component=output, coefficient=scoped.coefficient,
            simplify_each=False,
        )
        value = reduced(result["total"])
        if value != 0:
            betas[canonical(output)] = value
        if report_every and (index % report_every == 0 or index == len(components)):
            print(f"[{label}] {index}/{len(components)} components, "
                  f"{time.monotonic()-t0:.1f}s", file=sys.stderr, flush=True)
    return betas


def _text(expr, max_chars=1600):
    s = sp.sstr(sp.factor(expr))
    return s if len(s) <= max_chars else s[:max_chars] + "...[truncated]"


def _key(key: Key) -> str:
    return ",".join(map(str, key))


def residual_report(residual: Mapping[Key, sp.Expr], context, max_examples=16) -> dict:
    by_class = Counter(classify(k, context) for k in residual)
    examples = []
    for k in sorted(residual)[:max_examples]:
        examples.append({"component": list(k),
                         "class": classify(k, context),
                         "expression": _text(residual[k])})
    return {"component_count": len(residual),
            "component_counts_by_class": dict(sorted(by_class.items())),
            "examples": examples,
            "examples_truncated": len(residual) > max_examples}


def run_closure(context, *, scope: str = "llec", feedback: bool = False,
                max_examples: int = 16, report_every: int = 30) -> dict:
    """Check full psi2phi2 two-template closure, optionally 2nd-order return.

    ``ll`` is a faster explicitly PARTIAL scan; ``llec`` tests all one-gen
    external fermion pairs used by the repo's scalar-only closure framework.
    Feedback evaluates beta(residual) on LLSS and LLHH candidate components;
    its output is a second-order diagnostic, not a complete resummation.
    """
    if scope not in {"ll", "llec"}:
        raise ValueError("scope must be 'll' or 'llec'")
    from RGE.running.intermediate.AuditC12C5Mixing import make_seeds
    from RGE.running.intermediate.ScalarOnlyWilsonTensorRGE import closure_output_components
    c12, c5 = make_seeds(context)
    templates = {"C12": c12, "C5": c5}
    components = closure_output_components(context)
    if scope == "ll":
        components = [k for k in components
                      if classify(k, context) != "has_eC_external"]
    keys = [canonical(k) for k in components]
    retained = [k for k in keys if classify(k, context) in {"LL_HH", "LL_SS"}]

    results = {}
    for name, tensor in templates.items():
        betas = evaluated_betas(context, tensor, components,
                               label=f"beta({name})", report_every=report_every)
        factors, residual = project_onto_two_templates(betas, templates, components)
        entry = {
            "nonzero_beta_components": len(betas),
            "pivot_coefficients": {k: _text(v) for k, v in factors.items()},
            "residual": residual_report(residual, context, max_examples),
            "closed_in_scanned_sector": len(residual) == 0,
        }
        if feedback:
            # Only residuals outside the putative two-template span are used;
            # the retained pieces are deliberately NOT inserted again.
            if residual:
                second_betas = evaluated_betas(
                    context, residual, retained,
                    label=f"beta(residual({name}))",
                    report_every=report_every,
                )
                pivot, second_residual = project_onto_two_templates(
                    second_betas, templates, retained
                )
                entry["second_step_return"] = {
                    "description": "One further master-RGE action on unretained tensor; "
                                   "pivot is not a physical mixing coefficient unless its own residual vanishes",
                    "projected_pivot_coefficients": {k: _text(v) for k, v in pivot.items()},
                    "residual": residual_report(second_residual, context, max_examples),
                    "nonzero_beta_components_in_LLSS_LLHH": len(second_betas),
                    "nonzero_return_to_retained_support": any(
                        reduced(second_betas.get(k, 0)) != 0
                        for k in set(templates["C12"]) | set(templates["C5"])
                    ),
                }
            else:
                entry["second_step_return"] = {
                    "description": "First-step residual zero; no outside-sector feedback at second step",
                    "nonzero_return_to_retained_support": False,
                }
        results[name] = entry

    closed = all(x["closed_in_scanned_sector"] for x in results.values())
    status = "TwoTemplateClosedInScannedSector" if closed else "AdditionalTensorComponentsGenerated"
    return {
        "status": status,
        "scope": scope,
        "scope_complete_one_generation_L_eC_psi2phi2": scope == "llec",
        "full_three_generation_closure_verified": False,
        "other_EFT_operator_classes_checked": False,
        "one_loop_tensor_subspace_closed": closed if scope == "llec" else None,
        "scan_component_count": len(components),
        "templates": {name: len(canonical_values(t)) for name,t in templates.items()},
        "feedback_requested": feedback,
        "results": results,
        "interpretation": (
            "The retained two-template span fails under the tested one-loop tensor map; "
            "a two-operator-only evolution needs an extended basis or justified truncation."
            if not closed else
            "No extra components in the scanned psi2phi2 tensor sector at one loop. "
            "This is NOT proof of complete EFT or 3-generation closure."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, required=True,
                        help="Directory containing eft1_after_F_*_seed.json and eft1_rgbeta_rge.json")
    parser.add_argument("--scope", choices=("ll", "llec"), default="llec")
    parser.add_argument("--feedback", action="store_true",
                        help="Also apply the master beta to generated residuals in LLSS/LLHH (slower)")
    parser.add_argument("--output", type=Path, default=Path("output/t3b_operator_closure.json"))
    parser.add_argument("--max-examples", type=int, default=16)
    parser.add_argument("--report-every", type=int, default=30)
    args = parser.parse_args()
    if args.max_examples < 0 or args.report_every < 0:
        parser.error("--max-examples and --report-every must be nonnegative")
    d = args.data_dir
    from RGE.running.intermediate.ScalarOnlyWilsonTensorRGE import build_scalar_only_rge_context
    context = build_scalar_only_rge_context(
        d / "eft1_after_F_wilson_seed.json",
        d / "eft1_after_F_scalar_quartic_seed.json",
        d / "eft1_rgbeta_rge.json",
    )
    result = run_closure(
        context, scope=args.scope, feedback=args.feedback,
        max_examples=args.max_examples, report_every=args.report_every,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    summary = {
        "status": result["status"],
        "scope": result["scope"],
        "scanned_components": result["scan_component_count"],
        "operator_closure_verified": result["one_loop_tensor_subspace_closed"],
        "feedback_requested": args.feedback,
        "sources": {
            name: {"additional_components": entry["residual"]["component_count"],
                   "classes": entry["residual"]["component_counts_by_class"],
                   "second_step_return": entry.get("second_step_return", {}).get(
                       "nonzero_return_to_retained_support")}
            for name, entry in result["results"].items()
        },
        "full_output": str(args.output),
    }
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
