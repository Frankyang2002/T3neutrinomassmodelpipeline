"""Audit exact UV T3 CG component support; do not infer loop multiplicity.

Uses Matchete-exported T3Y1CG, T3Y2CG, T3MixCG. This is a prerequisite
for computing charge-compatible loop contractions, not itself an integer
multiplicity or a prescription for rescaling lambdaT3.
"""
from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

import sympy as sp
from RGE.running.intermediate.IntermediateMatcheteParsing import load_cg_registry

NAMES = ("T3Y1CG", "T3Y2CG", "T3MixCG")


def expand_tensor(tensor: sp.MutableDenseNDimArray, dimensions: tuple[int, ...]):
    """Restore singleton axes omitted from Matchete's sparse CG representation."""
    kept = tuple(i for i, d in enumerate(dimensions) if d > 1)
    shape = tuple(int(v) for v in tensor.shape)
    expected = tuple(dimensions[i] for i in kept)
    if shape != expected:
        raise ValueError(f"CG axes {shape} do not agree with expected {expected}")
    for indices in itertools.product(*(range(1, d + 1) for d in dimensions)):
        active = tuple(indices[i] - 1 for i in kept)
        weight = sp.simplify(tensor[active])
        if weight != 0:
            yield indices, weight


def audit_registry(path: Path) -> dict:
    payload = json.loads(Path(path).read_text(encoding="utf-8-sig"))
    if payload.get("status") != "Success":
        raise ValueError(f"Incomplete CG export: {path}")
    meta = payload["metadata"]
    d1, d2, df = (int(meta[k]) for k in ("dS1", "dS2", "dF"))
    dimensions = {
        "T3Y1CG": (2, df, d1),
        "T3Y2CG": (2, df, d2),
        "T3MixCG": (2, 2, d1, d2),
    }
    registry = load_cg_registry(payload)
    if any(k not in registry for k in NAMES):
        raise ValueError(f"Missing CG tensors: {sorted(set(NAMES) - set(registry))}")
    vertices = {}
    for name in NAMES:
        dims = dimensions[name]
        rows = [
            {"indices": list(indices), "cg": str(weight),
             "cg_abs_squared": str(sp.simplify(weight * sp.conjugate(weight)))}
            for indices, weight in expand_tensor(registry[name].tensor, dims)
        ]
        vertices[name] = {
            "dimensions": list(dims),
            "nonzero_component_count": len(rows),
            "sum_cg_abs_squared": str(sp.simplify(sum(
                (sp.sympify(row["cg_abs_squared"]) for row in rows), sp.S.Zero
            ))),
            "components": rows,
        }
    return {
        "source": str(path),
        "model": {k: meta[k] for k in ("dS1", "dS2", "dF", "alpha")},
        "vertices": vertices,
        "loop_multiplicity": None,
        "compensation_factor": None,
        "status": "component_support_only",
        "warning": ("Nonzero CG entries include charged states, different weights, "
                    "and field conjugation. Their count is NOT the count of "
                    "equivalent neutrino-mass loop contributions. The full "
                    "charge-selected topology contraction is still required."),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path("output/full/hypercharge"))
    parser.add_argument("--output", type=Path, default=Path("Reports/output/full/multiplicity/component_support_audit.json"))
    args = parser.parse_args()
    paths = sorted(args.root.rglob("uv_cg_registry.json"))
    if not paths:
        print("No uv_cg_registry.json found. Export the UV CG tensors first; no counts were invented.")
        return 2
    cases, errors = [], []
    for path in paths:
        try:
            cases.append(audit_registry(path))
        except (ValueError, KeyError, TypeError, IndexError) as exc:
            errors.append({"path": str(path), "error": str(exc)})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps({"cases": cases, "errors": errors}, indent=2), encoding="utf-8")
    print(f"Component CG audits: {len(cases)}; failures: {len(errors)}")
    print(f"Report: {args.output}")
    print("No multiplicity correction has been assigned or applied.")
    return 1 if errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
