"""Derive the representation-dependent C5 -> C12 mixing coefficient.

Works in the repository's one-generation, separately symmetric real-tensor
convention. A pure Weinberg seed C5=kappa is inserted in the actual one-loop
Wilson tensor master beta. The term proportional to conjugate(lambdaT3) is
compared with the Matchete-CG mixed C12 operator tensor normalized to
C12 = conjugate(y1)*conjugate(y2)/MF.

A group factor is marked Verified only if the full generated LLSS tensor in
this coupling channel is proportional to that matched tensor, including
components where the matched C12 seed vanishes. No extrapolation between
representations and no assumption of full three-flavour operator closure.
"""
from __future__ import annotations

import argparse
from dataclasses import replace
import hashlib
import json
from itertools import combinations_with_replacement
from pathlib import Path

import sympy as sp

SEED = "eft1_after_F_wilson_seed.json"
QUARTIC = "eft1_after_F_scalar_quartic_seed.json"
RGBETA = "eft1_rgbeta_rge.json"
OUTPUT = "c12_reverse_c5_group_factor.json"
SCHEMA = "t3_c12_reverse_mixing_projection_v1"


def _parse(value: object) -> sp.Expr:
    return sp.sympify(str(value), locals={"I": sp.I, "sqrt": sp.sqrt, "conjugate": sp.conjugate})


def _canonical(key: tuple[int, int, int, int]) -> tuple[int, int, int, int]:
    i, j, a, b = map(int, key)
    return (min(i, j), max(i, j), min(a, b), max(a, b))


def _normalised_c12_seed(entries: dict) -> dict[tuple[int, int, int, int], sp.Expr]:
    """Coefficient tensor of conjugate(y1)conjugate(y2)/MF, with CG factors."""
    y1, y2, mf = sp.symbols("y1 y2 MF")
    u, v = sp.conjugate(y1), sp.conjugate(y2)
    out = {}
    for raw_key, value in entries.items():
        key = _canonical(raw_key)
        coeff = sp.simplify(mf * sp.expand(_parse(value)).coeff(u).coeff(v))
        if coeff == 0:
            continue
        if key in out and sp.simplify(coeff - out[key]) != 0:
            raise ValueError(f"C12 tensor symmetry discrepancy at {key}")
        out[key] = coeff
    if not out:
        raise ValueError("No mixed conjugate(y1)*conjugate(y2) C12 operator in matching seed")
    if any(val.free_symbols & {y1, y2, mf} for val in out.values()):
        raise ValueError("The C12 tensor still depends on the matching Yukawa/mass seed")
    return out


def project_reverse_mixing(
    c12_seed: dict[tuple[int, int, int, int], sp.Expr],
    beta_linear_in_lamstar_over_kappa: dict[tuple[int, int, int, int], sp.Expr],
    *,
    all_llss_components: list[tuple[int, int, int, int]],
) -> dict:
    """Exact full-channel proportionality check, not just pivot division."""
    seed = {_canonical(k): sp.simplify(v) for k, v in c12_seed.items() if v != 0}
    if not seed:
        raise ValueError("No C12 template")
    generated = {_canonical(k): sp.simplify(v) for k, v in beta_linear_in_lamstar_over_kappa.items()
                 if sp.simplify(v) != 0}
    keys = {_canonical(k) for k in all_llss_components} | set(seed) | set(generated)
    pivot = min(seed)
    ratio = sp.factor(sp.simplify(generated.get(pivot, sp.S.Zero) / seed[pivot]))
    residuals = []
    for key in sorted(keys):
        difference = sp.simplify(generated.get(key, sp.S.Zero) - ratio * seed.get(key, sp.S.Zero))
        if difference != 0:
            residuals.append({"component": list(key), "residual": sp.sstr(sp.factor(difference))})
    verified = not residuals and bool(generated)
    return {
        "status": "Verified" if verified else "NotProportionalOrAbsent",
        "factor": sp.sstr(ratio) if verified else None,
        "candidate_factor": sp.sstr(ratio),
        "coupling": "conjugate(lambdaT3)",
        "reference_component": list(pivot),
        "template_component_count": len(seed),
        "generated_component_count": len(generated),
        "scanned_llss_component_count": len(keys),
        "residual_count": len(residuals),
        "residuals": residuals[:12],
        "scope": "one-generation, lambdaT3-star channel only; not full 3-generation closure",
    }


def _sha256(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def derive_reverse_mixing(data_dir: Path) -> dict:
    """Re-evaluate a *pure C5 seed* using the source master tensor beta."""
    d = Path(data_dir)
    input_paths = {name: d / name for name in (SEED, QUARTIC, RGBETA)}
    for path in input_paths.values():
        if not path.is_file():
            raise FileNotFoundError(path)

    from RGE.running.intermediate.ScalarOnlyWilsonTensorRGE import build_scalar_only_rge_context
    from RGE.matching.WeinbergTensorAdapter import SparseWilsonLookup, build_weinberg_wilson_tensor
    from RGE.general.WilsonTensorRGE import calculate_complete_wilson_tensor_rge

    context = build_scalar_only_rge_context(input_paths[SEED], input_paths[QUARTIC], input_paths[RGBETA])
    rep_meta = json.loads(input_paths[RGBETA].read_text(encoding="utf-8"))["metadata"]
    if rep_meta.get("SharedScalar", False):
        raise NotImplementedError("Shared-scalar reverse C5 to C12 operator requires a separate template")
    rep = {"dS1": int(rep_meta["dS1"]), "dS2": int(rep_meta["dS2"]),
           "dF": int(rep_meta["dF"]), "alpha": int(rep_meta["alpha"]),
           "shared_scalar": False}
    c12 = _normalised_c12_seed(context.coefficient.entries)

    kappa = sp.Symbol("kappa")
    lamstar = sp.conjugate(sp.Symbol("lambdaT3"))
    c5_entries = build_weinberg_wilson_tensor(context.model, context.fermion_basis, kappa)
    c5_context = replace(context, coefficient=SparseWilsonLookup(c5_entries))

    f = context.fermion_basis
    leptons = [f.global_index("L", 1), f.global_index("L", 2)]
    higgs = set(context.model.block("H").indices)
    bsm_scalars = [i for i in range(1, context.model.total_real_scalar_dimension + 1)
                   if i not in higgs]
    llss_keys = [(i, j, a, b)
                 for i, j in combinations_with_replacement(leptons, 2)
                 for a, b in combinations_with_replacement(bsm_scalars, 2)]
    generated = {}
    for key in llss_keys:
        tensor_beta = calculate_complete_wilson_tensor_rge(
            model=c5_context.model, inputs=c5_context.inputs,
            output_component=key, coefficient=c5_context.coefficient,
            simplify_each=False,
        )
        total = sp.expand(_parse(tensor_beta["total"]))
        lam_coefficient = sp.simplify(total.coeff(lamstar))
        if lam_coefficient == 0:
            continue
        reduced = sp.simplify(lam_coefficient / kappa)
        if reduced.has(kappa, sp.conjugate(kappa), lamstar):
            raise ValueError(f"Reverse C5 beta has unexpected kappa/lambda dependence at {key}")
        generated[key] = reduced

    result = project_reverse_mixing(c12, generated, all_llss_components=llss_keys)
    # The reference T3-B alpha=-1 coefficient fixes the ordered-tensor
    # normalization used in the numerical coupled benchmark.
    if (rep["dS1"], rep["dS2"], rep["dF"], rep["alpha"]) == (2, 2, 1, -1):
        expected = 8 * sp.sqrt(3) / 3
        candidate = _parse(result["factor"]) if result["factor"] is not None else None
        if candidate is None or sp.simplify(candidate - expected) != 0:
            result["reference_mismatch"] = {
                "expected": sp.sstr(expected), "found": result["candidate_factor"],
            }
            result["status"] = "ReferenceNormalizationMismatch"
            result["factor"] = None
    result.update({
        "schema": SCHEMA,
        "representation": rep,
        "source_files_sha256": {name: _sha256(path) for name, path in input_paths.items()},
        "source_paths": {name: str(path) for name, path in input_paths.items()},
    })
    return result


def validate_cached_projection(data_dir: Path, rep: dict) -> dict | None:
    """Return verified factor only for matching representation and current input files."""
    d = Path(data_dir)
    cache = d / OUTPUT
    if not cache.is_file():
        return None
    result = json.loads(cache.read_text(encoding="utf-8"))
    if result.get("schema") != SCHEMA or result.get("status") != "Verified":
        return None
    expected_rep = {k: rep.get(k) for k in ("dS1", "dS2", "dF", "alpha")}
    if any(int(result.get("representation", {}).get(k, -999)) != int(v) for k, v in expected_rep.items()):
        return None
    recorded = result.get("source_files_sha256", {})
    if any(not (d / name).is_file() or recorded.get(name) != _sha256(d / name)
           for name in (SEED, QUARTIC, RGBETA)):
        return None
    if result.get("coupling") != "conjugate(lambdaT3)" or int(result.get("residual_count", -1)) != 0:
        return None
    return result


def _run_one(data_dir: Path, output: Path | None) -> dict:
    result = derive_reverse_mixing(data_dir)
    dest = output if output is not None else data_dir / OUTPUT
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    return {"model": str(data_dir), "status": result["status"], "factor": result["factor"],
            "residual_count": result["residual_count"], "output": str(dest)}


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    group = p.add_mutually_exclusive_group(required=True)
    group.add_argument("--data-dir", type=Path, help="Single saved T3 model data directory")
    group.add_argument("--source", type=Path, help="Root containing saved models; scan data folders recursively")
    p.add_argument("--output", type=Path, help="Only with --data-dir: explicit JSON output")
    args = p.parse_args(argv)
    if args.output and not args.data_dir:
        p.error("--output requires --data-dir")
    if args.data_dir:
        directories = [args.data_dir]
    else:
        directories = [path.parent for path in sorted(args.source.rglob(SEED))]
        if not directories:
            raise FileNotFoundError(f"No {SEED} found under {args.source}")
    summaries = []
    failed = False
    for directory in directories:
        try:
            item = _run_one(directory, args.output)
            if item["status"] != "Verified":
                failed = True
        except Exception as exc:
            item = {"model": str(directory), "status": "Error", "error": str(exc)}
            failed = True
        summaries.append(item)
        print(json.dumps(item))
    print(json.dumps({"count": len(summaries), "verified": sum(s["status"] == "Verified" for s in summaries),
                      "failed_or_unverified": sum(s["status"] != "Verified" for s in summaries)}, indent=2))
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
