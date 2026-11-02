"""All-16 source-derived ordered Wilson real-component tensor stage.

This command is a physical-basis *precondition* for EFT1 gamma evolution,
NOT a claim of completed general Wilson running or scalar-threshold C5.

14 models: contract actual 3-heavy-flavour Matchete PL tree terms from the
previous full-general boundary and expand their CG/scalar factors into
ordered real components (Higgs real components 1..4 included in the index).
The six self-bilinear models retain distinct neutral real mass eigenfields.

B_-1/C_-1: validate the CURRENT mu_F mass-basis Matchete LLSS source,
without identifying historical or missing CG definitions as valid tensors.
Their native mixed-scalar tree exporter does not provide a CGRegistry in the
same format as RunMatching.wl, so source-to-real-tensor conversion is pending.

Optional numerical Wilson evolution requires an *external* full gamma input
with exactly matching ordered tensor+flavour coordinates. No anomalous
matrix is invented and no gamma source is shipped with this batch.

Usage:
 python -X utf8 -m scripts.RunT3GeneralOrderedWilsonTensor --all
 python -X utf8 -m scripts.RunT3GeneralOrderedWilsonTensor --model A_-4
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np

from Numerical.running.GeneralT3OrderedTensor import project_ordered_term
from Numerical.running.GeneralT3LinearWilsonFlow import evolve_explicit_gamma

ROOT = Path(__file__).resolve().parents[1]
MODELS = tuple(f"{family}_{alpha:+d}" for family, alphas in {
    "A": (-4, -2, 0), "B": (-3, -1, 1), "C": (-3, -1, 1),
    "D": (-2, 0, 2), "E": (-4, -2, 0, 2),
}.items() for alpha in alphas)
MIXED = frozenset(("B_-1", "C_-1"))
FORMAT = "t3-ordered-source-flavour-and-real-scalar-tensor-v1"
BOUNDARY = Path("output/general_observable/general_wilson_boundary_v1")
SYNC = Path("output/general_observable/mixed_f_synchronised_v1")
OUT = Path("output/general_observable/ordered_wilson_tensor_v1")


def tag(model: str) -> str:
    return model.replace("+", "p").replace("-", "m")


def absolute(path: str | Path) -> Path:
    p = Path(path)
    return p.resolve() if p.is_absolute() else (ROOT / p).resolve()


def load(path: Path) -> dict:
    obj = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(obj, dict):
        raise ValueError(f"Not an object: {path}")
    return obj


def save(path: Path, obj: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def flatten_for_gamma(operators: list[dict]) -> tuple[list[str], np.ndarray]:
    """Flatten ordered tensor components, without silent pair symmetrisation."""
    labels, data = [], []
    for operator in operators:
        source = operator["source_index"]
        for component in operator["components"]:
            i, j = component["lepton_su2"]
            a, b = component["scalar_real"]
            raw = component["ordered_flavour_coefficient_gev_inv"]
            matrix = np.asarray(raw["real"], dtype=float) + 1j * np.asarray(raw["imag"], dtype=float)
            if matrix.shape != (3, 3) or not np.isfinite(matrix).all():
                raise ValueError("Invalid source-derived 3x3 tensor component")
            for p in range(3):
                for q in range(3):
                    labels.append(f"term{source}:L{i}{j}:S{a},{b}:F{p+1}{q+1}")
                    data.append(matrix[p, q])
    if len(labels) != len(set(labels)):
        raise ValueError("Repeated ordered tensor coordinate")
    return labels, np.asarray(data, dtype=complex)


def nonmixed(model: str, boundary_root: Path, output_root: Path,
             gamma_root: Path | None, *, max_gamma_dim: int = 2500) -> dict:
    src = boundary_root / tag(model) / "wilson_boundary.json"
    b = load(src)
    if (b.get("format") != "t3-source-exact-F-tree-flavour-operator-handoff-v1" or
            b.get("model") != model or
            b.get("status") != "SourceFactoredOrderedFlavourMatricesAtF_NotCanonicalWilsonTensor" or
            b.get("numeric_flavour_contractions_computed") is not True or
            b.get("canonical_real_wilson_tensor_computed") is not False):
        raise ValueError(f"{model}: invalid numerical source-factorised Matchete boundary")
    source = b.get("source", {})
    for name, path_key in (("running_sha256", "running_json"),
                           ("matchete_operator_sha256", "matchete_operator_source"),
                           ("first_threshold_sha256", "first_threshold_certificate")):
        file = absolute(source[path_key])
        if not file.is_file() or sha(file) != source[name]:
            raise ValueError(f"{model}: stale or missing source: {path_key}")
    terms = b.get("operator_terms")
    cg = b.get("source_CG_registry")
    if not isinstance(terms, list) or not terms or not isinstance(cg, list):
        raise ValueError(f"{model}: missing exact Matchete operator terms/CG registry")
    operators = [project_ordered_term(t, cg, model) for t in terms]
    payload = {
        "format": FORMAT, "model": model,
        "status": "OrderedNumericalRealComponentTensorFromMatchete",
        "source_boundary": str(src), "source_boundary_sha256": sha(src),
        "mu_F_gev": b["mu_F_gev"],
        "scalar_basis_explicit": True, "matchete_cg_contractions_retained": True,
        "external_three_family_flavour_order_not_symmetrized": True,
        "tree_operator_count_PL": len(operators),
        "total_ordered_real_component_count": sum(x["component_count"] for x in operators),
        "operators": operators,
        "canonical_C_ijab_normalization_certified": False,
        "physical_matchete_rgbeta_coupling_dictionary_certified": False,
        "full_eft1_wilson_anomalous_dimension_implemented": False,
        "eft1_wilson_rge_computed": False,
        "scalar_threshold_c5_computed": False,
    }
    dest = output_root / tag(model) / "ordered_tensor.json"
    save(dest, payload)
    report = {"model": model, "status": payload["status"],
              "tree_operator_count_PL": len(operators),
              "real_component_count": payload["total_ordered_real_component_count"],
              "tensor_file": str(dest), "eft1_wilson_rge_computed": False}
    if gamma_root is not None:
        gamma_path = gamma_root / tag(model) / "gamma.json"
        if gamma_path.is_file():
            keys, initial = flatten_for_gamma(operators)
            if len(keys) > max_gamma_dim:
                raise ValueError(f"Gamma basis would have {len(keys)} coordinates; raise --max-gamma-dim explicitly")
            gamma = load(gamma_path)
            if gamma.get("matching_boundary_sha256") != sha(src):
                raise ValueError("Gamma source does not match the exact matched boundary")
            mu_s = float(gamma["mu_S_gev"])
            result = evolve_explicit_gamma(initial, keys, gamma, float(b["mu_F_gev"]), mu_s)
            save(output_root / tag(model) / "conditional_gamma_running.json", result)
            report.update(gamma_status=result["status"],
                          gamma_file=str(gamma_path),
                          conditional_wilson_running=True)
    return report


def mixed(model: str, sync_root: Path, output_root: Path) -> dict:
    """Validate current B/C tree source, fail closed on missing CG dictionary."""
    status_file = sync_root / tag(model) / "status.json"
    s = load(status_file)
    if (s.get("model") != model or
            s.get("status") != "CurrentRunningMassBasisFThresholdTreeLLSSExported" or
            s.get("source_operator_matching_at_running_mu_F") is not True or
            s.get("historical_mu_S_basis_reused") is not False):
        raise ValueError(f"{model}: no current-F synchronized Matchete tree certificate")
    src = absolute(s["matched_llss_file"])
    mass_file = absolute(s["mass_basis_file"])
    running_file = absolute(s["running_file"])
    for p in (src, mass_file, running_file):
        if not p.is_file():
            raise FileNotFoundError(f"Missing synchronized B/C source {p}")
    basis = load(mass_file)
    provenance = basis.get("running_boundary_provenance", {})
    if provenance.get("running_sha256") != sha(running_file):
        raise ValueError("B/C mixed mass eigenbasis no longer matches running.json")
    if float(provenance.get("boundary_scale_gev", -1)) != float(s["mu_F_gev"]):
        raise ValueError("B/C scalar rotation not evaluated at the F threshold")
    if "GammaCC" not in src.read_text(encoding="utf-8"):
        raise ValueError("B/C source lacks an LLSS Majorana operator")
    result = {
        "format": FORMAT, "model": model,
        "status": "SynchronizedMixedLLSSSource_CGAndFlavourTensorPending",
        "source_status_sha256": sha(status_file),
        "matched_llss_sha256": sha(src),
        "matched_llss_file": str(src),
        "mass_basis_file": str(mass_file),
        "mu_F_gev": s["mu_F_gev"],
        "scalar_fields": ["NewScalarX", "NewScalarY"],
        "scalar_masses_gev": s["scalar_eigenmasses_gev"],
        "raw_tree_llss_term_count": s["llss_term_count"],
        "current_running_basis_verified": True,
        "missing_source_item": "CGRegistry for X/Y native Matchete fields, plus exact source flavour-template factorisation",
        "canonical_C_ijab_normalization_certified": False,
        "eft1_wilson_rge_computed": False,
        "scalar_threshold_c5_computed": False,
    }
    save(output_root / tag(model) / "ordered_tensor_status.json", result)
    return {"model": model, "status": result["status"],
            "tree_operator_count_PL": None,
            "raw_llss_term_count": s["llss_term_count"],
            "eft1_wilson_rge_computed": False}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--all", action="store_true")
    group.add_argument("--model", choices=MODELS)
    group.add_argument("--models", help="Comma separated models")
    parser.add_argument("--boundary-root", type=Path, default=BOUNDARY)
    parser.add_argument("--synchronised-root", type=Path, default=SYNC)
    parser.add_argument("--output-root", type=Path, default=OUT)
    parser.add_argument("--gamma-root", type=Path,
                        help="Optional explicit exact-basis gamma matrix samples. No default gamma exists.")
    parser.add_argument("--max-gamma-dim", type=int, default=2500)
    args = parser.parse_args(argv)
    selection = MODELS if args.all else (args.model,) if args.model else tuple(
        x.strip() for x in args.models.split(",") if x.strip())
    if not selection or len(set(selection)) != len(selection) or not set(selection) <= set(MODELS):
        parser.error("Invalid/duplicate model selection")
    out = absolute(args.output_root)
    boundary = absolute(args.boundary_root)
    sync = absolute(args.synchronised_root)
    gamma = absolute(args.gamma_root) if args.gamma_root else None
    rows = []
    for i, model in enumerate(selection, 1):
        try:
            row = mixed(model, sync, out) if model in MIXED else nonmixed(
                model, boundary, out, gamma, max_gamma_dim=args.max_gamma_dim)
        except Exception as exc:
            row = {"model": model, "status": "Blocked", "error": f"{type(exc).__name__}: {exc}"}
            save(out / tag(model) / "error.json", row)
        rows.append(row)
        print(f"[{i}/{len(selection)}] {model}: {row['status']}"
              + ("; " + row["error"] if "error" in row else ""), flush=True)
        save(out / "summary.json", {
            "format": FORMAT, "results": rows,
            "ordered_real_tensor_count": sum(x["status"] == "OrderedNumericalRealComponentTensorFromMatchete" for x in rows),
            "synchronized_mixed_source_pending_count": sum(x["status"].startswith("SynchronizedMixedLLSS") for x in rows),
            "conditional_gamma_run_count": sum(x.get("conditional_wilson_running", False) for x in rows),
            "blocked_count": sum(x["status"] == "Blocked" for x in rows),
            "physical_sequential_matching_certified": False,
            "full_eft1_wilson_rge_implemented": False,
        })
    print("Ordered real tensors: " + str(sum(x["status"] == "OrderedNumericalRealComponentTensorFromMatchete" for x in rows))
          + "; mixed source pending CG registry: " + str(sum(x["status"].startswith("SynchronizedMixedLLSS") for x in rows))
          + "; physical EFT1 running: NOT CERTIFIED")
    return 1 if any(row["status"] == "Blocked" for row in rows) else 0


if __name__ == "__main__":
    raise SystemExit(main())
