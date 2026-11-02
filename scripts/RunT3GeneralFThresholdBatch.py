"""All-16 source-backed *first* threshold: integrate F, retain scalar EFT1.

For 14 models, invoke the existing universal RunSequentialT3Matching / tree
Wilson seed exporter through Lagrangian/RunT3GeneralFThreshold.wl. For mixed
B_-1/C_-1, use the repository's separately validated scalar eigen-doublet
Matchete tree matching. No inferred SU(2) Wilson tensors are used.

The output is a genuine symbolic tree-level F-threshold matching artifact,
NOT a numerically evolved EFT1 Wilson coefficient and NOT final physical C5.

  python -X utf8 -m scripts.RunT3GeneralFThresholdBatch --all --preflight
  python -X utf8 -m scripts.RunT3GeneralFThresholdBatch --all
  python -X utf8 -m scripts.RunT3GeneralFThresholdBatch --models A_-4,B_-1
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import traceback

ROOT = Path(__file__).resolve().parents[1]
NATIVE = ROOT / "Lagrangian" / "RunT3GeneralFThreshold.wl"
MATCHING = ROOT / "Lagrangian" / "RunMatching.wl"
DEFAULT_OUT = Path("output/general_observable/fermion_first_all16_v1")
ALPHAS = {"A": (-4, -2, 0), "B": (-3, -1, 1), "C": (-3, -1, 1),
          "D": (-2, 0, 2), "E": (-4, -2, 0, 2)}
MODELS = tuple(f"{c}_{a:+d}" for c, alphas in ALPHAS.items() for a in alphas)
MIXED = frozenset(("B_-1", "C_-1"))
SELF = frozenset(("A_-2", "A_+0", "D_-2", "D_+0", "E_-2", "E_+0"))
FORMAT = "t3-general-first-fermion-threshold-source-matched-v1"


def tag(model: str) -> str:
    if model not in MODELS:
        raise ValueError(f"Unknown T3 model {model!r}")
    return model.replace("+", "p").replace("-", "m")


def read_json(path: Path) -> dict:
    value = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(value, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return value


def save_json(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def sha(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def verify_native_tree(folder: Path, model: str, loop_order: int) -> dict:
    """Fail-closed on Matchete's exported *operator* and CG data."""
    certificate = folder / "native_status.json"
    report = read_json(certificate)
    expected = {
        "model": model, "status": "FThresholdTreeWilsonExported",
        "loop_order": loop_order, "integrated_fields": ["F"],
        "light_bsm_fields": ["S1", "S2"], "heavy_flavour_count": 3,
        "operator_basis_is_full_general_requested": True,
        "full_sequential_eft1_computed": False,
        "physical_c5_certified": False,
        "hard_one_loop_MSbar_certified": False,
        "eft1_wilson_anomalous_dimensions_implemented": False,
        "scalar_threshold_matching_implemented": False,
    }
    for k, value in expected.items():
        if report.get(k) != value or type(report.get(k)) is not type(value):
            raise ValueError(f"{model}: native certificate mismatch: {k}")
    files = list((folder / "data").glob("*wilson_seed.json"))
    if len(files) != 1 or files[0].name != report.get("wilson_seed_file"):
        raise ValueError(f"{model}: expected exactly one matching Wilson seed")
    seed = read_json(files[0])
    if seed.get("IntegratedFields") != ["F"] or seed.get("ActiveHeavyFields") != ["S1", "S2"]:
        raise ValueError(f"{model}: wrong threshold/heavy/light field inventory")
    terms = seed.get("TreeWilsonTerms")
    n = seed.get("TreeWilsonTermCount")
    if (not isinstance(n, int) or isinstance(n, bool) or n <= 0 or
            not isinstance(terms, list) or len(terms) != n or
            report.get("tree_wilson_term_count") != n):
        raise ValueError(f"{model}: missing or inconsistent tree Wilson operator terms")
    for term in terms:
        if not isinstance(term, dict) or "Field[NewFermion" in term.get("TermInputForm", ""):
            raise ValueError(f"{model}: invalid tree Wilson term containing integrated F")
        if not isinstance(term.get("TermInputForm"), str) or not term["TermInputForm"]:
            raise ValueError(f"{model}: missing exact Wilson operator expression")
    if model in SELF:
        basis = report.get("self_bilinear_mass_basis")
        if not isinstance(basis, dict):
            raise ValueError(f"{model}: missing self-bilinear mass-basis certificate")
        species = 2 if model.endswith("-2") else 1
        required = {
            "mass_basis_kind": "RealSelfBilinearSplitLightFields",
            "split_species": species,
            "scalar_mass_components_remain_light": True,
            "real_bilinear_mass_basis_candidate": True,
            "complex_self_bilinear_phase_supported": False,
            "holomorphic_bilinear_absorbed_into_split_free_masses": True,
            "field_basis_normalization_physically_certified": False,
        }
        for key, expected in required.items():
            if basis.get(key) != expected or type(basis.get(key)) is not type(expected):
                raise ValueError(f"{model}: invalid mass-basis certificate: {key}")
        if seed.get("SelfBilinearMassBasis") != basis or seed.get("WilsonSeedSource") != (
                "MatcheteFThresholdTreeReclassifiedForSplitLightScalarFields"):
            raise ValueError(f"{model}: split-field Wilson seed provenance mismatch")
        split_prefix = f"Field[T3Scalar{species}"
        old_prefix = f"Field[NewScalar{species}"
        if not any(split_prefix in item["TermInputForm"] for item in terms):
            raise ValueError(f"{model}: split-scalar fields absent from Wilson operators")
        if any(old_prefix in item["TermInputForm"] for item in terms):
            raise ValueError(f"{model}: unrotated original scalar found in Wilson operators")
        if not (folder / "self_mass_basis.json").is_file() or read_json(
                folder / "self_mass_basis.json") != basis:
            raise ValueError(f"{model}: missing independent mass-basis metadata")
    cg = seed.get("CGRegistry")
    if not isinstance(cg, (list, dict)):
        raise ValueError(f"{model}: missing CG tensor registry")
    return {
        "status": "FermionThresholdTreeWilsonExported",
        "tree_wilson_term_count": n,
        "wilson_seed": str(files[0]),
        "cg_registry_entries": len(cg),
        "raw_one_loop_F_matching_included": loop_order == 1,
        "mass_basis": "real-self-bilinear-split-light-scalar" if model in SELF
                      else "original-scalar-basis",
        "complex_self_bilinear_phase_supported": False if model in SELF else None,
        "fermion_first_matching_scheme_certified": False,
        "eft1_wilson_running_performed": False,
        "scalar_threshold_matching_performed": False,
    }


def _run_general(model: str, folder: Path, args: argparse.Namespace) -> dict:
    folder.mkdir(parents=True, exist_ok=True)
    cmd = [args.wolframscript, "-file", str(NATIVE), str(folder), model,
           str(args.loop_order)]
    env = dict(os.environ)
    for key in ("T3_EXPERIMENTAL_GENERAL_FERMIONS",
                "T3_EXPERIMENTAL_GENERAL_QUARTICS",
                "T3_EXPERIMENTAL_GENERAL_QUADRATICS",
                "T3_EXPERIMENTAL_FLAVOR_MASSES"):
        env[key] = "1"
    try:
        completed = subprocess.run(cmd, cwd=ROOT, env=env, capture_output=True,
                                   text=True, timeout=args.timeout_sec, check=False)
    except subprocess.TimeoutExpired as error:
        raise TimeoutError(f"{model}: Wolfram timed out after {args.timeout_sec}s") from error
    (folder / "native.log").write_text(
        "command: " + " ".join(cmd) + "\nexit: " + str(completed.returncode) +
        "\nstdout:\n" + completed.stdout + "\nstderr:\n" + completed.stderr,
        encoding="utf-8")
    if completed.returncode:
        native = folder / "native_status.json"
        detail = read_json(native) if native.is_file() else {}
        raise RuntimeError(f"{model}: Matchete F-threshold failed exit={completed.returncode}; "
                           f"{detail.get('error', '')}; see {folder / 'native.log'}; "
                           f"tail={(completed.stdout + completed.stderr)[-650:]}")
    return verify_native_tree(folder, model, args.loop_order)


def _run_mixed(model: str, folder: Path, args: argparse.Namespace) -> dict:
    # The only allowed source for this mass basis is the existing checked
    # experimental Matchete driver. It does NOT use the unrelated one-scale C5.
    from scripts import RunGeneralT3FermionFirstEFT1Tree as tree
    from scripts import RunGeneralT3FermionFirstEFT1OneLoop as one_loop
    folder.mkdir(parents=True, exist_ok=True)
    t = tree.run_model(model, which="final", outdir=folder / "tree",
                       wolframscript=args.wolframscript, timeout_sec=args.timeout_sec)
    if t.get("status") != tree.STATUS or t.get("physical_C5_computed") is not False:
        raise ValueError(f"{model}: wrong existing mass-eigenstate tree certificate")
    report = {
        "status": "FermionThresholdMixedMassTreeLLSSExported",
        "tree_wilson_term_count": t["llss_term_count"],
        "llss_operator_file": t["llss_file"],
        "matched_eft_file": t["matched_eft_file"],
        "mass_basis": "validated experimental mixed-doublet X,Y",
        "mass_basis_input_provenance": "existing experimental B/C mass-basis files; NOT automatically reconstructed from the current all16 numerical trajectory",
        "raw_one_loop_F_matching_included": False,
        "fermion_first_matching_scheme_certified": False,
        "eft1_wilson_running_performed": False,
        "scalar_threshold_matching_performed": False,
    }
    if args.loop_order == 1:
        x = one_loop.run_model(model, which="final", tree_dir=folder / "tree",
                               outdir=folder / "one_loop",
                               wolframscript=args.wolframscript,
                               timeout_sec=args.timeout_sec)
        if x.get("status") != one_loop.STATUS or x.get("physical_C5_computed") is not False:
            raise ValueError(f"{model}: wrong existing raw one-loop certificate")
        report["raw_one_loop_F_matching_included"] = True
        report["raw_one_loop_hard_candidate_terms"] = x["raw_hard_llhh_candidate_terms"]
    return report


def select_models(args: argparse.Namespace) -> tuple[str, ...]:
    if args.all:
        return MODELS
    if args.model:
        return (args.model,)
    chosen = [m.strip() for m in args.models.split(",")]
    if not chosen or any(m not in MODELS for m in chosen) or len(chosen) != len(set(chosen)):
        raise ValueError("--models requires distinct comma-separated identifiers from the 16-model catalogue")
    return tuple(chosen)


def run(args: argparse.Namespace) -> dict:
    selected = select_models(args)
    sources = [NATIVE, MATCHING,
               ROOT / "Lagrangian/model/LagrangianBuilder.wl",
               ROOT / "Lagrangian/model/T3ModelCatalog.wl",
               ROOT / "Lagrangian/model/T3Fields.wl",
               ROOT / "Lagrangian/model/SU2Invariants.wl",
               ROOT / "Lagrangian/interactions/ScalarPotential.wl",
               ROOT / "Lagrangian/interactions/T3Topology.wl"]
    missing = [str(path) for path in sources if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing required source file(s): " + ", ".join(missing))
    fingerprints = {str(s.relative_to(ROOT)): sha(s) for s in sources}
    root = args.output_root.resolve()
    if not args.preflight:
        root.mkdir(parents=True, exist_ok=True)
    entries = []
    for i, model in enumerate(selected, 1):
        dest = root / tag(model)
        if args.preflight:
            source_type = "existing_mixed_mass_eigenstate_matchete" if model in MIXED else "general_f_only_matchete"
            row = {"model": model, "status": "Preflight", "source_type": source_type,
                   "expected_output": str(dest)}
            if model in MIXED:
                for required in (
                    ROOT / "output/general_t3_mass_eigen_matchete" / (tag(model) + "_final_matchete.json"),
                    ROOT / "output/general_t3_mass_basis" / (tag(model) + "_final_mass_basis.json")):
                    if not required.is_file():
                        row.setdefault("missing_mixed_prerequisites", []).append(str(required))
                if row.get("missing_mixed_prerequisites"):
                    row["status"] = "PreflightMissingExistingMixedMassBasis"
        else:
            print(f"[{i}/{len(selected)}] {model}: F-first Matchete tree (loop={args.loop_order})", flush=True)
            try:
                if model not in MIXED and args.reuse and (dest / "native_status.json").is_file():
                    previous = read_json(dest / "status.json")
                    if previous.get("source_sha256") != fingerprints:
                        raise ValueError("Cannot reuse F-threshold export: native source hashes changed")
                    evidence = verify_native_tree(dest, model, args.loop_order)
                    evidence["cache"] = "ReusedSourceHashVerified"
                elif model in MIXED:
                    evidence = _run_mixed(model, dest, args)
                else:
                    evidence = _run_general(model, dest, args)
                row = {"model": model, **evidence, "source_sha256": fingerprints}
                print(f"  {row['status']}; LLSS/Wilson terms={row['tree_wilson_term_count']}", flush=True)
            except (Exception,) as error:
                row = {"model": model, "status": "Blocked", "error": f"{type(error).__name__}: {error}",
                       "traceback": traceback.format_exc(limit=7)}
                print("  BLOCKED: " + row["error"], file=sys.stderr, flush=True)
            save_json(dest / "status.json", row)
        entries.append(row)
        if not args.preflight:
            save_json(root / "summary.json", {
                "format": FORMAT, "status": "CompletedWithExplicitPerModelFailures",
                "loop_order": args.loop_order,
                "source_sha256": fingerprints,
                "model_count": len(selected),
                "tree_success_count": sum(x["status"].startswith("FermionThreshold") for x in entries),
                "full_sequential_wilson_matching_performed": False,
                "eft1_dim5_anomalous_dimension_running_performed": False,
                "scalar_threshold_renormalized_matching_performed": False,
                "matchete_rgbeta_physical_dictionary_certified": False,
                "physical_c5_prediction_certified": False,
                "models": entries,
            })
    return {"models": entries, "output_root": str(root),
            "successes": sum(x["status"].startswith("FermionThreshold") for x in entries),
            "count": len(selected)}


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    choices = p.add_mutually_exclusive_group(required=True)
    choices.add_argument("--all", action="store_true")
    choices.add_argument("--model", choices=MODELS)
    choices.add_argument("--models", type=str)
    p.add_argument("--loop-order", type=int, default=0, choices=(0, 1),
                   help="0: exact tree Wilson; 1: additionally export RAW one-loop hard correction")
    p.add_argument("--preflight", action="store_true", help="Check sources without Matchete calculation")
    p.add_argument("--reuse", action="store_true", help="Reuse previously verified native exports only when source hashes match")
    p.add_argument("--wolframscript", default="wolframscript")
    p.add_argument("--timeout-sec", type=int, default=7200)
    p.add_argument("--output-root", type=Path, default=DEFAULT_OUT)
    return p


def main(argv: list[str] | None = None) -> int:
    args = parser().parse_args(argv)
    if args.timeout_sec <= 0:
        raise SystemExit("--timeout-sec must be positive")
    if not args.preflight and not shutil.which(args.wolframscript) and not Path(args.wolframscript).is_file():
        raise SystemExit("wolframscript not found")
    report = run(args)
    print(f"F-first tree: {report['successes']}/{report['count']}; "
          "EFT1 Wilson RGEs / scalar matching: NOT performed", flush=True)
    return 0 if args.preflight or report["successes"] == report["count"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
