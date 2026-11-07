"""Batch 5: native scalar-only hard one-loop matching for six ordinary T3 models.

This is the NEXT stage after the saved F-only Matchete threshold, not another
fermion matching. It does not add the inherited one-loop EFT1 operator or
EFT1 direct Weinberg running; it therefore DOES NOT output complete C5.

Examples:
  python -X utf8 -m scripts.RunT3GeneralOrdinaryScalarHardBatch5 --preflight-only
  python -X utf8 -m scripts.RunT3GeneralOrdinaryScalarHardBatch5 --run-id ordinary6_hard_01
  python -X utf8 -m scripts.RunT3GeneralOrdinaryScalarHardBatch5 --run-id ordinary6_hard_01 --resume
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import traceback
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
SOURCE_ROOT = ROOT / "output/general_observable/fermion_first_all16_v1"
OUTPUT_ROOT = ROOT / "output/general_sequential_pilot/ordinary_scalar_hard_batch5"
WL = ROOT / "Lagrangian/RunT3GeneralOrdinaryScalarHardBatch5.wl"
SUPPORTED = ("B_-3", "B_+1", "C_-3", "C_+1", "E_-4", "E_+2")
EXPECTED_WILSON_TERMS = {"B": 2, "C": 4, "E": 4}
EXPECTED_EXTRA_QUARTICS = {"B": 1, "C": 1, "E": 0}
SOURCES = (
    "data/eft1_after_F_tree_inputform.txt",
    "data/eft1_after_F_wilson_seed.json",
    "data/eft1_after_F_scalar_quartic_seed.json",
    "data/eft1_after_F_rge_export.json",
    "native_status.json",
    "status.json",
)
SUCCESS = "HardScalarCandidateExported"


def alpha_argument(model: str) -> str:
    """Encode signed alpha without a leading option-like '-' for wolframscript."""
    alpha = int(model.split("_", 1)[1])
    return ("m" if alpha < 0 else "p") + str(abs(alpha))



class PreflightError(ValueError):
    pass


def slug(model: str) -> str:
    return model.replace("+", "p").replace("-", "m")


def read(path: Path) -> dict[str, Any]:
    data = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(data, dict):
        raise PreflightError(f"Expected JSON object: {path}")
    return data


def write(path: Path, data: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    tmp.replace(path)


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def source_fingerprints(model: str, source_root: Path = SOURCE_ROOT,
                        wl: Path = WL) -> dict[str, str]:
    base = source_root / slug(model)
    names = {str(n): base / n for n in SOURCES}
    names["matching_wl"] = wl
    missing = [f"{k}: {v}" for k, v in names.items() if not v.is_file()]
    if missing:
        raise PreflightError("Missing first-stage source inputs:\n" + "\n".join(missing))
    return {k: sha(v) for k, v in names.items()}


def preflight(model: str, source_root: Path = SOURCE_ROOT,
              wl: Path = WL) -> dict[str, Any]:
    if model not in SUPPORTED:
        raise PreflightError(f"{model}: ordinary hard source-reuse pilot supports only {SUPPORTED}")
    hashes = source_fingerprints(model, source_root, wl)
    folder = source_root / slug(model)
    native = read(folder / "native_status.json")
    status = read(folder / "status.json")
    wilson = read(folder / "data/eft1_after_F_wilson_seed.json")
    quartic = read(folder / "data/eft1_after_F_scalar_quartic_seed.json")
    metadata = read(folder / "data/eft1_after_F_rge_export.json")
    if status.get("status") != "FermionThresholdTreeWilsonExported" or \
       native.get("status") != "FThresholdTreeWilsonExported" or \
       native.get("integrated_fields") != ["F"] or \
       native.get("light_bsm_fields") != ["S1", "S2"] or \
       native.get("heavy_flavour_count") != 3:
        raise PreflightError("Saved F-only matching provenance is not verified")
    if native.get("loop_order") != 0:
        raise PreflightError("F-only reference must use tree matching")
    if not native.get("operator_basis_is_full_general_requested"):
        raise PreflightError("F-only source did not request full general operators")
    family = model[0]
    if (native.get("full_general_extra_yukawa_count") != 0 or
        native.get("full_general_quadratic_count") != 0 or
        native.get("full_general_quartic_count") != EXPECTED_EXTRA_QUARTICS[family]):
        raise PreflightError("Full-general additional operator inventory differs from verified model class")
    if any(doc.get("IntegratedFields") != ["F"] or
           doc.get("ActiveHeavyFields") != ["S1", "S2"]
           for doc in (wilson, quartic, metadata)):
        raise PreflightError("Tree, seeds, and EFT1 exporter disagree on active fields")
    expected_terms = EXPECTED_WILSON_TERMS[family]
    if (wilson.get("TreeWilsonTermCount") != expected_terms or
        len(wilson.get("TreeWilsonTerms", [])) != expected_terms or
        not {"PL", "PR"} <= {term.get("Chirality")
                             for term in wilson.get("TreeWilsonTerms", [])}):
        raise PreflightError(f"Expected exactly {expected_terms} native LLSS terms with both chiralities")
    if quartic.get("ScalarQuarticTermCount", 0) <= 0:
        raise PreflightError("No native scalar quartics found")
    tree = (folder / "data/eft1_after_F_tree_inputform.txt").read_text(
        encoding="utf-8-sig")
    if "Field[NewFermion," in tree or \
       not all(f"Field[{k}," in tree for k in ("NewScalar1", "NewScalar2")):
        raise PreflightError("First-threshold tree does not preserve exactly both scalar species")
    if re.search(r"Coupling\[(MS1|MS2),\s*\{\},\s*1\]", tree) is None:
        raise PreflightError("First-threshold light scalar mass counting is absent")
    diagonal_masses = re.findall(
        r"Coupling\[MF,\s*\{Index\[[^]\n]+,\s*NFlavor\]\},\s*0\]", tree
    )
    if len(diagonal_masses) != expected_terms:
        raise PreflightError(
            f"Saved LLSS tree has {len(diagonal_masses)} indexed MF factors, "
            f"expected {expected_terms}; unsafe to derive one-eigenvalue hard template")
    model_name = native.get("model")
    if model_name != model or status.get("model") != model:
        raise PreflightError("Source model mismatch")
    return {"model": model, "status": "Ready", "source_fingerprints": hashes,
            "tree_bytes": len(tree.encode("utf-8")),
            "tree_wilson_terms": wilson["TreeWilsonTermCount"],
            "additional_quartics_retained": EXPECTED_EXTRA_QUARTICS[family],
            "scalar_quartic_terms": quartic["ScalarQuarticTermCount"],
            "candidate_backend": "independent_native_hard_scalar_per_mass_eigenvalue_template",
            "indexed_fermion_mass_factors": len(diagonal_masses),
            "three_heavy_mass_numerical_sum_required": True}


def select_models(csv: str) -> tuple[str, ...]:
    models = tuple(x.strip() for x in csv.split(",") if x.strip())
    if not models or len(set(models)) != len(models) or any(x not in SUPPORTED for x in models):
        raise PreflightError(f"Choose unique comma-separated names from {SUPPORTED}")
    return models


def run(args: argparse.Namespace) -> int:
    models = select_models(args.models)
    preflights = []
    for model in models:
        try:
            item = preflight(model, args.source_root, args.wolfram_script)
            print(f"{model}: Ready; F tree reused, no fermion rematching", flush=True)
        except (OSError, ValueError, KeyError, TypeError) as exc:
            item = {"model": model, "status": "Blocked", "error": str(exc)}
            print(f"{model}: BLOCKED: {exc}", flush=True)
        preflights.append(item)
    if args.preflight_only or any(x["status"] != "Ready" for x in preflights):
        return 0 if all(x["status"] == "Ready" for x in preflights) else 2
    if not re.fullmatch(r"[A-Za-z0-9_-]+", args.run_id):
        raise PreflightError("--run-id must be a simple directory name")
    out = OUTPUT_ROOT / args.run_id
    settings = {"models": models, "source_root": str(args.source_root.resolve()),
                "wolfram_script": str(args.wolfram_script.resolve()),
                "timeout_seconds": args.timeout_seconds}
    cfg = out / "settings.json"
    if out.exists():
        if not args.resume or not cfg.is_file() or read(cfg) != dict(settings, models=list(models)):
            raise PreflightError("Output already exists or resume settings do not match")
    else:
        if args.resume:
            raise PreflightError("--resume requested but output directory does not exist")
        out.mkdir(parents=True)
        write(cfg, dict(settings, models=list(models)))
    successes = 0
    for model, info in zip(models, preflights):
        target = out / slug(model)
        report_path = target / "hard_scalar_report.json"
        prior = read(report_path) if args.resume and report_path.exists() else {}
        previous_meta = target / "provenance.json"
        if prior.get("Status") == SUCCESS and prior.get("PerHeavyMassEigenvalueTemplate") is True and previous_meta.is_file() and \
           read(previous_meta).get("source_fingerprints") == info["source_fingerprints"] and \
           (target / "hard_hbar_sector.wxf").is_file() and \
           (target / "hard_c5_candidate.txt").is_file():
            successes += 1
            print(f"{model}: reused completed hard stage (no Wolfram rerun)")
            continue
        if args.resume and report_path.exists() and prior.get("Status") == SUCCESS:
            raise PreflightError(f"{model}: completed source hashes or saved outputs changed; refusing overwrite")
        target.mkdir(parents=True, exist_ok=True)
        write(previous_meta, dict(info, generated_utc=datetime.now(timezone.utc).isoformat(),
                                  complete_C5=False))
        cmd = [args.wolframscript, "-file", str(args.wolfram_script),
               model.split("_")[0], alpha_argument(model),
               str(args.source_root / slug(model) / "data"), str(target),
               str(args.timeout_seconds)]
        env = os.environ.copy()
        env["T3_EXPERIMENTAL_GENERAL_QUADRATICS"] = "1"
        env["T3_EXPERIMENTAL_GENERAL_FERMIONS"] = "1"
        env["T3_EXPERIMENTAL_GENERAL_QUARTICS"] = "1"
        try:
            with (target / "wolfram_stdout.log").open("w", encoding="utf-8") as stdout, \
                 (target / "wolfram_stderr.log").open("w", encoding="utf-8") as stderr:
                proc = subprocess.run(cmd, stdout=stdout, stderr=stderr, env=env,
                                      cwd=ROOT, timeout=args.timeout_seconds * 6,
                                      check=False)
            result = read(report_path) if report_path.is_file() else {}
            if proc.returncode != 0 or result.get("Status") != SUCCESS or \
               result.get("PerHeavyMassEigenvalueTemplate") is not True:
                raise RuntimeError(f"exit={proc.returncode}, report_status={result.get('Status')}; "
                                   f"inspect {target / 'wolfram_stdout.log'}")
            if not (target / "hard_hbar_sector.wxf").is_file() or \
               not (target / "hard_c5_candidate.txt").is_file():
                raise RuntimeError("Missing postprocessed hard EFT outputs")
            successes += 1
            print(f"{model}: native one-eigenvalue hard scalar template exported")
        except (OSError, subprocess.TimeoutExpired, RuntimeError) as exc:
            write(target / "error.json", {"model": model, "status": "Failed",
                    "error": str(exc), "traceback": traceback.format_exc(limit=5)})
            print(f"{model}: hard scalar stage BLOCKED: {exc}", flush=True)
    write(out / "summary.json", {
        "status": "HardScalarCandidatesExported" if successes == len(models) else "PartialOrBlocked",
        "requested_models": list(models), "hard_candidates_exported": successes,
        "full_fixed_order_C5_assembled": False,
        "nf1_exact_pole_cancellation_derived": False,
        "independently_verified_native_normalization": False,
        "fermion_matching_recomputed": False,
        "hard_scalar_per_eigenvalue_template_only": True,
        "three_heavy_mass_sum_assembled": False,
        "note": "Independent one-eigenvalue scalar hard template; 3-mass flavour lift, inherited EFT1 and direct running still separate",
    })
    return 0 if successes == len(models) else 2


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--models", default=",".join(SUPPORTED))
    p.add_argument("--run-id", default="ordinary6_hard_01")
    p.add_argument("--source-root", type=Path, default=SOURCE_ROOT)
    p.add_argument("--wolfram-script", type=Path, default=WL)
    p.add_argument("--wolframscript", default="wolframscript")
    p.add_argument("--timeout-seconds", type=int, default=900)
    p.add_argument("--preflight-only", action="store_true")
    p.add_argument("--resume", action="store_true")
    a = p.parse_args(argv)
    if not 30 <= a.timeout_seconds <= 7200:
        p.error("timeout must be in 30..7200 seconds")
    try:
        return run(a)
    except (PreflightError, OSError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
