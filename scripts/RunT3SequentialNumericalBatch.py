"""Resumable all-16 (default three) conditional LLSS sequential numerical campaign.

IMPORTANT: hard one-loop matching from the verified D(-2) Wolfram pipeline is
NOT generalized here. Results are LLSS-induced conditional C5 contributions,
not physically complete C5. Does not invoke or overwrite D(-2) matching.

Examples from repository root:
 python -X utf8 -m scripts.RunT3SequentialNumericalBatch --check-only
 python -X utf8 -m scripts.RunT3SequentialNumericalBatch --run-id first3_v1
 python -X utf8 -m scripts.RunT3SequentialNumericalBatch --run-id first3_v1 --resume
 python -X utf8 -m scripts.RunT3SequentialNumericalBatch --run-id all16_v1 --all16
 python -X utf8 -m scripts.RunT3SequentialNumericalBatch --backend archived-v6 --run-id archive3_v1
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import sys
import traceback
from typing import Any

from Numerical.running.T3SequentialBatchNumerics import (
    ARCHIVED_REL, DEFAULT_MODELS, FORMAT, MODELS, PhysicsInputError,
    dminus2_existing_reference, fingerprint, from_archived_v6,
    from_dynamic_case, mass_class, model_tag, normalize_model, read_json,
)

ROOT = Path(__file__).resolve().parents[1]
OUTPUT_ROOT = ROOT / "output/general_sequential_pilot/all16_numerical_batch"
DEFAULT_RELEASE = ROOT / "output/general_observable/lo_fermion_boundary/releases/source_021dcc1e8fe8ba0f"
DEFAULT_CHANNELS = ROOT / "output/general_observable/lo_joint_scalar_matching/source_quartic_v4/native_quartic_channels"
DEFAULT_SEEDS = ROOT / "output/full16_symbolic"
DEFAULT_RUNNING = ROOT / "output/general_observable/all16_v1"
DEFAULT_REGISTRY = ROOT / "RGE/general/generated/general_rgbeta_quartic_registry.json"
DEFAULT_D2 = ROOT / "output/general_sequential_pilot/D_m2_stage1/data/general_d_m2_scalar_threshold/sequential_run/sequential_matching_manifest.json"


def write_json(path: Path, obj: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(json.dumps(obj, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
                    encoding="utf-8")
    temp.replace(path)


def select_models(args: argparse.Namespace) -> tuple[str, ...]:
    raw = MODELS if args.all16 else (
        tuple(part.strip() for part in args.models.split(",") if part.strip())
        if args.models else DEFAULT_MODELS)
    selected = tuple(normalize_model(x) for x in raw)
    if not selected or len(set(selected)) != len(selected):
        raise PhysicsInputError("Model selection must be nonempty and unique")
    return selected


def settings(args: argparse.Namespace, models: tuple[str, ...]) -> dict[str, Any]:
    return {
        "format": FORMAT,
        "models": list(models), "backend": args.backend,
        "spinor_factor": args.spinor_factor, "vev_gev": args.vev,
        "with_smeft": args.with_smeft,
        "release_dir": str(args.release_dir.resolve()),
        "v4_dir": str(args.v4_channel_dir.resolve()),
        "seeds_dir": str(args.quartic_seed_dir.resolve()),
        "running_dir": str(args.running_dir.resolve()),
        "registry": str(args.registry.resolve()),
        "archive_dir": str(args.archived_dir.resolve()),
    }


def dynamic_paths(args: argparse.Namespace) -> dict[str, Path]:
    return {"release": args.release_dir, "v4": args.v4_channel_dir,
            "seeds": args.quartic_seed_dir, "running": args.running_dir,
            "registry": args.registry}


def input_paths(model: str, args: argparse.Namespace) -> dict[str, Path]:
    if args.backend == "archived-v6":
        return {"archived_c5": args.archived_dir / f"{model_tag(model)}.json"}
    from scripts.RunT3GeneralSequentialLeadingLog import source_paths
    return source_paths(dynamic_paths(args), model)


def input_hashes(model: str, args: argparse.Namespace) -> dict[str, str]:
    paths = input_paths(model, args)
    missing = {k: str(v) for k, v in paths.items() if not v.is_file()}
    if args.backend == "dynamic":
        if not args.registry.is_file():
            missing["general_quartic_registry"] = str(args.registry)
        else:
            paths["general_quartic_registry"] = args.registry
    if missing:
        raise PhysicsInputError(f"{model}: missing source inputs: {missing}")
    return {key: fingerprint(path) for key, path in paths.items()}


def preflight(model: str, args: argparse.Namespace) -> dict[str, Any]:
    try:
        hashes = input_hashes(model, args)
        if args.backend == "dynamic":
            from Numerical.joint_scalar_matching.fermion_boundary_bridge import load_release_manifest
            release = load_release_manifest(args.release_dir)
            if set(release) != set(MODELS) or model not in release:
                raise PhysicsInputError("Full 16-model source-boundary manifest missing")
            reg = read_json(args.registry)
            if model not in reg.get("models", {}):
                raise PhysicsInputError(f"{model}: missing quartic registry row")
        else:
            # Model-specific schema/physics checks are repeated during computation.
            raw = read_json(args.archived_dir / f"{model_tag(model)}.json")
            if raw.get("model") != model:
                raise PhysicsInputError(f"{model}: archived model does not match")
        return {"model": model, "mass_class": mass_class(model),
                "status": "Ready", "sha256": hashes}
    except (OSError, PhysicsInputError, ValueError, ImportError) as exc:
        return {"model": model, "mass_class": mass_class(model),
                "status": "Blocked", "error": f"{type(exc).__name__}: {exc}"}


def evaluate_model(model: str, args: argparse.Namespace) -> dict[str, Any]:
    if args.backend == "archived-v6":
        path = args.archived_dir / f"{model_tag(model)}.json"
        return from_archived_v6(path, expected_model=model,
                                vev_gev=args.vev, normalization_factor=args.spinor_factor)
    from scripts.RunT3GeneralSequentialLeadingLog import process, load
    registry = load(args.registry)
    if model not in registry.get("models", {}):
        raise PhysicsInputError(f"{model}: model absent from source quartic registry")
    target = None
    if args.with_smeft:
        from Numerical.fitting.OscillationFit import OscillationFitTarget
        from Numerical.orchestration.RunSmallLambdaRefit import TARGET
        target = OscillationFitTarget.from_json(TARGET)
    raw = process(model, dynamic_paths(args), registry,
                  spinor_factor=complex(1),  # explicit factor applied *once* in assembler
                  include_extra=True, include_portal=True,
                  with_smeft=args.with_smeft, vev=args.vev, target=target)
    result = from_dynamic_case(raw, model=model, vev_gev=args.vev,
                               normalization_factor=args.spinor_factor,
                               with_smeft=args.with_smeft)
    # Preserve a few key scope markers; avoid storing unnecessary ablation cases.
    result["general_extra_hhss_enabled"] = raw.get("general_extra_hhss_enabled")
    result["full_rgbeta_portal_completion_enabled"] = raw.get(
        "full_rgbeta_portal_completion_enabled")
    return result


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    columns = ("model", "mass_class", "status", "m1_eV_per_abs_N", "m2_eV_per_abs_N",
               "m3_eV_per_abs_N", "C5_LLSS_norm_GeV_inverse", "projection_residual",
               "matching_scale_residual", "complete_C5", "error")
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        for item in rows:
            masses = item.get("masses_eV_per_abs_unknown_NGamma") or [None] * 3
            writer.writerow(dict(model=item["model"], mass_class=item["mass_class"],
                status=item["status"], m1_eV_per_abs_N=masses[0],
                m2_eV_per_abs_N=masses[1], m3_eV_per_abs_N=masses[2],
                C5_LLSS_norm_GeV_inverse=item.get("C5_LLSS_frobenius_GeV_inverse"),
                projection_residual=item.get("SU2_Weinberg_projection_residual"),
                matching_scale_residual=item.get("matching_scale_cancellation_residual"),
                complete_C5=item.get("complete_C5", False), error=item.get("error")))


def run(args: argparse.Namespace) -> int:
    models = select_models(args)
    info = settings(args, models)
    if args.check_only:
        preflights = [preflight(m, args) for m in models]
        for row in preflights:
            print(f"{row['model']} ({row['mass_class']}): {row['status']}", flush=True)
            if "error" in row:
                print("  " + row["error"], flush=True)
        print("Local D(-2): " + dminus2_existing_reference(args.dminus2_manifest)["status"])
        return 0 if all(x["status"] == "Ready" for x in preflights) else 2

    if args.run_id in ("", ".", "..") or any(x in args.run_id for x in "/\\:"):
        raise PhysicsInputError("run-id must be a single folder name")
    out = OUTPUT_ROOT / args.run_id
    config = out / "settings.json"
    if out.exists():
        if not args.resume:
            raise PhysicsInputError(f"Output exists: {out}. Use --resume or a new --run-id")
        if read_json(config) != info:
            raise PhysicsInputError("Cannot resume: settings changed")
    else:
        if args.resume:
            raise PhysicsInputError(f"Cannot resume: run directory absent: {out}")
        out.mkdir(parents=True)
        write_json(config, info)
    rows = []
    for model in models:
        path = out / model_tag(model) / "result.json"
        old: dict[str, Any] | None = None
        try:
            hashes = input_hashes(model, args)
            if args.resume and path.is_file():
                candidate = read_json(path)
                if (candidate.get("status") == "ConditionalLLSSContributionComputed" and
                        candidate.get("input_sha256") == hashes and
                        candidate.get("complete_C5") is False):
                    old = candidate
            if old is not None:
                print(f"{model}: reused validated previous output", flush=True)
                row = old
            else:
                row = evaluate_model(model, args)
                row["input_sha256"] = hashes
                row["generated_utc"] = datetime.now(timezone.utc).isoformat()
                write_json(path, row)
                print(f"{model}: computed conditional LLSS C5", flush=True)
        except (Exception) as exc:
            # An incomplete model must be visible in the summary, never zeroed.
            row = {"model": model, "mass_class": mass_class(model),
                   "status": "Blocked", "complete_C5": False,
                   "error": f"{type(exc).__name__}: {exc}",
                   "traceback": traceback.format_exc(limit=8)}
            write_json(out / model_tag(model) / "error.json", row)
            print(f"{model}: BLOCKED: {row['error']}", flush=True)
        rows.append(row)
    write_csv(out / "comparison.csv", rows)
    summary = {
        "format": FORMAT, "models_requested": list(models),
        "conditional_LLSS_completed": sum(r["status"] == "ConditionalLLSSContributionComputed"
                                          for r in rows),
        "model_count": len(rows),
        "complete_sequential_models_newly_certified": 0,
        "source_backend": args.backend,
        "D_minus2_verified_reference": dminus2_existing_reference(args.dminus2_manifest),
        "interpretation": (
            "Conditional, source-backed LLSS scalar-pair O(1/16pi^2) contribution "
            "after integrating all fermions, running EFT1 and integrating both "
            "scalars. Independent hard one-loop matching, full dim-5 closure, "
            "absolute GammaCC normalization and full SMEFT C5 NOT certified. "
            "Values are eV per absolute unknown spinor convention factor."
        ),
        "missing_for_complete_Dplus2": [
            "independent hard one-loop matching from physical EFT1 boundary",
            "full operator-inventory closure and finite-part scheme alignment",
            "Matchete/RGBeta quartic and spinor/CG normalization dictionary",
            "exact one-loop pole-cancellation test in target model",
        ],
        "rows": [{k: v for k, v in r.items() if k not in (
            "C5_LLSS_total_GeV_inverse", "C5_LLSS_direct_log_GeV_inverse",
            "C5_LLSS_joint_scalar_bubble_GeV_inverse", "traceback")}
            for r in rows],
    }
    write_json(out / "summary.json", summary)
    print(f"Conditional {summary['conditional_LLSS_completed']}/{len(rows)}; "
          f"full sequential independently certified 0 new models.", flush=True)
    print("Summary: " + str(out / "summary.json"), flush=True)
    return 0 if summary["conditional_LLSS_completed"] == len(rows) else 2


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--run-id", default="first3_v1")
    p.add_argument("--models", default=None,
                   help="Comma-separated models, default D_+2,A_-4,B_-1")
    p.add_argument("--all16", action="store_true")
    p.add_argument("--backend", choices=("dynamic", "archived-v6"), default="dynamic")
    p.add_argument("--check-only", action="store_true")
    p.add_argument("--resume", action="store_true")
    p.add_argument("--with-smeft", action="store_true",
                   help="Run existing conditional SMEFT C5 RGE, dynamic backend only")
    p.add_argument("--spinor-factor", type=float, default=1.0,
                   help="Explicit provisional signed normalization multiplier (default per N=1)")
    p.add_argument("--vev", type=float, default=246.22)
    p.add_argument("--release-dir", type=Path, default=DEFAULT_RELEASE)
    p.add_argument("--v4-channel-dir", type=Path, default=DEFAULT_CHANNELS)
    p.add_argument("--quartic-seed-dir", type=Path, default=DEFAULT_SEEDS)
    p.add_argument("--running-dir", type=Path, default=DEFAULT_RUNNING)
    p.add_argument("--registry", type=Path, default=DEFAULT_REGISTRY)
    p.add_argument("--archived-dir", type=Path, default=ROOT / ARCHIVED_REL)
    p.add_argument("--dminus2-manifest", type=Path, default=DEFAULT_D2)
    args = p.parse_args(argv)
    if not math.isfinite(args.spinor_factor) or args.spinor_factor == 0 or (
            not math.isfinite(args.vev) or args.vev <= 0):
        p.error("Invalid normalization factor or vev")
    if args.backend == "archived-v6" and args.with_smeft:
        p.error("--with-smeft requires --backend dynamic")
    try:
        return run(args)
    except (PhysicsInputError, OSError, ValueError) as exc:
        print(f"[T3 sequential batch] ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
