#!/usr/bin/env python3
"""Full-inventory *raw* one-loop fermion-first Matchete stage for A(-4), D(-2).

Consumes the checksummed numerical RG endpoint from thesis Batch 1.  Calls
Lagrangian/RunT3GeneralFThreshold.wl, which ALREADY uses
RunSequentialT3Matching[{F}, EFTOrder->5, LoopOrder->1].  Records every
Matchete-matched EFT1 one-loop term.  Does not substitute an analytic triangle
or claim UV-minus-EFT1-soft renormalization, a complete numerical Wilson
coefficient basis, or the subsequently required EFT1 Wilson anomalous dimensions.

A successful raw run is a SOURCE checkpoint, not a completed physical Batch 2.
The return code is 3 unless --accept-raw-source-only is explicitly selected.

PowerShell:
 python -X utf8 -m scripts.RunT3ThesisGeneralBatch2 --all --preflight
 python -X utf8 -m scripts.RunT3ThesisGeneralBatch2 --all \
   --run-id thesis_full_general_batch2_raw_v1 --accept-raw-source-only
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import traceback
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from Numerical.parallel_t3.background.general_batch1_api import load_batch1
from Numerical.parallel_t3.wilson.general_batch2_source import (
    MODELS, SCHEMA, BLOCKERS, MatchingIncomplete, batch2_scope, load_json,
    native_numbers, require, sha, verify_ledger, verify_native_wolfram, write_json,
)

BATCH1 = ROOT / "output/parallel_t3/full_general_batch1/thesis_full_general_batch1_v1"
OUTPUT = ROOT / "output/parallel_t3/full_general_batch2"
NATIVE_MATCHER = ROOT / "Lagrangian/RunT3GeneralFThreshold.wl"
LEDGER = ROOT / "Lagrangian/RunT3Batch2OneLoopLedger.wl"
REQUIRED_SOURCE = (
    NATIVE_MATCHER, LEDGER, ROOT / "Lagrangian/RunMatching.wl",
    ROOT / "Lagrangian/model/LagrangianBuilder.wl",
    ROOT / "Lagrangian/interactions/GeneralScalarOperators.wl",
    ROOT / "Lagrangian/interactions/ScalarPotential.wl",
    ROOT / "Lagrangian/model/T3Fields.wl",
    ROOT / "scripts/RunT3GeneralFThresholdBatch.py",
    ROOT / "RGE/general/generated/general_rgbeta_quartic_registry.json",
    ROOT / "Numerical/parallel_t3/wilson/general_batch2_source.py",
)
ENV_GENERAL = (
    "T3_EXPERIMENTAL_GENERAL_FERMIONS",
    "T3_EXPERIMENTAL_GENERAL_QUARTICS",
    "T3_EXPERIMENTAL_GENERAL_QUADRATICS",
    "T3_EXPERIMENTAL_FLAVOR_MASSES",
)


def safe_id(value: str) -> str:
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]{0,63}", value):
        raise argparse.ArgumentTypeError("ID must be alphanumeric/underscore/hyphen (<=64 chars)")
    return value


def check_b1(model: str, base: Path) -> tuple[Any, dict, dict]:
    run = load_batch1(base / model)
    require(run.model == model, "wrong Batch1 trajectory model")
    require(run.manifest.get("source_sha256") is not None,
            "Batch1 lacks symbolic source provenance")
    manifest = run.manifest
    require(manifest.get("fermion_threshold_wilson_matching_performed") is False
            and manifest.get("EFT1_wilson_running_performed") is False,
            "input mislabelled as already having Wilson physics")
    require(manifest.get("scales_gev") == {"UV": 1e7, "F": 1e5, "S": 1e3, "low": 100},
            "unexpected thesis threshold scales")
    src = run.fermion_matching_input()
    require(src.get("model") == model and src.get("scale_GeV") == 1e5,
            "wrong F-matching JSON")
    require(src.get("all_heavy_fermions_to_integrate_out_together") is True
            and src.get("both_heavy_scalar_multiplets_remain_light") is True,
            "Batch1 threshold topology does not retain both scalars")
    values = src.get("native_uv_rgbeta_couplings_at_muF")
    require(isinstance(values, dict), "missing native RG at F")
    # Bitwise equality not required: native JSON stores machine precision
    # whereas NPZ is complex128. Compare every coupling with exact saved state.
    from Numerical.parallel_t3.background.general_batch1_api import decode_native
    require(set(values) == set(run.at_fermion_threshold()),
            "fermion input JSON differs from the complete native trajectory")
    import numpy as np
    for key, value in values.items():
        require(np.allclose(decode_native(value), run.at_fermion_threshold()[key],
                            rtol=1e-13, atol=1e-14),
                "F-threshold input/trajectory mismatch: " + key)
    return run, src, native_numbers(values, model)


def source_hashes() -> dict[str, str]:
    return {str(p.relative_to(ROOT)).replace("\\", "/"): sha(p) for p in REQUIRED_SOURCE}


def launch(command: list[str], workdir: Path, log: Path,
           env: dict[str, str], timeout: int) -> None:
    try:
        result = subprocess.run(command, cwd=ROOT, capture_output=True, text=True,
                                check=False, env=env, timeout=timeout)
    except subprocess.TimeoutExpired as exc:
        raise RuntimeError(f"Wolfram stage exceeded configured timeout ({timeout}s)") from exc
    log.parent.mkdir(parents=True, exist_ok=True)
    log.write_text("command: " + " ".join(command) + "\nexit: "
                   + str(result.returncode) + "\nstdout:\n" + result.stdout
                   + "\nstderr:\n" + result.stderr, encoding="utf-8")
    if result.returncode:
        raise RuntimeError(f"Wolfram exited {result.returncode}; see {log}; "
                           f"tail={(result.stdout + result.stderr)[-800:]}")


def build_model(model: str, args: argparse.Namespace, hashes: dict) -> dict:
    run, fsrc, numeric = check_b1(model, args.batch1_root)
    b1_dir = args.batch1_root / model
    inputs = {
        "batch1_manifest": sha(b1_dir / "manifest.json"),
        "batch1_fermion_matching_input": sha(b1_dir / "fermion_matching_input.json"),
        "batch1_scalar_matching_input": sha(b1_dir / "scalar_matching_input.json"),
        "batch1_trajectory": sha(b1_dir / "trajectories.npz"),
    }
    record = {
        "format": SCHEMA,
        "model": model,
        "status": "PreflightPassed_SourceMatchedNotExecuted",
        "batch1_run_id": run.manifest["run_id"],
        "input_sha256": inputs,
        "matcher_source_sha256": hashes,
        "native_fermion_threshold_values": numeric,
        "all_heavy_fermions_one_threshold": True,
        "both_scalar_multiplets_light": True,
        "matching_scale_GeV": 1e5,
        "matching_source": "Lagrangian/RunT3GeneralFThreshold.wl",
        "matching_function": "RunSequentialT3Matching",
        "full_symmetry_allowed_source_interactions_requested": True,
        "native_running_values_inserted_into_Matchete_Lagrangian": False,
        **batch2_scope(),
    }
    if args.preflight:
        return record
    out = args.output_root / args.run_id / model
    require(not out.exists(), "refusing to overwrite or mix outputs: " + str(out))
    out.mkdir(parents=True)
    write_json(out / "numerical_fermion_boundary.json", numeric)
    write_json(out / "manifest.json", record)
    match_dir = out / "native_matching"
    match_dir.mkdir()
    environment = dict(os.environ)
    environment.update({k: "1" for k in ENV_GENERAL})
    try:
        launch([args.wolframscript, "-file", str(NATIVE_MATCHER),
                str(match_dir), model, "1"], ROOT, out / "native_matching.log",
               environment, args.timeout_sec)
        native = verify_native_wolfram(match_dir, model)
        raw = match_dir / "data/eft1_after_F_one_loop_inputform.txt"
        ledger_path = out / "one_loop_operator_ledger.json"
        launch([args.wolframscript, "-file", str(LEDGER), str(raw),
                str(ledger_path)], ROOT, out / "operator_ledger.log",
               environment, args.timeout_sec)
        ledger = verify_ledger(load_json(ledger_path), raw)
        record.update({
            "status": "RawFullGeneralOneLoopFThresholdSourceExported_NotFiniteWilsons",
            "matching_source_validated": True,
            "native_wolfram": native,
            "operator_ledger": ledger,
            "match_outputs_directory": "native_matching",
            "numerical_fermion_boundary_file": "numerical_fermion_boundary.json",
            "output_files_sha256": {
                "operator_ledger": sha(ledger_path),
                "native_status": sha(match_dir / "native_status.json"),
            },
        })
    except Exception as exc:
        record.update({"status": "Blocked_" + type(exc).__name__,
                       "reason": str(exc),
                       "traceback_tail": traceback.format_exc(limit=5)})
    write_json(out / "manifest.json", record)
    return record


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--all", action="store_true")
    p.add_argument("--model", choices=MODELS)
    p.add_argument("--batch1-root", type=Path, default=BATCH1)
    p.add_argument("--output-root", type=Path, default=OUTPUT)
    p.add_argument("--run-id", type=safe_id, default="thesis_full_general_batch2_raw_v1")
    p.add_argument("--wolframscript", default="wolframscript")
    p.add_argument("--timeout-sec", type=int, default=7200)
    p.add_argument("--preflight", action="store_true")
    p.add_argument("--accept-raw-source-only", action="store_true",
                   help="Return 0 on complete raw one-loop Matchete export, NOT finite Wilsons")
    args = p.parse_args(argv)
    if not args.all and args.model is None:
        p.error("provide --all or --model")
    if args.timeout_sec < 1:
        p.error("--timeout-sec must be positive")
    if not args.preflight and not (shutil.which(args.wolframscript)
                                   or Path(args.wolframscript).is_file()):
        p.error("wolframscript unavailable; cannot perform the real one-loop matching")
    if not args.preflight and (args.output_root / args.run_id).exists():
        p.error("output run-id already exists; choose a new --run-id to protect old results")
    selected = MODELS if args.all else (args.model,)
    hashes = source_hashes()
    if not args.preflight:
        (args.output_root / args.run_id).mkdir(parents=True, exist_ok=False)
    outputs = []
    for model in selected:
        try:
            rec = build_model(model, args, hashes)
        except Exception as exc:
            rec = {"model": model, "status": "Blocked_" + type(exc).__name__,
                   "reason": str(exc), "scope": batch2_scope()}
            if not args.preflight:
                folder = args.output_root / args.run_id / model
                folder.mkdir(exist_ok=True)
                write_json(folder / "manifest.json", rec)
        outputs.append(rec)
        print(f"{model}: {rec['status']}", flush=True)
        if "reason" in rec: print(f"  Reason: {rec['reason']}", flush=True)
    if args.preflight:
        ok = all(x["status"] == "PreflightPassed_SourceMatchedNotExecuted" for x in outputs)
        print("Source and Batch1 numerical preflight:", "PASS" if ok else "BLOCKED")
        return 0 if ok else 2
    successful_raw = all(x["status"] ==
                         "RawFullGeneralOneLoopFThresholdSourceExported_NotFiniteWilsons"
                         for x in outputs)
    summary = {
        "format": SCHEMA,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "models": [{"model": x["model"], "status": x["status"],
                    "reason": x.get("reason")} for x in outputs],
        "models_with_complete_raw_matchete_export": sum(
            x["status"] == "RawFullGeneralOneLoopFThresholdSourceExported_NotFiniteWilsons"
            for x in outputs),
        **batch2_scope(),
    }
    write_json(args.output_root / args.run_id / "summary.json", summary)
    print("Complete finite one-loop F Wilsons and EFT1 Wilson running: NOT CALCULATED")
    print("One-loop native Matchete raw source:", "OK" if successful_raw else "BLOCKED")
    if not successful_raw:
        return 2
    return 0 if args.accept_raw_source_only else 3


if __name__ == "__main__":
    raise SystemExit(main())
