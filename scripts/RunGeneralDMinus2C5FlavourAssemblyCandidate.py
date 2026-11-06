"""Provisional fixed-order C5 flavour-index map; reuse saved EFTs, NO Match.

From the repository root:
    python -X utf8 -m scripts.RunGeneralDMinus2C5FlavourAssemblyCandidate
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import signal
import subprocess
import time
from typing import Any

from scripts.RunGeneralDMinus2ScalarThreshold import preflight

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "Lagrangian/RunGeneralDMinus2C5FlavourAssemblyCandidate.wl"
DEFAULT_MATCHING = ROOT / "output/general_sequential_pilot/D_m2_stage1"
PROFILE_REL = Path("data/general_d_m2_scalar_threshold/profile_tree")
REPORT = "profile_c5_flavour_candidate_result.json"
SUMMARY = "profile_c5_flavour_candidate_summary.json"
LOG = "profile_c5_flavour_candidate.log"


def classify(code: int | None, timed_out: bool,
             report: dict[str, Any] | None) -> str:
    if timed_out:
        return "ExternalTimeout"
    if code != 0:
        return "WolframFailed"
    if not isinstance(report, dict) or report.get("Status") != "C5FlavourCandidateExported":
        return "NoValidCandidate"
    checks = {
        "NoNewMatchPerformed": True,
        "EachHolomorphicTermHasTwoDifferentLeptonFlavourLabels": True,
        "HardSymmetricByConstruction": True,
        "InheritedSymmetricByConstruction": True,
        "RunningLeftUnsymmetrized": True,
        "HardCoefficientHasNativeHbar": True,
        "InheritedCoefficientHasNativeHbar": True,
        "RunningInsertedWithHbar": True,
        "TotalCandidateLinearInHbar": True,
        "FixedOrderC5CandidateAssembled": True,
        "OperatorNormalizationIndependentlyCertified": False,
        "InheritedSMCanonicalizationCertified": False,
        "PhysicalScalarMassRelationsImposed": False,
        "FixedOrderSMEFTCertified": False,
    }
    if any(report.get(k) is not v for k, v in checks.items()):
        return "InconsistentCandidateMetadata"
    if report.get("HardHolomorphicTermCount") != 34:
        return "HardTermInventoryChanged"
    if report.get("InheritedHolomorphicTermCount") != 6:
        return "InheritedTermInventoryChanged"
    return "C5FlavourCandidateExported"


def kill_process_tree(p: subprocess.Popen[bytes]) -> None:
    if p.poll() is not None:
        return
    if os.name == "nt":
        try:
            subprocess.run(["taskkill", "/PID", str(p.pid), "/T", "/F"],
                           check=False, capture_output=True, timeout=20)
        except (OSError, subprocess.TimeoutExpired):
            p.kill()
    else:
        try:
            os.killpg(p.pid, signal.SIGTERM)
        except ProcessLookupError:
            pass
    try:
        p.wait(timeout=5)
    except subprocess.TimeoutExpired:
        p.kill()
        p.wait()


def run(matching_dir: Path = DEFAULT_MATCHING,
        timeout_seconds: int = 180,
        startup_allowance_seconds: int = 120,
        wolframscript: str = "wolframscript") -> dict[str, Any]:
    if not 15 <= timeout_seconds <= 3600:
        raise ValueError("--timeout-seconds must be 15..3600")
    if not 30 <= startup_allowance_seconds <= 3600:
        raise ValueError("--startup-allowance-seconds must be 30..3600")
    matching_dir = matching_dir.resolve()
    ready = preflight(matching_dir, script=SCRIPT)
    profile = matching_dir / PROFILE_REL
    profile.mkdir(parents=True, exist_ok=True)
    (profile / REPORT).unlink(missing_ok=True)
    env = os.environ.copy()
    env.update(T3_EXPERIMENTAL_GENERAL_QUADRATICS="1",
               T3_EXPERIMENTAL_GENERAL_FERMIONS="1")
    cmd = [wolframscript, "-file", str(SCRIPT.resolve()),
           str(matching_dir), str(Path(ready["insertion_file"]).resolve())]
    cap = timeout_seconds + startup_allowance_seconds
    started = time.monotonic()
    print(f"[C5 flavour candidate] Saved EFTs only; no Match; external limit {cap}s", flush=True)
    timed_out = False
    with (profile / LOG).open("wb") as stream:
        flags = getattr(subprocess, "CREATE_NEW_PROCESS_GROUP", 0) if os.name == "nt" else 0
        proc = subprocess.Popen(cmd, cwd=ROOT, env=env,
                                stdout=stream, stderr=subprocess.STDOUT,
                                creationflags=flags,
                                start_new_session=(os.name != "nt"))
        heartbeat = started
        while proc.poll() is None:
            now = time.monotonic()
            if now - started >= cap:
                timed_out = True
                kill_process_tree(proc)
                break
            if now - heartbeat >= 20:
                print(f"[C5 flavour candidate] Active for {int(now-started)}s", flush=True)
                heartbeat = now
            time.sleep(0.5)
        code = proc.poll()
    p = profile / REPORT
    report = json.loads(p.read_text(encoding="utf-8-sig")) if p.is_file() else None
    status = classify(code, timed_out, report)
    summary = {
        "status": status, "exit_code": code,
        "elapsed_seconds": round(time.monotonic() - started, 2),
        "candidate_not_certified": True,
        "no_new_matching": True,
        "result": str(p), "log": str(profile / LOG),
        "output_roots": (report or {}).get("OutputRoots"),
        "reason": (report or {}).get("Reason"),
    }
    (profile / SUMMARY).write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"[C5 flavour candidate] {status} (exit {code})\n"
          f"  Result: {p}\n  Log: {profile / LOG}", flush=True)
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matching-dir", type=Path, default=DEFAULT_MATCHING)
    parser.add_argument("--timeout-seconds", type=int, default=180)
    parser.add_argument("--startup-allowance-seconds", type=int, default=120)
    parser.add_argument("--wolframscript", default="wolframscript")
    args = parser.parse_args()
    try:
        summary = run(args.matching_dir, args.timeout_seconds,
                      args.startup_allowance_seconds, args.wolframscript)
    except (OSError, ValueError, RuntimeError, json.JSONDecodeError) as exc:
        print(f"[C5 flavour candidate] ERROR: {exc}", flush=True)
        return 1
    return 0 if summary["status"] == "C5FlavourCandidateExported" else 1


if __name__ == "__main__":
    raise SystemExit(main())
