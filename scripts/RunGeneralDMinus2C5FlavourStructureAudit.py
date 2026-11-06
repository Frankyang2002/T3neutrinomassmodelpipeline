"""Read-only C5 field-content and flavour-dummy inventory for general D(-2)."""
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
SCRIPT = ROOT / "Lagrangian" / "RunGeneralDMinus2C5FlavourStructureAudit.wl"
DEFAULT_MATCHING = ROOT / "output/general_sequential_pilot/D_m2_stage1"
PROFILE_REL = Path("data/general_d_m2_scalar_threshold/profile_tree")
REPORT = "profile_c5_flavour_structure_result.json"
SUMMARY = "profile_c5_flavour_structure_summary.json"
LOG = "profile_c5_flavour_structure.log"


def classify(code: int | None, timed_out: bool, report: dict[str, Any] | None) -> str:
    if timed_out:
        return "ExternalTimeout"
    if code != 0:
        return "WolframFailed"
    if not isinstance(report, dict):
        return "MissingResult"
    if report.get("Status") != "WeinbergStructureAudited":
        return "UncertifiedStructure"
    if not all(report.get(x) is True for x in (
        "NoNewMatching", "AllCandidateTermsHaveTwoLeptonsTwoHiggs",
        "AllCandidateTermsHaveUnambiguousChirality",
        "HardC5IsLinearInHbar", "InheritedC5IsLinearInHbar")):
        return "InvalidStructureChecks"
    if not all(report.get(x) is False for x in (
        "FlavourIndexMappingCertified", "C5SymmetricNormalizationCertified",
        "FixedOrderSMEFTAssembled")):
        return "PrematureC5Certification"
    return "WeinbergStructureAudited"


def terminate_tree(p: subprocess.Popen[bytes]) -> None:
    if p.poll() is not None:
        return
    if os.name == "nt":
        try:
            subprocess.run(["taskkill", "/PID", str(p.pid), "/T", "/F"],
                           capture_output=True, timeout=20, check=False)
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


def run(matching_dir: Path = DEFAULT_MATCHING, timeout_seconds: int = 180,
        startup_allowance_seconds: int = 120,
        wolframscript: str = "wolframscript") -> dict[str, Any]:
    if not 15 <= timeout_seconds <= 3600:
        raise ValueError("--timeout-seconds must be 15..3600")
    if not 30 <= startup_allowance_seconds <= 3600:
        raise ValueError("--startup-allowance-seconds must be 30..3600")
    matching_dir = matching_dir.resolve()
    preflight(matching_dir, script=SCRIPT)
    profile = matching_dir / PROFILE_REL
    profile.mkdir(parents=True, exist_ok=True)
    (profile / REPORT).unlink(missing_ok=True)
    cmd = [wolframscript, "-file", str(SCRIPT.resolve()), str(matching_dir)]
    env = os.environ.copy()
    env.update(T3_EXPERIMENTAL_GENERAL_QUADRATICS="1",
               T3_EXPERIMENTAL_GENERAL_FERMIONS="1")
    cap = timeout_seconds + startup_allowance_seconds
    started = time.monotonic()
    print(f"[C5 structure] Reading saved WXF; no Match; external limit {cap}s", flush=True)
    timed_out = False
    with (profile / LOG).open("wb") as stream:
        flags = getattr(subprocess, "CREATE_NEW_PROCESS_GROUP", 0) if os.name == "nt" else 0
        p = subprocess.Popen(cmd, cwd=ROOT, env=env,
                             stdout=stream, stderr=subprocess.STDOUT,
                             creationflags=flags,
                             start_new_session=(os.name != "nt"))
        heartbeat = started
        while p.poll() is None:
            now = time.monotonic()
            if now - started >= cap:
                timed_out = True
                terminate_tree(p)
                break
            if now - heartbeat >= 20:
                print(f"[C5 structure] Active for {int(now-started)}s", flush=True)
                heartbeat = now
            time.sleep(0.5)
        exit_code = p.poll()
    path = profile / REPORT
    report = json.loads(path.read_text(encoding="utf-8-sig")) if path.is_file() else None
    status = classify(exit_code, timed_out, report)
    summary = {
        "status": status, "exit_code": exit_code,
        "elapsed_seconds": round(time.monotonic() - started, 2),
        "report": str(path), "log": str(profile / LOG),
        "hard": (report or {}).get("HardWeinberg"),
        "inherited": (report or {}).get("InheritedWeinberg"),
    }
    (profile / SUMMARY).write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"[C5 structure] {status} (exit {exit_code})\n"
          f"  Result: {path}\n  Log: {profile / LOG}", flush=True)
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matching-dir", type=Path, default=DEFAULT_MATCHING)
    parser.add_argument("--timeout-seconds", type=int, default=180)
    parser.add_argument("--startup-allowance-seconds", type=int, default=120)
    parser.add_argument("--wolframscript", default="wolframscript")
    args = parser.parse_args()
    try:
        r = run(args.matching_dir, args.timeout_seconds,
                args.startup_allowance_seconds, args.wolframscript)
    except (OSError, ValueError, RuntimeError, json.JSONDecodeError) as exc:
        print(f"[C5 structure] ERROR: {exc}", flush=True)
        return 1
    return 0 if r["status"] == "WeinbergStructureAudited" else 1


if __name__ == "__main__":
    raise SystemExit(main())
