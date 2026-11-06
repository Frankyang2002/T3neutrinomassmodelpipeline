"""Audit existing D(-2) hard, inherited and running C5 sources without new matching.

Run from the project root:
    python -X utf8 -m scripts.RunGeneralDMinus2FixedOrderComponents

This is a diagnostic; it does not assemble or certify SMEFT C5.
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
SCRIPT = ROOT / "Lagrangian" / "RunGeneralDMinus2FixedOrderComponents.wl"
DEFAULT_MATCHING = ROOT / "output/general_sequential_pilot/D_m2_stage1"
PROFILE_REL = Path("data/general_d_m2_scalar_threshold/profile_tree")
REPORT = "profile_fixed_order_components_result.json"
SUMMARY = "profile_fixed_order_components_summary.json"
LOG = "profile_fixed_order_components.log"


def classify(code: int | None, timed_out: bool,
             report: dict[str, Any] | None) -> str:
    if timed_out:
        return "ExternalTimeout"
    if code != 0:
        return "WolframFailed"
    if report is None or report.get("Status") != "ComponentsAudited":
        return "InvalidComponentsReport"
    expected = {
        "NoNewMatchPerformed": True,
        "TreeHasLinearHeavySource": False,
        "SavedHardSectorMatchesProcessedO1": True,
        "SavedHardSectorContainsHeavy": False,
        "HardC5LinearInHbar": True,
        "HardCandidateMatchesSavedText": True,
        "DirectRunningHasLog": True,
        "DirectRunningHasExplicitHbar": False,
        "InheritedPureSMHasZeroHbarPart": True,
        "InheritedPureSMHasHbarSquared": False,
        "FixedOrderSMEFTAssembled": False,
        "FullGeneralScalarThresholdCertified": False,
    }
    if any(report.get(k) is not v for k, v in expected.items()):
        return "InconsistentComponentsReport"
    if report.get("HardWeinbergExtractionStatus") != "Success":
        return "HardC5ExtractionNotVerified"
    if report.get("InheritedPureSMWeinbergExtractionStatus") not in {"Success", "NotFound"}:
        return "InheritedC5ExtractionPending"
    if report.get("FlavourAssemblyCertified") is not False:
        return "InvalidFlavourCertification"
    return "ComponentsAudited"


def read_json(path: Path) -> dict[str, Any] | None:
    if not path.is_file():
        return None
    value = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(value, dict):
        raise ValueError(f"Expected JSON object in {path}")
    return value


def kill_process_tree(p: subprocess.Popen[bytes]) -> None:
    if p.poll() is not None:
        return
    if os.name == "nt":
        try:
            subprocess.run(["taskkill", "/PID", str(p.pid), "/T", "/F"],
                           capture_output=True, check=False, timeout=20)
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
            try:
                os.killpg(p.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
    try:
        p.wait(timeout=10)
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
    print(f"[fixed-order] Auditing saved EFTs; no new Match; external limit {cap}s", flush=True)
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
                print("[fixed-order] External timeout; terminating process tree", flush=True)
                kill_process_tree(p)
                break
            if now - heartbeat >= 20:
                print(f"[fixed-order] Active for {int(now - started)}s", flush=True)
                heartbeat = now
            time.sleep(0.5)
        code = p.poll()
    report = read_json(profile / REPORT)
    status = classify(code, timed_out, report)
    summary = {
        "status": status,
        "exit_code": code,
        "elapsed_seconds": round(time.monotonic() - started, 2),
        "diagnostic_only": True,
        "new_match_performed": False,
        "fixed_order_smeft_assembled": False,
        "hard_weinberg_term_count": (report or {}).get("HardWeinbergTermCount"),
        "inherited_weinberg_status": (report or {}).get("InheritedPureSMWeinbergExtractionStatus"),
        "result_json": str(profile / REPORT),
        "wolfram_log": str(profile / LOG),
        "reason": (report or {}).get("Reason")
    }
    (profile / SUMMARY).write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"[fixed-order] {status} (exit {code})\n"
          f"  Result: {profile / REPORT}\n  Log: {profile / LOG}", flush=True)
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matching-dir", type=Path, default=DEFAULT_MATCHING)
    parser.add_argument("--timeout-seconds", type=int, default=180)
    parser.add_argument("--startup-allowance-seconds", type=int, default=120)
    parser.add_argument("--wolframscript", default="wolframscript")
    args = parser.parse_args()
    try:
        summary = run(matching_dir=args.matching_dir,
                      timeout_seconds=args.timeout_seconds,
                      startup_allowance_seconds=args.startup_allowance_seconds,
                      wolframscript=args.wolframscript)
    except (OSError, ValueError, RuntimeError, json.JSONDecodeError) as exc:
        print(f"[fixed-order] ERROR: {exc}", flush=True)
        return 1
    return 0 if summary["status"] == "ComponentsAudited" else 1


if __name__ == "__main__":
    raise SystemExit(main())
