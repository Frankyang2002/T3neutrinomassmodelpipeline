"""Audit masses, match full D(-2) EFT1 at one loop, reduce and inspect C5.

Diagnostic only: does not assemble fixed-order SMEFT or change the scalar driver.
Run from the project root:
    python -X utf8 -m scripts.RunGeneralDMinus2FullScalarLoopProbe \
      --match-timeout-seconds 90 --startup-allowance-seconds 180
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
SCRIPT = ROOT / "Lagrangian" / "RunGeneralDMinus2FullScalarLoopProbe.wl"
DEFAULT_MATCHING = ROOT / "output/general_sequential_pilot/D_m2_stage1"
PROFILE_REL = Path("data/general_d_m2_scalar_threshold/profile_tree")
AUDIT = "profile_full_scalar_one_loop_audit.json"
RESULT = "profile_full_scalar_one_loop_result.json"
SUMMARY = "profile_full_scalar_one_loop_summary.json"
LOG = "profile_full_scalar_one_loop.log"


def classify(*, code: int | None, timed_out: bool,
             audit: dict[str, Any] | None,
             result: dict[str, Any] | None) -> str:
    if timed_out:
        return "ExternalTimeout"
    if code != 0:
        return "WolframFailed"
    if not audit or audit.get("Status") != "Passed":
        return "MassAuditFailed"
    if (audit.get("MassPowerCountingPromotionOnly") is not True or
            audit.get("SMAlphaComparisonCertified") is not True or
            audit.get("AllPhysicalScalarsRegisteredHeavy") is not True or
            audit.get("PureHeavyQuadraticAlphaEqualToFreeLag") is not True or
            audit.get("CheckLagrangianPassed") is not True or
            audit.get("FullPhysicalSMIncluded") is not True or
            audit.get("NoLinearHeavyScalarSources") is not True):
        return "MissingMassOrFullEFT1Validation"
    if (not result or result.get("Status") != "FullEFT1OneLoopPostprocessed" or
            result.get("OneLoopRawMatchReturned") is not True or
            result.get("FullPhysicalSMIncluded") is not True or
            result.get("LoopOrder") != 1 or result.get("EFTOrder") != 5 or
            result.get("RemainingHeavyFields") != [] or
            result.get("MatchResultContainsHbar") is not True or
            result.get("ProcessedEFTHasHbar") is not True or
            result.get("HardHbarSectorExtracted") is not True or
            result.get("SimplificationPerformed") is not True or
            len(result.get("ProcessingStages", [])) != 4 or
            any(s.get("Status") != "Returned" for s in result.get("ProcessingStages", []))):
        return "InvalidOneLoopMatchResult"
    return "FullEFT1OneLoopPostprocessed"


def _json(path: Path) -> dict[str, Any] | None:
    if not path.is_file():
        return None
    value = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(value, dict):
        raise ValueError(f"Expected JSON object in {path}")
    return value


def _kill_process_tree(p: subprocess.Popen[bytes]) -> None:
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


def run(*, matching_dir: Path = DEFAULT_MATCHING,
        timeout_seconds: int = 90, startup_allowance_seconds: int = 180,
        wolframscript: str = "wolframscript") -> dict[str, Any]:
    if not 5 <= timeout_seconds <= 3600:
        raise ValueError("match timeout must be between 5 and 3600 seconds")
    if not 30 <= startup_allowance_seconds <= 3600:
        raise ValueError("startup allowance must be between 30 and 3600 seconds")
    matching_dir = matching_dir.resolve()
    ready = preflight(matching_dir, script=SCRIPT)
    profile = matching_dir / PROFILE_REL
    profile.mkdir(parents=True, exist_ok=True)
    for name in (AUDIT, RESULT):
        (profile / name).unlink(missing_ok=True)
    env = os.environ.copy()
    env.update(T3_EXPERIMENTAL_GENERAL_QUADRATICS="1",
               T3_EXPERIMENTAL_GENERAL_FERMIONS="1")
    cmd = [wolframscript, "-file", str(SCRIPT.resolve()),
           str(matching_dir), str(Path(ready["insertion_file"]).resolve()),
           str(timeout_seconds)]
    cap = timeout_seconds * 7 + startup_allowance_seconds
    started = time.monotonic()
    print(f"[full scalar 1L] Mass audit + full-EFT1 1L Match + reduction; "
          f"per stage {timeout_seconds}s; external limit {cap}s", flush=True)
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
                print("[full scalar 1L] External timeout; terminating process tree", flush=True)
                _kill_process_tree(p)
                break
            if now - heartbeat >= 20:
                print(f"[full scalar 1L] Active for {int(now - started)}s", flush=True)
                heartbeat = now
            time.sleep(0.5)
        code = p.poll()
    audit, result = _json(profile / AUDIT), _json(profile / RESULT)
    status = classify(code=code, timed_out=timed_out, audit=audit, result=result)
    summary = {
        "status": status,
        "exit_code": code,
        "elapsed_seconds": round(time.monotonic() - started, 2),
        "diagnostic_only": True,
        "full_scalar_threshold_certified": False,
        "raw_one_loop_match_returned": bool(result and result.get("OneLoopRawMatchReturned")),
        "one_loop_postprocessed": status == "FullEFT1OneLoopPostprocessed",
        "weinberg_extraction_status": result.get("WeinbergExtractionStatus") if result else None,
        "mass_audit_json": str(profile / AUDIT),
        "result_json": str(profile / RESULT),
        "wolfram_log": str(profile / LOG),
        "reason": result.get("Reason") if result else None,
    }
    (profile / SUMMARY).write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"[full scalar 1L] {status} (exit {code})\n"
          f"  Audit: {profile / AUDIT}\n"
          f"  Result: {profile / RESULT}\n"
          f"  Log: {profile / LOG}", flush=True)
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matching-dir", type=Path, default=DEFAULT_MATCHING)
    parser.add_argument("--match-timeout-seconds", type=int, default=90)
    parser.add_argument("--startup-allowance-seconds", type=int, default=180)
    parser.add_argument("--wolframscript", default="wolframscript")
    args = parser.parse_args()
    try:
        result = run(matching_dir=args.matching_dir,
                     timeout_seconds=args.match_timeout_seconds,
                     startup_allowance_seconds=args.startup_allowance_seconds,
                     wolframscript=args.wolframscript)
    except (OSError, ValueError, RuntimeError, json.JSONDecodeError) as exc:
        print(f"[full scalar 1L] ERROR: {exc}", flush=True)
        return 1
    return 0 if result["status"] == "FullEFT1OneLoopPostprocessed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
