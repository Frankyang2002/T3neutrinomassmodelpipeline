"""Diagnose D(-2) scalar Match with exact heavy sector + canonical SM kinetic terms.

The earlier EFT1-minus-SM input is not a canonically normalized Matchete
Lagrangian and yielded integer zero, NOT a matched EFT. This diagnostic
restores the SM kinetic core from the saved EFT1, verifies the exact
partition, and refuses to simplify numeric/fieldless Match results.

This is NOT the physical scalar threshold: selected pure-SM interactions
are omitted solely to localize the slowdown. No one-loop C5 is certified.

Usage from repository root (PowerShell):
  python -X utf8 -m scripts.RunGeneralDMinus2CanonicalCoreProbe \
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
SCRIPT = ROOT / "Lagrangian" / "RunGeneralDMinus2CanonicalCoreProbe.wl"
DEFAULT_MATCHING = ROOT / "output/general_sequential_pilot/D_m2_stage1"
PROFILE_REL = Path("data/general_d_m2_scalar_threshold/profile_tree")
AUDIT_NAME = "profile_canonical_core_audit.json"
RESULT_NAME = "profile_canonical_core_result.json"
LOG_NAME = "profile_canonical_core.log"
SUMMARY_NAME = "profile_canonical_core_summary.json"


def make_command(matching_dir: Path, insertion_file: Path, *,
                 timeout_seconds: int = 90,
                 wolframscript: str = "wolframscript",
                 script: Path = SCRIPT) -> list[str]:
    if not 5 <= timeout_seconds <= 3600:
        raise ValueError("--match-timeout-seconds must be between 5 and 3600")
    return [wolframscript, "-file", str(script.resolve()),
            str(matching_dir.resolve()), str(insertion_file.resolve()),
            str(timeout_seconds)]


def _stop_tree(process: subprocess.Popen[bytes]) -> None:
    if process.poll() is not None:
        return
    if os.name == "nt":
        try:
            subprocess.run(["taskkill", "/PID", str(process.pid), "/T", "/F"],
                           check=False, capture_output=True, timeout=20)
        except (OSError, subprocess.TimeoutExpired):
            process.kill()
    else:
        try:
            os.killpg(process.pid, signal.SIGTERM)
        except ProcessLookupError:
            return
        try:
            process.wait(timeout=5)
        except subprocess.TimeoutExpired:
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
    try:
        process.wait(timeout=10)
    except subprocess.TimeoutExpired:
        process.kill()
        process.wait()


def _json(path: Path) -> dict[str, Any]:
    content = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(content, dict):
        raise ValueError(f"JSON output is not an object: {path}")
    return content


def classify_result(*, exit_code: int | None, external_timeout: bool,
                    audit: dict[str, Any] | None,
                    result: dict[str, Any] | None) -> str:
    """Never infer the physical result from an exit code or Match return alone."""
    if external_timeout:
        return "ExternalTimeout"
    if exit_code != 0:
        return "WolframFailed"
    if not audit or audit.get("Status") != "Passed":
        return "AuditMissingOrFailed"
    if (audit.get("SMAlphaComparisonCertified") is not True or
            not (audit.get("SavedSMAlphaEquivalentToOriginalLSM") is True or
                 audit.get("SavedSMAlphaEquivalentToRelabelledLSM") is True)):
        return "SMAlphaComparisonNotCertified"
    if (audit.get("StandaloneSMTermsRemaining") != 0 or
            audit.get("AllHeavyTermsPreservedExactly") is not True or
            audit.get("ReconstructionResidualExactlyZero") is not True or
            audit.get("SMSubtractandContainsHeavyFields") is not False or
            audit.get("SMSubtractandExactlyEqualsSavedPureSM") is not True):
        return "SMSubtractionAuditIncomplete"
    if (audit.get("InputMode") != "ExactHeavyPlusSavedSMKineticCore" or
            audit.get("CorePlusDroppedSMReconstructsFullEFT1") is not True or
            audit.get("CoreHeavySectorExactlyPreserved") is not True or
            audit.get("UnclassifiedCouplinglessSMTerms") != 0 or
            audit.get("CheckLagrangianPassed") is not True or
            audit.get("FullPhysicalSMIncluded") is not False or
            audit.get("SMKineticTermCount", 0) <= 0 or
            audit.get("DroppedSMInteractionTermCount", 0) <= 0):
        return "CanonicalCoreNotCertified"
    if not result or result.get("Status") != "CompleteCanonicalCoreTreeDiagnostic":
        return "TreeProcessingIncomplete"
    if (result.get("InputMode") != "ExactHeavyPlusSavedSMKineticCore" or
            result.get("MatchInputCheckLagrangianPassed") is not True or
            result.get("FullPhysicalSMIncluded") is not False or
            result.get("SMKineticCoreRestored") is not True):
        return "InvalidDiagnosticInputContract"
    if result.get("OneLoopScalarHardMatchingPerformed") is not False:
        return "UnexpectedResultContract"
    if result.get("RemainingHeavyFields") != []:
        return "HeavyFieldsRemain"
    stages = result.get("Stages")
    if not isinstance(stages, list) or [s.get("Stage") for s in stages] != [
        "Match", "GreensSimplify", "EOMSimplify", "EvaluateLoopFunctions",
        "ReplaceEffectiveCouplings",
    ]:
        return "MissingTreeProcessingStage"
    if any(s.get("Status") != "Returned" for s in stages):
        return "FailedTreeProcessingStage"
    if any(s.get("Head") in {"Integer", "Real", "Rational", "Complex", "Symbol"}
           for s in stages):
        return "InvalidEFTStageReturn"
    return "CompleteCanonicalCoreTreeDiagnostic"


def run(*, matching_dir: Path = DEFAULT_MATCHING,
        timeout_seconds: int = 90,
        startup_allowance_seconds: int = 180,
        wolframscript: str = "wolframscript",
        script: Path = SCRIPT) -> dict[str, Any]:
    matching_dir = matching_dir.resolve()
    if not 5 <= timeout_seconds <= 3600:
        raise ValueError("--match-timeout-seconds must be between 5 and 3600")
    if not 30 <= startup_allowance_seconds <= 3600:
        raise ValueError("--startup-allowance-seconds must be between 30 and 3600")
    # Import the installed local preflight rather than inventing or rebuilding
    # missing first-threshold data.  The new WL file is a separate new entry.
    ready = preflight(matching_dir, script=script)
    profile_dir = matching_dir / PROFILE_REL
    profile_dir.mkdir(parents=True, exist_ok=True)
    audit_path = profile_dir / AUDIT_NAME
    result_path = profile_dir / RESULT_NAME
    log_path = profile_dir / LOG_NAME
    summary_path = profile_dir / SUMMARY_NAME
    for path in (audit_path, result_path):
        path.unlink(missing_ok=True)
    command = make_command(matching_dir, Path(ready["insertion_file"]),
                           timeout_seconds=timeout_seconds,
                           wolframscript=wolframscript, script=script)
    env = os.environ.copy()
    env["T3_EXPERIMENTAL_GENERAL_QUADRATICS"] = "1"
    env["T3_EXPERIMENTAL_GENERAL_FERMIONS"] = "1"
    hard_cap = timeout_seconds * 5 + startup_allowance_seconds
    # Five bounded processing steps, plus setup.  Native/front-end calls may
    # ignore Mathematica TimeConstrained, so the external cap remains necessary.
    started = time.monotonic()
    timed_out = False
    print(f"[canonical core] Fresh Wolfram kernel; per-stage limit {timeout_seconds}s; "
          f"external process limit {hard_cap}s", flush=True)
    with log_path.open("wb") as output:
        flags = getattr(subprocess, "CREATE_NEW_PROCESS_GROUP", 0) if os.name == "nt" else 0
        process = subprocess.Popen(command, cwd=ROOT, env=env,
                                   stdout=output, stderr=subprocess.STDOUT,
                                   creationflags=flags,
                                   start_new_session=(os.name != "nt"))
        heartbeat = started
        while process.poll() is None:
            now = time.monotonic()
            if now - started >= hard_cap:
                timed_out = True
                print("[canonical core] External timeout; terminating Wolfram process tree",
                      flush=True)
                _stop_tree(process)
                break
            if now - heartbeat >= 20:
                print(f"[canonical core] Active for {int(now - started)}s", flush=True)
                heartbeat = now
            time.sleep(0.5)
        exit_code = process.poll()
    audit = _json(audit_path) if audit_path.is_file() else None
    result = _json(result_path) if result_path.is_file() else None
    status = classify_result(exit_code=exit_code, external_timeout=timed_out,
                             audit=audit, result=result)
    summary: dict[str, Any] = {
        "status": status,
        "diagnostic_only": True,
        "input_mode": "ExactHeavyPlusSavedSMKineticCore",
        "full_physical_sm_included": False,
        "full_general_scalar_threshold_certified": False,
        "one_loop_scalar_matching_performed": False,
        "elapsed_seconds": round(time.monotonic() - started, 2),
        "exit_code": exit_code,
        "failure_reason": result.get("Reason") if result else None,
        "failed_stages": result.get("Stages", []) if result and exit_code != 0 else [],
        "external_timeout_seconds": hard_cap,
        "wolfram_log": str(log_path),
        "audit_json": str(audit_path) if audit is not None else None,
        "result_json": str(result_path) if result is not None else None,
    }
    summary_path.write_text(json.dumps(summary, indent=2) + "\n",
                            encoding="utf-8")
    print(f"[canonical core] {status} (exit {exit_code})\n"
          f"  Summary: {summary_path}\n  Log: {log_path}", flush=True)
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
    except (OSError, RuntimeError, ValueError, KeyError) as error:
        parser.exit(1, f"[canonical core] ERROR: {error}\n")
    return 0 if result["status"] == "CompleteCanonicalCoreTreeDiagnostic" else 1


if __name__ == "__main__":
    raise SystemExit(main())
