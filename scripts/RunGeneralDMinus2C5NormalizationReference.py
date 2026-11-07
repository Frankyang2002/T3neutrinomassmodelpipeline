"""Reference-check the D(-2) provisional Weinberg coefficient normalization.

No matching, WXF rewriting, or parameter substitution is performed.

PowerShell:
    python -X utf8 -m scripts.RunGeneralDMinus2C5NormalizationReference
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

ROOT = Path(__file__).resolve().parents[1]
WL = ROOT / "Lagrangian" / "RunGeneralDMinus2C5NormalizationReference.wl"
DEFAULT_OUTPUT = ROOT / "output/general_sequential_pilot/D_m2_stage1"
PROFILE_REL = Path("data/general_d_m2_scalar_threshold/profile_tree")
REPORT = "profile_c5_normalization_reference_result.json"
SUMMARY = "profile_c5_normalization_reference_summary.json"
LOG = "profile_c5_normalization_reference.log"
CANDIDATE_REPORT = "profile_c5_flavour_candidate_result.json"
CANDIDATE_WXF = "profile_c5_flavour_candidate_total_provisional.wxf"


def read_object(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise FileNotFoundError(path)
    try:
        value = json.loads(path.read_text(encoding="utf-8-sig"))
    except json.JSONDecodeError as exc:
        raise ValueError(f"Invalid JSON in {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return value


def preflight(output_dir: Path) -> Path:
    if not WL.is_file():
        raise FileNotFoundError(WL)
    profile = output_dir.resolve() / PROFILE_REL
    prior = read_object(profile / CANDIDATE_REPORT)
    required = {
        "Status": "C5FlavourCandidateExported",
        "FixedOrderC5CandidateAssembled": True,
        "RunningSymmetricModuloSummationDummies": True,
        "TotalCandidateLinearInHbar": True,
    }
    for key, expected in required.items():
        if prior.get(key) != expected:
            raise ValueError(f"Prerequisite {key} is {prior.get(key)!r}; need {expected!r}")
    if not (profile / CANDIDATE_WXF).is_file():
        raise FileNotFoundError(profile / CANDIDATE_WXF)
    return profile


def kill_tree(proc: subprocess.Popen[bytes]) -> None:
    if proc.poll() is not None:
        return
    if os.name == "nt":
        try:
            subprocess.run(["taskkill", "/PID", str(proc.pid), "/T", "/F"],
                           capture_output=True, timeout=20, check=False)
        except (OSError, subprocess.TimeoutExpired):
            proc.kill()
    else:
        try:
            os.killpg(proc.pid, signal.SIGTERM)
        except ProcessLookupError:
            pass
        try:
            proc.wait(timeout=5)
        except subprocess.TimeoutExpired:
            try:
                os.killpg(proc.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
    try:
        proc.wait(timeout=10)
    except subprocess.TimeoutExpired:
        proc.kill()
        proc.wait()


def classify(exit_code: int | None, timed_out: bool,
             report: dict[str, Any] | None) -> str:
    if timed_out:
        return "ExternalTimeout"
    if exit_code != 0:
        return "WolframFailed"
    if report is None:
        return "MissingOrInvalidResult"
    if report.get("Status") != "WarsawNormalizationReferenceAudited":
        return "ReferenceNotAudited"
    required = {
        "NoNewMatching": True,
        "ProvisionalCandidateUnchanged": True,
        "ProjectExtractorRecognizesWarsawWeinbergOperator": True,
        "ProjectExtractorLinearOnWarsawReference": True,
        "WarsawOperatorNeutrinoVEVComponentIsPlusVSquaredOverTwo": True,
        "FullPhysicalC5NormalizationCertified": False,
        "PhysicalMassRelationsImposed": False,
        "NumericalBenchmarkPerformed": False,
    }
    if any(report.get(key) is not expected for key, expected in required.items()):
        return "InconsistentReferenceReport"
    return "WarsawNormalizationReferenceAudited"


def run(output_dir: Path = DEFAULT_OUTPUT, timeout_seconds: int = 180,
        startup_allowance_seconds: int = 120,
        wolframscript: str = "wolframscript") -> dict[str, Any]:
    if not (30 <= timeout_seconds <= 3600):
        raise ValueError("--timeout-seconds must be 30..3600")
    if not (30 <= startup_allowance_seconds <= 3600):
        raise ValueError("--startup-allowance-seconds must be 30..3600")
    profile = preflight(output_dir)
    profile.mkdir(parents=True, exist_ok=True)
    report_file = profile / REPORT
    report_file.unlink(missing_ok=True)
    env = os.environ.copy()
    env.update(T3_EXPERIMENTAL_GENERAL_QUADRATICS="1",
               T3_EXPERIMENTAL_GENERAL_FERMIONS="1")
    command = [wolframscript, "-file", str(WL.resolve()), str(output_dir.resolve())]
    limit = timeout_seconds + startup_allowance_seconds
    started = time.monotonic()
    print(f"[C5 normalization] Warsaw-basis reference; no Match; external limit {limit}s",
          flush=True)
    timed_out = False
    with (profile / LOG).open("wb") as stream:
        kwargs: dict[str, Any] = {"cwd": ROOT, "env": env,
            "stdout": stream, "stderr": subprocess.STDOUT}
        if os.name == "nt":
            kwargs["creationflags"] = getattr(subprocess, "CREATE_NEW_PROCESS_GROUP", 0)
        else:
            kwargs["start_new_session"] = True
        proc = subprocess.Popen(command, **kwargs)
        heartbeat = started
        while proc.poll() is None:
            now = time.monotonic()
            if now - started >= limit:
                timed_out = True
                kill_tree(proc)
                break
            if now - heartbeat >= 20:
                print(f"[C5 normalization] Active for {int(now - started)}s", flush=True)
                heartbeat = now
            time.sleep(0.5)
        exit_code = proc.poll()
    bad_json = None
    report = None
    if report_file.exists():
        try:
            report = read_object(report_file)
        except (OSError, ValueError) as exc:
            bad_json = str(exc)
    state = classify(exit_code, timed_out, report)
    summary = {
        "status": state,
        "exit_code": exit_code,
        "elapsed_seconds": round(time.monotonic() - started, 2),
        "no_new_matching": True,
        "full_physical_normalization_certified": False,
        "project_extractor_linear_on_warsaw_reference":
            (report or {}).get("ProjectExtractorLinearOnWarsawReference"),
        "reference_factor": (report or {}).get(
            "WarsawOrderedCoefficientWithReferenceCouplingsSetToOne"),
        "report": str(report_file),
        "log": str(profile / LOG),
        "report_parse_error": bad_json,
        "reason": (report or {}).get("Reason"),
    }
    (profile / SUMMARY).write_text(json.dumps(summary, indent=2) + "\n",
                                    encoding="utf-8")
    print(f"[C5 normalization] {state} (exit {exit_code})\n"
          f"  Result: {report_file}\n  Log: {profile / LOG}", flush=True)
    if state != "WarsawNormalizationReferenceAudited":
        if bad_json:
            print(f"  JSON: {bad_json}", flush=True)
        log_file = profile / LOG
        if log_file.is_file():
            lines = log_file.read_text(encoding="utf-8", errors="replace").splitlines()
            print("  Wolfram log (last 30 lines):", flush=True)
            for line in lines[-30:]:
                print(f"    {line[:300]}", flush=True)
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--timeout-seconds", type=int, default=180)
    parser.add_argument("--startup-allowance-seconds", type=int, default=120)
    parser.add_argument("--wolframscript", default="wolframscript")
    args = parser.parse_args()
    try:
        summary = run(output_dir=args.output_dir,
                      timeout_seconds=args.timeout_seconds,
                      startup_allowance_seconds=args.startup_allowance_seconds,
                      wolframscript=args.wolframscript)
    except (OSError, ValueError, RuntimeError) as exc:
        print(f"[C5 normalization] ERROR: {exc}", flush=True)
        return 1
    return 0 if summary["status"] == "WarsawNormalizationReferenceAudited" else 1


if __name__ == "__main__":
    raise SystemExit(main())
