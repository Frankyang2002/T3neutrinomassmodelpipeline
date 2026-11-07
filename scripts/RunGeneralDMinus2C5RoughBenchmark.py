"""Run a provisional D(-2) neutrino mass benchmark from the saved C5 WXF.

No matching or modification of previously generated EFTs is performed.

From repository root:
    python -X utf8 -m scripts.RunGeneralDMinus2C5RoughBenchmark
"""
from __future__ import annotations

import argparse
from collections import deque
import json
import os
from pathlib import Path
import signal
import subprocess
import time
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "Lagrangian/RunGeneralDMinus2C5RoughBenchmark.wl"
DEFAULT_OUTPUT = ROOT / "output/general_sequential_pilot/D_m2_stage1"
DEFAULT_CONFIG = ROOT / "benchmarks/general_d_m2_rough_benchmark.json"
PROFILE_REL = Path("data/general_d_m2_scalar_threshold/profile_tree")
REPORT = "profile_c5_rough_benchmark_result.json"
SUMMARY = "profile_c5_rough_benchmark_summary.json"
LOG = "profile_c5_rough_benchmark.log"


def inspect_config(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise FileNotFoundError(f"Missing benchmark configuration: {path}")
    with path.open("r", encoding="utf-8-sig") as file:
        config = json.load(file)
    if not isinstance(config, dict):
        raise ValueError("Benchmark JSON must be an object")
    if config.get("NF") not in range(1, 7):
        raise ValueError("NF must be an integer 1..6")
    if config.get("PoleHandling", "Reject") not in ("Reject", "FinitePartDiagnostic"):
        raise ValueError("PoleHandling must be Reject or FinitePartDiagnostic")
    for name in ("y1", "y2", "y2Conj"):
        matrix = config.get("Yukawas", {}).get(name)
        if not isinstance(matrix, list) or len(matrix) != 3 or any(
            not isinstance(row, list) or len(row) != config["NF"] or any(
                not isinstance(x, (float, int)) or isinstance(x, bool) for x in row)
            for row in matrix
        ):
            raise ValueError(f"{name} must be a numeric 3 x NF matrix")
    return config


def read_report(path: Path) -> tuple[dict[str, Any] | None, str | None]:
    if not path.is_file():
        return None, "Wolfram did not create a report"
    try:
        raw = path.read_text(encoding="utf-8-sig")
    except (OSError, UnicodeError) as exc:
        return None, f"Cannot read report: {exc}"
    if not raw.strip():
        return None, "Empty Wolfram report"
    try:
        value = json.loads(raw)
    except json.JSONDecodeError as exc:
        return None, f"Invalid Wolfram JSON at line {exc.lineno}: {exc.msg}"
    if not isinstance(value, dict):
        return None, "Wolfram JSON is not an object"
    return value, None


def log_tail(path: Path, lines: int = 35) -> str:
    if not path.is_file():
        return "(Log missing)"
    try:
        with path.open(encoding="utf-8", errors="replace") as stream:
            return "".join(deque(stream, maxlen=lines)).rstrip()
    except OSError as exc:
        return f"(Log unreadable: {exc})"


def terminate_process(proc: subprocess.Popen[bytes]) -> None:
    if proc.poll() is not None:
        return
    if os.name == "nt":
        try:
            subprocess.run(["taskkill", "/PID", str(proc.pid), "/T", "/F"],
                           capture_output=True, check=False, timeout=15)
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
        proc.kill()
        proc.wait()


def classify(returncode: int | None, timed_out: bool,
             report: dict[str, Any] | None) -> str:
    if timed_out:
        return "ExternalTimeout"
    if returncode != 0:
        return "WolframFailed"
    if not isinstance(report, dict):
        return "InvalidResultReport"
    valid_statuses = {
        "ProvisionalNumericalBenchmarkEvaluated",
        "ProvisionalFinitePartBenchmarkEvaluated",
    }
    if report.get("Status") not in valid_statuses:
        return "NoNumericalBenchmark"
    if not (report.get("NoNewMatchPerformed") is True
            and report.get("NumericalBenchmarkCertifiedPhysical") is False
            and report.get("PhysicalC5NormalizationIndependentlyCertified") is False
            and isinstance(report.get("NeutrinoMasses_eV"), list)
            and len(report["NeutrinoMasses_eV"]) == 3):
        return "InconsistentBenchmarkMetadata"
    if report.get("Status") == "ProvisionalFinitePartBenchmarkEvaluated":
        if not (report.get("FinitePartOnly") is True
                and report.get("RenormalizationCertified") is False
                and report.get("PoleHandling") == "FinitePartDiagnostic"):
            return "InconsistentBenchmarkMetadata"
    return report["Status"]


def run(output_dir: Path = DEFAULT_OUTPUT, config: Path = DEFAULT_CONFIG,
        timeout_seconds: int = 360,
        wolframscript: str = "wolframscript") -> dict[str, Any]:
    output_dir = output_dir.resolve()
    config = config.resolve()
    inspect_config(config)
    if not SCRIPT.is_file():
        raise FileNotFoundError(f"Wolfram script missing: {SCRIPT}")
    profile = output_dir / PROFILE_REL
    required = [
        profile / "profile_c5_flavour_candidate_total_provisional.wxf",
        profile / "profile_c5_flavour_candidate_result.json",
        profile / "profile_c5_normalization_reference_result.json",
    ]
    missing = [str(f) for f in required if not f.is_file()]
    if missing:
        raise FileNotFoundError("Missing saved inputs:\n" + "\n".join(missing))
    if not 30 <= timeout_seconds <= 3600:
        raise ValueError("--timeout-seconds must be 30..3600")
    profile.mkdir(parents=True, exist_ok=True)
    (profile / REPORT).unlink(missing_ok=True)
    cmd = [wolframscript, "-file", str(SCRIPT.resolve()),
           str(output_dir), str(config)]
    env = os.environ.copy()
    env.update(T3_EXPERIMENTAL_GENERAL_QUADRATICS="1",
               T3_EXPERIMENTAL_GENERAL_FERMIONS="1")
    print(f"[C5 benchmark] Saved C5 only; no Match; external limit {timeout_seconds}s", flush=True)
    started = time.monotonic()
    timed_out = False
    with (profile / LOG).open("wb") as out:
        flags = getattr(subprocess, "CREATE_NEW_PROCESS_GROUP", 0) if os.name == "nt" else 0
        proc = subprocess.Popen(cmd, cwd=ROOT, stdout=out,
                                stderr=subprocess.STDOUT, env=env,
                                creationflags=flags,
                                start_new_session=(os.name != "nt"))
        heartbeat = started
        while proc.poll() is None:
            now = time.monotonic()
            if now - started >= timeout_seconds:
                timed_out = True
                terminate_process(proc)
                break
            if now - heartbeat >= 20:
                print(f"[C5 benchmark] Active for {int(now-started)}s", flush=True)
                heartbeat = now
            time.sleep(0.5)
        code = proc.poll()
    report, error = read_report(profile / REPORT)
    status = classify(code, timed_out, report)
    summary = {
        "status": status,
        "exit_code": code,
        "elapsed_seconds": round(time.monotonic() - started, 2),
        "result": str(profile / REPORT),
        "log": str(profile / LOG),
        "report_error": error,
        "reason": (report or {}).get("Reason"),
        "candidate_not_certified": True,
        "no_new_matching": True,
    }
    (profile / SUMMARY).write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"[C5 benchmark] {status} (exit {code})\n"
          f"  Result: {profile / REPORT}\n  Log: {profile / LOG}", flush=True)
    if error:
        print(f"  JSON error: {error}", flush=True)
    if status in ("ProvisionalNumericalBenchmarkEvaluated", "ProvisionalFinitePartBenchmarkEvaluated"):
        print(f"  Neutrino masses (eV; finite part only if pole present): {report['NeutrinoMasses_eV']}", flush=True)
        print(f"  Max UV-pole residue of C5 (GeV^-1): {report.get('C5PoleResidueMaxAbs_GeVinv')}; renormalization certified: {report.get('RenormalizationCertified')}", flush=True)
    else:
        print(f"[C5 benchmark] {log_tail(profile / LOG, 35)}", flush=True)
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matching-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--timeout-seconds", type=int, default=360)
    parser.add_argument("--wolframscript", default="wolframscript")
    args = parser.parse_args()
    try:
        summary = run(args.matching_dir, args.config,
                      args.timeout_seconds, args.wolframscript)
    except (OSError, ValueError, RuntimeError, json.JSONDecodeError) as exc:
        print(f"[C5 benchmark] ERROR: {exc}", flush=True)
        return 1
    return 0 if summary["status"] in {
        "ProvisionalNumericalBenchmarkEvaluated", "ProvisionalFinitePartBenchmarkEvaluated"
    } else 1


if __name__ == "__main__":
    raise SystemExit(main())
