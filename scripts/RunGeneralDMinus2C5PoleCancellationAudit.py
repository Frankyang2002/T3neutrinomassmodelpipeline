"""Read-only exact dimensional-pole audit for the saved D(-2) C5 candidate.

From repository root:
    python -X utf8 -m scripts.RunGeneralDMinus2C5PoleCancellationAudit
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
SCRIPT = ROOT / "Lagrangian/RunGeneralDMinus2C5PoleCancellationAudit.wl"
DEFAULT_OUTPUT = ROOT / "output/general_sequential_pilot/D_m2_stage1"
PROFILE_REL = Path("data/general_d_m2_scalar_threshold/profile_tree")
RESULT = "profile_c5_pole_cancellation_result.json"
LOG = "profile_c5_pole_cancellation.log"
SUMMARY = "profile_c5_pole_cancellation_summary.json"


def read_result(path: Path) -> tuple[dict[str, Any] | None, str | None]:
    if not path.exists():
        return None, "Missing Wolfram result JSON"
    try:
        raw = path.read_text(encoding="utf-8-sig")
        if not raw.strip():
            return None, "Empty Wolfram result JSON"
        obj = json.loads(raw)
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        return None, f"Invalid Wolfram result JSON: {exc}"
    if not isinstance(obj, dict):
        return None, "Wolfram result must be a JSON object"
    return obj, None


def tail(path: Path, lines: int = 30) -> str:
    if not path.exists():
        return "(Wolfram log missing)"
    with path.open(encoding="utf-8", errors="replace") as f:
        return "".join(deque(f, maxlen=lines)).rstrip()


def terminate(proc: subprocess.Popen[bytes]) -> None:
    if proc.poll() is not None:
        return
    if os.name == "nt":
        try:
            subprocess.run(["taskkill", "/PID", str(proc.pid), "/T", "/F"],
                           check=False, capture_output=True, timeout=15)
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


def classify(code: int | None, timed_out: bool,
             result: dict[str, Any] | None) -> str:
    if timed_out:
        return "ExternalTimeout"
    if code != 0:
        return "WolframFailed"
    if not isinstance(result, dict):
        return "InvalidResultReport"
    status = result.get("Status")
    if status not in {"ExactNF1PoleCancellationVerified",
                      "ExactNF1PoleCancellationNotEstablished"}:
        return "InvalidResultStatus"
    if (result.get("NoNewMatchPerformed") is not True
            or result.get("RenormalizationCertified") is not False
            or result.get("NF") != 1):
        return "InconsistentMetadata"
    return status


def run(output: Path = DEFAULT_OUTPUT, timeout: int = 360,
        wolframscript: str = "wolframscript") -> dict[str, Any]:
    output = output.resolve()
    profile = output / PROFILE_REL
    required = [profile / "profile_c5_flavour_candidate_total_provisional.wxf",
                profile / "profile_c5_flavour_candidate_result.json"]
    missing = [str(p) for p in required if not p.is_file()]
    if missing:
        raise FileNotFoundError("Missing saved candidate inputs: " + ", ".join(missing))
    if not SCRIPT.is_file():
        raise FileNotFoundError(f"Missing Wolfram source: {SCRIPT}")
    if not 180 <= timeout <= 3600:
        raise ValueError("--timeout-seconds must be 180..3600")
    (profile / RESULT).unlink(missing_ok=True)
    cmd = [wolframscript, "-file", str(SCRIPT), str(output), "1"]
    env = os.environ.copy()
    env.update(T3_EXPERIMENTAL_GENERAL_QUADRATICS="1",
               T3_EXPERIMENTAL_GENERAL_FERMIONS="1")
    print(f"[C5 pole] Exact NF=1 residue audit; no Match; external limit {timeout}s",
          flush=True)
    started = time.monotonic()
    timed_out = False
    with (profile / LOG).open("wb") as log:
        flags = getattr(subprocess, "CREATE_NEW_PROCESS_GROUP", 0) if os.name == "nt" else 0
        proc = subprocess.Popen(cmd, cwd=ROOT, env=env, stdout=log,
                                stderr=subprocess.STDOUT, creationflags=flags,
                                start_new_session=(os.name != "nt"))
        heartbeat = started
        while proc.poll() is None:
            now = time.monotonic()
            if now - started > timeout:
                timed_out = True
                terminate(proc)
                break
            if now - heartbeat >= 20:
                print(f"[C5 pole] Active for {int(now-started)}s", flush=True)
                heartbeat = now
            time.sleep(0.5)
        code = proc.poll()
    result, error = read_result(profile / RESULT)
    status = classify(code, timed_out, result)
    summary = {"status": status, "exit_code": code,
               "elapsed_seconds": round(time.monotonic()-started, 2),
               "result": str(profile / RESULT), "log": str(profile / LOG),
               "json_error": error, "reason": (result or {}).get("Reason"),
               "no_new_matching": True}
    (profile / SUMMARY).write_text(json.dumps(summary, indent=2)+"\n", encoding="utf-8")
    print(f"[C5 pole] {status} (exit {code})\n"
          f"  Result: {profile / RESULT}\n  Log: {profile / LOG}", flush=True)
    if error:
        print(f"  JSON: {error}", flush=True)
    if status not in {"ExactNF1PoleCancellationVerified",
                      "ExactNF1PoleCancellationNotEstablished"}:
        print(tail(profile / LOG), flush=True)
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matching-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--timeout-seconds", type=int, default=360)
    parser.add_argument("--wolframscript", default="wolframscript")
    args = parser.parse_args()
    try:
        summary = run(args.matching_dir, args.timeout_seconds, args.wolframscript)
    except (OSError, ValueError) as exc:
        print(f"[C5 pole] ERROR: {exc}", flush=True)
        return 1
    return int(summary["status"] not in {
        "ExactNF1PoleCancellationVerified",
        "ExactNF1PoleCancellationNotEstablished"})


if __name__ == "__main__":
    raise SystemExit(main())
