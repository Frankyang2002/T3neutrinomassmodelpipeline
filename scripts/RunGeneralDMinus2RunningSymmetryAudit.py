"""Read-only direct-running p<->q symmetry diagnostic for general D(-2).

From the repository root:
  python -X utf8 -m scripts.RunGeneralDMinus2RunningSymmetryAudit
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
SCRIPT = ROOT / "Lagrangian/RunGeneralDMinus2RunningSymmetryAudit.wl"
DEFAULT_MATCHING = ROOT / "output/general_sequential_pilot/D_m2_stage1"
PROFILE_REL = Path("data/general_d_m2_scalar_threshold/profile_tree")
REPORT = "profile_c5_running_symmetry_audit_result.json"
SUMMARY = "profile_c5_running_symmetry_audit_summary.json"
LOG = "profile_c5_running_symmetry_audit.log"


def read_report(path: Path) -> tuple[dict[str, Any] | None, str | None]:
    if not path.is_file():
        return None, "Wolfram result JSON does not exist"
    try:
        raw = path.read_text(encoding="utf-8-sig")
        if not raw.strip():
            return None, "Wolfram result JSON is empty"
        result = json.loads(raw)
        if not isinstance(result, dict):
            return None, "Wolfram result JSON must be an object"
        return result, None
    except (UnicodeError, OSError, json.JSONDecodeError) as exc:
        return None, f"Wolfram result JSON invalid: {exc}"


def log_tail(path: Path, lines: int = 35) -> str:
    if not path.is_file():
        return "(no Wolfram log)"
    from collections import deque
    with path.open(encoding="utf-8", errors="replace") as stream:
        return "".join(deque(stream, maxlen=lines)).strip()


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
        p.kill()
        p.wait()


def run(matching_dir: Path = DEFAULT_MATCHING,
        timeout_seconds: int = 180, wolframscript: str = "wolframscript") -> dict[str, Any]:
    matching_dir = matching_dir.resolve()
    if timeout_seconds < 30 or timeout_seconds > 3600:
        raise ValueError("--timeout-seconds must be within 30..3600")
    ready = preflight(matching_dir, script=SCRIPT)
    profile = matching_dir / PROFILE_REL
    profile.mkdir(parents=True, exist_ok=True)
    (profile / REPORT).unlink(missing_ok=True)
    cmd = [wolframscript, "-file", str(SCRIPT.resolve()),
           str(matching_dir), str(Path(ready["insertion_file"]).resolve())]
    env = os.environ.copy()
    env.update(T3_EXPERIMENTAL_GENERAL_QUADRATICS="1",
               T3_EXPERIMENTAL_GENERAL_FERMIONS="1")
    print(f"[running symmetry] Direct-running definition only; no Match; external limit {timeout_seconds}s", flush=True)
    start = time.monotonic()
    with (profile / LOG).open("wb") as stream:
        p = subprocess.Popen(
            cmd, cwd=ROOT, env=env, stdout=stream, stderr=subprocess.STDOUT,
            creationflags=getattr(subprocess, "CREATE_NEW_PROCESS_GROUP", 0) if os.name == "nt" else 0,
            start_new_session=(os.name != "nt"))
        timed_out = False
        last_heartbeat = start
        while p.poll() is None:
            now = time.monotonic()
            if now - start >= timeout_seconds:
                timed_out = True
                kill_process_tree(p)
                break
            if now - last_heartbeat >= 20:
                print(f"[running symmetry] Active for {int(now-start)}s", flush=True)
                last_heartbeat = now
            time.sleep(0.5)
        code = p.poll()
    report, report_error = read_report(profile / REPORT)
    status = (
        "ExternalTimeout" if timed_out else
        "WolframFailed" if code != 0 else
        "InvalidWolframResultJSON" if report_error else
        "RunningSymmetryAudited" if report.get("Status") == "RunningSymmetryAudited" else
        "UnrecognizedWolframStatus"
    )
    summary = {
        "status": status,
        "exit_code": code,
        "elapsed_seconds": round(time.monotonic() - start, 2),
        "no_matching": True,
        "report": str(profile / REPORT),
        "log": str(profile / LOG),
        "report_error": report_error,
    }
    (profile / SUMMARY).write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"[running symmetry] {status} (exit {code})\n  Result: {profile / REPORT}\n  Log: {profile / LOG}", flush=True)
    if status != "RunningSymmetryAudited":
        if report_error:
            print(f"  JSON: {report_error}", flush=True)
        print(log_tail(profile / LOG), flush=True)
    return summary


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--matching-dir", type=Path, default=DEFAULT_MATCHING)
    ap.add_argument("--timeout-seconds", type=int, default=180)
    ap.add_argument("--wolframscript", default="wolframscript")
    args = ap.parse_args()
    try:
        result = run(args.matching_dir, args.timeout_seconds, args.wolframscript)
    except (OSError, ValueError, RuntimeError) as exc:
        print(f"[running symmetry] ERROR: {exc}", flush=True)
        return 1
    return 0 if result["status"] == "RunningSymmetryAudited" else 1


if __name__ == "__main__":
    raise SystemExit(main())
