"""Run the existing Matchete matching engine with the full-general interactions.

Does not create a new matching formula. RunModel.wl executes actual Match[].
B_-1 and C_-1 require a scalar mass-basis transformation and are marked blocked.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
RUNNER = ROOT / "Lagrangian" / "RunModel.wl"
MODELS = {
    "A": (-4, -2, 0),
    "B": (-3, -1, 1),
    "C": (-3, -1, 1),
    "D": (-2, 0, 2),
    "E": (-4, -2, 0, 2),
}
BLOCKED = {"B_-1", "C_-1"}
ALL = tuple(f"{family}_{alpha:+d}" for family, alphas in MODELS.items()
            for alpha in alphas)


def tag(model: str) -> str:
    return model.replace("+", "p").replace("-", "m")


def run_model(model: str, root: Path, executable: str, *, reuse: bool) -> dict:
    output = root / tag(model)
    output.mkdir(parents=True, exist_ok=True)
    log = output / "general_matching.log"
    result_file = output / "general_matching_status.json"
    if model in BLOCKED:
        result = {
            "model": model, "status": "BlockedMixedScalarMassBasis",
            "reason": "Source LagrangianBuilder.wl does not authorize Match[] "
                      "before off-diagonal scalar mass diagonalization",
        }
    elif reuse and result_file.exists():
        previous = json.loads(result_file.read_text(encoding="utf-8"))
        if previous.get("status") == "Matched":
            return previous
        result = None
    else:
        result = None

    if result is None:
        family, alpha = model.split("_")
        # wolframscript can consume a leading minus as an option; RunModel.wl
        # explicitly accepts m2/p2 notation for signed model parameters.
        alpha_cli = ("m" + alpha[1:]) if alpha.startswith("-") else ("p" + alpha[1:]) if alpha.startswith("+") else alpha
        candidate = output / "c5_coefficient.txt"
        if candidate.exists():
            candidate.unlink()  # Never report a stale coefficient as new matching.
        command = [executable, "-file", str(RUNNER), str(output),
                   "5", "1", family, alpha_cli]
        env = os.environ.copy()
        env["T3_EXPERIMENTAL_GENERAL_FERMIONS"] = "1"
        env["T3_EXPERIMENTAL_GENERAL_QUADRATICS"] = "1"
        env["T3_EXPERIMENTAL_GENERAL_QUARTICS"] = "1"
        try:
            process = subprocess.run(command, cwd=ROOT, env=env,
                                     capture_output=True, text=True, check=False)
            log.write_text(process.stdout + "\n" + process.stderr,
                           encoding="utf-8", errors="replace")
            expected_alpha = int(alpha)
            expected_model_line = f"Model: T3-{family}, alpha = {expected_alpha}"
            model_correct = expected_model_line in process.stdout
            result = {
                "model": model,
                "status": "Matched" if process.returncode == 0
                          and model_correct and candidate.is_file() and candidate.stat().st_size > 0
                          else "Failed",
                "exit_code": process.returncode,
                "model_selection_verified": model_correct,
                "wolfram_alpha_argument": alpha_cli,
                "c5_coefficient_file": str(candidate) if candidate.is_file() else None,
                "log": str(log),
            }
        except OSError as exc:
            result = {"model": model, "status": "Failed",
                      "error": str(exc), "log": str(log)}
    result_file.write_text(json.dumps(result, indent=2) + "\n",
                           encoding="utf-8")
    return result


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", choices=ALL)
    parser.add_argument("--output-dir", type=Path,
                        default=ROOT / "output" / "general_matchete_matching")
    parser.add_argument("--wolframscript", default="wolframscript")
    parser.add_argument("--reuse-success", action="store_true")
    parser.add_argument("--list", action="store_true")
    args = parser.parse_args(argv)
    if args.list:
        for model in ALL:
            print(model, "(mass-basis blocked)" if model in BLOCKED else "")
        return 0
    if not RUNNER.is_file():
        parser.error("Missing existing Matchete runner: " + str(RUNNER))
    models = (args.model,) if args.model else ALL
    args.output_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for index, model in enumerate(models, 1):
        row = run_model(model, args.output_dir, args.wolframscript,
                        reuse=args.reuse_success)
        rows.append(row)
        print(f"[{index}/{len(models)}] {model}: {row['status']}", flush=True)
        (args.output_dir / "summary.json").write_text(
            json.dumps({"models": rows,
                        "matched_count": sum(r["status"] == "Matched" for r in rows),
                        "blocked_count": sum(r["status"].startswith("Blocked") for r in rows),
                        "general_numerical_c5_implemented": False}, indent=2) + "\n",
            encoding="utf-8")
    return int(any(r["status"] == "Failed" for r in rows))


if __name__ == "__main__":
    sys.exit(main())
