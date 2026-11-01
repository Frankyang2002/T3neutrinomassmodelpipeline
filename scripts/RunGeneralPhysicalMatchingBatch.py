"""Run the existing experimental full-general physical matching for all 16 T3 models.

This is a batch orchestrator, not a new matching algorithm. It never silently
substitutes restricted-basis C5 for missing general matching. Mixed-bilinear
B_-1/C_-1 are expected unsupported until their mass basis is implemented.

PowerShell:
 python -X utf8 -m scripts.RunGeneralPhysicalMatchingBatch --all
 python -X utf8 -m scripts.RunGeneralPhysicalMatchingBatch --model D_-2
 python -X utf8 -m scripts.RunGeneralPhysicalMatchingBatch --all --run-running
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import re
import subprocess
import sys
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
REPRESENTATIONS = {
    "A": (-4, -2, 0),
    "B": (-3, -1, 1),
    "C": (-3, -1, 1),
    "D": (-2, 0, 2),
    "E": (-4, -2, 0, 2),
}
MODELS = [f"{name}_{alpha:+d}" for name, alphas in REPRESENTATIONS.items() for alpha in alphas]
# This manifest is authoritative only for known obstacles, not for model success.
MIXED_SCALAR_MASS_MODELS = {"B_-1", "C_-1"}


def wolfram_alpha_token(alpha: int) -> str:
    """Match RunGeneralPhysicalMatching.wl CLI: m2, p0, p2, etc."""
    return f"m{abs(alpha)}" if alpha < 0 else f"p{alpha}"


def model_tag(model: str) -> str:
    return model.replace("+", "p").replace("-", "m")


def atomically_write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    tmp.replace(path)


def run_command(command: list[str], log_path: Path, env: dict[str, str]) -> int:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w", encoding="utf-8", errors="replace") as log:
        log.write("Command: " + subprocess.list2cmdline(command) + "\n\n")
        log.flush()
        try:
            result = subprocess.run(command, cwd=ROOT, env=env, stdout=log,
                                    stderr=subprocess.STDOUT, check=False)
            return result.returncode
        except OSError as exc:
            log.write("Execution error: " + str(exc) + "\n")
            return 127


def classify_matching(exit_code: int, coefficient: Path) -> str:
    if exit_code == 0 and coefficient.is_file() and coefficient.stat().st_size > 0:
        return "GeneralSymbolicC5Exported"
    if exit_code == 18:
        return "UnsupportedScalarMassDiagonalization"
    if exit_code == 14:
        return "MatchingSucceeded_C5ExtractionUnavailable"
    return "Failed" if exit_code else "MissingC5DespiteZeroExit"


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    models_group = p.add_mutually_exclusive_group(required=True)
    models_group.add_argument("--all", action="store_true", help="Run the full 16-model catalog")
    models_group.add_argument("--model", choices=MODELS, help="Run one model")
    p.add_argument("--output-root", type=Path, default=Path("output/general_physical_matching_batch"))
    p.add_argument("--wolframscript", default="wolframscript")
    p.add_argument("--skip-existing", action="store_true")
    p.add_argument("--run-running", action="store_true", help="Also run canonical UV/EFT1 numerical RGEs")
    p.add_argument("--canonical-run-id", default="full16_general_v1",
                   help="Cached general symbolic RGE run ID; generate cache separately if needed")
    p.add_argument("--mu-uv", type=float, default=1e7)
    p.add_argument("--mu-f", type=float, default=1e5)
    p.add_argument("--mu-s", type=float, default=1e3)
    p.add_argument("--initial-json", type=Path,
                   help="Optional shared initial couplings; unsupported overrides are rejected per model")
    p.add_argument("--skip-known-blocked", action="store_true",
                   help="Do not call Matchete for known B_-1/C_-1 mixed-bilinear cases")
    args = p.parse_args(argv)
    selected = MODELS if args.all else [args.model]
    root = args.output_root.resolve()
    if not (ROOT / "Lagrangian/RunGeneralPhysicalMatching.wl").is_file():
        p.error("Missing Lagrangian/RunGeneralPhysicalMatching.wl; install working A(-2) matching driver first")
    if args.run_running and not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]{0,63}", args.canonical_run_id):
        p.error("--canonical-run-id must be a simple directory name")
    env = os.environ.copy()
    env.update({
        "T3_EXPERIMENTAL_GENERAL_FERMIONS": "1",
        "T3_EXPERIMENTAL_GENERAL_QUADRATICS": "1",
        "T3_EXPERIMENTAL_GENERAL_QUARTICS": "1",
    })
    rows: list[dict[str, Any]] = []
    root.mkdir(parents=True, exist_ok=True)
    for i, model in enumerate(selected, 1):
        cls, alpha_text = model.split("_")
        alpha = int(alpha_text)
        tag = model_tag(model)
        dest = root / tag
        coeff = dest / "c5_coefficient.txt"
        row: dict[str, Any] = {"model": model,
            "known_mixed_mass_basis_obstacle": model in MIXED_SCALAR_MASS_MODELS,
            "matching_output_directory": str(dest),
            "coefficient_path": str(coeff),
            "matching_approach": "experimental_general_mass_basis_driver",
            "matching_physical_normalization_certified": False,
            "higher_dimensional_eft1_wilson_running": False}
        if model in MIXED_SCALAR_MASS_MODELS and args.skip_known_blocked:
            row.update(matching_status="UnsupportedScalarMassDiagonalization", matching_exit_code=18)
        elif args.skip_existing and coeff.is_file() and coeff.stat().st_size > 0:
            row.update(matching_status="ExistingNonemptyC5_UnverifiedCache", matching_exit_code=None)
        else:
            # Clear stale coefficient, otherwise a failed subprocess could appear successful.
            coeff.unlink(missing_ok=True)
            command = [args.wolframscript, "-file", "Lagrangian/RunGeneralPhysicalMatching.wl",
                       str(dest), cls, wolfram_alpha_token(alpha)]
            returncode = run_command(command, dest / "matching.log", env)
            row.update(matching_status=classify_matching(returncode, coeff),
                       matching_exit_code=returncode,
                       matching_log=str(dest / "matching.log"))
            if row["matching_status"] != "GeneralSymbolicC5Exported":
                log_text = (dest / "matching.log").read_text(encoding="utf-8", errors="replace")
                row["matching_log_tail"] = "\n".join(log_text.splitlines()[-24:])
        if args.run_running:
            running_path = dest / "general_running.json"
            command = [sys.executable, "-X", "utf8", "-m",
                       "Numerical.orchestration.GeneralRunningBenchmark",
                       "--canonical-run-id", args.canonical_run_id,
                       "--model", model, "--mu-uv", str(args.mu_uv),
                       "--mu-f", str(args.mu_f), "--mu-s", str(args.mu_s),
                       "--output-json", str(running_path)]
            if args.initial_json:
                command += ["--initial-json", str(args.initial_json.resolve())]
            rc = run_command(command, dest / "running.log", env)
            row.update(running_exit_code=rc,
                       running_status="Success" if rc == 0 and running_path.is_file() else "Failed",
                       running_path=str(running_path), running_log=str(dest / "running.log"))
        rows.append(row)
        atomically_write_json(root / "batch_summary.json", {
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "requested": len(selected), "processed": len(rows),
            "full_general_symbolic_c5_exports": sum(r["matching_status"] == "GeneralSymbolicC5Exported" for r in rows),
            "models": rows,
        })
        print(f"[{i}/{len(selected)}] {model}: {row['matching_status']}" +
              (f"; running={row['running_status']}" if args.run_running else ""), flush=True)
    successes = {"GeneralSymbolicC5Exported", "ExistingNonemptyC5_UnverifiedCache"}
    print("Summary:", root / "batch_summary.json")
    return 0 if all(row["matching_status"] in successes for row in rows) else 1


if __name__ == "__main__":
    raise SystemExit(main())
