"""Audit stored T3 benchmark Weinberg continuations without rerunning physics.

Run: python -m Numerical.diagnostics.VerifyFullContinuation
Reads output/full/comparison/<scenario>/**/running_diagnostics.json.
Writes CSV and JSON audit to Reports/output/full/continuation_verification.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from Numerical.plotting.WeinbergMatchedContinuation import (
    reconstruct_matched_continuation,
)

PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT = PROJECT_ROOT / "output" / "full" / "comparison"
DEFAULT_OUTPUT = PROJECT_ROOT / "Reports" / "output" / "full" / "continuation_verification"
SCENARIOS = ("smallY_smallL", "smallY_largeL", "largeY_smallL", "largeY_largeL")


def _model(path: Path, payload: dict[str, Any]) -> str:
    search = payload.get("benchmark_search")
    model = search.get("model") if isinstance(search, dict) else None
    if isinstance(model, dict) and model.get("model_key"):
        return str(model["model_key"])
    return path.parent.parent.name


def audit(input_root: Path = DEFAULT_INPUT, output_root: Path = DEFAULT_OUTPUT) -> dict[str, Any]:
    """Check matching, direct-curve reconstruction and all 64 expected cases."""
    root, output = Path(input_root), Path(output_root)
    rows: list[dict[str, Any]] = []
    seen: set[tuple[str, str]] = set()
    for scenario in SCENARIOS:
        for path in sorted((root / scenario).rglob("running_diagnostics.json")):
            row: dict[str, Any] = {
                "scenario": scenario, "model": path.parent.parent.name,
                "status": "FAIL", "method": "", "points": "", "max_relative_residual": "",
                "hard_norm_gev_inv": "", "direct_norm_gev_inv": "",
                "combined_norm_gev_inv": "", "explanation": "", "source": str(path),
            }
            try:
                payload = json.loads(path.read_text(encoding="utf-8"))
                row["model"] = _model(path, payload)
                if payload.get("status") != "Success":
                    raise ValueError("Saved diagnostic has non-success status")
                key = (scenario, row["model"])
                if key in seen:
                    raise ValueError("Duplicate diagnostics for scenario/model")
                seen.add(key)
                result = reconstruct_matched_continuation(payload)
                row.update(
                    status="PASS", method=result.method, points=int(result.mu_gev.size),
                    max_relative_residual=float(result.maximum_relative_magnitude_residual),
                    hard_norm_gev_inv=float(np.linalg.norm(result.hard)),
                    direct_norm_gev_inv=float(np.linalg.norm(result.direct[-1])),
                    combined_norm_gev_inv=float(np.linalg.norm(result.combined[-1])),
                )
            except (KeyError, ValueError, TypeError, OSError, json.JSONDecodeError,
                    IndexError, FloatingPointError) as exc:
                row["explanation"] = f"{type(exc).__name__}: {exc}"
            rows.append(row)
    models = sorted({row["model"] for row in rows})
    pass_count = sum(row["status"] == "PASS" for row in rows)
    result = {
        "status": "PASS" if len(rows) == 64 and len(models) == 16 and pass_count == 64 else "INCOMPLETE_OR_FAILED",
        "expected_cases": 64, "observed_cases": len(rows),
        "unique_model_count": len(models), "pass_count": pass_count,
        "failed_count": len(rows) - pass_count,
        "cases_per_scenario": {s: sum(row["scenario"] == s for row in rows) for s in SCENARIOS},
        "rows": rows,
    }
    output.mkdir(parents=True, exist_ok=True)
    csv_path = output / "continuation_verification.csv"
    with csv_path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]) if rows else [
            "scenario", "model", "status", "method", "points", "max_relative_residual",
            "hard_norm_gev_inv", "direct_norm_gev_inv", "combined_norm_gev_inv", "explanation", "source",
        ])
        writer.writeheader()
        writer.writerows(rows)
    (output / "continuation_verification.json").write_text(
        json.dumps(result, indent=2, allow_nan=False), encoding="utf-8"
    )
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    result = audit(args.input, args.output)
    print(f"Observed: {result['observed_cases']}/64; models: {result['unique_model_count']}/16; "
          f"passed: {result['pass_count']}; failed: {result['failed_count']}")
    for scenario, count in result["cases_per_scenario"].items():
        print(f"  {scenario}: {count}/16")
    for row in result["rows"]:
        if row["status"] != "PASS":
            print(f"  FAIL {row['scenario']} / {row['model']}: {row['explanation']}")
    print(f"CSV: {args.output / 'continuation_verification.csv'}")
    print(f"JSON: {args.output / 'continuation_verification.json'}")
    return 0 if result["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
