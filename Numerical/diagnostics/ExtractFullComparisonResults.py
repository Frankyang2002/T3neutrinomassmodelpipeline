"""Extract fixed-comparison low-energy T3 observables from existing JSON outputs.

Read-only with respect to pipeline results. No matching, RG evolution, or plots.
Run: python -m Numerical.diagnostics.ExtractFullComparisonResults
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import re
from collections import Counter
from pathlib import Path
from typing import Any

PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT = PROJECT_ROOT / "output" / "full" / "comparison"
DEFAULT_OUTPUT = PROJECT_ROOT / "Reports" / "output" / "full"
SCENARIOS = ("smallY_smallL", "smallY_largeL", "largeY_smallL", "largeY_largeL")
MODEL_PATTERN = re.compile(r"^T3_dS1_(\d+)_dS2_(\d+)_dF_(\d+)_alpha_([mp]\d+)$")
CLASS_PATTERN = re.compile(r"^T3_([ABCDE])_alpha_([mp]\d+)$")
COLUMNS = (
    "scenario", "model_class", "alpha", "d_s1", "d_s2", "d_f", "model_key",
    "ordering", "m1_ev", "m2_ev", "m3_ev", "delta_m21_sq_ev2",
    "delta_m31_sq_ev2", "delta_m32_sq_ev2", "sin2_theta12", "sin2_theta13",
    "sin2_theta23", "takagi_residual", "diagnostics_path",
)


def _alpha(value: str) -> int:
    return (-1 if value[0] == "m" else 1) * int(value[1:])


def _finite(value: Any, label: str) -> float:
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"{label} is non-finite")
    return result


def extract_case(path: Path, input_root: Path) -> dict[str, Any]:
    """Read and validate one diagnostic; raise ValueError on malformed input."""
    relative = path.relative_to(input_root)
    if len(relative.parts) != 5 or relative.parts[-1] != "running_diagnostics.json" or relative.parts[-2] != "data":
        raise ValueError(f"Unexpected comparison output layout: {relative}")
    scenario, model_key, class_key = relative.parts[:3]
    if scenario not in SCENARIOS:
        raise ValueError(f"Unknown comparison scenario: {scenario}")
    match = MODEL_PATTERN.fullmatch(model_key)
    class_match = CLASS_PATTERN.fullmatch(class_key)
    if not match or not class_match or match.group(4) != class_match.group(2):
        raise ValueError(f"Invalid or mismatched model keys: {model_key}/{class_key}")
    dims = tuple(map(int, match.groups()[:3]))
    cls = class_match.group(1)
    if {"A": (1, 3, 2), "B": (2, 2, 1), "C": (2, 2, 3), "D": (3, 1, 2), "E": (3, 3, 2)}[cls] != dims:
        raise ValueError("Class/representation mismatch")

    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload.get("status") != "Success":
        raise ValueError(f"Diagnostics status is not Success: {payload.get('status')!r}")
    low = payload["low_energy"]
    masses = low["masses_ev"]
    if len(masses) != 3:
        raise ValueError("Expected three neutrino masses")
    pmns = low["pmns_abs"]
    if len(pmns) != 3 or any(len(row) != 3 for row in pmns):
        raise ValueError("Expected 3x3 PMNS absolute-value matrix")
    ue3 = _finite(pmns[0][2], "Ue3")
    denominator = 1.0 - ue3**2
    if denominator <= 0.0:
        raise ValueError("PMNS angle denominator 1-|Ue3|^2 is non-positive")
    row = {
        "scenario": scenario,
        "model_class": cls,
        "alpha": _alpha(match.group(4)),
        "d_s1": dims[0], "d_s2": dims[1], "d_f": dims[2],
        "model_key": model_key,
        "ordering": payload["ordering"],
        "m1_ev": _finite(masses[0], "m1"),
        "m2_ev": _finite(masses[1], "m2"),
        "m3_ev": _finite(masses[2], "m3"),
        "delta_m21_sq_ev2": _finite(low["delta_m21_sq_ev2"], "dm21"),
        "delta_m31_sq_ev2": _finite(low["delta_m31_sq_ev2"], "dm31"),
        "delta_m32_sq_ev2": _finite(low["delta_m32_sq_ev2"], "dm32"),
        "sin2_theta12": _finite(pmns[0][1], "Ue2")**2 / denominator,
        "sin2_theta13": ue3**2,
        "sin2_theta23": _finite(pmns[1][2], "Umu3")**2 / denominator,
        "takagi_residual": _finite(low["takagi_residual"], "Takagi residual"),
        "diagnostics_path": relative.as_posix(),
    }
    if row["ordering"] not in {"NO", "IO"}:
        raise ValueError(f"Unknown neutrino ordering: {row['ordering']}")
    return row


def collect(input_root: Path) -> tuple[list[dict[str, Any]], list[dict[str, str]]]:
    rows: list[dict[str, Any]] = []
    errors: list[dict[str, str]] = []
    for path in sorted(input_root.glob("*/*/T3_*/data/running_diagnostics.json")):
        try:
            rows.append(extract_case(path, input_root))
        except (OSError, ValueError, KeyError, TypeError, IndexError, json.JSONDecodeError) as exc:
            errors.append({"path": str(path), "error": str(exc)})
    rows.sort(key=lambda row: (SCENARIOS.index(row["scenario"]), row["model_class"], row["alpha"]))
    return rows, errors


def write_results(input_root: Path, output_root: Path) -> dict[str, Any]:
    rows, errors = collect(input_root)
    output_root.mkdir(parents=True, exist_ok=True)
    csv_path = output_root / "t3_full_comparison_observables.csv"
    json_path = output_root / "t3_full_comparison_observables.json"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS)
        writer.writeheader()
        writer.writerows(rows)
    counts = Counter(row["scenario"] for row in rows)
    summary = {
        "status": "Complete" if len(rows) == 64 and not errors and all(counts[s] == 16 for s in SCENARIOS) else "Incomplete",
        "expected_cases": 64,
        "extracted_cases": len(rows),
        "per_scenario": {scenario: counts[scenario] for scenario in SCENARIOS},
        "errors": errors,
        "note": "Observables are predictions at fixed benchmarks, not experimental-fit claims; angles derived from PMNS magnitudes with the production RunningDiagnostics convention.",
        "rows": rows,
    }
    json_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"Extracted {len(rows)}/64 cases; status={summary['status']}; malformed={len(errors)}")
    print(f"CSV: {csv_path}\nJSON: {json_path}")
    return summary


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args(argv)
    if not args.input_root.is_dir():
        parser.error(f"Input directory does not exist: {args.input_root}")
    result = write_results(args.input_root, args.output_root)
    return 0 if result["status"] == "Complete" else 1


if __name__ == "__main__":
    raise SystemExit(main())
