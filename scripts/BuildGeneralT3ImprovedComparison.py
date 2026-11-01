"""Build CSV and Markdown comparison tables from existing improved 14-model neutrino JSONs.

Read-only with respect to all physics outputs; creates files only in --output-dir.
Run from repository root:
    python -X utf8 -m scripts.BuildGeneralT3ImprovedComparison
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

MODELS = (
    "A_-4", "A_-2", "A_+0", "B_-3", "B_+1", "C_-3", "C_+1",
    "D_-2", "D_+0", "D_+2", "E_-4", "E_-2", "E_+0", "E_+2",
)


def tag(model: str) -> str:
    return model.replace("-", "m").replace("+", "p")


def finite_number(value: object, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (float, int)):
        raise ValueError(f"{label}: expected numeric value")
    value = float(value)
    if not math.isfinite(value):
        raise ValueError(f"{label}: non-finite value")
    return value


def load_json(path: Path) -> dict:
    """Read standard JSON and the literal-backslash-n suffix from the prior batch writer."""
    content = path.read_text(encoding="utf-8-sig")
    try:
        result = json.loads(content)
    except json.JSONDecodeError:
        # The previous low-energy driver mistakenly appended the two characters
        # backslash and n, instead of an actual newline, to summary.json.
        if not content.endswith("\\n"):
            raise
        result = json.loads(content[:-2])
    if not isinstance(result, dict):
        raise ValueError(f"{path}: expected JSON object")
    return result


def extract(model: str, path: Path) -> dict:
    data = load_json(path)
    if data.get("model") != model:
        raise ValueError(f"{path}: model mismatch: {data.get('model')!r} != {model!r}")
    masses = data["neutrino_masses_ev"]
    pmns = data["pmns_abs"]
    if not isinstance(masses, list) or len(masses) != 3:
        raise ValueError(f"{path}: expected 3 neutrino masses")
    if not isinstance(pmns, list) or len(pmns) != 3 or any(
        not isinstance(row, list) or len(row) != 3 for row in pmns
    ):
        raise ValueError(f"{path}: expected 3x3 PMNS magnitudes")
    result = {"model": model, "low_energy_rge_status": str(data["low_energy_rge_status"]),
              "observables_scale_gev": finite_number(data["observables_scale_gev"], "scale")}
    for i, mass in enumerate(masses, 1):
        result[f"m{i}_ev"] = finite_number(mass, f"m{i}")
    for key in ("delta_m21_sq_ev2", "delta_m31_sq_ev2"):
        result[key] = finite_number(data[key], key)
    for i in range(3):
        for j in range(3):
            result[f"Uabs_{i+1}{j+1}"] = finite_number(pmns[i][j], f"PMNS[{i},{j}]")
    result["takagi_residual"] = finite_number(data["takagi_residual"], "Takagi residual")
    return result


def fmt(x: float) -> str:
    return f"{x:.6g}"


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input-dir", type=Path,
                   default=Path("output/general_numerical_all14_improved/neutrinos"))
    p.add_argument("--output-dir", type=Path,
                   default=Path("output/general_numerical_all14_improved/comparison"))
    args = p.parse_args(argv)

    summary_path = args.input_dir / "summary.json"
    if not summary_path.is_file():
        p.error(f"Missing summary: {summary_path}")
    summary = load_json(summary_path)
    status_by_model = {r["model"]: r for r in summary.get("models", [])}
    rows = []
    for model in MODELS:
        record = status_by_model.get(model)
        if record is None or record.get("status") != "Success":
            raise ValueError(f"Batch summary does not record success for {model}: {record}")
        path = args.input_dir / tag(model) / "neutrinos.json"
        if not path.is_file():
            raise FileNotFoundError(f"Missing neutrino result for {model}: {path}")
        row = extract(model, path)
        if row["low_energy_rge_status"] != "Success":
            raise ValueError(f"{model}: low-energy RGE failed or used threshold-only fallback")
        rows.append(row)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = args.output_dir / "improved_14_model_comparison.csv"
    md_path = args.output_dir / "improved_14_model_comparison.md"
    fields = list(rows[0])
    with csv_path.open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    cols = ("model", "m1_ev", "m2_ev", "m3_ev", "delta_m21_sq_ev2",
            "delta_m31_sq_ev2", "low_energy_rge_status")
    hdrs = ("Model", "m1 (eV)", "m2 (eV)", "m3 (eV)",
            "Delta m21^2 (eV^2)", "Delta m31^2 (eV^2)", "RGE status")
    lines = ["# Improved T3 14-model numerical comparison", "",
             "Exploratory direct matching, not a certified physical prediction.",
             "Input: existing per-model neutrino JSONs; no physics recalculated.", "",
             "| " + " | ".join(hdrs) + " |",
             "| " + " | ".join("---" for _ in cols) + " |"]
    for row in rows:
        lines.append("| " + " | ".join(
            str(row[col]) if col in ("model", "low_energy_rge_status") else fmt(row[col])
            for col in cols
        ) + " |")
    lines += ["", "## Absolute PMNS matrices", "",
              "Each entry is |U_ij|, not U_ij with its complex phase.", ""]
    for row in rows:
        lines += [f"### {row['model']}", "",
                  "| | 1 | 2 | 3 |", "|---|---:|---:|---:|"]
        for i in (1, 2, 3):
            lines.append("| " + str(i) + " | " + " | ".join(
                fmt(row[f"Uabs_{i}{j}"]) for j in (1, 2, 3)
            ) + " |")
        lines.append("")
    md_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Validated {len(rows)}/{len(MODELS)} successful model results")
    print(f"CSV: {csv_path}")
    print(f"Markdown: {md_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
