"""Build a compact supervisor-facing copy of the comparison-only T3 study.

Presentation-only: no physics is recomputed and no source results are modified.

Inputs:
    output/full/comparison/
    Reports/output/full/comparison/

Default output:
    Reports/output/display/

The display directory is recreated on each run.  Legacy ``optimal`` output and
optimizer benchmark caches are deliberately ignored.
"""
from __future__ import annotations

import argparse
import csv
import json
import shutil
from pathlib import Path
from typing import Any

PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_DATA_INPUT = PROJECT_ROOT / "output" / "full"
DEFAULT_REPORT_INPUT = PROJECT_ROOT / "Reports" / "output" / "full"
DEFAULT_OUTPUT = PROJECT_ROOT / "Reports" / "output" / "display"

KEEP_FIGURES = (
    "uv_coupling_running.png",
    "intermediate_direct_weinberg_running.png",
    "c5_threshold_contributions.png",
    "c5_running.png",
    "neutrino_mass_splitting_running.png",
    "neutrino_mixing_running.png",
)

SCENARIOS = {
    "smallY_smallL": (0.005, 0.1),
    "smallY_largeL": (0.005, 1.0),
    "largeY_smallL": (0.5, 0.1),
    "largeY_largeL": (0.5, 1.0),
}
SCENARIO_ORDER = tuple(SCENARIOS)


def _load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"{path} must contain a JSON object.")
    return value


def _model_name(path: Path, payload: dict[str, Any]) -> str:
    model = payload.get("benchmark_search", {}).get("model")
    if isinstance(model, dict) and model.get("model_key"):
        return str(model["model_key"])
    for part in reversed(path.parts):
        if part.startswith("T3_dS1_"):
            return part
    return path.parent.parent.name


def _last(value: Any) -> Any:
    return value[-1] if isinstance(value, list) and value else ""


def _summary_row(
    scenario: str,
    path: Path,
    payload: dict[str, Any],
) -> dict[str, Any]:
    low = payload.get("low_energy", {})
    final = payload.get("final_running", {})
    masses = low.get("masses_ev", [])
    model = payload.get("benchmark_search", {}).get("model", {})
    if not isinstance(model, dict):
        model = {}

    return {
        "scenario": scenario,
        "model": _model_name(path, payload),
        "dS1": model.get("d_s1", ""),
        "dS2": model.get("d_s2", ""),
        "dF": model.get("d_f", ""),
        "alpha": model.get("alpha", ""),
        "m1_eV": masses[0] if len(masses) > 0 else "",
        "m2_eV": masses[1] if len(masses) > 1 else "",
        "m3_eV": masses[2] if len(masses) > 2 else "",
        "delta_m21_sq_eV2": low.get(
            "delta_m21_sq_ev2", _last(final.get("delta_m21_sq_ev2"))
        ),
        "delta_m3l_sq_eV2": low.get(
            "delta_m31_sq_ev2", _last(final.get("delta_m3l_sq_ev2"))
        ),
        "sin2_theta12": _last(final.get("sin2_theta12")),
        "sin2_theta13": _last(final.get("sin2_theta13")),
        "sin2_theta23": _last(final.get("sin2_theta23")),
        "takagi_residual": low.get("takagi_residual", ""),
    }


def _copy_dashboard(source: Path, destination: Path) -> bool:
    dashboard = source / "interactive_comparison.html"
    if not dashboard.is_file():
        return False
    destination.mkdir(parents=True, exist_ok=True)
    shutil.copy2(dashboard, destination / dashboard.name)
    return True


def _report_model_dirs(root: Path) -> list[Path]:
    if not root.is_dir():
        return []
    return sorted(
        path for path in root.iterdir()
        if path.is_dir() and path.name.startswith("T3_dS1_")
    )


def _find_figures_dir(model_dir: Path) -> Path | None:
    direct = model_dir / "figures"
    if direct.is_dir():
        return direct
    candidates = sorted(model_dir.glob("T3_*/figures"))
    return candidates[0] if candidates else None


def _copy_selected_figures(model_dir: Path, destination: Path) -> int:
    figures = _find_figures_dir(model_dir)
    if figures is None:
        return 0
    copied = 0
    for name in KEEP_FIGURES:
        source = figures / name
        if source.is_file():
            destination.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, destination / name)
            copied += 1
    return copied


def _copy_report_section(source: Path, destination: Path) -> tuple[int, bool]:
    dashboard = _copy_dashboard(source, destination)
    count = 0
    for model_dir in _report_model_dirs(source):
        count += _copy_selected_figures(
            model_dir,
            destination / "models" / model_dir.name,
        )
    return count, dashboard


def _collect_rows(data_root: Path) -> list[dict[str, Any]]:
    """Collect only the four comparison studies.

    Any legacy ``output/full/optimal`` tree is intentionally ignored so stale
    optimizer output cannot enter a newly generated standard display.
    """
    rows: list[dict[str, Any]] = []
    comparison = data_root / "comparison"
    if not comparison.is_dir():
        return rows

    for scenario_dir in sorted(comparison.iterdir()):
        if not scenario_dir.is_dir():
            continue
        for path in sorted(scenario_dir.rglob("running_diagnostics.json")):
            try:
                payload = _load_json(path)
            except (OSError, json.JSONDecodeError, ValueError):
                continue
            if payload.get("status") == "Success":
                rows.append(_summary_row(scenario_dir.name, path, payload))

    rows.sort(key=lambda row: (row["scenario"], row["model"]))
    return rows


def _write_csv(rows: list[dict[str, Any]], path: Path) -> None:
    if not rows:
        return
    fields = tuple(rows[0])
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def _write_readme(
    path: Path,
    *,
    rows: list[dict[str, Any]],
    report_model_count: int,
) -> None:
    lines = [
        "# T3 numerical study — comparison display",
        "",
        "This is a compact presentation copy of the standard comparison-only T3 study.",
        "No numerical calculation has been rerun and the full output is unchanged.",
        "Legacy optimal/optimizer outputs are not read or copied.",
        "",
        "## Included",
        "",
        "- Interactive cross-model comparison dashboards.",
        "- Selected figures for each model, including:",
        "  - UV Yukawa / lambda_T3 running",
        "  - one-loop operator-running overview across the EFT thresholds",
        "  - direct intermediate LLSS -> Weinberg generation",
        "  - hard/direct/combined C5 at the scalar threshold",
        "  - final Weinberg coefficient C5 running",
        "  - neutrino mass-splitting and mixing-angle running",
    ]
    if rows:
        lines.append(
            "- `model_summary.csv` with low-energy observables from comparison diagnostics."
        )

    lines.extend([
        "",
        "## Operator-running convention",
        "",
        "The intermediate curve is the implemented direct O(hbar) LLSS -> Weinberg",
        "contribution Delta C5^direct(mu). The final curve is C5(mu) in the",
        "SM+Weinberg EFT. A numerical LLSS self-running trajectory is not inferred",
        "or fabricated by the plotting layer. Feeding one-loop LLSS self-running",
        "through the scalar loop would first modify C5 at O(hbar^2).",
        "",
        "## Common comparison benchmarks",
        "",
    ])
    for name in SCENARIO_ORDER:
        y, lam = SCENARIOS[name]
        lines.append(f"- `{name}`: Y = {y:g}, lambda_T3 = {lam:g}")

    lines.extend([
        "",
        "The comparison mode uses the same fixed non-diagonal Yukawa textures for",
        "every model, scaled by the common Y value. These are controlled comparison",
        "benchmarks rather than oscillation-data fits.",
        "",
        f"Report model directories copied: {report_model_count}.",
    ])
    if rows:
        lines.append(f"Successful numerical rows summarized: {len(rows)}.")
    else:
        lines.append(
            "No comparison numerical JSON directory was supplied/found, so "
            "model_summary.csv was not generated."
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_display(
    data_input: Path,
    report_input: Path,
    output: Path,
) -> Path:
    data_input = Path(data_input).resolve()
    report_input = Path(report_input).resolve()
    output = Path(output).resolve()

    if not report_input.is_dir():
        raise FileNotFoundError(f"Report directory not found: {report_input}")
    if output in (data_input, report_input):
        raise ValueError("Display output must differ from source directories.")

    if output.exists():
        shutil.rmtree(output)
    output.mkdir(parents=True, exist_ok=True)

    report_model_count = 0
    comparison = report_input / "comparison"
    if comparison.is_dir():
        scenario_dirs = sorted(
            (path for path in comparison.iterdir() if path.is_dir()),
            key=lambda path: (
                SCENARIO_ORDER.index(path.name)
                if path.name in SCENARIO_ORDER else len(SCENARIO_ORDER),
                path.name,
            ),
        )
        for source in scenario_dirs:
            destination = output / "comparison" / source.name
            _copy_report_section(source, destination)
            report_model_count += len(_report_model_dirs(source))

    rows = _collect_rows(data_input)
    _write_csv(rows, output / "model_summary.csv")
    _write_readme(
        output / "README.md",
        rows=rows,
        report_model_count=report_model_count,
    )

    print(f"Display results written to: {output}")
    print(f"Comparison report model directories: {report_model_count}")
    print(f"Comparison numerical summary rows: {len(rows)}")
    return output


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--data-input", type=Path, default=DEFAULT_DATA_INPUT,
        help="numerical data root (default: output/full)",
    )
    parser.add_argument(
        "--report-input", type=Path, default=DEFAULT_REPORT_INPUT,
        help="report root (default: Reports/output/full)",
    )
    parser.add_argument(
        "--output", type=Path, default=DEFAULT_OUTPUT,
        help="display root to recreate (default: Reports/output/display)",
    )
    args = parser.parse_args()
    build_display(args.data_input, args.report_input, args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
