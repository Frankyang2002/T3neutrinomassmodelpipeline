"""Build a compact supervisor-facing copy of the full T3 study.

Presentation-only: no physics is recomputed and no source results are modified.

Inputs used by the current project layout:
    output/full/          numerical JSON data
    Reports/output/full/  generated reports, figures and dashboards

Default output:
    Reports/output/display/

The display directory is recreated on each run.
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
    "c5_running.png",
    "neutrino_mass_splitting_running.png",
    "neutrino_mixing_running.png",
)

SCENARIOS = {
    "smallY_smallL": (0.005, 0.01),
    "smallY_largeL": (0.005, 0.25),
    "largeY_smallL": (0.5, 0.01),
    "largeY_largeL": (0.5, 0.25),
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
    mode: str,
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
        "mode": mode,
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
    return sorted(
        path
        for path in root.iterdir()
        if path.is_dir() and path.name.startswith("T3_dS1_")
    ) if root.is_dir() else []


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
    rows: list[dict[str, Any]] = []
    if not data_root.is_dir():
        return rows

    comparison = data_root / "comparison"
    if comparison.is_dir():
        for scenario_dir in sorted(comparison.iterdir()):
            if not scenario_dir.is_dir():
                continue
            for path in sorted(scenario_dir.rglob("running_diagnostics.json")):
                try:
                    payload = _load_json(path)
                except (OSError, json.JSONDecodeError, ValueError):
                    continue
                if payload.get("status") == "Success":
                    rows.append(
                        _summary_row("comparison", scenario_dir.name, path, payload)
                    )

    optimal = data_root / "optimal"
    if optimal.is_dir():
        for path in sorted(optimal.rglob("running_diagnostics.json")):
            try:
                payload = _load_json(path)
            except (OSError, json.JSONDecodeError, ValueError):
                continue
            if payload.get("status") == "Success":
                rows.append(_summary_row("optimal", "", path, payload))

    rows.sort(key=lambda row: (row["mode"], row["scenario"], row["model"]))
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
        "# T3 numerical study — supervisor display",
        "",
        "This is a compact presentation copy of the full T3 study.",
        "No numerical calculation has been rerun and the full output is unchanged.",
        "",
        "## Included",
        "",
        "- Interactive cross-model comparison dashboards.",
        "- Four selected figures for each model:",
        "  - UV Yukawa / lambda_T3 running",
        "  - Weinberg coefficient C5 running",
        "  - neutrino mass-splitting running",
        "  - neutrino mixing-angle running",
    ]
    if rows:
        lines.append("- `model_summary.csv` with low-energy observables from the stored diagnostics.")

    lines.extend(
        [
            "",
            "## Omitted from this display copy",
            "",
            "Group-factor derivations, full Lagrangian/RGE PDFs and TeX files,",
            "LaTeX auxiliary/log files, and secondary debugging/intermediate plots.",
            "These remain in the complete project output for reproducibility.",
            "",
            "## Common comparison benchmarks",
            "",
        ]
    )
    for name in SCENARIO_ORDER:
        y, lam = SCENARIOS[name]
        lines.append(f"- `{name}`: Y = {y:g}, lambda_T3 = {lam:g}")

    lines.extend(
        [
            "",
            "The comparison mode uses the same normalized non-diagonal Yukawa",
            "textures for every model, scaled by the common Y value. These are",
            "controlled comparison benchmarks rather than oscillation-data fits.",
            "",
            f"Report model directories copied: {report_model_count}.",
        ]
    )
    if rows:
        lines.append(f"Successful numerical rows summarized: {len(rows)}.")
    else:
        lines.append(
            "No numerical JSON directory was supplied/found, so model_summary.csv "
            "was not generated."
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _benchmark_model_key(path: Path, payload: dict[str, Any]) -> str:
    compatibility = payload.get("compatibility", {})
    model = compatibility.get("model", {}) if isinstance(compatibility, dict) else {}
    if isinstance(model, dict):
        try:
            ds1 = int(model["d_s1"])
            ds2 = int(model["d_s2"])
            df = int(model["d_f"])
            alpha = int(model["alpha"])
            a = f"p{alpha}" if alpha >= 0 else f"m{abs(alpha)}"
            return f"T3_dS1_{ds1}_dS2_{ds2}_dF_{df}_alpha_{a}"
        except (KeyError, TypeError, ValueError):
            pass
    return path.parent.name


def _benchmark_parameter_rows(benchmark_root: Path) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Read the authoritative cached optimal benchmarks.

    The search varies Re(lambdaT3) plus the 18 real y1/y2 entries.  Complete
    matrices are retained in JSON; the CSV uses compact Frobenius/max norms.
    """
    rows: list[dict[str, Any]] = []
    full: dict[str, Any] = {}

    if not benchmark_root.is_dir():
        return rows, full

    for path in sorted(benchmark_root.rglob("sobol_benchmark.json")):
        try:
            payload = _load_json(path)
        except (OSError, json.JSONDecodeError, ValueError):
            continue
        best = payload.get("best", {})
        params = best.get("parameters", {}) if isinstance(best, dict) else {}
        if not isinstance(params, dict):
            continue

        def matrix(prefix: str) -> list[list[float]]:
            return [
                [float(params.get(f"{prefix}_{i}{j}", 0.0)) for j in range(1, 4)]
                for i in range(1, 4)
            ]

        y1 = matrix("y1")
        y2 = matrix("y2")

        def frob(m: list[list[float]]) -> float:
            return sum(x*x for row in m for x in row) ** 0.5

        def maxabs(m: list[list[float]]) -> float:
            return max(abs(x) for row in m for x in row)

        key = _benchmark_model_key(path, payload)
        compatibility = payload.get("compatibility", {})
        model = compatibility.get("model", {}) if isinstance(compatibility, dict) else {}
        prediction = best.get("prediction", {}) if isinstance(best, dict) else {}
        if not isinstance(prediction, dict):
            prediction = {}

        row = {
            "model": key,
            "dS1": model.get("d_s1", "") if isinstance(model, dict) else "",
            "dS2": model.get("d_s2", "") if isinstance(model, dict) else "",
            "dF": model.get("d_f", "") if isinstance(model, dict) else "",
            "alpha": model.get("alpha", "") if isinstance(model, dict) else "",
            "chi2": best.get("chi2", ""),
            "benchmark_status": payload.get("status", ""),
            "lambdaT3_real": params.get("lambdaT3_real", ""),
            "y1_frobenius": frob(y1),
            "y1_max_abs": maxabs(y1),
            "y2_frobenius": frob(y2),
            "y2_max_abs": maxabs(y2),
            "delta_m21_sq_eV2": prediction.get("delta_m21_sq_ev2", ""),
            "delta_m3l_sq_eV2": prediction.get("delta_m3l_sq_ev2", ""),
            "sin2_theta12": prediction.get("sin2_theta12", ""),
            "sin2_theta13": prediction.get("sin2_theta13", ""),
            "sin2_theta23": prediction.get("sin2_theta23", ""),
        }
        rows.append(row)
        full[key] = {
            "chi2": best.get("chi2"),
            "benchmark_status": payload.get("status"),
            "method": payload.get("method"),
            "parameters": {
                "lambdaT3_real": params.get("lambdaT3_real"),
                "y1_real": y1,
                "y2_real": y2,
            },
            "prediction": prediction,
            "ordering": best.get("ordering"),
            "source_cache": str(path),
        }

    rows.sort(key=lambda row: row["model"])
    return rows, full


def _write_rows_csv(rows: list[dict[str, Any]], path: Path) -> None:
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=tuple(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _write_optimal_parameter_outputs(
    benchmark_root: Path,
    output: Path,
) -> int:
    rows, full = _benchmark_parameter_rows(benchmark_root)
    if not rows:
        return 0
    optimal = output / "optimal"
    optimal.mkdir(parents=True, exist_ok=True)
    _write_rows_csv(rows, optimal / "optimal_parameters.csv")
    (optimal / "optimal_parameters.json").write_text(
        json.dumps(
            {
                "description": (
                    "Per-model optimal benchmark selected by the configured "
                    "oscillation chi-square search. CSV contains compact Yukawa "
                    "norms; this JSON retains the full selected real y1/y2 matrices."
                ),
                "chi2_definition": (
                    "Multivariate Gaussian chi-square of delta_m21_sq_ev2, "
                    "delta_m3l_sq_ev2, sin2_theta12, sin2_theta13 and "
                    "sin2_theta23 against the configured oscillation target."
                ),
                "models": full,
            },
            indent=2,
        ) + "\n",
        encoding="utf-8",
    )
    return len(rows)


def build_display(
    data_input: Path,
    report_input: Path,
    output: Path,
    benchmark_input: Path | None = None,
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
            (p for p in comparison.iterdir() if p.is_dir()),
            key=lambda p: (
                SCENARIO_ORDER.index(p.name)
                if p.name in SCENARIO_ORDER else len(SCENARIO_ORDER),
                p.name,
            ),
        )
        for source in scenario_dirs:
            destination = output / "comparison" / source.name
            _copy_report_section(source, destination)
            report_model_count += len(_report_model_dirs(source))

    optimal = report_input / "optimal"
    if optimal.is_dir():
        _copy_report_section(optimal, output / "optimal")
        report_model_count += len(_report_model_dirs(optimal))

    rows = _collect_rows(data_input)
    _write_csv(rows, output / "model_summary.csv")

    if benchmark_input is None:
        benchmark_input = PROJECT_ROOT / "output" / "benchmarks"
    optimal_parameter_count = _write_optimal_parameter_outputs(
        Path(benchmark_input).resolve(), output
    )

    _write_readme(
        output / "README.md",
        rows=rows,
        report_model_count=report_model_count,
    )

    print(f"Display results written to: {output}")
    print(f"Report model directories: {report_model_count}")
    print(f"Numerical summary rows: {len(rows)}")
    print(f"Optimal benchmark parameter rows: {optimal_parameter_count}")
    return output


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--data-input",
        type=Path,
        default=DEFAULT_DATA_INPUT,
        help="numerical data root (default: output/full)",
    )
    parser.add_argument(
        "--report-input",
        type=Path,
        default=DEFAULT_REPORT_INPUT,
        help="report root (default: Reports/output/full)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help="display root to recreate (default: Reports/output/display)",
    )
    parser.add_argument(
        "--benchmark-input",
        type=Path,
        default=PROJECT_ROOT / "output" / "benchmarks",
        help="cached optimal benchmark root (default: output/benchmarks)",
    )
    args = parser.parse_args()
    build_display(
        args.data_input,
        args.report_input,
        args.output,
        args.benchmark_input,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
