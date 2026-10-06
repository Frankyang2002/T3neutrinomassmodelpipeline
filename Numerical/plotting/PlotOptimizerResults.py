"""
Plot and summarize a T3 OscillationOptimizer result JSON.

Usage
-----
python -m Numerical.PlotOptimizerResults OPTIMIZER_JSON

Optional:
python -m Numerical.PlotOptimizerResults OPTIMIZER_JSON --output-dir PATH

Outputs
-------
- chi2_trajectory.png
- best_chi2_trajectory.png
- final_pulls.png
- sensitivity_ranking.png
- active_parameter_trajectory.png
- summary.txt

The script does not recompute the T3 physics. It visualizes quantities already
stored by Numerical.OscillationOptimizer.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np


OBSERVABLE_LABELS = {
    "delta_m21_sq_ev2": r"$\Delta m_{21}^2$",
    "delta_m3l_sq_ev2": r"$\Delta m_{3\ell}^2$",
    "sin2_theta12": r"$\sin^2\theta_{12}$",
    "sin2_theta13": r"$\sin^2\theta_{13}$",
    "sin2_theta23": r"$\sin^2\theta_{23}$",
}


def _load_json(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"{path} must contain a JSON object.")
    return payload


def _successful_history(payload: dict[str, Any]) -> list[dict[str, Any]]:
    history = payload.get("history")
    if not isinstance(history, list):
        raise ValueError("Optimizer JSON must contain a 'history' list.")

    success = [
        item
        for item in history
        if (
            isinstance(item, dict)
            and item.get("status") == "Success"
            and item.get("chi2") is not None
        )
    ]
    if not success:
        raise ValueError("Optimizer result has no successful history entries.")
    return success


def _save(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(path, dpi=180, bbox_inches="tight")
    plt.close()


def plot_chi2_trajectory(
    history: list[dict[str, Any]],
    output_path: Path,
) -> None:
    x = np.asarray([int(item["evaluation"]) for item in history], dtype=int)
    y = np.asarray([float(item["chi2"]) for item in history], dtype=float)

    plt.figure(figsize=(7.2, 4.8))
    plt.plot(x, y, marker="o", markersize=2.5, linewidth=1.0)
    plt.yscale("log")
    plt.xlabel("Pipeline evaluation")
    plt.ylabel(r"$\chi^2$")
    plt.title(r"Optimizer $\chi^2$ trajectory")
    plt.grid(True, alpha=0.25)
    _save(output_path)


def plot_best_chi2_trajectory(
    history: list[dict[str, Any]],
    output_path: Path,
) -> None:
    x = np.asarray([int(item["evaluation"]) for item in history], dtype=int)
    y = np.asarray([float(item["chi2"]) for item in history], dtype=float)
    best = np.minimum.accumulate(y)

    plt.figure(figsize=(7.2, 4.8))
    plt.plot(x, best, linewidth=1.4)
    plt.yscale("log")
    plt.xlabel("Pipeline evaluation")
    plt.ylabel(r"Best $\chi^2$ so far")
    plt.title(r"Best-fit convergence")
    plt.grid(True, alpha=0.25)
    _save(output_path)


def plot_final_pulls(
    payload: dict[str, Any],
    output_path: Path,
) -> None:
    best = payload.get("best")
    if not isinstance(best, dict):
        raise ValueError("Optimizer JSON has no 'best' object.")

    names = best.get("observable_names")
    residual = best.get("whitened_residual")
    if not isinstance(names, list) or not isinstance(residual, list):
        raise ValueError(
            "Best point must contain observable_names and whitened_residual."
        )

    values = np.asarray([float(x) for x in residual], dtype=float)
    labels = [OBSERVABLE_LABELS.get(str(name), str(name)) for name in names]

    plt.figure(figsize=(8.0, 4.8))
    x = np.arange(len(values))
    plt.bar(x, values)
    plt.axhline(0.0, linewidth=1.0)
    plt.axhline(1.0, linestyle="--", linewidth=0.8)
    plt.axhline(-1.0, linestyle="--", linewidth=0.8)
    plt.xticks(x, labels)
    plt.ylabel("Whitened residual / pull")
    plt.title("Final fitted pulls")
    plt.grid(True, axis="y", alpha=0.25)
    _save(output_path)


def plot_sensitivity_ranking(
    payload: dict[str, Any],
    output_path: Path,
) -> None:
    ranking = payload.get("sensitivity_ranking")
    if not isinstance(ranking, list) or not ranking:
        raise ValueError("Optimizer JSON has no sensitivity_ranking list.")

    names = [str(item["name"]) for item in ranking]
    scores = np.asarray([float(item["score"]) for item in ranking], dtype=float)

    order = np.argsort(scores)
    names = [names[i] for i in order]
    scores = scores[order]

    plt.figure(figsize=(8.0, 6.2))
    y = np.arange(len(names))
    plt.barh(y, scores)
    plt.yticks(y, names)
    plt.xlabel("Local residual sensitivity score")
    plt.title("Parameter sensitivity ranking at warm start")
    plt.grid(True, axis="x", alpha=0.25)
    _save(output_path)


def plot_active_parameter_trajectory(
    payload: dict[str, Any],
    history: list[dict[str, Any]],
    output_path: Path,
) -> None:
    active = payload.get("active_parameters")
    if not isinstance(active, list) or not active:
        raise ValueError("Optimizer JSON has no active_parameters list.")

    ls_entries = [
        item for item in history
        if str(item.get("phase", "")).startswith("least_squares")
    ]
    if not ls_entries:
        raise ValueError("No least_squares history entries were found.")

    x = np.asarray([int(item["evaluation"]) for item in ls_entries], dtype=int)

    plt.figure(figsize=(8.2, 5.6))
    for name in active:
        values = []
        for item in ls_entries:
            parameters = item.get("parameters")
            if not isinstance(parameters, dict) or name not in parameters:
                raise ValueError(
                    f"History entry is missing active parameter {name!r}."
                )
            values.append(float(parameters[name]))
        plt.plot(x, values, label=str(name), linewidth=1.1)

    plt.xlabel("Pipeline evaluation")
    plt.ylabel("Parameter value")
    plt.title("Active-parameter evolution during least squares")
    plt.grid(True, alpha=0.25)
    plt.legend(fontsize=8, ncol=2)
    _save(output_path)


def write_summary(
    payload: dict[str, Any],
    output_path: Path,
) -> None:
    best = payload.get("best")
    if not isinstance(best, dict):
        raise ValueError("Optimizer JSON has no 'best' object.")

    lines = [
        "T3 OscillationOptimizer summary",
        "",
        f"status: {payload.get('status')}",
        f"seed_source_kind: {payload.get('seed_source_kind')}",
        f"seed_stored_chi2: {payload.get('seed_stored_chi2')}",
        f"seed_recomputed_chi2: {payload.get('seed_recomputed_chi2')}",
        f"active_parameter_count: {payload.get('active_parameter_count')}",
        f"total_pipeline_evaluations: {payload.get('total_pipeline_evaluations')}",
        f"successful_evaluations: {payload.get('successful_evaluations')}",
        f"failed_evaluations: {payload.get('failed_evaluations')}",
        f"best_chi2: {best.get('chi2')}",
        "",
        "Active parameters:",
    ]

    for name in payload.get("active_parameters", []):
        lines.append(f"  {name}")

    lines.extend(["", "Best-point predictions:"])
    prediction = best.get("prediction", {})
    for name, value in prediction.items():
        lines.append(f"  {name}: {float(value):.12g}")

    lines.extend(["", "Best-point pulls:"])
    names = best.get("observable_names", [])
    residual = best.get("whitened_residual", [])
    for name, value in zip(names, residual):
        lines.append(f"  {name}: {float(value):+.8g}")

    lines.extend(["", "Best-point parameters:"])
    parameters = best.get("parameters", {})
    for name, value in sorted(parameters.items()):
        lines.append(f"  {name}: {float(value):.12g}")

    lines.extend(["", "Sensitivity ranking:"])
    for item in payload.get("sensitivity_ranking", []):
        lines.append(
            f"  {item['name']}: score={float(item['score']):.8g}"
        )

    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def generate_plots(
    optimizer_json: Path,
    output_dir: Path,
) -> Path:
    payload = _load_json(optimizer_json)
    history = _successful_history(payload)

    output_dir.mkdir(parents=True, exist_ok=True)

    plot_chi2_trajectory(
        history,
        output_dir / "chi2_trajectory.png",
    )
    plot_best_chi2_trajectory(
        history,
        output_dir / "best_chi2_trajectory.png",
    )
    plot_final_pulls(
        payload,
        output_dir / "final_pulls.png",
    )
    plot_sensitivity_ranking(
        payload,
        output_dir / "sensitivity_ranking.png",
    )
    plot_active_parameter_trajectory(
        payload,
        history,
        output_dir / "active_parameter_trajectory.png",
    )
    write_summary(
        payload,
        output_dir / "summary.txt",
    )

    return output_dir


def _default_output_dir(path: Path) -> Path:
    return Path("output") / "plots" / path.stem


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Plot a T3 OscillationOptimizer result JSON."
    )
    parser.add_argument("optimizer_json", type=Path)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Default: output/plots/<optimizer json stem>/",
    )
    args = parser.parse_args()

    output_dir = (
        args.output_dir
        if args.output_dir is not None
        else _default_output_dir(args.optimizer_json)
    )

    generated = generate_plots(
        args.optimizer_json,
        output_dir,
    )
    print(f"Plots written to: {generated}")


if __name__ == "__main__":
    main()
