"""
Standard plotting/diagnostics for T3 numerical parameter-scan JSON output.

Usage
-----
python -m Numerical.PlotScanResults SCAN_JSON --target TARGET_JSON

The target JSON is the Gaussianized oscillation-fit target used by
Numerical.OscillationFit, with ``central_values`` and ``one_sigma_errors``.

Outputs
-------
The script creates one output directory containing:
- best_chi2_progress.png
- chi2_distribution.png
- best_point_pulls.png
- best_point_observable_ratios.png
- parameter_vs_chi2/*.png
- summary.txt

No physics is recomputed here: the script visualizes the predictions and chi^2
values already stored in the scan result.
"""

from __future__ import annotations

import argparse
import json
import math
import re
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
    with path.open("r", encoding="utf-8") as handle:
        data = json.load(handle)
    if not isinstance(data, dict):
        raise ValueError(f"{path} must contain a JSON object.")
    return data


def _successful_points(scan: dict[str, Any]) -> list[dict[str, Any]]:
    points = scan.get("points")
    if not isinstance(points, list):
        raise ValueError("Scan JSON must contain a 'points' list.")

    successful: list[dict[str, Any]] = []
    for point in points:
        if (
            isinstance(point, dict)
            and point.get("status") == "Success"
            and point.get("chi2") is not None
            and point.get("prediction") is not None
        ):
            successful.append(point)

    if not successful:
        raise ValueError("Scan contains no successful points with chi2/prediction.")
    return successful


def _target_data(
    target: dict[str, Any],
) -> tuple[dict[str, float], dict[str, float]]:
    central = target.get("central_values")
    sigma = target.get("one_sigma_errors")
    if not isinstance(central, dict) or not isinstance(sigma, dict):
        raise ValueError(
            "Target JSON must contain 'central_values' and 'one_sigma_errors'."
        )

    central_out = {str(k): float(v) for k, v in central.items()}
    sigma_out = {str(k): float(v) for k, v in sigma.items()}

    for name in OBSERVABLE_LABELS:
        if name not in central_out or name not in sigma_out:
            raise ValueError(f"Target JSON is missing observable {name!r}.")
        if sigma_out[name] <= 0.0:
            raise ValueError(f"Target uncertainty for {name!r} must be positive.")

    return central_out, sigma_out


def _safe_filename(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", name)


def _save_close(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(path, dpi=180, bbox_inches="tight")
    plt.close()


def plot_best_chi2_progress(
    points: list[dict[str, Any]],
    output_path: Path,
) -> None:
    ordered = sorted(points, key=lambda point: int(point["index"]))
    indices = np.asarray([int(point["index"]) for point in ordered], dtype=int)
    chi2 = np.asarray([float(point["chi2"]) for point in ordered], dtype=float)
    running_best = np.minimum.accumulate(chi2)

    plt.figure(figsize=(7.2, 4.8))
    plt.plot(indices, running_best)
    plt.yscale("log")
    plt.xlabel("Completed scan point index")
    plt.ylabel(r"Best $\chi^2$ so far")
    plt.title(r"Best $\chi^2$ versus scan progress")
    plt.grid(True, alpha=0.25)
    _save_close(output_path)


def plot_chi2_distribution(
    points: list[dict[str, Any]],
    output_path: Path,
) -> None:
    chi2 = np.asarray([float(point["chi2"]) for point in points], dtype=float)
    positive = chi2[chi2 > 0.0]
    if positive.size == 0:
        raise ValueError("Cannot plot log10 chi2 distribution: no positive values.")

    plt.figure(figsize=(7.2, 4.8))
    plt.hist(np.log10(positive), bins=50)
    plt.xlabel(r"$\log_{10}\chi^2$")
    plt.ylabel("Number of successful points")
    plt.title(r"Distribution of scan $\chi^2$")
    plt.grid(True, alpha=0.25)
    _save_close(output_path)


def _best_point(points: list[dict[str, Any]]) -> dict[str, Any]:
    return min(points, key=lambda point: float(point["chi2"]))


def _pulls(
    point: dict[str, Any],
    central: dict[str, float],
    sigma: dict[str, float],
) -> dict[str, float]:
    prediction = point["prediction"]
    return {
        name: (float(prediction[name]) - central[name]) / sigma[name]
        for name in OBSERVABLE_LABELS
    }


def plot_best_point_pulls(
    point: dict[str, Any],
    central: dict[str, float],
    sigma: dict[str, float],
    output_path: Path,
) -> None:
    pulls = _pulls(point, central, sigma)
    names = list(OBSERVABLE_LABELS)
    values = np.asarray([pulls[name] for name in names], dtype=float)
    labels = [OBSERVABLE_LABELS[name] for name in names]

    plt.figure(figsize=(8.0, 4.8))
    x = np.arange(len(names))
    plt.bar(x, values)
    plt.axhline(0.0, linewidth=1.0)
    plt.axhline(1.0, linestyle="--", linewidth=0.8)
    plt.axhline(-1.0, linestyle="--", linewidth=0.8)
    plt.xticks(x, labels)
    plt.ylabel(r"Pull $(O_{\rm pred}-O_{\rm target})/\sigma$")
    plt.title(
        rf"Best-point pulls: index {int(point['index'])}, "
        rf"$\chi^2={float(point['chi2']):.3g}$"
    )
    plt.grid(True, axis="y", alpha=0.25)
    _save_close(output_path)


def plot_best_point_observable_ratios(
    point: dict[str, Any],
    central: dict[str, float],
    output_path: Path,
) -> None:
    prediction = point["prediction"]
    names = list(OBSERVABLE_LABELS)
    ratios = np.asarray(
        [float(prediction[name]) / central[name] for name in names],
        dtype=float,
    )
    labels = [OBSERVABLE_LABELS[name] for name in names]

    plt.figure(figsize=(8.0, 4.8))
    x = np.arange(len(names))
    plt.bar(x, ratios)
    plt.axhline(1.0, linestyle="--", linewidth=1.0)
    plt.xticks(x, labels)
    plt.ylabel("Prediction / target central value")
    plt.title("Best-point observable ratios")
    plt.grid(True, axis="y", alpha=0.25)
    _save_close(output_path)


def plot_parameter_vs_chi2(
    points: list[dict[str, Any]],
    output_dir: Path,
) -> list[tuple[str, float]]:
    parameter_names = sorted(
        {
            str(name)
            for point in points
            for name in point.get("parameters", {}).keys()
        }
    )

    chi2 = np.asarray([float(point["chi2"]) for point in points], dtype=float)
    if np.any(chi2 <= 0.0):
        mask = chi2 > 0.0
    else:
        mask = np.ones_like(chi2, dtype=bool)
    log_chi2 = np.log10(chi2[mask])

    correlations: list[tuple[str, float]] = []

    for name in parameter_names:
        values = []
        valid_mask = []
        for point in points:
            parameters = point.get("parameters", {})
            value = parameters.get(name)
            if value is None:
                values.append(np.nan)
                valid_mask.append(False)
            else:
                values.append(float(value))
                valid_mask.append(True)

        values_array = np.asarray(values, dtype=float)
        combined = mask & np.asarray(valid_mask, dtype=bool)
        x = values_array[combined]
        y = np.log10(chi2[combined])

        if x.size < 2:
            continue

        if np.std(x) > 0.0 and np.std(y) > 0.0:
            correlation = float(np.corrcoef(x, y)[0, 1])
        else:
            correlation = float("nan")
        correlations.append((name, correlation))

        plt.figure(figsize=(6.6, 4.8))
        plt.scatter(x, 10.0**y, s=12, alpha=0.6)
        plt.yscale("log")
        plt.xlabel(name)
        plt.ylabel(r"$\chi^2$")
        plt.title(rf"{name} versus $\chi^2$")
        plt.grid(True, alpha=0.25)
        _save_close(output_dir / f"{_safe_filename(name)}_vs_chi2.png")

    correlations.sort(
        key=lambda item: abs(item[1]) if math.isfinite(item[1]) else -1.0,
        reverse=True,
    )
    return correlations


def write_summary(
    scan: dict[str, Any],
    point: dict[str, Any],
    central: dict[str, float],
    sigma: dict[str, float],
    correlations: list[tuple[str, float]],
    output_path: Path,
) -> None:
    pulls = _pulls(point, central, sigma)
    prediction = point["prediction"]

    lines = [
        "T3 parameter-scan plotting summary",
        "",
        f"status: {scan.get('status')}",
        f"point_count: {scan.get('point_count')}",
        f"successful_point_count: {scan.get('successful_point_count')}",
        f"failed_point_count: {scan.get('failed_point_count')}",
        f"best_point_index: {point['index']}",
        f"best_chi2: {float(point['chi2']):.12g}",
        "",
        "Best-point observables:",
    ]

    for name in OBSERVABLE_LABELS:
        contribution = pulls[name] ** 2
        lines.append(
            f"  {name}: prediction={float(prediction[name]):.12g}, "
            f"target={central[name]:.12g}, sigma={sigma[name]:.12g}, "
            f"pull={pulls[name]:+.6g}, pull^2={contribution:.6g}"
        )

    lines.extend(["", "Best-point parameters:"])
    for name, value in sorted(point.get("parameters", {}).items()):
        lines.append(f"  {name}: {float(value):.12g}")

    lines.extend(
        [
            "",
            "Linear correlation with log10(chi2) across successful scan points:",
            "(diagnostic only; a small linear correlation does not imply irrelevance)",
        ]
    )
    for name, correlation in correlations:
        lines.append(f"  {name}: {correlation:+.6f}")

    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def generate_plots(
    scan_path: Path,
    target_path: Path,
    output_dir: Path,
) -> Path:
    scan = _load_json(scan_path)
    target = _load_json(target_path)
    central, sigma = _target_data(target)
    points = _successful_points(scan)
    best = _best_point(points)

    output_dir.mkdir(parents=True, exist_ok=True)

    plot_best_chi2_progress(
        points,
        output_dir / "best_chi2_progress.png",
    )
    plot_chi2_distribution(
        points,
        output_dir / "chi2_distribution.png",
    )
    plot_best_point_pulls(
        best,
        central,
        sigma,
        output_dir / "best_point_pulls.png",
    )
    plot_best_point_observable_ratios(
        best,
        central,
        output_dir / "best_point_observable_ratios.png",
    )
    correlations = plot_parameter_vs_chi2(
        points,
        output_dir / "parameter_vs_chi2",
    )
    write_summary(
        scan,
        best,
        central,
        sigma,
        correlations,
        output_dir / "summary.txt",
    )

    return output_dir


def _default_output_dir(scan_path: Path) -> Path:
    return Path("output") / "plots" / scan_path.stem


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate standard plots from a T3 scan result JSON."
    )
    parser.add_argument("scan_json", type=Path)
    parser.add_argument(
        "--target",
        type=Path,
        required=True,
        help="Gaussian oscillation-fit target JSON.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Output directory. Default: output/plots/<scan stem>/",
    )

    args = parser.parse_args()
    output_dir = (
        args.output_dir
        if args.output_dir is not None
        else _default_output_dir(args.scan_json)
    )

    generated = generate_plots(
        args.scan_json,
        args.target,
        output_dir,
    )
    print(f"Plots written to: {generated}")


if __name__ == "__main__":
    main()
