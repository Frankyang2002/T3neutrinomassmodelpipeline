"""Four-scenario, within-model plots from saved T3 numerical diagnostics only.

The dashed intermediate 'hard + direct' C5 is a threshold-anchored display
reconstruction, NOT an intermediate-EFT Weinberg Wilson coefficient.
"""
from __future__ import annotations

from pathlib import Path
from typing import Mapping

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from Numerical.plotting.RunningTrajectoryData import (
    QUANTITIES, SCENARIOS, SavedRun, extract_trajectory,
)
from Numerical.plotting.WeinbergMatchedContinuation import reconstruct_matched_continuation

COLORS = {"smallY": "tab:blue", "largeY": "tab:orange"}
STYLE = {"smallL": "-", "largeL": "--"}
LABELS = {
    "smallY_smallL": "Y=0.005, lambda(ref)=0.1",
    "smallY_largeL": "Y=0.005, lambda(ref)=1.0",
    "largeY_smallL": "Y=0.5, lambda(ref)=0.1",
    "largeY_largeL": "Y=0.5, lambda(ref)=1.0",
}


def _style(scenario: str) -> dict[str, object]:
    y, lam = scenario.split("_")
    return {"color": COLORS[y], "linestyle": STYLE[lam], "linewidth": 1.8,
            "label": LABELS[scenario]}


def _finish(ax, output: Path, *, log_y: bool = True, threshold: float | None = None) -> None:
    ax.set_xscale("log")
    if log_y:
        ax.set_yscale("log")
    if threshold is not None:
        ax.axvline(threshold, color="0.4", alpha=0.5, lw=1, ls=":")
    ax.set_xlabel(r"Renormalisation scale $\mu$ [GeV]")
    ax.grid(True, alpha=0.2)
    ax.legend(fontsize=8, ncol=2)
    fig = ax.figure
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=210, bbox_inches="tight")
    plt.close(fig)


def plot_comparison(runs: Mapping[str, SavedRun], quantity: str, output: Path) -> None:
    fig, ax = plt.subplots(figsize=(8, 4.8))
    positives = []
    for scenario in SCENARIOS:
        if scenario not in runs:
            continue
        mu, value = extract_trajectory(runs[scenario].payload, quantity)
        # A zero is meaningful but cannot be put on a logarithmic ordinate.
        positive = value > 0
        positives.append(bool(np.all(positive)))
        ax.plot(mu, np.where(positive, value, np.nan), **_style(scenario))
    ax.set_ylabel(QUANTITIES[quantity][3])
    ax.set_title(f"{quantity}: four benchmarks for {next(iter(runs.values())).model_key}")
    _finish(ax, output, log_y=all(positives))


def plot_weinberg_eft(runs: Mapping[str, SavedRun], output: Path) -> list[str]:
    """Join checked display reconstructions with stored final C5, when possible."""
    fig, ax = plt.subplots(figsize=(8.6, 5.1))
    issues: list[str] = []
    scalar_thresholds = []
    for scenario in SCENARIOS:
        if scenario not in runs:
            continue
        run = runs[scenario]
        mu, final = extract_trajectory(run.payload, "c5")
        style = _style(scenario)
        ax.plot(mu, np.where(final > 0, final, np.nan), **style)
        scalar_thresholds.append(float(run.payload["scales_gev"]["mu_scalar_threshold"]))
        try:
            matched = reconstruct_matched_continuation(run.payload)
        except (KeyError, ValueError, TypeError, IndexError) as exc:
            issues.append(f"{scenario}: no checked intermediate continuation ({exc})")
            continue
        intermediate = np.linalg.norm(matched.combined.reshape(len(matched.mu_gev), 9), axis=1)
        ax.plot(matched.mu_gev, np.where(intermediate > 0, intermediate, np.nan),
                color=style["color"], linewidth=1.25, linestyle=":", alpha=0.8)
    if scalar_thresholds and np.allclose(scalar_thresholds, scalar_thresholds[0]):
        ax.axvline(scalar_thresholds[0], lw=1, ls="-.", alpha=0.5, color="0.5")
    ax.set_ylabel(r"$\Vert C_5\Vert_F$ [GeV$^{-1}$]")
    ax.set_title("Solid/dashed: final SM+Weinberg; dotted: checked display continuation")
    _finish(ax, output)
    return issues


def plot_threshold_components(runs: Mapping[str, SavedRun], output: Path) -> None:
    fig, ax = plt.subplots(figsize=(8, 4.6))
    names = ("hard", "direct_running", "combined")
    x = np.arange(len(names))
    for index, scenario in enumerate(SCENARIOS):
        if scenario not in runs:
            continue
        data = runs[scenario].payload["c5_threshold_contributions"]
        magnitudes = [float(np.linalg.norm(np.asarray(data[name]["abs"], dtype=float)))
                      for name in names]
        ax.bar(x + (index - 1.5) * 0.19, magnitudes, width=0.18,
               color=_style(scenario)["color"], alpha=0.5 if "largeL" in scenario else 0.95,
               label=LABELS[scenario])
    ax.set_xticks(x, ["Hard at $M_S$", "Direct at $M_S$", "Combined at $M_S$"])
    ax.set_ylabel(r"Frobenius norm [GeV$^{-1}$]")
    ax.set_title("Scalar-threshold matching contributions (magnitudes)")
    ax.grid(axis="y", alpha=0.2)
    ax.legend(fontsize=8, ncol=2)
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=210, bbox_inches="tight")
    plt.close(fig)


def generate_model_figures(runs: Mapping[str, SavedRun], figure_dir: Path) -> tuple[list[str], list[str]]:
    figures: list[str] = []
    for name in ("y1", "y2", "lambda", "direct", "c5", "dm21", "dm3l", "m1", "m2", "m3"):
        target = figure_dir / f"{name}_comparison.png"
        plot_comparison(runs, name, target)
        figures.append(target.name)
    target = figure_dir / "weinberg_eft_comparison.png"
    issues = plot_weinberg_eft(runs, target)
    figures.append(target.name)
    target = figure_dir / "c5_threshold_components.png"
    plot_threshold_components(runs, target)
    figures.append(target.name)
    return figures, issues
