"""Four-scenario, within-model plots from saved T3 numerical diagnostics only.

The goal is to make the within-model static figures visually consistent with the
existing comparison figures: clean single-panel plots, consistent typography,
clear legends, threshold markers where relevant, and restrained line styling.

The dotted intermediate 'hard + direct' C5 curve is a threshold-anchored display
reconstruction, NOT an intermediate-EFT Weinberg Wilson coefficient.
"""
from __future__ import annotations

from pathlib import Path
from typing import Mapping

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, NullFormatter, ScalarFormatter
import numpy as np

from Numerical.plotting.RunningTrajectoryData import (
    QUANTITIES,
    SCENARIOS,
    SavedRun,
    extract_trajectory,
)
from Numerical.plotting.WeinbergMatchedContinuation import (
    reconstruct_matched_continuation,
)

# Use four distinct scenario colours to improve readability while keeping the
# labels explicit.  The output remains recognisably matplotlib-default and sits
# closer to the existing comparison figures than the earlier custom two-colour
# encoding.
SCENARIO_STYLES: dict[str, dict[str, object]] = {
    "smallY_smallL": {
        "color": "tab:blue",
        "linestyle": "-",
        "linewidth": 2.0,
        "label": r"smallY\_smallL  $(Y=0.005,\ \lambda_{\rm ref}=0.1)$",
    },
    "smallY_largeL": {
        "color": "tab:orange",
        "linestyle": "-",
        "linewidth": 2.0,
        "label": r"smallY\_largeL  $(Y=0.005,\ \lambda_{\rm ref}=1.0)$",
    },
    "largeY_smallL": {
        "color": "tab:green",
        "linestyle": "-",
        "linewidth": 2.0,
        "label": r"largeY\_smallL  $(Y=0.5,\ \lambda_{\rm ref}=0.1)$",
    },
    "largeY_largeL": {
        "color": "tab:red",
        "linestyle": "-",
        "linewidth": 2.0,
        "label": r"largeY\_largeL  $(Y=0.5,\ \lambda_{\rm ref}=1.0)$",
    },
}

# Short plain-text labels remain useful in the HTML report table.
LABELS = {
    key: str(style["label"]).replace("$", "").replace("\\", "")
    for key, style in SCENARIO_STYLES.items()
}


def _style(scenario: str) -> dict[str, object]:
    return dict(SCENARIO_STYLES[scenario])


# Layout constants used across all static figures.
_FIGSIZE = (8.4, 5.0)
_TITLE_SIZE = 11
_LABEL_SIZE = 10
_TICK_SIZE = 9
_LEGEND_SIZE = 8


def _apply_common_style(ax, *, ylabel: str, title: str) -> None:
    ax.set_xscale("log")
    ax.set_xlabel(r"Renormalisation scale $\mu$ [GeV]", fontsize=_LABEL_SIZE)
    ax.set_ylabel(ylabel, fontsize=_LABEL_SIZE)
    ax.set_title(title, fontsize=_TITLE_SIZE)
    ax.grid(True, which="major", alpha=0.28)
    ax.grid(True, which="minor", alpha=0.10)
    ax.tick_params(axis="both", which="major", labelsize=_TICK_SIZE)
    ax.tick_params(axis="both", which="minor", labelsize=_TICK_SIZE - 1)
    ax.xaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1))
    ax.xaxis.set_minor_formatter(NullFormatter())



def _maybe_set_log_y(ax, values: list[np.ndarray]) -> None:
    positive = True
    spread_values: list[float] = []
    for arr in values:
        if arr.size == 0:
            continue
        positive = positive and bool(np.all(arr > 0.0))
        spread_values.extend(arr[np.isfinite(arr)].tolist())

    if positive and spread_values:
        finite = np.asarray(spread_values, dtype=float)
        if finite.min() > 0.0 and finite.max() / finite.min() > 20.0:
            ax.set_yscale("log")



def _format_linear_y(ax) -> None:
    if ax.get_yscale() == "linear":
        formatter = ScalarFormatter(useMathText=True)
        formatter.set_powerlimits((-2, 3))
        ax.yaxis.set_major_formatter(formatter)



def _finish(ax, output: Path, *, thresholds: list[tuple[float, str]] | None = None) -> None:
    if thresholds:
        ymin, ymax = ax.get_ylim()
        for value, label in thresholds:
            ax.axvline(value, color="0.45", alpha=0.75, lw=1.1, ls="--")
            if ymin > 0 and ymax > ymin:
                y_text = ymax / (10 ** 0.06) if ax.get_yscale() == "log" else ymax - 0.06 * (ymax - ymin)
            else:
                y_text = ymax
            ax.text(
                value,
                y_text,
                label,
                rotation=90,
                va="top",
                ha="right",
                fontsize=8,
                color="0.35",
                bbox={"boxstyle": "round,pad=0.12", "facecolor": "white", "edgecolor": "none", "alpha": 0.7},
            )

    ax.legend(fontsize=_LEGEND_SIZE, ncol=2, frameon=True, loc="best")
    _format_linear_y(ax)

    fig = ax.figure
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(fig)



def _extract_thresholds(runs: Mapping[str, SavedRun]) -> tuple[float | None, float | None]:
    fermion_values: list[float] = []
    scalar_values: list[float] = []
    for run in runs.values():
        scales = run.payload.get("scales_gev", {})
        mf = scales.get("mu_fermion_threshold")
        ms = scales.get("mu_scalar_threshold")
        if mf is not None:
            fermion_values.append(float(mf))
        if ms is not None:
            scalar_values.append(float(ms))

    def _unique_or_none(values: list[float]) -> float | None:
        if not values:
            return None
        first = values[0]
        return first if np.allclose(values, first) else None

    return _unique_or_none(fermion_values), _unique_or_none(scalar_values)



def _title(model_key: str, quantity_label: str) -> str:
    return f"{quantity_label} — four benchmark scenarios for {model_key}"



def plot_comparison(runs: Mapping[str, SavedRun], quantity: str, output: Path) -> None:
    fig, ax = plt.subplots(figsize=_FIGSIZE)
    all_values: list[np.ndarray] = []
    model_key = next(iter(runs.values())).model_key

    for scenario in SCENARIOS:
        if scenario not in runs:
            continue
        mu, value = extract_trajectory(runs[scenario].payload, quantity)
        all_values.append(value)
        ax.plot(mu, value, **_style(scenario))

    ylabel = QUANTITIES[quantity][3]
    _apply_common_style(ax, ylabel=ylabel, title=_title(model_key, quantity))
    _maybe_set_log_y(ax, all_values)

    mf, ms = _extract_thresholds(runs)
    thresholds: list[tuple[float, str]] = []
    if quantity in {"y1", "y2", "lambda"} and mf is not None:
        thresholds.append((mf, r"$M_F$"))
    if quantity == "direct":
        if mf is not None:
            thresholds.append((mf, r"$M_F$"))
        if ms is not None:
            thresholds.append((ms, r"$M_S$"))
    if quantity in {"c5", "dm21", "dm3l", "m1", "m2", "m3"} and ms is not None:
        thresholds.append((ms, r"$M_S$"))

    _finish(ax, output, thresholds=thresholds)



def plot_weinberg_eft(runs: Mapping[str, SavedRun], output: Path) -> list[str]:
    """Join checked display reconstructions with stored final C5, when possible."""
    fig, ax = plt.subplots(figsize=(8.6, 5.2))
    issues: list[str] = []
    all_values: list[np.ndarray] = []
    model_key = next(iter(runs.values())).model_key
    mf, ms = _extract_thresholds(runs)

    for scenario in SCENARIOS:
        if scenario not in runs:
            continue
        run = runs[scenario]
        mu, final = extract_trajectory(run.payload, "c5")
        style = _style(scenario)
        ax.plot(mu, final, **style)
        all_values.append(final)

        try:
            matched = reconstruct_matched_continuation(run.payload)
        except (KeyError, ValueError, TypeError, IndexError) as exc:
            issues.append(f"{scenario}: no checked intermediate continuation ({exc})")
            continue

        intermediate = np.linalg.norm(
            matched.combined.reshape(len(matched.mu_gev), 9),
            axis=1,
        )
        ax.plot(
            matched.mu_gev,
            intermediate,
            color=style["color"],
            linewidth=1.6,
            linestyle=":",
            alpha=0.9,
            label=f"{style['label']}  [matched display]",
        )
        all_values.append(intermediate)

    _apply_common_style(
        ax,
        ylabel=r"$\Vert C_5\Vert_F$ [GeV$^{-1}$]",
        title=_title(model_key, "weinberg_eft"),
    )
    _maybe_set_log_y(ax, all_values)

    subtitle = (
        "solid: stored final SM+Weinberg running; dotted: matched hard+direct display continuation"
    )
    ax.text(
        0.01,
        1.01,
        subtitle,
        transform=ax.transAxes,
        fontsize=8.2,
        color="0.35",
        ha="left",
        va="bottom",
    )

    thresholds: list[tuple[float, str]] = []
    if mf is not None:
        thresholds.append((mf, r"$M_F$"))
    if ms is not None:
        thresholds.append((ms, r"$M_S$"))
    _finish(ax, output, thresholds=thresholds)
    return issues



def plot_threshold_components(runs: Mapping[str, SavedRun], output: Path) -> None:
    fig, ax = plt.subplots(figsize=_FIGSIZE)
    model_key = next(iter(runs.values())).model_key
    names = ("hard", "direct_running", "combined")
    labels = [r"Hard at $M_S$", r"Direct at $M_S$", r"Combined at $M_S$"]
    x = np.arange(len(names), dtype=float)
    bar_width = 0.18

    all_magnitudes: list[np.ndarray] = []
    for index, scenario in enumerate(SCENARIOS):
        if scenario not in runs:
            continue
        data = runs[scenario].payload["c5_threshold_contributions"]
        magnitudes = np.asarray(
            [
                float(np.linalg.norm(np.asarray(data[name]["abs"], dtype=float)))
                for name in names
            ],
            dtype=float,
        )
        all_magnitudes.append(magnitudes)
        offset = (index - 1.5) * bar_width
        ax.bar(
            x + offset,
            magnitudes,
            width=bar_width,
            edgecolor="black",
            linewidth=0.5,
            alpha=0.82,
            label=LABELS[scenario],
            color=str(_style(scenario)["color"]),
        )

    ax.set_xticks(x, labels)
    _apply_common_style(
        ax,
        ylabel=r"Frobenius norm [GeV$^{-1}$]",
        title=_title(model_key, "c5_threshold_components"),
    )
    ax.set_xscale("linear")
    _maybe_set_log_y(ax, all_magnitudes)
    _finish(ax, output)



def generate_model_figures(
    runs: Mapping[str, SavedRun],
    figure_dir: Path,
) -> tuple[list[str], list[str]]:
    """Generate the full static within-model figure set.

    Returns
    -------
    figures:
        Filenames written beneath ``figure_dir``.
    issues:
        Non-fatal reconstruction warnings from the matched Weinberg view.
    """
    figures: list[str] = []
    for name in (
        "y1",
        "y2",
        "lambda",
        "direct",
        "c5",
        "dm21",
        "dm3l",
        "m1",
        "m2",
        "m3",
    ):
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
