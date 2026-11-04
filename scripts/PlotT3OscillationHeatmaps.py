from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np


ROOT = Path(__file__).resolve().parents[1]

DEFAULT_INPUT = ROOT / "Reports" / "output" / "thesis_qualitative" / "representation_running" / "csv" / "frozen_stage_observable_changes.csv"
DEFAULT_PLOT_DIR = ROOT / "Reports" / "output" / "thesis_qualitative" / "representation_running" / "plots"
DEFAULT_CSV_DIR = ROOT / "Reports" / "output" / "thesis_qualitative" / "representation_running" / "csv"

MODEL_ORDER = [
    "A_-4", "A_-2", "A_+0",
    "B_-3", "B_-1", "B_+1",
    "C_-3", "C_-1", "C_+1",
    "D_-2", "D_+0", "D_+2",
    "E_-4", "E_-2", "E_+0", "E_+2",
]

STAGE_ORDER = ["UV", "EFT1", "SMEFT"]

MASSSQ_OBSERVABLES = [
    ("delta_m21_sq_ev2", r"$\Delta m_{21}^2$"),
    ("delta_m3l_sq_ev2", r"$\Delta m_{3\ell}^2$"),
]

MIXING_OBSERVABLES = [
    ("sin2_theta12", r"$\sin^2\theta_{12}$"),
    ("sin2_theta13", r"$\sin^2\theta_{13}$"),
    ("sin2_theta23", r"$\sin^2\theta_{23}$"),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build annotated heatmaps of signed percentage changes in neutrino "
            "oscillation observables from the saved representation-running CSV."
        )
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=DEFAULT_INPUT,
        help="Input CSV containing frozen-stage oscillation-observable changes.",
    )
    parser.add_argument(
        "--plot-dir",
        type=Path,
        default=DEFAULT_PLOT_DIR,
        help="Directory for output PNG plots.",
    )
    parser.add_argument(
        "--csv-dir",
        type=Path,
        default=DEFAULT_CSV_DIR,
        help="Directory for output matrix CSV tables.",
    )
    return parser.parse_args()


def read_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(f"Input CSV not found: {path}")

    rows: list[dict[str, str]] = []
    with path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {
            "model",
            "switch",
            "observable",
            "relative_change_percent",
            "comparison_kind",
        }
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Missing required CSV columns: {sorted(missing)}")

        for row in reader:
            if row["comparison_kind"] != "independently_evolved_running_stage":
                continue
            if row["switch"] not in STAGE_ORDER:
                continue
            rows.append(row)

    if not rows:
        raise ValueError("No usable rows found in input CSV.")
    return rows


def build_matrix(
    rows: list[dict[str, str]],
    observable_defs: list[tuple[str, str]],
) -> tuple[np.ndarray, list[str], list[str]]:
    value_lookup: dict[tuple[str, str, str], float] = {}
    for row in rows:
        key = (row["model"], row["switch"], row["observable"])
        value_lookup[key] = float(row["relative_change_percent"])

    column_keys: list[tuple[str, str]] = []
    column_labels: list[str] = []
    for stage in STAGE_ORDER:
        for observable_key, observable_label in observable_defs:
            column_keys.append((stage, observable_key))
            column_labels.append(f"{stage}\n{observable_label}")

    matrix = np.full((len(MODEL_ORDER), len(column_keys)), np.nan, dtype=float)

    for i, model in enumerate(MODEL_ORDER):
        for j, (stage, observable_key) in enumerate(column_keys):
            value = value_lookup.get((model, stage, observable_key))
            if value is not None:
                matrix[i, j] = value

    return matrix, MODEL_ORDER[:], column_labels


def write_matrix_csv(
    output_path: Path,
    matrix: np.ndarray,
    row_labels: list[str],
    column_labels: list[str],
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["model", *column_labels])
        for model, row in zip(row_labels, matrix):
            writer.writerow(
                [model] + ["" if np.isnan(x) else f"{x:.10g}" for x in row]
            )


def format_cell_value(value: float) -> str:
    if np.isnan(value):
        return ""
    abs_value = abs(value)
    if abs_value >= 10:
        return f"{value:.1f}"
    if abs_value >= 1:
        return f"{value:.2f}"
    if abs_value >= 0.1:
        return f"{value:.3f}"
    if abs_value >= 0.01:
        return f"{value:.4f}"
    return f"{value:.5f}"


def plot_heatmap(
    matrix: np.ndarray,
    row_labels: list[str],
    column_labels: list[str],
    output_path: Path,
    title: str,
    stage_group_size: int,
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)

    finite_values = matrix[np.isfinite(matrix)]
    if finite_values.size == 0:
        raise ValueError(f"No finite values available for plot: {output_path}")

    max_abs = float(np.max(np.abs(finite_values)))
    # Protect against the degenerate all-zero case.
    if max_abs == 0.0:
        max_abs = 1.0

    norm = mcolors.TwoSlopeNorm(vmin=-max_abs, vcenter=0.0, vmax=max_abs)

    fig_width = max(8.5, 0.9 * matrix.shape[1] + 4.0)
    fig_height = max(8.0, 0.38 * matrix.shape[0] + 2.8)

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    image = ax.imshow(matrix, cmap="RdBu_r", norm=norm, aspect="auto")

    ax.set_title(title, pad=16)
    ax.set_xticks(np.arange(len(column_labels)))
    ax.set_xticklabels(column_labels, rotation=0, ha="center")
    ax.set_yticks(np.arange(len(row_labels)))
    ax.set_yticklabels(row_labels)

    ax.set_xlabel("Running stage and oscillation observable")
    ax.set_ylabel("T3 model")

    # Cell annotations.
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            value = matrix[i, j]
            if np.isnan(value):
                continue
            text_color = "white" if abs(value) > 0.55 * max_abs else "black"
            ax.text(
                j,
                i,
                format_cell_value(value),
                ha="center",
                va="center",
                fontsize=8.5,
                color=text_color,
            )

    # Minor gridlines around all cells.
    ax.set_xticks(np.arange(-0.5, matrix.shape[1], 1), minor=True)
    ax.set_yticks(np.arange(-0.5, matrix.shape[0], 1), minor=True)
    ax.grid(which="minor", color="black", linestyle="-", linewidth=0.6)
    ax.tick_params(which="minor", bottom=False, left=False)

    # Horizontal separators between model families.
    family_breaks = [3, 6, 9, 12]
    for y in family_breaks:
        ax.axhline(y - 0.5, color="black", linewidth=1.8)

    # Vertical separators between stage groups.
    if stage_group_size > 0:
        for x in range(stage_group_size, matrix.shape[1], stage_group_size):
            ax.axvline(x - 0.5, color="black", linewidth=1.8)

    cbar = fig.colorbar(image, ax=ax, fraction=0.035, pad=0.03)
    cbar.set_label(r"Signed percentage change: $100\times (O_{\mathrm{full}}-O_{X\mathrm{-off}})/O_{X\mathrm{-off}}$")

    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()

    rows = read_rows(args.input)

    masssq_matrix, models, masssq_columns = build_matrix(rows, MASSSQ_OBSERVABLES)
    mixing_matrix, _, mixing_columns = build_matrix(rows, MIXING_OBSERVABLES)

    masssq_csv = args.csv_dir / "oscillation_heatmap_masssq_table.csv"
    mixing_csv = args.csv_dir / "oscillation_heatmap_mixing_table.csv"

    write_matrix_csv(masssq_csv, masssq_matrix, models, masssq_columns)
    write_matrix_csv(mixing_csv, mixing_matrix, models, mixing_columns)

    masssq_plot = args.plot_dir / "oscillation_masssq_heatmap_percent.png"
    mixing_plot = args.plot_dir / "oscillation_mixing_heatmap_percent.png"

    plot_heatmap(
        masssq_matrix,
        models,
        masssq_columns,
        masssq_plot,
        "Signed percentage changes in mass-squared differences",
        stage_group_size=len(MASSSQ_OBSERVABLES),
    )

    plot_heatmap(
        mixing_matrix,
        models,
        mixing_columns,
        mixing_plot,
        "Signed percentage changes in mixing-angle observables",
        stage_group_size=len(MIXING_OBSERVABLES),
    )

    print("Created:")
    print(f"  {masssq_plot}")
    print(f"  {mixing_plot}")
    print(f"  {masssq_csv}")
    print(f"  {mixing_csv}")


if __name__ == "__main__":
    main()