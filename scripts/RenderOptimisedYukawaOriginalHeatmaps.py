"""Render optimised Yukawa heatmaps using the EXISTING repository report creator.

This adapter only reads saved 16-case results and calls the unchanged
scripts.RunT3SimultaneousEFT1StageFreezeHeatmaps functions. It suppresses
plot titles and exposes font sizes through a temporary Matplotlib style hook;
no independent plotting routine or physics calculation is introduced.

PowerShell from repository root:
  python -X utf8 -m scripts.RenderOptimisedYukawaOriginalHeatmaps `
      --run-id optimal_yukawa_thesis_heatmaps_3stage_v1 `
      --font-size 11 --cell-font-size 9 --stage-font-size 11 `
      --axis-font-size 12 --tick-font-size 10 --overwrite

Case definitions: 1111 full; 0111 UV off; 1001 EFT1+Wilson off; 1110 SMEFT off.
"""
from __future__ import annotations

import argparse
from contextlib import contextmanager
import json
import math
from pathlib import Path
from unittest.mock import patch

from scripts import RunT3SimultaneousEFT1StageFreezeHeatmaps as original

ROOT = Path(__file__).resolve().parents[1]
DATA_ROOT = (ROOT / "output/studies/general_scalar_extended"
             / "optimal_yukawa_thesis_refresh")
REFS = {"D_-2": "Dm2", "A_-4": "Am4"}


@contextmanager
def _original_plot_style(*, font_size: float, cell_font_size: float,
                         stage_font_size: float, axis_font_size: float,
                         tick_font_size: float, colorbar_font_size: float):
    """Apply typography to original renderer for this call only.

    The original heatmap function sets title and cell/stage annotation font
    sizes explicitly. Temporarily intercept those calls; everything else
    (matrices, Symlog colours, annotations, percentage formula, output names)
    remains the original code. Restore Matplotlib methods on exit.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.axes import Axes
    from matplotlib.colorbar import Colorbar

    set_title = Axes.set_title
    text = Axes.text
    set_cbar_label = Colorbar.set_label
    update_cbar_ticks = Colorbar.update_ticks

    def without_title(self, _label, *args, **kwargs):
        kwargs.pop("pad", None)
        return set_title(self, "", *args, pad=0, **kwargs)

    def styled_text(self, x, y, s, *args, **kwargs):
        if s == "EFT1 off":
            s = "EFT1"
        specified = kwargs.get("fontsize", kwargs.get("size"))
        if specified is not None:
            try:
                original_size = float(specified)
            except (TypeError, ValueError):
                original_size = None
            if original_size is not None:
                if math.isclose(original_size, 8.0):
                    new_size = cell_font_size
                elif math.isclose(original_size, 10.0):
                    new_size = stage_font_size
                else:
                    new_size = original_size
                kwargs.pop("size", None)
                kwargs["fontsize"] = new_size
        return text(self, x, y, s, *args, **kwargs)

    def styled_colorbar_label(self, label, *args, **kwargs):
        kwargs["fontsize"] = axis_font_size
        return set_cbar_label(self, label, *args, **kwargs)

    def styled_colorbar_ticks(self, *args, **kwargs):
        value = update_cbar_ticks(self, *args, **kwargs)
        self.ax.tick_params(labelsize=colorbar_font_size)
        return value

    style = {
        "font.size": font_size,
        "axes.titlesize": font_size,
        "axes.labelsize": axis_font_size,
        "xtick.labelsize": tick_font_size,
        "ytick.labelsize": tick_font_size,
        "legend.fontsize": font_size,
    }
    with plt.rc_context(style), patch.object(Axes, "set_title", without_title), \
            patch.object(Axes, "text", styled_text), \
            patch.object(Colorbar, "set_label", styled_colorbar_label), \
            patch.object(Colorbar, "update_ticks", styled_colorbar_ticks):
        yield


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--font-size", type=float, default=11.0,
                        help="Base font size, points (default: 11)")
    parser.add_argument("--cell-font-size", type=float, default=8.0,
                        help="Numbers printed inside cells, points (default: 8)")
    parser.add_argument("--stage-font-size", type=float, default=10.0,
                        help="UV off / EFT1 off / SMEFT off headings, points (default: 10)")
    parser.add_argument("--axis-font-size", type=float, default=11.0,
                        help="Axis titles and colourbar title, points (default: 11)")
    parser.add_argument("--tick-font-size", type=float, default=10.0,
                        help="Row/column/colourbar tick labels, points (default: 10)")
    parser.add_argument("--colorbar-font-size", type=float, default=10.0,
                        help="Colourbar tick labels, points (default: 10)")
    args = parser.parse_args(argv)
    size_names = ("font_size", "cell_font_size", "stage_font_size",
                  "axis_font_size", "tick_font_size", "colorbar_font_size")
    for name in size_names:
        value = getattr(args, name)
        if not math.isfinite(value) or not 4 <= value <= 36:
            parser.error(f"--{name.replace('_', '-')} must be between 4 and 36 pt")
    if not args.run_id or not all(c.isalnum() or c in "_-" for c in args.run_id):
        parser.error("Unsafe run ID")

    source = DATA_ROOT / args.run_id
    result_path = source / "results.json"
    manifest_path = source / "manifest.json"
    if not (result_path.is_file() and manifest_path.is_file()):
        parser.error("Missing completed physics results or run manifest: " + str(source))
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("full_cases") is not True:
        parser.error("Source run has only 1111; rerun physics with --full-cases")
    results = json.loads(result_path.read_text(encoding="utf-8"))
    output = args.output or (source / "original_heatmap_reports")
    if output.exists() and any(output.iterdir()) and not args.overwrite:
        parser.error("Output directory exists; use --overwrite or a fresh --output")
    if not original.ORIGINAL_PLOTTER.is_file():
        parser.error("Original plotting module is missing: " + str(original.ORIGINAL_PLOTTER))

    # Validate both donors and all 16 models before writing any report.
    source_data = {}
    for ref in REFS:
        models = results.get(ref, {})
        ref_cases = {}
        for model in original.MODELS:
            item = models.get(model)
            if not isinstance(item, dict) or item.get("status") != "Computed":
                parser.error(f"{ref}/{model} has no completed transfer result")
            cases = item.get("cases")
            if not isinstance(cases, dict) or len(cases) != 16:
                parser.error(f"{ref}/{model} does not have 16 U/E/W/S cases")
            if any(case not in cases for case in ("1111", "0111", "1001", "1110")):
                parser.error(f"{ref}/{model} missing required three-stage case")
            ref_cases[model] = cases
        source_data[ref] = ref_cases

    font_options = {name: getattr(args, name) for name in size_names}
    for ref, short in REFS.items():
        dest = output / short
        all_long, all_wide, missing = [], [], []
        for model in original.MODELS:
            long_rows, wide_rows, errors = original._records_for_model(
                model, source_data[ref][model])
            all_long.extend(long_rows)
            all_wide.extend(wide_rows)
            missing.extend(errors)
        dest.mkdir(parents=True, exist_ok=True)
        csv_path = dest / "frozen_stage_observable_changes.csv"
        original.write_csv(csv_path, all_long, original.LONG_FIELDS)
        original.write_csv(dest / "running_comparisons.csv", all_wide, original.WIDE_FIELDS)
        original.write_json(dest / "missing_cells.json", missing)
        original.write_json(dest / "provenance.json", {
            "source_physics_results": str(result_path),
            "fixed_yukawa_reference": ref,
            "original_report_generator": "scripts/RunT3SimultaneousEFT1StageFreezeHeatmaps.py",
            "original_matrix_plotter": "scripts/PlotT3OscillationHeatmaps.py",
            "reference_case": "1111",
            "stages": original.STAGES,
            "percentage_formula": "100*(O_1111-O_off)/O_off",
            "matching_certified": False,
            "missing_cells": len(missing),
            "no_new_physics_evaluation": True,
            "plot_title_removed": True,
            "font_sizes_pt": font_options,
        })
        with _original_plot_style(**font_options):
            original.render_original_heatmap_layout(csv_path, dest)
        print(ref, "original heatmaps (title removed):", dest / "plots", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
