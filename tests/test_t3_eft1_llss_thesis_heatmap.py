"""Tests for scripts/RunT3EFT1LLSSThesisHeatmap.py (no expensive physics)."""
from __future__ import annotations
import csv
import importlib.util
import math
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

FILE = Path(__file__).resolve().parents[1] / "scripts/RunT3EFT1LLSSThesisHeatmap.py"
spec = importlib.util.spec_from_file_location("t3_thesis_llss", FILE)
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)


class ThesisConsistent(unittest.TestCase):
    def test_frozen_is_common_denominator(self):
        baseline = {o: (10.0, 2.0, 400.0) for o in module.OBSERVABLES}
        off = {o: 3.0 for o in module.OBSERVABLES}
        row = module.make_row("A_-4", "A_-4", baseline, off, source="thesis.csv")
        for o in module.OBSERVABLES:
            self.assertEqual(row[f"{o}_Running_pct"], 400.0)
            self.assertEqual(row[f"{o}_Mixing_off_pct"], 50.0)
            self.assertEqual(row[f"{o}_LLSS_effect_pct"], 350.0)

    def test_zero_denominator_rejected(self):
        self.assertIsNone(module.percentage(3.0, 0.0))
        self.assertIsNone(module.percentage(float("nan"), 2.0))

    def test_thesis_percentage_file_is_verified(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "stage.csv"
            row = {"reference": "A_-4", "model": "D_-2",
                   "stage_frozen": "EFT1", "status": "VALID_EXPLORATORY"}
            for obs in module.OBSERVABLES:
                row[f"{obs}_full"] = "10"
                row[f"{obs}_frozen"] = "2"
                row[f"{obs}_pct"] = "400"
            with path.open("w", encoding="utf-8", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=row)
                writer.writeheader(); writer.writerow(row)
            got = module.stagefreeze_lookup(path)
            self.assertEqual(got[("A_-4", "D_-2")][module.OBSERVABLES[0]], (10., 2., 400.))
            row[f"{module.OBSERVABLES[0]}_pct"] = "401"
            with path.open("w", encoding="utf-8", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=row)
                writer.writeheader(); writer.writerow(row)
            with self.assertRaisesRegex(ValueError, "Corrupt thesis percentage"):
                module.stagefreeze_lookup(path)

    def test_font_controls(self):
        names = ("axis_label_fontsize", "bottom_tick_fontsize", "top_tick_fontsize",
                 "y_tick_fontsize", "cell_fontsize", "colorbar_fontsize")
        args = {name: None for name in names}
        args.update(font_size=16, figure_width=23, row_height=.6, dpi=300, cell_fontsize=12)
        got = module.style_from_args(SimpleNamespace(**args))
        self.assertEqual(got["cell_fontsize"], 12)
        self.assertEqual(got["top_tick_fontsize"], 18)
        self.assertEqual(got["bottom_tick_fontsize"], 16)

    def test_plot_16_by_10_and_report(self):
        from PIL import Image
        baseline = {o: (10.0, 2.0, 400.0) for o in module.OBSERVABLES}
        rows = [module.make_row("A_-4", f"M_{i}", baseline,
                                {o: 3.0 for o in module.OBSERVABLES}, source="thesis.csv")
                for i in range(16)]
        with tempfile.TemporaryDirectory() as d:
            out = Path(d)
            # Import available production plot renderer through standard 'scripts' package.
            style = {"axis_label_fontsize": 11, "top_tick_fontsize": 11,
                     "bottom_tick_fontsize": 10, "y_tick_fontsize": 10,
                     "cell_fontsize": 8, "colorbar_fontsize": 10,
                     "figure_width": 16, "row_height": .5, "dpi": 90}
            module.plot_three_way(rows, tuple(f"M_{i}" for i in range(16)), out, style)
            module.write_reports(out, "A_-4", rows, tuple(f"M_{i}" for i in range(16)))
            image = out / "plots/all16_heatmap_percent.png"
            self.assertTrue(image.is_file())
            with Image.open(image) as im:
                self.assertGreater(im.width, 1000)
            self.assertTrue(image.with_suffix(".pdf").is_file())
            self.assertTrue((out / "reports/thesis_consistent_LLSS_report.md").is_file())


if __name__ == "__main__":
    unittest.main()
