"""Unit checks for the two-donor plot-only wrapper (no physics reruns)."""
from __future__ import annotations

import csv
import json
import tempfile
import unittest
from argparse import Namespace
from pathlib import Path
from unittest.mock import patch

from scripts import RunT3EFT1LLSSBothDonors as both


class TestBothDonorReplot(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def make_case(self, donor: str):
        source = self.root / donor
        source.mkdir()
        fieldnames = ["reference", "model", "status"] + [
            f"{obs}_{case}"
            for obs in ("dm21_sq_eV2", "dm31_sq_eV2", "sin2_theta12", "sin2_theta13", "sin2_theta23")
            for case in ("W_FULL", "W_FROZEN", "LLSS_OFF")
        ]
        with (source / "all16_three_way.csv").open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for model in both.MODELS:
                row = {name: "0.5" for name in fieldnames}
                row.update(reference=donor, model=model, status="VALID_EXPLORATORY")
                writer.writerow(row)
        (source / "provenance.json").write_text(json.dumps({"reference": donor}))
        return source

    def test_validate_both_donors(self):
        for donor in ("A_-4", "D_-2"):
            source = self.make_case(donor)
            self.assertEqual(both._validate_saved(source, donor)["models_valid"], 16)

    def test_invalid_donor_blocked(self):
        source = self.make_case("A_-4")
        with self.assertRaises(ValueError):
            both._validate_saved(source, "D_-2")

    def test_default_font_sizes_are_larger(self):
        args = Namespace(font_size=14.0, figure_width=21.0, row_height=0.67, dpi=300,
                         **{name: None for name in both.FONT_FACTORS})
        items = both._style_options(args)
        flags = dict(zip(items[::2], items[1::2]))
        self.assertEqual(float(flags["--cell-fontsize"]), 12.6)
        self.assertEqual(float(flags["--bottom-tick-fontsize"]), 14.0)
        self.assertEqual(float(flags["--figure-width"]), 21.0)
        args.cell_fontsize = 17.0
        flags = dict(zip(both._style_options(args)[::2], both._style_options(args)[1::2]))
        self.assertEqual(float(flags["--cell-fontsize"]), 17.0)

    def test_default_results_path(self):
        self.assertIn("eft1_llss_both_donors_end", str(both.PUBLIC_RESULTS))
        self.assertEqual([x[0] for x in both.DONORS], ["A_-4", "D_-2"])


if __name__ == "__main__":
    unittest.main()
