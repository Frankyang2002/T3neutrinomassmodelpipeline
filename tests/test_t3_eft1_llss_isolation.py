"""Algebraic unit tests for scripts/RunT3EFT1LLSSIsolation.py.

Synthetic tensor tests only; not end-to-end physical 16-model integration.
"""
import unittest
from types import SimpleNamespace
import numpy as np

from scripts.RunT3EFT1LLSSIsolation import (
    llss_to_llhh_quartic_source, make_llss_off_beta, summarize_model, CASES,
    OBSERVABLES, _build_heatmap_matrix, _rederive_heatmap_percentages, _render_heatmaps,
    make_cached_batch1_verifier, ALLOWED_CACHED_BATCH1_SOURCE_MISMATCH,
    _write_donor_report, _plot_style)


class TestSourceIsolation(unittest.TestCase):
    def setUp(self):
        self.basis = SimpleNamespace(nf=7, ns=6, higgs=(0, 1, 2, 3),
                                     leptons=((0, 1), (2, 3), (4, 5)))
        self.c = np.zeros((7, 7, 6, 6), complex)
        self.lam = np.zeros((6, 6, 6, 6), float)
        self.c[0, 2, 4, 5] = 3e-12
        self.c[2, 0, 4, 5] = 3e-12
        self.c[0, 2, 5, 4] = 3e-12
        self.c[2, 0, 5, 4] = 3e-12
        # lambda_[H0 H2 S4 S5] fully Bose symmetric
        from itertools import permutations
        for key in set(permutations((0, 2, 4, 5))):
            self.lam[key] = 0.1

    def test_source_contraction_and_support(self):
        result = llss_to_llhh_quartic_source(self.c, self.lam, self.basis)
        expected = 2 * 0.1 * 3e-12
        self.assertAlmostEqual(result[0, 2, 0, 2].real / expected, 1.0)
        self.assertAlmostEqual(result[2, 0, 2, 0].real / expected, 1.0)
        self.assertEqual(result[0, 2, 4, 5], 0)
        self.assertEqual(result[6, 0, 0, 2], 0)
        self.assertTrue(np.allclose(result, result.swapaxes(0, 1)))
        self.assertTrue(np.allclose(result, result.swapaxes(2, 3)))
        # Boundary LLHH alone must not generate an LLSS -> LLHH term.
        only_h = np.zeros_like(self.c)
        only_h[0, 0, 2, 2] = 1e-10
        self.assertEqual(np.max(abs(llss_to_llhh_quartic_source(only_h, self.lam, self.basis))), 0)

    def test_beta_isolation_leaves_every_other_term_intact(self):
        baseline = np.ones_like(self.c) * (2 + 1j)
        counter = {"calls": 0}

        def original(c, y, lam, gauge=(), *, variant=None):
            counter["calls"] += 1
            return {"scalar_pair": baseline.copy(),
                    "fermion_anomalous_dimension": 3 * baseline.copy()}, "fixture"

        capture = {}
        isolated = make_llss_off_beta(original, self.basis, capture=capture)
        result, tag = isolated(self.c, None, self.lam, variant="fixture")
        source = llss_to_llhh_quartic_source(self.c, self.lam, self.basis)
        self.assertTrue(np.allclose(result["scalar_pair"], baseline - source))
        self.assertTrue(np.array_equal(result["fermion_anomalous_dimension"], 3 * baseline))
        self.assertEqual(counter["calls"], 1)
        self.assertEqual(capture["rhs_calls"], 1)
        self.assertIn("LLSS_to_LLHH", tag)

    def test_three_way_attribution(self):
        # A clean synthetic example: all running is LLSS mixing, and the
        # finite scalar threshold is kept identical in all three cases.
        zz = np.eye(3, dtype=complex)
        data = {}
        for case in CASES:
            scale = {"W_FROZEN": 0, "LLSS_OFF": 0, "W_FULL": 4}[case]
            data[case] = {
                "matrices": {"eft1": scale * zz, "match": zz,
                             "smeft_ms": (scale + 1) * zz, "low": (scale + 1) * zz},
                "observables": {"dm21_sq_eV2": float((scale + 1)**2),
                                "dm31_sq_eV2": float((scale + 1)**2),
                                "sin2_theta12": 0.3, "sin2_theta13": 0.02,
                                "sin2_theta23": 0.5},
                "chi2": 0.0,
            }
        result, arrays = summarize_model(data)
        self.assertAlmostEqual(result["LLSS_fraction_norm"], 1)
        self.assertAlmostEqual(result["LLSS_all_Wilson_alignment_cosine"], 1)
        self.assertAlmostEqual(result["C5_LLSS_over_scalar_matching"], 4)
        self.assertAlmostEqual(result["dm21_sq_eV2_pct_vs_W_FROZEN"], 2400)
        self.assertAlmostEqual(result["dm21_sq_eV2_pct_LLSS_OFF_vs_W_FROZEN"], 0)
        self.assertAlmostEqual(result["dm21_sq_eV2_pct_LLSS_effect_vs_W_FROZEN"], 2400)
        self.assertTrue(np.array_equal(arrays["C5_effect_LLSS"], 4 * zz))


class TestHeatmapComparisons(unittest.TestCase):
    def _source_row(self, model):
        row = {"model": model, "status": "VALID_EXPLORATORY"}
        for obs in OBSERVABLES:
            row[f"{obs}_W_FROZEN"] = "2.0"
            row[f"{obs}_LLSS_OFF"] = "4.0"
            row[f"{obs}_W_FULL"] = "10.0"
            # Obsolete earlier-version percentage intentionally inconsistent:
            row[f"{obs}_pct_vs_LLSS_OFF"] = "9999"
        return row

    def test_common_frozen_denominator_and_columns(self):
        row = _rederive_heatmap_percentages([self._source_row("A_+0")])[0]
        obs = OBSERVABLES[0]
        self.assertEqual(row[f"{obs}_pct_vs_W_FROZEN"], 400.0)
        self.assertEqual(row[f"{obs}_pct_LLSS_OFF_vs_W_FROZEN"], 100.0)
        self.assertEqual(row[f"{obs}_pct_LLSS_effect_vs_W_FROZEN"], 300.0)
        a, models, labels = _build_heatmap_matrix([row], ("A_+0", "E_-4"))
        self.assertEqual(a.shape, (2, 10))
        self.assertEqual(models, ["A_+0", "E_-4"])
        self.assertAlmostEqual(a[0, 0], 400)
        self.assertAlmostEqual(a[0, 1], 100)
        self.assertTrue(np.isnan(a[1]).all())
        self.assertEqual(labels[0], "Running")
        self.assertEqual(labels[1], "Mixing off")

    def test_a4_and_dm2_reports_from_existing_observables(self):
        from pathlib import Path
        import tempfile
        import json
        converted = _rederive_heatmap_percentages([self._source_row("A_+0")])
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            _write_donor_report(out, converted, ("A_+0", "D_-2"), "A_-4", {})
            _write_donor_report(out, converted, ("A_+0", "D_-2"), "D_-2", {})
            a4 = out / "reports/eft1_llss_A_m4_report.md"
            d2 = out / "reports/eft1_llss_D_m2_report.md"
            self.assertTrue(a4.exists())
            self.assertTrue(d2.exists())
            self.assertIn("A(-4) donor", a4.read_text(encoding="utf-8"))
            self.assertIn("D(-2) donor", d2.read_text(encoding="utf-8"))
            self.assertIn("+400.000", a4.read_text(encoding="utf-8"))
            self.assertTrue((out / "reports/eft1_llss_D_m2_report.tex").exists())
            info = json.loads((out / "reports/eft1_llss_D_m2_summary.json").read_text())
            self.assertEqual(info["completed"], 1)
            self.assertEqual(info["blocked"], ["D_-2"])

    def test_default_matplotlib_font_and_style_validation(self):
        from types import SimpleNamespace
        opts = SimpleNamespace(axis_label_fontsize=13,
            bottom_tick_fontsize=9, top_tick_fontsize=12, y_tick_fontsize=11,
            cell_fontsize=7, colorbar_fontsize=10, figure_width=15,
            row_height=0.5, dpi=300)
        style = _plot_style(opts)
        self.assertNotIn("font_family", style)
        self.assertNotIn("math_fontset", style)
        self.assertEqual(style["cell_fontsize"], 7)
        opts.row_height = 0
        with self.assertRaisesRegex(ValueError, "row_height"):
            _plot_style(opts)

    def test_zero_denominator_blanks_both_columns(self):
        row = self._source_row("A_+0")
        row[f"{OBSERVABLES[0]}_W_FROZEN"] = "0.0"
        converted = _rederive_heatmap_percentages([row])[0]
        self.assertIsNone(converted[f"{OBSERVABLES[0]}_pct_vs_W_FROZEN"])
        self.assertIsNone(converted[f"{OBSERVABLES[0]}_pct_LLSS_OFF_vs_W_FROZEN"])
        array, _, _ = _build_heatmap_matrix([converted], ("A_+0",))
        self.assertTrue(np.isnan(array[0, :2]).all())

    def test_generates_heatmap_files(self):
        from pathlib import Path
        import tempfile
        converted = _rederive_heatmap_percentages([self._source_row("A_+0")])
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            _render_heatmaps(out, converted, ("A_+0", "E_-4"), "all16")
            self.assertTrue((out / "plots/all16_heatmap_percent.png").is_file())
            self.assertTrue((out / "plots/all16_heatmap_percent.pdf").is_file())
            self.assertTrue((out / "plots/all16_masssq_heatmap_percent.png").is_file())
            self.assertTrue((out / "plots/all16_mixing_heatmap_percent.png").is_file())
            self.assertTrue((out / "csv/all16_heatmap_table.csv").is_file())




class TestCachedBatch1Reuse(unittest.TestCase):
    def test_known_source_mismatch_accepted_and_recorded(self):
        calls = []
        known = ALLOWED_CACHED_BATCH1_SOURCE_MISMATCH[0]
        def fake_verifier(*args, strict_sources=True):
            calls.append(strict_sources)
            if strict_sources:
                raise ValueError("old source hash mismatch")
            # stand-in for repository's retained cached beta export validation
            return {"local_sources": {
                known: {"sha256": "newhash", "batch1_sha256": "oldhash",
                        "same_as_batch1": False, "required_for_B_reuse": True},
            }, "uv_export_sha256": "uv_validated", "eft1_export_sha256": "eft_validated"}
        audit = {}
        verifier = make_cached_batch1_verifier(fake_verifier, audit)
        response = verifier(None, None, "A_-4", {}, None, None)
        self.assertEqual(calls, [False])
        self.assertTrue(audit["cached_batch1_export_hashes_verified"])
        self.assertTrue(audit["batch1_source_mismatch_accepted_for_exploratory_use"])
        self.assertEqual(audit["batch1_source_mismatches"][known]["current_sha256"], "newhash")
        self.assertEqual(response["uv_export_sha256"], "uv_validated")

    def test_other_required_source_mismatch_blocked(self):
        def fake_verifier(*args, strict_sources=True):
            return {"local_sources": {
                "Lagrangian/model/T3ModelCatalog.wl": {"sha256": "a", "batch1_sha256": "b",
                     "same_as_batch1": False, "required_for_B_reuse": True}
            }}
        verifier = make_cached_batch1_verifier(fake_verifier, {})
        with self.assertRaisesRegex(ValueError, "unapproved Batch 1 source mismatches"):
            verifier(None, None, "A_-4", {}, None, None)

    def test_cached_beta_verification_failure_not_suppressed(self):
        def fake_verifier(*args, strict_sources=True):
            raise ValueError("Batch 1 cached uv JSON has changed")
        verifier = make_cached_batch1_verifier(fake_verifier, {})
        with self.assertRaisesRegex(ValueError, "cached uv JSON has changed"):
            verifier(None, None, "A_-4", {}, None, None)

if __name__ == "__main__":
    unittest.main()
