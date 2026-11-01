"""Offline regression tests for the read-only full-general symbolic physics audit."""
from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest

from scripts import AuditFull16GeneralSymbolicPhysics as a


class GeneralSymbolicPhysicsAuditTests(unittest.TestCase):
    def test_selection_rules_and_group_factors(self) -> None:
        a.self_test()
        self.assertEqual(len(a.general_models()), 16)
        self.assertEqual(a.exact_group_invariants(3), (2, 2))
        with self.assertRaises(ValueError):
            a.check_monomial([0, 0, 1, 0, 0, 0], -2, 1)

    def test_matchete_rgbeta_structure_multiplicities(self) -> None:
        # Two independent contractions of the same four-field monomial:
        # one new RGBeta coupling, one structure already present in its base.
        counts = [1, 1, 0, 1, 0, 1]
        manifest = {"extra_quartic_structures": [
            {"counts": counts, "invariants": 2}],
            "extra_quartic_complex_couplings": 2}
        entry = {"total_additional_quartics": 2,
                 "new_rgbeta_symbols": [
                     {"name": "lambdaT3GeneralQ110101I1",
                      "counts": counts, "mapping_status": "NotNormalised"}],
                 "preexisting_rgbeta_symbols": [
                     {"existing_name": "lambdaHHdagS1barS2barCross", "counts": counts}]}
        interactions = ["T3ExtraQ_1_1_0_1_0_1_1", "T3ExtraQ_1_1_0_1_0_1_2"]
        rows, new, reused = a.declared_quartics("B_-1", entry, manifest, interactions)
        self.assertEqual((new, reused), (1, 1))
        self.assertEqual(len(rows), 1)
        entry["preexisting_rgbeta_symbols"] = []
        with self.assertRaises(ValueError):
            a.declared_quartics("B_-1", entry, manifest, interactions)

    def test_quadratic_mixed_mass_inventory(self) -> None:
        manifest = {"quadratic_structures_not_enabled": [
            {"counts": [0, 0, 0, 1, 0, 1], "invariants": 1}],
            "quadratic_complex_couplings_not_enabled": 1}
        interactions = ["T3ExtraM_0_0_0_1_0_1_1"]
        self.assertEqual(a.check_quadratics("B_-1", manifest, interactions, "mixed"), 1)
        with self.assertRaises(ValueError):
            a.check_quadratics("B_-1", manifest, [], "mixed")

    def test_beta_key_and_fermion_decoupling(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            entry = {"new_rgbeta_symbols": [{"name": "lambdaT3GeneralQ000004I1"}],
                     "preexisting_rgbeta_symbols": []}
            common = dict.fromkeys(a.EXPECTED_RGE_BASE, "1")
            uv = {**common, "y1": "1", "y2": "1", "MF": "1",
                  "y2Conj": "1", "lambdaT3GeneralQ000004I1": "1"}
            eft1 = {**common, "lambdaT3GeneralQ000004I1": "1"}
            for label, betas in (("uv", uv), ("eft1", eft1)):
                (folder / f"{label}_rgbeta.json").write_text(
                    json.dumps({"report_betas": betas}), encoding="utf-8")
            result = a.check_native_rges("A_-2", folder, entry, "y2Conj")
            self.assertEqual(result["UV"]["extra_yukawa_beta"], "y2Conj")
            eft1["MF"] = "1"
            (folder / "eft1_rgbeta.json").write_text(
                json.dumps({"report_betas": eft1}), encoding="utf-8")
            with self.assertRaises(ValueError):
                a.check_native_rges("A_-2", folder, entry, "y2Conj")

    def test_failed_is_never_certified(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            a.save_report(Path(tmp), [{"model": "A_-2", "status": a.FAILURE,
                                      "error": "missing data"}], [], [], ["A_-2"], Path(tmp))
            info = a.read(Path(tmp) / "audit.json")
            self.assertFalse(info["physical_matchete_rgbeta_coupling_dictionary_certified"])
            self.assertFalse(info["per_model_structural_checks_passed"])


if __name__ == "__main__":
    unittest.main()
