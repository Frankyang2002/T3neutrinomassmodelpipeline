"""Unit tests for the all-16 orchestration's interface and serialization.

Run: python -X utf8 -m unittest discover -s tests -p test_run_general_t3_observable.py
These tests do not claim to validate RGBeta, Matchete or the physical coupling map.
"""
import argparse
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from scripts import RunGeneralT3Observable as m


class GeneralObservableHelperTests(unittest.TestCase):
    def test_model_catalogue_and_mixed_exceptions(self):
        self.assertEqual(len(m.MODELS), 16)
        self.assertEqual(set(m.MIXED), {"B_-1", "C_-1"})
        self.assertEqual(m.select_models(argparse.Namespace(model=None,
                              models="A_-2,B_-1,A_-2", all=False)), ("A_-2", "B_-1"))
        self.assertEqual(m.select_models(argparse.Namespace(model=None,
                              models=None, all=True)), tuple(m.MODELS))

    def test_component_csv_flatten_and_complex_encoding(self):
        arr = np.array([[1+2j, 0, 0], [0, 3, 0], [0, 0, 4]], dtype=complex)
        vals = list(m.flatten_parameter("MF", m.encode(arr)))
        self.assertEqual(len(vals), 9)
        self.assertEqual(vals[0], ("MF[1][1]", 1., 2.))

    def test_flatten_entries(self):
        rows = list(m.flatten_parameter("y", [[{"real": 1.0, "imag": 2.0}, 0], [0, 3]]))
        self.assertEqual(rows, [("y[1][1]", 1.0, 2.0),
                                ("y[1][2]", 0.0, 0.0),
                                ("y[2][1]", 0.0, 0.0),
                                ("y[2][2]", 3.0, 0.0)])

    def test_per_model_override_precedence(self):
        with tempfile.TemporaryDirectory() as temp:
            p = Path(temp) / "values.json"
            p.write_text(json.dumps({"shared": {"g3": 0.75, "y1": 1},
                "models": {"B_+1": {"g3": 0.76}}}), encoding="utf-8")
            self.assertEqual(m.user_overrides(p, "B_+1"), {"g3": 0.76, "y1": 1})
            self.assertEqual(m.user_overrides(p, "A_-2"), {"g3": 0.75, "y1": 1})

    def test_hierarchy_checks_running_not_initial_masses(self):
        class UV: pass
        class EFT: pass
        uv, eft = UV(), EFT()
        uv.MF = np.diag([1.5e5, 1.8e5, 2e5])
        eft.mS1Sq, eft.mS2Sq = (1300.0)**2, (1500.0)**2
        d = m.hierarchy_diagnostics(uv, eft, 1e5, 1e3)
        self.assertTrue(d["all_fermion_masses_above_mu_F"])
        self.assertTrue(d["all_unsplit_scalar_masses_above_mu_S"])
        uv.MF = np.diag([9e4, 1.8e5, 2e5])
        self.assertFalse(m.hierarchy_diagnostics(uv, eft, 1e5, 1e3)["all_fermion_masses_above_mu_F"])

    def test_json_record_preserves_explicit_false_certificates(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "status.json"
            m.json_save(path, {"model": "B_+1", "C5": m.encode(np.eye(3)*1j),
                "physical_matching_certified": False})
            result = m.json_load(path)
            self.assertIs(result["physical_matching_certified"], False)
            self.assertEqual(result["C5"][1][1], {"real": 0., "imag": 1.})


if __name__ == "__main__":
    unittest.main()
