"""Integration contracts, not a physical validation of the B/C sign dictionary."""
from __future__ import annotations

import json
import math
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest

import numpy as np

from Numerical.matching.GeneralT3MixedFlavourC5 import (
    DIRECT_RGBETA_TO_MAT, GeneralT3MixedFlavourC5, _source_files,
    evaluate_ordered_coefficient, map_direct_quartics, validate_source,
)


class TestMixedConditionalC5(unittest.TestCase):
    def test_model_sources(self):
        self.assertEqual(_source_files("B_-1", Path("source"))[0].name,
                         "B_m1_final_one_loop_probe_weinberg_ordered_coefficient.inputform.txt")
        self.assertEqual(_source_files("C_-1", Path("source"))[0].name,
                         "C_m1_final_one_loop_probe_weinberg_ordered_coefficient.inputform.txt")
        with self.assertRaises(ValueError):
            _source_files("B_+1", Path("source"))

    def test_explicit_uncertified_signed_map(self):
        state = SimpleNamespace(lambdaHHdagS2S2=0.1+0.2j,
                                lambdaHHdagS1barS1bar=-0.03)
        out = map_direct_quartics(state)
        self.assertEqual(out[DIRECT_RGBETA_TO_MAT["lambdaHHdagS2S2"]], -0.1-0.2j)
        self.assertEqual(out[DIRECT_RGBETA_TO_MAT["lambdaHHdagS1barS1bar"]], 0.03)
        with self.assertRaises(ValueError):
            map_direct_quartics(SimpleNamespace(lambdaHHdagS2S2=0.0))

    def test_provenance_gate_and_false_certificates(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            coefficient, probe = _source_files("B_-1", root)
            coefficient.write_text("test-only", encoding="utf-8")
            report = {
                "model":"B_-1",
                "format":"t3-matchete-native-heavy-mass-eigenbasis-v1",
                "status":"OneLoopDirectSymbolicWeinbergProbePassedNumericalC5Pending",
                "one_loop_match_stage_completed":True,
                "weinberg_extraction_status":"Success",
                "weinberg_holomorphic_term_count":118,
                "ordered_coefficient_is_canonical_flavour_symmetric_c5":False,
                "physical_matching_performed":False,
                "production_sequential_matching_performed":False,
                "component_basis_map_certified":False,
                "matchete_rgbeta_quartic_normalization_certified":False,
                "weinberg_ordered_coefficient_inputform_file":coefficient.name,
                "mass_basis_ordinary_check_passed":True,
                "heavy_mass_basis_detailed_check_passed":True,
            }
            probe.write_text(json.dumps(report), encoding="utf-8")
            self.assertEqual(len(validate_source("B_-1", coefficient, probe)["coefficient_sha256"]),64)
            obj=GeneralT3MixedFlavourC5("B_-1", root)
            self.assertEqual(obj.coefficient_path,coefficient)
            self.assertAlmostEqual(obj.hbar_value, 1/(16*math.pi**2))
            report["physical_matching_performed"]=True
            probe.write_text(json.dumps(report), encoding="utf-8")
            with self.assertRaises(ValueError):
                validate_source("B_-1", coefficient, probe)

    def test_flavour_and_majorana_symmetrisation(self):
        # Synthetic expression tests A5_ij+A5_ji; it is NOT a Matchete source.
        formula=compile("hbar*lambdaT3*y1a*y2b/MF", "<synthetic>", "eval")
        basis=SimpleNamespace(U=np.eye(2,dtype=complex),
                              eigenvalues_sq=np.array([1.7e6,2.0e6]),
                              tachyonic=False,degenerate=False)
        mass=np.array([150000.,180000.,200000.])
        y1=np.diag([4e-4,3e-4,2e-4]).astype(complex)
        y2=np.diag([1e-4,2e-4,5e-4]).astype(complex)
        q={key:0j for key in DIRECT_RGBETA_TO_MAT.values()}
        out=evaluate_ordered_coefficient(formula,mass,y1,y2,basis,0.02,q)
        expected=2*0.02/(16*math.pi**2)*np.diag(y1).real*np.diag(y2).real/mass
        np.testing.assert_allclose(out,np.diag(expected),rtol=2e-14,atol=0)
        self.assertTrue(np.allclose(out,out.T,atol=0))

    def test_archived_sources_if_present(self):
        from Numerical.matching.GeneralT3MixedFlavourC5 import DEFAULT_SOURCE
        for model in ("B_-1", "C_-1"):
            coefficient, probe = _source_files(model, DEFAULT_SOURCE)
            if not coefficient.is_file() or not probe.is_file():
                continue  # Standalone unit-test fixture does not bundle the repository data.
            validated = validate_source(model, coefficient, probe)
            self.assertEqual(len(validated["probe_sha256"]), 64)

    def test_scalar_mass_guard(self):
        basis=SimpleNamespace(U=np.eye(2,dtype=complex),
                              eigenvalues_sq=np.array([1.0,1.0]),
                              tachyonic=False,degenerate=True)
        with self.assertRaises(ValueError):
            evaluate_ordered_coefficient(compile("0","<test>","eval"),
                np.array([1.,2.,3.]),np.eye(3),np.eye(3),basis,0.02,
                {x:0j for x in DIRECT_RGBETA_TO_MAT.values()})


if __name__ == "__main__":
    unittest.main()
