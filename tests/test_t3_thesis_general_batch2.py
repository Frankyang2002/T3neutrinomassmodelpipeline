"""Focused Batch-2 source/interface tests; no Wolfram or physics certification."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from Numerical.parallel_t3.background.general_batch1_api import encode_native
from Numerical.parallel_t3.wilson.general_batch2_source import (
    BLOCKERS, D_QUARTICS, MatchingIncomplete, native_numbers, verify_ledger,
)


def fixture(model):
    out = {
        "MF": encode_native(np.diag([105426., 115968., 126511.]).astype(complex)),
        "y1": encode_native(np.eye(3) * 1e-4),
        "y2": encode_native(np.eye(3) * 2e-4),
        "mS1Sq": 1.8e6,
        "mS2Sq": 1.0e6,
        "lambdaT3": 0.1,
    }
    if model == "D_-2":
        out.update({"y2Conj": encode_native(np.zeros((3, 3), complex)),
                    "muT3GeneralS2Sq": encode_native(-4338.3 + 1e-23j)})
        out.update({n: 0.01 for n in D_QUARTICS})
    return out


class InterfaceTests(unittest.TestCase):
    def test_a_has_eight_components_no_extra_terms(self):
        got = native_numbers(fixture("A_-4"), "A_-4")
        self.assertEqual(len(got["candidate_scalar_quadratic_mass_eigenvalues_GeV2"]), 8)
        self.assertEqual(got["neutral_scalar_self_bilinear_GeV2"], {"real": 0.0, "imag": 0.0})
        self.assertFalse(got["matchete_to_rgbeta_dictionary_applied"])

    def test_d_keeps_all_five_and_generated_bilinear(self):
        got = native_numbers(fixture("D_-2"), "D_-2")
        self.assertEqual(set(D_QUARTICS).issubset(got["native_parameter_names"]), True)
        eigen = got["candidate_scalar_quadratic_mass_eigenvalues_GeV2"]
        self.assertAlmostEqual(eigen[0], 1e6 - 8676.6)
        self.assertAlmostEqual(eigen[1], 1e6 + 8676.6)
        self.assertFalse(got["source_exact_one_loop_wilson_numbers_produced"])

    def test_d_missing_extra_quartic_blocks(self):
        bad = fixture("D_-2")
        del bad[next(iter(D_QUARTICS))]
        with self.assertRaises(MatchingIncomplete):
            native_numbers(bad, "D_-2")

    def test_d_complex_mass_phase_blocks(self):
        bad = fixture("D_-2")
        bad["muT3GeneralS2Sq"] = {"real": -4338.3, "imag": 50.0}
        with self.assertRaisesRegex(MatchingIncomplete, "complex phase"):
            native_numbers(bad, "D_-2")

    def test_d_y2conj_cannot_disappear(self):
        bad = fixture("D_-2")
        del bad["y2Conj"]
        with self.assertRaises(MatchingIncomplete):
            native_numbers(bad, "D_-2")

    def test_full_raw_ledger_requires_exact_source_hash(self):
        with tempfile.TemporaryDirectory() as d:
            source = Path(d) / "one_loop.inputform.txt"
            source.write_text("hbar * Coupling[y1]", encoding="utf-8")
            right_hash = hashlib.sha256(source.read_bytes()).hexdigest()
            result = {"format": "t3-batch2-complete-loop-operator-ledger-v1",
                      "status": "UnrenormalisedOneLoopTermsEnumerated",
                      "source_sha256": right_hash,
                      "term_count": 1,
                      "terms": [{"raw_inputform": "hbar * Coupling[y1]"}],
                      "sector_counts": {"OtherEFT1Operator": 1},
                      "exact_reconstruction_by_expansion": True,
                      "only_fermions_integrated_out": True,
                      "finite_msbar_wilson_coefficients_computed": False}
            self.assertEqual(verify_ledger(result, source)["term_count"], 1)
            result["source_sha256"] = "a" * 64
            with self.assertRaisesRegex(MatchingIncomplete, "SHA"):
                verify_ledger(result, source)

    def test_missing_complete_numerical_physics_is_explicit(self):
        self.assertGreaterEqual(len(BLOCKERS), 5)
        self.assertTrue(any("operator mixing" in b for b in BLOCKERS))


if __name__ == "__main__":
    unittest.main()
