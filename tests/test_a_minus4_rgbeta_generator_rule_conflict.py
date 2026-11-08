"""Offline validator tests; real native Wolfram execution is Windows-only."""
from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "InspectAminus4RGBetaGeneratorRuleConflict.py"
SPEC = importlib.util.spec_from_file_location("a4_rule_validator", SCRIPT)
assert SPEC and SPEC.loader
v = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(v)


def fixture(status="NATIVE_SAME_ORIENTATION_FIERZ_PATH_DISAGREEMENT"):
    conflict = status == "NATIVE_SAME_ORIENTATION_FIERZ_PATH_DISAGREEMENT"
    return {
        "format": v.FORMAT, "model": "A_-4", "stage": "EFT1",
        "test_kind": "same-orientation generator contraction; no beta functions",
        "status": status,
        "fund_direct_reduction_inputform": "-3*del[...]/4",
        "fund_four_distinct_reduction_inputform": "(...)",
        "fund_fierz_then_contract_reduction_inputform": "del[...]/4",
        "fund_direct_minus_fierz_path_inputform": "-del[...]" if conflict else "0",
        "native_direct_matches_generic": conflict,
        "native_direct_matches_fierz": not conflict,
        "native_fierz_path_matches_reference": True,
        "native_two_paths_disagree": conflict,
        "physical_scalar_beta_error_proven": False,
        "physical_matchete_coupling_map_certified": False,
        "source_modified": False,
    }


class TestNativeRuleGate(unittest.TestCase):
    def test_competing_source_rule_hash_is_pinned(self):
        self.assertEqual(len(v.EXPECTED_GROUPS_SHA256), 64)

    def test_detected_conflict_still_not_beta_certified(self):
        r = v.classify_native(fixture())
        self.assertEqual(r["native_algebra_status"],
                         "NATIVE_GROUP_REDUCTION_CONFLICT_OBSERVED")
        self.assertFalse(r["beta_function_causality_established"])

    def test_agree_paths(self):
        r = v.classify_native(fixture("NATIVE_SAME_ORIENTATION_PATHS_AGREE"))
        self.assertEqual(r["native_algebra_status"],
                         "NATIVE_PATHS_AGREE_FOR_TESTED_CONTRACTION")

    def test_unresolved_paths(self):
        d = fixture("NATIVE_SAME_ORIENTATION_PATH_COMPARISON_UNRESOLVED")
        d["native_direct_matches_fierz"] = False
        self.assertEqual(v.classify_native(d)["native_algebra_status"],
                         "INCONCLUSIVE_NATIVE_REDUCTION")

    def test_reject_fabricated_conflict(self):
        d = fixture()
        d["native_fierz_path_matches_reference"] = False
        with self.assertRaises(v.AuditBlocked):
            v.classify_native(d)

    def test_reject_source_modified(self):
        d = fixture()
        d["source_modified"] = True
        with self.assertRaises(v.AuditBlocked):
            v.classify_native(d)

    def test_reject_fabricated_certification(self):
        d = fixture()
        d["physical_scalar_beta_error_proven"] = True
        with self.assertRaises(v.AuditBlocked):
            v.classify_native(d)

    def test_reject_wrong_model(self):
        d = fixture()
        d["model"] = "A_-2"
        with self.assertRaises(v.AuditBlocked):
            v.classify_native(d)

    def test_local_source_and_native_snapshot_must_match(self):
        with tempfile.TemporaryDirectory() as t:
            p = Path(t)
            model = p / "T3RGBetaModel.wl"
            groups = p / "GroupsAndIndices.m"
            native = p / "witness.json"
            model.write_text("physics untouched", encoding="utf-8")
            groups.write_text("installed package", encoding="utf-8")
            d = fixture()
            d["model_source_sha256"] = hashlib.sha256(model.read_bytes()).hexdigest()
            d["installed_groups_source_sha256"] = hashlib.sha256(groups.read_bytes()).hexdigest()
            native.write_text(json.dumps(d), encoding="utf-8")
            with self.assertRaisesRegex(v.AuditBlocked, "differs from uploaded"):
                v.inspect(native, model, groups)

    def test_packaged_wolfram_script_is_read_only(self):
        wl = (ROOT / "scripts" / "RunAminus4RGBetaGeneratorRuleConflict.wl").read_text()
        for token in ("T3RGBetaBuildEFT1[1,3,2,-4]", "RefineGroupStructures",
                      "native_two_paths_disagree", "physical_scalar_beta_error_proven"):
            self.assertIn(token, wl)
        self.assertNotIn("SetOptions[", wl)
        self.assertNotIn("Unprotect[", wl)


if __name__ == "__main__":
    unittest.main()
