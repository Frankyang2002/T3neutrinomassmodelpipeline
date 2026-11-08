"""No-Wolfram regression tests for A(-4) native operator capture gate.

The tests use synthetic fixtures. They do not claim that the Windows native
RGBeta runtime has been executed or its component basis certified.
"""
from __future__ import annotations

import importlib.util
from hashlib import sha256
from pathlib import Path
import tempfile
import unittest

SCRIPT = (Path(__file__).resolve().parents[1] / "scripts" /
          "CheckAminus4RGBetaOperatorWitness.py")
spec = importlib.util.spec_from_file_location("opcheck", SCRIPT)
mod = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(mod)


INVARIANTS = {
    "lambdaH": "del[SU2L[fund], opA, opB]*del[SU2L[fund], opC, opD]/2",
    "lambdaS1": "1/2",
    "lambdaS2": "del[SU2L[S2], opA, opB]*del[SU2L[S2], opC, opD]/2",
    "lambdaH1": "del[SU2L[fund], opA, opB]",
    "lambdaH2": "del[SU2L[fund], opA, opB]*del[SU2L[S2], opC, opD]",
    "lambda12": "del[SU2L[S2], opC, opD]",
    "lambdaT3": "delS2[SU2L, opD, opA, opB]",
    "lambdaH2Adj": "tGen[SU2L[fund], A, opA, opB]*tGen[SU2L[S2], A, opC, opD]",
    "lambdaS2Adj": "tGen[SU2L[S2], A, opA, opB]*tGen[SU2L[S2], A, opC, opD]/2",
}


class OperatorWitnessTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        source = self.root / mod.MODEL_SOURCE
        source.parent.mkdir(parents=True)
        source.write_bytes(b"(* tracked native model source *)")
        self.hash = sha256(source.read_bytes()).hexdigest()

    def fixture(self):
        rows = []
        for name, fields in mod.QUARTIC_FIELDS.items():
            rows.append({
                "name": name,
                "registration_complete": True,
                "fields": fields,
                "invariant_inputform": INVARIANTS[name],
                "coupling_inputform": f"6*{name}",
                "coupling_bar_inputform":
                    "6*Bar[lambdaT3]" if name == "lambdaT3" else f"6*{name}",
                "self_conjugate": name != "lambdaT3",
                "native_unique_arrangements": 12,
            })
        return {
            "format": mod.FORMAT, "status": mod.STATUS,
            "model": "A_-4", "stage": "EFT1",
            "dimensions_alpha": [1, 3, 2, -4],
            "source_model_path": mod.MODEL_SOURCE.as_posix(),
            "source_model_sha256": self.hash,
            "rgbeta_package_file": "C:/Wolfram/Applications/RGBeta/Kernel/init.m",
            "native_quartic_symbols": sorted(mod.QUARTIC_FIELDS),
            "quartic_operators": rows,
            "scalar_species_indices": {
                "H": 1, "Bar[H]": 2,
                "S1": 3, "Bar[S1]": 4,
                "S2": 5, "Bar[S2]": 6,
            },
            "native_Lam_symbol_names": ["RGBeta`PackageScope`Lam"],
            "native_Lam_formal_tensor_available": True,
            "native_Lam_visible_quartic_symbols": sorted(mod.QUARTIC_FIELDS),
            "native_Lam_formal_tensor_inputform":
                "TStructure[$scalar[opA], $scalar[opB]][SparseArray[Automatic,{8,8,8,8},0,{" +
                ",".join(mod.QUARTIC_FIELDS) + "}]]",
            "native_Lam_string_length": 1000,
            "native_Lam_truncated": False,
            "native_fund_generator_components_numeric": False,
            "native_triplet_generator_components_numeric": False,
            "native_triplet_CG_components_numeric": False,
            "native_fund_generator_components_inputform": "unavailable",
            "native_triplet_generator_components_inputform": "unavailable",
            "native_triplet_CG_components_inputform": "unavailable",
            "native_operator_registry_complete": True,
            "matchete_physical_coupling_map_certified": False,
            "native_to_reference_component_axes_certified": False,
            "rgbeta_quartic_betas_certified": False,
            "production_modified": False,
        }

    def test_minimal_synthetic_registry_accepted_as_symbolic_only(self):
        result = mod.analyze(self.fixture(), self.root)
        self.assertEqual(result["registered_operators_verified"], 9)
        self.assertFalse(result["physical_matchete_quartic_map_certified"])
        self.assertFalse(result["scalar_beta_polynomials_certified"])

    def test_stale_source_rejected(self):
        native = self.fixture()
        native["source_model_sha256"] = "0" * 64
        with self.assertRaisesRegex(mod.OperatorWitnessError, "stale"):
            mod.analyze(native, self.root)

    def test_wrong_stage_rejected(self):
        native = self.fixture()
        native["stage"] = "UV"
        with self.assertRaises(mod.OperatorWitnessError):
            mod.analyze(native, self.root)

    def test_missing_adj_quartic_rejected(self):
        native = self.fixture()
        native["quartic_operators"].pop()
        with self.assertRaises(mod.OperatorWitnessError):
            mod.analyze(native, self.root)

    def test_missing_triplet_generator_rejected(self):
        native = self.fixture()
        native["quartic_operators"][-1]["invariant_inputform"] = "1/2"
        with self.assertRaisesRegex(mod.OperatorWitnessError, "tGen"):
            mod.analyze(native, self.root)

    def test_portal_field_order_rejected(self):
        native = self.fixture()
        native["quartic_operators"][4]["fields"] = "{Bar[H], H, S2, Bar[S2]}"
        with self.assertRaisesRegex(mod.OperatorWitnessError, "field order"):
            mod.analyze(native, self.root)

    def test_native_internal_prefactor_change_rejected(self):
        native = self.fixture()
        native["quartic_operators"][8]["coupling_inputform"] = "12*lambdaS2Adj"
        with self.assertRaisesRegex(mod.OperatorWitnessError, "ratio"):
            mod.analyze(native, self.root)

    def test_lam_tensor_missing_quartic_rejected(self):
        native = self.fixture()
        native["native_Lam_visible_quartic_symbols"].remove("lambdaS2Adj")
        with self.assertRaisesRegex(mod.OperatorWitnessError, "every quartic"):
            mod.analyze(native, self.root)

    def test_tensor_unevaluated_rejected(self):
        native = self.fixture()
        native["native_Lam_formal_tensor_available"] = False
        with self.assertRaises(mod.OperatorWitnessError):
            mod.analyze(native, self.root)

    def test_lambda_symbol_inside_other_symbol_is_not_enough(self):
        native = self.fixture()
        native["native_Lam_formal_tensor_inputform"] = native[
            "native_Lam_formal_tensor_inputform"].replace(
                "{lambdaH,", "{lambdaH2Adj,")
        with self.assertRaisesRegex(mod.OperatorWitnessError,
                                    "independent symbol lambdaH"):
            mod.analyze(native, self.root)

    def test_native_formal_tensor_truncation_rejected(self):
        native = self.fixture()
        native["native_Lam_truncated"] = True
        with self.assertRaisesRegex(mod.OperatorWitnessError, "truncated"):
            mod.analyze(native, self.root)

    def test_false_physical_certification_cannot_be_promoted(self):
        native = self.fixture()
        native["matchete_physical_coupling_map_certified"] = True
        with self.assertRaisesRegex(mod.OperatorWitnessError, "Unsafe"):
            mod.analyze(native, self.root)

    def test_scalar_species_collision_rejected(self):
        native = self.fixture()
        native["scalar_species_indices"]["S2"] = 2
        with self.assertRaisesRegex(mod.OperatorWitnessError, "ambiguous"):
            mod.analyze(native, self.root)

    def test_native_component_exposure_does_not_certify_matchete_map(self):
        native = self.fixture()
        for which in ("fund_generator", "triplet_generator", "triplet_CG"):
            native[f"native_{which}_components_numeric"] = True
            native[f"native_{which}_components_inputform"] = "{{1,0},{0,1}}"
        result = mod.analyze(native, self.root)
        self.assertTrue(result["native_numeric_triplet_generators"])
        self.assertFalse(result["physical_matchete_quartic_map_certified"])

    def test_wolfram_script_references_native_lam_and_registry(self):
        script = SCRIPT.with_name("RunAminus4RGBetaOperatorWitness.wl").read_text(
            encoding="utf-8")
        for required in ("T3RGBetaBuildEFT1[1, 3, 2, -4]", "$quartics",
                         'Names["RGBeta`PackageScope`Lam"]',
                         'Names["RGBeta`PackageScope`$fieldIndexMap"]',
                         "tGen[SU2L[S2]", "delS2[SU2L",
                         "matchete_physical_coupling_map_certified"):
            self.assertIn(required, script)
        self.assertNotIn("BetaTerm[", script)
        self.assertNotIn("Import[\"https://", script)


if __name__ == "__main__":
    unittest.main()
