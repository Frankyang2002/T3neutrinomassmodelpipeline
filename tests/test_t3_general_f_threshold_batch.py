"""Checks source-backed, fail-closed first-threshold orchestration data."""
from __future__ import annotations

import argparse
import tempfile
import unittest
from pathlib import Path

from scripts import RunT3GeneralFThresholdBatch as mod


class GeneralFThresholdTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.folder = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def _valid(self, model="A_-4"):
        native = {
            "model": model, "status": "FThresholdTreeWilsonExported",
            "loop_order": 0, "integrated_fields": ["F"],
            "light_bsm_fields": ["S1", "S2"], "heavy_flavour_count": 3,
            "operator_basis_is_full_general_requested": True,
            "full_sequential_eft1_computed": False, "physical_c5_certified": False,
            "hard_one_loop_MSbar_certified": False,
            "eft1_wilson_anomalous_dimensions_implemented": False,
            "scalar_threshold_matching_implemented": False,
            "tree_wilson_term_count": 1, "wilson_seed_file": "stage1_wilson_seed.json",
        }
        seed = {
            "IntegratedFields": ["F"], "ActiveHeavyFields": ["S1", "S2"],
            "TreeWilsonTermCount": 1,
            "TreeWilsonTerms": [{"TermInputForm": "GammaCC[3]*Field[NewScalar1]*Field[NewScalar2]"}],
            "CGRegistry": ["CG1"],
        }
        mod.save_json(self.folder / "native_status.json", native)
        mod.save_json(self.folder / "data/stage1_wilson_seed.json", seed)
        return native, seed


    def _valid_self(self, model="A_-2"):
        native, seed = self._valid(model)
        species = 2 if model.endswith("-2") else 1
        basis = {
            "mass_basis_kind": "RealSelfBilinearSplitLightFields",
            "split_species": species,
            "scalar_mass_components_remain_light": True,
            "real_bilinear_mass_basis_candidate": True,
            "complex_self_bilinear_phase_supported": False,
            "holomorphic_bilinear_absorbed_into_split_free_masses": True,
            "field_basis_normalization_physically_certified": False,
        }
        native["self_bilinear_mass_basis"] = basis
        seed["SelfBilinearMassBasis"] = basis
        seed["WilsonSeedSource"] = "MatcheteFThresholdTreeReclassifiedForSplitLightScalarFields"
        seed["TreeWilsonTerms"][0]["TermInputForm"] = (
            f"GammaCC[3]*Field[T3Scalar{species}Minus]*Field[T3Scalar{species}Plus]")
        mod.save_json(self.folder / "native_status.json", native)
        mod.save_json(self.folder / "data/stage1_wilson_seed.json", seed)
        mod.save_json(self.folder / "self_mass_basis.json", basis)
        return native, seed

    def test_wolfram_split_fields_remain_light_and_seed_selector_recognizes_them(self):
        source = mod.NATIVE.read_text(encoding="utf-8")
        for name in ("T3Scalar1Minus", "T3Scalar1Plus", "T3Scalar2Minus", "T3Scalar2Plus"):
            self.assertIn(f"Field[{name}", source)
        for name in ("T31MinusMass", "T31PlusMass", "T32MinusMass", "T32PlusMass"):
            self.assertIn(f"Mass->{{Light,{name}}}", source)
        self.assertIn('"CGRegistry"->splitCGRegistry', source)
        self.assertIn('"SelfBilinearMassBasis"->selfMetadata', source)
        self.assertIn('!hasExtraMassQ[#] &', source)
        self.assertNotIn('Mass->{Heavy,T31MinusMass}', source)

    def test_six_self_models_declared(self):
        self.assertEqual(len(mod.SELF), 6)
        self.assertEqual(mod.SELF & mod.MIXED, frozenset())

    def test_self_bilinear_source_contract(self):
        self._valid_self("A_-2")
        result = mod.verify_native_tree(self.folder, "A_-2", 0)
        self.assertEqual(result["mass_basis"], "real-self-bilinear-split-light-scalar")

    def test_self_bilinear_missing_basis_blocked(self):
        self._valid("A_+0")
        with self.assertRaisesRegex(ValueError, "mass-basis certificate"):
            mod.verify_native_tree(self.folder, "A_+0", 0)

    def test_self_bilinear_unrotated_scalar_blocked(self):
        native, seed = self._valid_self("D_+0")
        seed["TreeWilsonTerms"][0]["TermInputForm"] += "*Field[NewScalar1]"
        mod.save_json(self.folder / "data/stage1_wilson_seed.json", seed)
        with self.assertRaisesRegex(ValueError, "unrotated original scalar"):
            mod.verify_native_tree(self.folder, "D_+0", 0)

    def test_self_bilinear_metadata_mismatch_blocked(self):
        native, _ = self._valid_self("E_-2")
        native["self_bilinear_mass_basis"]["scalar_mass_components_remain_light"] = False
        mod.save_json(self.folder / "native_status.json", native)
        with self.assertRaisesRegex(ValueError, "scalar_mass_components_remain_light"):
            mod.verify_native_tree(self.folder, "E_-2", 0)

    def test_catalogue_full16(self):
        self.assertEqual(len(mod.MODELS), 16)
        self.assertEqual(len(set(mod.MODELS)), 16)
        self.assertIn("B_-1", mod.MODELS)
        self.assertIn("C_-1", mod.MODELS)

    def test_model_selection(self):
        x = argparse.Namespace(all=False, model=None, models="A_-4,B_-1,C_-1")
        self.assertEqual(mod.select_models(x), ("A_-4", "B_-1", "C_-1"))
        x.models = "A_-4,A_-4"
        with self.assertRaises(ValueError):
            mod.select_models(x)

    def test_valid_seed(self):
        self._valid()
        result = mod.verify_native_tree(self.folder, "A_-4", 0)
        self.assertEqual(result["tree_wilson_term_count"], 1)
        self.assertFalse(result["scalar_threshold_matching_performed"])

    def test_wrong_loop_order_rejected(self):
        self._valid()
        with self.assertRaisesRegex(ValueError, "loop_order"):
            mod.verify_native_tree(self.folder, "A_-4", 1)

    def test_heavy_particle_contamination_rejected(self):
        _, seed = self._valid()
        seed["TreeWilsonTerms"][0]["TermInputForm"] += "*Field[NewFermion]"
        mod.save_json(self.folder / "data/stage1_wilson_seed.json", seed)
        with self.assertRaisesRegex(ValueError, "integrated F"):
            mod.verify_native_tree(self.folder, "A_-4", 0)

    def test_seed_count_mismatch_rejected(self):
        _, seed = self._valid()
        seed["TreeWilsonTermCount"] = 2
        mod.save_json(self.folder / "data/stage1_wilson_seed.json", seed)
        with self.assertRaisesRegex(ValueError, "inconsistent tree Wilson"):
            mod.verify_native_tree(self.folder, "A_-4", 0)


if __name__ == "__main__":
    unittest.main()
