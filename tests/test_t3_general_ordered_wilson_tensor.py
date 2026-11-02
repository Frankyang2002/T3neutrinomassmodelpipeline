"""Pure Python/fixture tests; live Matchete integration is a separate check."""
from __future__ import annotations

import json
import math
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
from types import ModuleType

import numpy as np

# Works in both the normal repository and the standalone ZIP test tree.
ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from Numerical.running.GeneralT3OrderedTensor import (
    scalar_basis, ScalarLeg, scalar_legs, REPRESENTATIONS, project_ordered_term,
)
from Numerical.running.GeneralT3LinearWilsonFlow import (
    validate_gamma, evolve_explicit_gamma, LOOP,
)
from scripts import RunT3GeneralOrderedWilsonTensor as run


def example_gamma(coords, value=2.0):
    g = np.eye(len(coords)) * value
    return {
        "format": "t3-explicit-ordered-wilson-gamma-v1",
        "equation": "16pi2*dC/dlnmu=Gamma(mu)@C",
        "basis": "ordered_Matchete_real_components_with_flavour",
        "coordinates": coords,
        "samples": [
            {"mu_gev": scale, "gamma": {"real": g.tolist(), "imag": (g * 0).tolist()}}
            for scale in (1e3, 1e5)
        ],
    }


class BasisTests(unittest.TestCase):
    def test_all_16_catalogue(self):
        self.assertEqual(len(run.MODELS), 16)
        for model in run.MODELS:
            b = scalar_basis(model)
            d1, d2, _ = REPRESENTATIONS[model[0]]
            self.assertEqual(b.total_dimension, 4 + 2 * (d1 + d2))

    def test_six_real_split_models(self):
        for model in ("A_-2", "A_+0", "D_-2", "D_+0", "E_-2", "E_+0"):
            b = scalar_basis(model)
            self.assertEqual(len(b.dimensions), 3)
            self.assertEqual(len(b.real_fields), 2)

    def test_mixed_both_complex_doublets(self):
        for model in ("B_-1", "C_-1"):
            b = scalar_basis(model)
            self.assertEqual(b.dimensions, {"NewScalarX": 2, "NewScalarY": 2})
            self.assertEqual(b.total_dimension, 12)
            self.assertFalse(b.real_fields)

    def test_complex_field_normalization_and_bar(self):
        b = scalar_basis("A_-4")
        x = b.expansion(ScalarLeg("NewScalar1", None, None, False), 1)
        y = b.expansion(ScalarLeg("NewScalar1", None, None, True), 1)
        self.assertEqual([k for k, _ in x], [5, 6])
        self.assertAlmostEqual(x[0][1].real, 2**(-0.5))
        self.assertAlmostEqual(x[1][1].imag, 2**(-0.5))
        self.assertAlmostEqual(y[1][1].imag, -2**(-0.5))

    def test_real_split_not_sqrt_two_twice(self):
        b = scalar_basis("D_-2")
        leg = ScalarLeg("T3Scalar2Minus", None, None, True)
        self.assertEqual(b.expansion(leg, 1), ((b.offsets[leg.name], 1 + 0j),))

    def test_explicit_leg_parsing(self):
        b = scalar_basis("D_-2")
        raw = ("Bar[Field[T3Scalar2Minus, Scalar, {}, {}]] * "
               "Field[NewScalar1, Scalar, {Index[d$$1, T3BSMd3]}, {}]")
        legs = scalar_legs(raw, b)
        self.assertEqual(legs[0].name, "T3Scalar2Minus")
        self.assertTrue(legs[0].barred)
        self.assertEqual(legs[1].rep, "T3BSMd3")

    def test_unsupported_field_fails_closed(self):
        with self.assertRaises(ValueError):
            scalar_legs("Field[InventedScalar, Scalar, {}, {}]", scalar_basis("A_-4"))

    def test_bad_component_raises(self):
        b = scalar_basis("C_-3")
        with self.assertRaises(ValueError):
            b.expansion(ScalarLeg("NewScalar1", None, None, False), 99)


    def test_full_ordered_component_projection_with_stubbed_cg_adapter(self):
        """Numerical scalar + SU2 + flavour multiplication, with empty CG product."""
        import sys
        fake = ModuleType("RGE.running.intermediate.ScalarOnlyTensorAdapters")
        fake.load_wilson_cg_registry = lambda seed: {}
        fake._find_wilson_cg_calls = lambda term: []
        fake._find_lepton_dummies = lambda term: ("li", "lj")
        fake._wilson_cg_value = lambda *a: 1
        fake._term_scalar_expression = lambda term: "1"
        fakeparse = ModuleType("RGE.running.intermediate.IntermediateMatcheteParsing")
        fakeparse.parse_mathematica_scalar = lambda term: 1
        term = {
            "source_index": 7,
            "operator_template_unit_flavor_inputform":
                "Field[NewScalar1, Scalar, {}, {}]*"
                "Field[NewScalar2, Scalar, {Index[s, T3BSMd3]}, {}]*GammaCC",
            "ordered_flavour_matrix_gev_inv": {
                "real": np.eye(3).tolist(), "imag": np.zeros((3,3)).tolist()},
        }
        fakepkg = ModuleType("RGE")
        fakepkg.__path__ = []
        fakesub = ModuleType("RGE.running")
        fakesub.__path__ = []
        fakesub2 = ModuleType("RGE.running.intermediate")
        fakesub2.__path__ = []
        with patch.dict(sys.modules, {"RGE": fakepkg, "RGE.running": fakesub,
                                      "RGE.running.intermediate": fakesub2,
                                      "RGE.running.intermediate.ScalarOnlyTensorAdapters": fake,
                                      "RGE.running.intermediate.IntermediateMatcheteParsing": fakeparse}):
            data = project_ordered_term(term, [], "A_-4")
        self.assertEqual(data["component_count"], 48)
        first = data["components"][0]
        self.assertAlmostEqual(first["ordered_flavour_coefficient_gev_inv"]["real"][0][0], .5)
        self.assertTrue(data["ordered_not_canonical_symmetrized"])


class RepeatedDummyTests(unittest.TestCase):
    """Einstein contraction: repeated scalar/lepton indices are one sum."""

    @staticmethod
    def _project(model, scalar_fields):
        fake = ModuleType("RGE.running.intermediate.ScalarOnlyTensorAdapters")
        fake.load_wilson_cg_registry = lambda seed: {}
        fake._find_wilson_cg_calls = lambda source: []
        fake._find_lepton_dummies = lambda source: ("li", "lj")
        fake._wilson_cg_value = lambda *args: 1
        fake._term_scalar_expression = lambda term: "1"
        fakeparse = ModuleType("RGE.running.intermediate.IntermediateMatcheteParsing")
        fakeparse.parse_mathematica_scalar = lambda raw: 1
        fakepkg = ModuleType("RGE")
        fakepkg.__path__ = []
        fakesub = ModuleType("RGE.running")
        fakesub.__path__ = []
        fakesub2 = ModuleType("RGE.running.intermediate")
        fakesub2.__path__ = []
        term = {
            "source_index": 9,
            "operator_template_unit_flavor_inputform": "*".join(scalar_fields) + "*GammaCC",
            "ordered_flavour_matrix_gev_inv": {
                "real": np.eye(3).tolist(),
                "imag": np.zeros((3, 3)).tolist(),
            },
        }
        with patch.dict(sys.modules, {
            "RGE": fakepkg,
            "RGE.running": fakesub,
            "RGE.running.intermediate": fakesub2,
            "RGE.running.intermediate.ScalarOnlyTensorAdapters": fake,
            "RGE.running.intermediate.IntermediateMatcheteParsing": fakeparse,
        }):
            return project_ordered_term(term, [], model)

    def test_repeated_scalar_scalar_dummy_contracts_once(self):
        for model, rep, dimension in (
            ("C_-3", "SU2L[fund]", 2),
            ("E_-4", "T3BSMd3", 3),
        ):
            with self.subTest(model=model):
                same_index = (
                    f"Field[NewScalar1, Scalar, {{Index[s, {rep}]}}, {{}}]",
                    f"Field[NewScalar2, Scalar, {{Index[s, {rep}]}}, {{}}]",
                )
                result = self._project(model, same_index)
                # 2*2 external SU(2) labels, one scalar sum, 2*2 real parts.
                self.assertEqual(result["component_count"], 4 * dimension * 4)
                first_offset = scalar_basis(model).offsets["NewScalar1"]
                second_offset = scalar_basis(model).offsets["NewScalar2"]
                for entry in result["components"]:
                    aa = (entry["scalar_real"][0] - first_offset) // 2
                    bb = (entry["scalar_real"][1] - second_offset) // 2
                    self.assertEqual(aa, bb)

    def test_repeated_dummy_with_self_bilinear_real_fields(self):
        cases = (
            ("E_-2", (
                "Field[NewScalar1, Scalar, {Index[s, T3BSMd3]}, {}]",
                "Field[T3Scalar2Minus, Scalar, {Index[s, T3BSMd3]}, {}]",
            ), 24),
            ("E_+0", (
                "Field[T3Scalar1Minus, Scalar, {Index[s, T3BSMd3]}, {}]",
                "Field[T3Scalar1Plus, Scalar, {Index[s, T3BSMd3]}, {}]",
            ), 12),
        )
        for model, fields, expected_count in cases:
            with self.subTest(model=model):
                result = self._project(model, fields)
                self.assertEqual(result["component_count"], expected_count)

    def test_repeated_lepton_scalar_dummy_enforces_equal_components(self):
        scalar_fields = (
            "Field[NewScalar1, Scalar, {Index[li, SU2L[fund]]}, {}]",
            "Field[NewScalar2, Scalar, {Index[s, SU2L[fund]]}, {}]",
        )
        result = self._project("C_-3", scalar_fields)
        self.assertEqual(result["component_count"], 32)
        offset = scalar_basis("C_-3").offsets["NewScalar1"]
        for entry in result["components"]:
            i = entry["lepton_su2"][0]
            aa = (entry["scalar_real"][0] - offset) // 2 + 1
            self.assertEqual(i, aa)

    def test_mismatched_dimensions_remain_blocked(self):
        scalar_fields = (
            "Field[NewScalar1, Scalar, {Index[li, SU2L[fund]]}, {}]",
            "Field[NewScalar2, Scalar, {Index[s, T3BSMd3]}, {}]",
        )
        with self.assertRaisesRegex(ValueError, "Inconsistent dimensions"):
            self._project("E_-4", scalar_fields)


class GammaTests(unittest.TestCase):
    def test_constant_diagonal_gamma_exact_exp(self):
        out = evolve_explicit_gamma(np.array([1 + 0j, 1j]), ["a", "b"],
                                    example_gamma(["a", "b"]), 1e5, 1e3)
        factor = (1e3 / 1e5) ** (2 * LOOP)
        self.assertTrue(np.allclose(np.array(out["end_real"]) + 1j * np.array(out["end_imag"]),
                                    np.array([1 + 0j, 1j]) * factor, rtol=1e-7))
        self.assertFalse(out["physical_matching_certified"])

    def test_gamma_wrong_basis_rejected(self):
        bad = example_gamma(["wrong"])
        with self.assertRaises(ValueError):
            validate_gamma(bad, ["a"], 1e5, 1e3)

    def test_gamma_missing_scales_rejected(self):
        bad = example_gamma(["a"])
        bad["samples"][0]["mu_gev"] = 1e4
        with self.assertRaises(ValueError):
            validate_gamma(bad, ["a"], 1e5, 1e3)

    def test_gamma_no_implicit_zero(self):
        bad = example_gamma(["a", "b"])
        bad["samples"][1]["gamma"]["real"] = [[1]]
        with self.assertRaises(ValueError):
            validate_gamma(bad, ["a", "b"], 1e5, 1e3)


class RunnerTests(unittest.TestCase):
    def test_flatten_retains_term_and_flavour_orientation(self):
        mat = {"real": np.arange(9).reshape(3, 3).tolist(),
               "imag": np.zeros((3, 3)).tolist()}
        terms = [{"source_index": 4, "components": [{
            "lepton_su2": [2, 1], "scalar_real": [5, 8],
            "ordered_flavour_coefficient_gev_inv": mat}]}]
        keys, y = run.flatten_for_gamma(terms)
        self.assertEqual(len(keys), 9)
        self.assertEqual(keys[0], "term4:L21:S5,8:F11")
        self.assertEqual(keys[-1], "term4:L21:S5,8:F33")
        self.assertEqual(y[5], 5)

    def test_mixed_source_validator_fails_if_old_basis(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            d = root / "B_m1"
            d.mkdir()
            (d / "status.json").write_text(json.dumps({"model": "B_-1",
                "status": "CurrentRunningMassBasisFThresholdTreeLLSSExported",
                "historical_mu_S_basis_reused": True,
                "source_operator_matching_at_running_mu_F": True}))
            with self.assertRaisesRegex(ValueError, "current-F"):
                run.mixed("B_-1", root, root / "out")

    def test_cli_summary_with_mocked_sources(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            with patch.object(run, "nonmixed", return_value={"model": "A_-4",
                    "status": "OrderedNumericalRealComponentTensorFromMatchete"}):
                self.assertEqual(run.main(["--model", "A_-4", "--output-root", str(root)]), 0)
            s = json.loads((root / "summary.json").read_text())
            self.assertEqual(s["ordered_real_tensor_count"], 1)
            self.assertFalse(s["physical_sequential_matching_certified"])


if __name__ == "__main__":
    unittest.main()
