"""No Wolfram dependency: physical-field grammar and fixed-basis regressions."""
from __future__ import annotations

import importlib.util
from pathlib import Path
import sys
import unittest

import sympy as sp

ROOT = Path(__file__).resolve().parents[2]
MODULE = ROOT / "RGE/running/intermediate/GeneralDMinus2RealTensor.py"
spec = importlib.util.spec_from_file_location("D2RealAdapterTest", MODULE)
assert spec and spec.loader
adapter = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = adapter
spec.loader.exec_module(adapter)


class D2RealTensorGrammarTests(unittest.TestCase):
    def test_singlet_real_coordinates_have_no_second_sqrt2(self):
        for field, index in (("T3S2R", 11), ("T3S2I", 12)):
            leg = adapter.FieldLeg(field, None, None, False)
            self.assertEqual(adapter.scalar_coordinates(leg,1),((index,sp.S.One),))

    def test_complex_triplet_has_interleaved_real_components(self):
        leg = adapter.FieldLeg("NewScalar1","i","T3BSMd3",False)
        expansion = adapter.scalar_coordinates(leg,2)
        self.assertEqual([x[0] for x in expansion],[7,8])
        self.assertEqual([sp.simplify(x[1]) for x in expansion],
                         [1/sp.sqrt(2),sp.I/sp.sqrt(2)])

    def test_conjugated_complex_field_changes_imaginary_sign(self):
        leg = adapter.FieldLeg("H","i","SU2L[fund]",True)
        self.assertEqual(adapter.scalar_coordinates(leg,1)[1][1],-sp.I/sp.sqrt(2))

    def test_scalar_parser_treats_barred_real_as_real(self):
        text = "Bar[Field[T3S2R, Scalar, {}, {}]]*Field[T3S2I, Scalar, {}, {}]"
        legs = adapter.scalar_legs(text,2)
        self.assertEqual([x.name for x in legs],["T3S2R","T3S2I"])
        self.assertEqual([len(adapter.scalar_coordinates(x,1)) for x in legs],[1,1])

    def test_triplet_dummy_is_recovered(self):
        txt = "Field[NewScalar1, Scalar, {Index[d$$1, T3BSMd3]}, {}]*Field[T3S2R, Scalar, {}, {}]"
        legs = adapter.scalar_legs(txt,2)
        self.assertEqual((legs[0].dummy,legs[0].rep),("d$$1","T3BSMd3"))

    def test_wilson_prefactor_recovers_extra_yukawa_and_phase(self):
        text = ("-2*Bar[Coupling[y1, {Index[p, Flavor]}, 0]]*"
                "Coupling[y2Conj, {Index[q, Flavor]}, 0]*"
                "Bar[Coupling[T3S2Phase, {}, 0]]/Coupling[MF, {}, 0]")
        value,names=adapter.scalar_prefactor(text)
        self.assertEqual(set(names),{"y1","y2Conj","T3S2Phase","MF"})
        target=-2*sp.conjugate(sp.Symbol('y1'))*sp.Symbol('y2Conj')*sp.conjugate(sp.Symbol('T3S2Phase'))/sp.Symbol('MF')
        self.assertEqual(sp.simplify(value-target),0)

    def test_physical_real_quartic_prefactor(self):
        text=("-Coupling[lambdaT3ExtraQ000004I1, {}, 0]*"
              "Field[T3S2R, Scalar, {}, {}]^4/4")
        expr,names=adapter.scalar_prefactor(text)
        self.assertEqual(names,("lambdaT3ExtraQ000004I1",))
        self.assertEqual(sp.simplify(expr+sp.Symbol(names[0])/4),0)

    def test_real_quartic_phase_does_not_vanish(self):
        text=("Coupling[lambdaT3ExtraQ000013I1, {}, 0]*"
              "Coupling[T3S2Phase, {}, 0]^3*"
              "Bar[Coupling[T3S2Phase, {}, 0]]")
        expr,_=adapter.scalar_prefactor(text)
        self.assertTrue(expr.has(sp.Symbol('T3S2Phase')))
        self.assertTrue(expr.has(sp.conjugate(sp.Symbol('T3S2Phase'))))

    def test_barred_cg_is_not_treated_as_coupling(self):
        expr,names=adapter.scalar_prefactor(
            "Bar[CG[T3MyCG, {Index[i, SU2L[fund]]}]]*Coupling[lambdaH1, {}, 0]"
        )
        self.assertEqual(expr,sp.Symbol('lambdaH1'))
        self.assertEqual(names,('lambdaH1',))

    def test_missing_external_su2_labels_is_error(self):
        with self.assertRaises(ValueError):
            adapter._l_indices("Field[l, Fermion, {}, {}]")

    def test_unregistered_mass_basis_phase_is_error(self):
        with self.assertRaises(ValueError):
            adapter.scalar_prefactor('T3S2Phase[]')

    def test_extra_quartic_names_inventory(self):
        runner=ROOT/'scripts/RunGeneralDMinus2RealEFT1RGE.py'
        source=runner.read_text(encoding='utf-8')
        for name in ('lambdaT3ExtraQ000004I1', 'lambdaT3ExtraQ000013I1',
                     'lambdaT3ExtraQ001102I1', 'lambdaT3ExtraQ020101I1',
                     'lambdaT3ExtraQ110002I1'):
            self.assertIn(name,source)


if __name__ == '__main__':
    unittest.main()


class D2RealPowerParsingRegressionTests(unittest.TestCase):
    """Protect against Field[...]^n being counted as just one scalar."""

    def test_four_identical_real_singlet_fields_are_four_legs(self):
        term = "-Coupling[lambdaT3ExtraQ000004I1, {}, 0]*Field[T3S2R, Scalar, {}, {}]^4/4"
        legs = adapter.scalar_legs(term, 4)
        self.assertEqual([leg.name for leg in legs], ["T3S2R"] * 4)

    def test_power_outside_bar_wrapper_is_not_lost(self):
        term = "Bar[Field[T3S2I, Scalar, {}, {}]]^3*Field[T3S2R, Scalar, {}, {}]"
        legs = adapter.scalar_legs(term, 4)
        self.assertEqual([leg.name for leg in legs], ["T3S2I"] * 3 + ["T3S2R"])
        self.assertTrue(all(leg.barred for leg in legs[:3]))

    def test_complex_triplet_power_represents_repeated_gauge_label(self):
        term = "Field[NewScalar1, Scalar, {Index[d$$1, T3BSMd3]}, {}]^2"
        legs = adapter.scalar_legs(term, 2)
        self.assertEqual(len(legs), 2)
        self.assertEqual(legs[0].dummy, legs[1].dummy)

    def test_exponent_four_reduces_correctly_in_scalar_prefactor(self):
        term = "-Coupling[lambdaT3ExtraQ000004I1, {}, 0]*Field[T3S2R, Scalar, {}, {}]^4/4"
        val, names = adapter.scalar_prefactor(term)
        self.assertEqual(names, ("lambdaT3ExtraQ000004I1",))
        self.assertEqual(sp.simplify(val + sp.Symbol(names[0])/4), 0)

    def test_cg_powers_are_expanded_and_outer_bar_conjugates(self):
        from dataclasses import dataclass

        @dataclass(frozen=True)
        class CG:
            name: str
            conjugated: bool = False

        calls = adapter._expanded_cg_calls("Bar[CG[T3X, {}]]^2", (CG("T3X"),))
        self.assertEqual(len(calls), 2)
        self.assertTrue(all(call.conjugated for call in calls))
