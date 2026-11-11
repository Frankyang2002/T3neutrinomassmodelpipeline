"""Regression checks for the Batch 15 renormalisable scalar-threshold audit."""
import unittest
from dataclasses import dataclass
from types import SimpleNamespace

from Numerical.experimental.Batch15ScalarThresholdHandoff import (
    _decode_value, scalar_threshold_audit, GENERAL_MASS_KEYS,
)
from scripts.RunBatch15ScalarThresholdAll16 import MODELS


@dataclass(frozen=True)
class FakeBlock:
    name: str
    kind: str
    shape: tuple = ()


class Batch15ScalarThresholdTests(unittest.TestCase):
    def test_all_16_catalogue(self):
        self.assertEqual(len(MODELS),16)
        self.assertEqual(len(set(MODELS)),16)

    def test_decode_complex_quadratic(self):
        b=FakeBlock('muT3GeneralS1S2','complex_scalar')
        self.assertEqual(_decode_value({'real':3.0,'imag':-4.0},b),3-4j)

    def test_reject_nonzero_imaginary_real(self):
        with self.assertRaisesRegex(ValueError,'complex value'):
            _decode_value({'real':1,'imag':1e-6},FakeBlock('mS1Sq','real_scalar'))

    def test_reject_incorrect_shape(self):
        with self.assertRaisesRegex(ValueError,'incorrect value shape'):
            _decode_value({'real':[1,2],'imag':[0,0]},FakeBlock('mS1Sq','real_scalar'))

    def test_no_extra_bilinear_reports_diagnostic_not_certified(self):
        s=SimpleNamespace(mS1Sq=9e6,mS2Sq=16e6,general_couplings={})
        r=scalar_threshold_audit(s,mu_f_gev=1e5,mu_s_gev=1e3)
        self.assertTrue(r['diagonal_proxy_threshold_hierarchy_passed'])
        self.assertFalse(r['physical_scalar_mass_eigenvalues_computed'])
        self.assertFalse(r['mass_mixing_requires_resolution'])

    def test_extra_bilinear_must_not_be_silently_ignored(self):
        s=SimpleNamespace(mS1Sq=9e6,mS2Sq=16e6,
                          general_couplings={'muT3GeneralS1S2':1e5+3e4j})
        r=scalar_threshold_audit(s,mu_f_gev=1e5,mu_s_gev=1e3)
        self.assertTrue(r['mass_mixing_requires_resolution'])
        self.assertEqual(r['additional_quadratic_symbols'], ['muT3GeneralS1S2'])
        self.assertFalse(r['physical_scalar_mass_basis_certified'])

    def test_reject_multiple_bilinears(self):
        s=SimpleNamespace(mS1Sq=9e6,mS2Sq=16e6,
                          general_couplings={k:0j for k in list(GENERAL_MASS_KEYS)[:2]})
        with self.assertRaisesRegex(ValueError,'multiple'):
            scalar_threshold_audit(s,mu_f_gev=1e5,mu_s_gev=1e3)

    def test_reject_tachyonic_diagonal(self):
        s=SimpleNamespace(mS1Sq=-1,mS2Sq=16e6,general_couplings={})
        with self.assertRaisesRegex(ValueError,'Nonpositive'):
            scalar_threshold_audit(s,mu_f_gev=1e5,mu_s_gev=1e3)


if __name__ == '__main__':
    unittest.main()
