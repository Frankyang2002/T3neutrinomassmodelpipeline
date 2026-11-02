"""Focused Batch 1 tests: source contract, Weyl sign ambiguity, fail-closed gates."""
from __future__ import annotations
import copy
import json
from pathlib import Path
import tempfile
import unittest

from Numerical.lo_certification.benchmark import batch1 as b

ROOT = Path(__file__).resolve().parents[1]


def fake_frozen(contract: dict) -> dict:
    v = contract['frozen_values_at_mu_UV']
    s = contract['matching']
    mf = v['M_F_diagonal_gev']
    return {
        'representation': {'d_s1': 1, 'd_s2': 3, 'd_f': 2, 'alpha': -4},
        'scales': {'mu_fermion_threshold_gev': s['mu_F_gev'],
                   'mu_scalar_threshold_gev': s['mu_S_gev'],
                   'mu_low_gev': s['mu_low_gev']},
        'numerical': {'vev_gev': v['vev_gev']},
        'base_state': {'mu_gev': s['mu_UV_gev'], 'ordinary': {
            'y1': {'real': copy.deepcopy(v['y1_real']), 'imag': [[0]*3 for _ in range(3)]},
            'y2': {'real': copy.deepcopy(v['y2_real']), 'imag': [[0]*3 for _ in range(3)]},
            'MF': {'real': [[mf[i] if i==j else 0 for j in range(3)] for i in range(3)],
                   'imag': [[0]*3 for _ in range(3)]},
            'mS1Sq': v['m_S1_squared_gev2'], 'mS2Sq': v['m_S2_squared_gev2'],
            'lambdaS1': v['lambda_S1'], 'lambdaS2': v['lambda_S2'],
            'lambdaH1': v['lambda_H1'], 'lambdaH2': v['lambda_H2'],
            'lambda12': v['lambda_12'], 'lambdaT3': v['lambda_T3'],
            'lambdaH2Adj': v['lambda_H2Adj'], 'lambdaS2Adj': v['lambda_S2Adj'],
        }},
    }


class TestRestrictedBenchmark(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.contract = b.read_json(ROOT / b.INPUT_REL)

    def test_source_derived_numerical_contract(self):
        r = b.compare_frozen(self.contract, fake_frozen(self.contract))
        self.assertAlmostEqual(r['scalar_unbroken_mass_parameters_gev'][0], 1142.857142857143)
        self.assertAlmostEqual(r['mF_min_over_mS_max'], 87.5)
        self.assertFalse(r['couplings_run_to_mu_F'])
        self.assertFalse(r['scalar_broken_phase_masses_certified'])

    def test_frozen_yukawa_drift_fails(self):
        frozen = fake_frozen(self.contract)
        frozen['base_state']['ordinary']['y1']['real'][0][0] += 1e-6
        with self.assertRaisesRegex(b.AuditError, 'y1 real'):
            b.compare_frozen(self.contract, frozen)

    def test_incorrect_hypercharge_fails(self):
        c = copy.deepcopy(self.contract)
        c['representation']['scalar2']['hypercharge'] = +1
        with self.assertRaisesRegex(b.AuditError, 'hypercharge'):
            b.compare_frozen(c, fake_frozen(self.contract))

    def test_fermion_must_be_heavier_than_both_scalars(self):
        c = copy.deepcopy(self.contract)
        f = fake_frozen(c)
        c['frozen_values_at_mu_UV']['m_S1_squared_gev2'] = 101000.0**2
        f['base_state']['ordinary']['mS1Sq'] = 101000.0**2
        with self.assertRaisesRegex(b.AuditError, 'Fermion not heavier'):
            b.compare_frozen(c, f)

    def test_threshold_order_fails_closed(self):
        c = copy.deepcopy(self.contract)
        c['matching']['mu_S_gev'] = 100000.0
        f = fake_frozen(c)
        with self.assertRaisesRegex(b.AuditError, 'two-threshold ordering'):
            b.compare_frozen(c, f)

    def test_exact_signed_gamma_ambiguity(self):
        d = b.signed_ambiguity()
        self.assertEqual([x['PL_over_Weyl_epsilon'] for x in
                          d['two_indistinguishable_abstract_choices']], [1, -1])
        self.assertFalse(d['abstract_GammaCC_identity_sufficient_to_fix_sign'])
        self.assertIsNone(d['native_GammaCC_to_master_Weyl_factor'])
        self.assertTrue(all(x['GammaCC_square_minus_I'] and x['GammaCC_antisymmetric']
                            for x in d['two_indistinguishable_abstract_choices']))

    def test_normalization_gate_null_by_construction(self):
        self.assertIsNone(self.contract['source_verification']['gamma_cc_to_weyl_factor'])
        self.assertIsNone(self.contract['source_verification']['physical_LLSS_boundary'])
        self.assertFalse(self.contract['physical_certification'])

    def test_git_blob_sha_algorithm(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / 'x'
            p.write_bytes(b'test content\n')
            self.assertEqual(b.git_blob_sha(p), 'd670460b4b4aece5915caf5c68d12f560a9fe3e4')
            self.assertIsNone(b.git_blob_sha(Path(d) / 'absent'))

    def test_immutable_output(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / 'report.json'
            b.write_immutable(path, b'first')
            b.write_immutable(path, b'first')
            with self.assertRaisesRegex(b.AuditError, 'Refusing to overwrite'):
                b.write_immutable(path, b'different')
            self.assertEqual(path.read_bytes(), b'first')

    def test_path_is_unique_owned_workstream(self):
        with self.assertRaisesRegex(b.AuditError, 'Unsafe run ID'):
            b.run(ROOT, '../other', check_only=True)

    def test_real_source_integration_when_local_checkout_present(self):
        source = ROOT / self.contract['source']['provisional_boundary_path']
        if not source.exists():
            self.skipTest('Actual original source release is only in the Windows/GitHub repository')
        result = b.audit(ROOT)
        self.assertEqual(result['source_validation']['A_m4_independent_LLSS_rows'], 60)
        self.assertEqual(result['source_validation']['A_m4_native_PL_sources'], 1)
        self.assertIsNone(result['fermion_boundary']['numeric_physical_C_LLSS'])
        self.assertFalse(result['physical_certification']['fermion_boundary_certified'])


if __name__ == '__main__':
    unittest.main()
