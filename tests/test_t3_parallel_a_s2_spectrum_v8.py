"""Workstream A v8: source-independent SU(2) spectrum and guarded native witness.

Mathematica live-native execution must be performed on the user's Windows
checkout. The all-16 Matchete fixture check runs only with repository tests.
"""
from __future__ import annotations

from copy import deepcopy
from hashlib import sha256
from pathlib import Path
import tempfile
import unittest

import numpy as np

from Numerical.parallel_t3.wilson.s2_intertwiner_v7 import (
    S2EvidenceError, canonical_checks, adjoint_tensor, rel,
)
from Numerical.parallel_t3.wilson.s2_spectral_v8 import (
    FORMAT, STATUS_REDUCED, STATUS_UNEXPOSED,
    NATIVE_FORMAT, NATIVE_REDUCED, NATIVE_UNEXPOSED, SOURCE_FILES,
    EXPECTED_EIGENVALUES, _identity, _spectral, select_spectral_sign,
    inspect_witness, _v6_adjoint,
)


def native_fixture(reduced=True):
    return {
        'format':NATIVE_FORMAT,
        'status':NATIVE_REDUCED if reduced else NATIVE_UNEXPOSED,
        'fundamental_dimension_and_trace_valid':True,
        'S2_dimension_trace_and_casimir_valid':True,
        'native_s2_generator_still_symbolic':True,
        'native_delS2_still_symbolic':True,
        'native_S2_symmetric_isometry_valid':True,
        'native_s2_generator_reduction_observed':reduced,
        'native_delS2_reduction_observed':reduced,
        'native_s2_generator_refined_inputform':'tGen[SU2L[fund], A, a[1], b[1]]',
        'native_delS2_refined_inputform':'del[SU2L[fund], a[1], i]',
        'native_s2_numeric_components_exposed':False,
        'native_s2_numeric_frame_certified':False,
        'native_s2_to_matchete_intertwiner_certified':False,
        'physical_coupling_dictionary_certified':False,
        'physical_wilson_evolution_authorized':False,
    }


class TestBasisInvariantSpectrum(unittest.TestCase):
    def setUp(self):
        _, g, _=canonical_checks()
        self.k=adjoint_tensor(g)

    def test_total_casimir_spectrum_polynomial(self):
        check=_spectral(self.k)
        np.testing.assert_allclose(check['eigenvalues'],EXPECTED_EIGENVALUES,atol=1.e-13)
        self.assertLess(check['minimal_polynomial_residual'],1.e-13)
        self.assertAlmostEqual(check['operator_norm'],np.sqrt(3.))

    def test_opposite_sign_rejected_even_with_identical_norm(self):
        plus=_spectral(self.k)
        minus=_spectral(-self.k)
        self.assertAlmostEqual(plus['operator_norm'],minus['operator_norm'])
        self.assertLess(plus['max_residual'],1.e-12)
        self.assertGreater(minus['minimal_polynomial_residual'],.1)
        self.assertGreater(minus['spectrum_residual'],.1)

    def test_native_direction_chooses_sign_without_v6(self):
        # First nonzero native traceless direction has opposite sign.
        raw=[1.2*_identity(3)-.6*self.k,
             -0.3*_identity(3)+2.5*self.k]
        k,check=select_spectral_sign(raw)
        self.assertEqual(check['sign_relative_to_first_native_traceless_CG'],-1)
        self.assertLess(rel(k,self.k),1.e-13)
        self.assertLess(check['native_Matchete_portal_span_residual'],1.e-13)
        self.assertGreater(check['opposite_sign_spectrum_residual'],.1)

    def test_triplet_unitary_frame_rotations_preserve_spectrum(self):
        rng=np.random.default_rng(908211)
        eye2=np.eye(2,dtype=complex)
        for _ in range(15):
            a=rng.normal(size=(3,3))+1j*rng.normal(size=(3,3))
            q,r=np.linalg.qr(a)
            q=q@np.diag(np.exp(-1j*np.angle(np.diag(r))))
            w=np.kron(eye2,q)
            km=self.k.transpose(0,2,1,3).reshape(6,6)
            transformed=w@km@w.conj().T
            tensor=transformed.reshape(2,3,2,3).transpose(0,2,1,3)
            self.assertLess(_spectral(tensor)['max_residual'],1.e-12)

    def test_wrong_tensor_dimension_fails_closed(self):
        with self.assertRaisesRegex(S2EvidenceError,'must be'):
            _spectral(np.zeros((2,2,2,2)))

    def test_wrong_hermitian_structure_fails_selection(self):
        changed=self.k.copy()
        changed[0,0,0,1]+=0.15
        with self.assertRaises(S2EvidenceError):
            select_spectral_sign([_identity(3),changed])


class TestNativeWitness(unittest.TestCase):
    def test_native_reduction_is_symbolic_only(self):
        self.assertTrue(inspect_witness(native_fixture(True))['symbolic_reduction_observed'])
        self.assertFalse(inspect_witness(native_fixture(False))['symbolic_reduction_observed'])
        self.assertFalse(inspect_witness(native_fixture(True))['numeric_s2_frame_exposed'])

    def test_claiming_numeric_certification_rejected(self):
        w=native_fixture()
        for key in ('native_s2_numeric_frame_certified',
                    'native_s2_to_matchete_intertwiner_certified',
                    'physical_wilson_evolution_authorized'):
            changed=deepcopy(w);changed[key]=True
            with self.subTest(key=key), self.assertRaisesRegex(
                    S2EvidenceError,'unsupported numeric/physical'):
                inspect_witness(changed)

    def test_inconsistent_reduction_status_rejected(self):
        w=native_fixture(True)
        w['native_delS2_reduction_observed']=False
        with self.assertRaisesRegex(S2EvidenceError,'contradicts'):
            inspect_witness(w)

    def test_source_hash_guard_rejects_stale_native(self):
        w=native_fixture()
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            for n in SOURCE_FILES:
                p=root/n;p.parent.mkdir(parents=True,exist_ok=True)
                p.write_text(n,encoding='utf-8')
            w['source_sha256']={name:sha256((root/name).read_bytes()).hexdigest()
                                for name in SOURCE_FILES}
            inspect_witness(w,root)
            (root/SOURCE_FILES[0]).write_text('changed',encoding='utf-8')
            with self.assertRaisesRegex(S2EvidenceError,'Stale'):
                inspect_witness(w,root)


class TestAll16NativeFixture(unittest.TestCase):
    @staticmethod
    def prepared():
        try:
            from tests.test_t3_parallel_a_s2_intertwiner_v7 import TestV7Fixture
        except ImportError as exc:
            raise unittest.SkipTest('Full repository v5-v7 fixture unavailable') from exc
        return TestV7Fixture.prepared_v6()

    def test_all_14_sectors_spectral_sign_independent_of_v6_selection(self):
        from Numerical.parallel_t3.wilson.s2_spectral_v8 import audit
        mat,v6=self.prepared()
        for reduction in (True,False):
            result=audit(mat,v6,native_fixture(reduction))
            self.assertEqual(result['triplet_sectors_basis_invariant_spectrum_verified'],14)
            self.assertEqual(result['status'],STATUS_REDUCED if reduction else STATUS_UNEXPOSED)
            self.assertLess(result['maximum_spectral_and_v6_comparison_residual'],1.e-12)
            self.assertFalse(result['certification']['physical_Wilson_evolution_authorized'])
            for sectors in result['model_sectors'].values():
                for row in sectors:
                    self.assertTrue(row['source_CG_sign_selected_without_using_v6_map'])
                    self.assertGreater(row['opposite_sign_polynomial_residual'],.1)

    def test_v6_tampering_fails_spectral_sign_comparison(self):
        from Numerical.parallel_t3.wilson.s2_spectral_v8 import audit
        mat,v6=self.prepared()
        changed=deepcopy(v6)
        sector=next(x for y in changed['models'].values() for x in y['sectors']
                    if x['sector']=='triplet')
        sector['candidate_rgbeta_from_matchete_real_matrix'][1]=[
            -x for x in sector['candidate_rgbeta_from_matchete_real_matrix'][1]]
        with self.assertRaisesRegex(S2EvidenceError,'spectrum contradicts v6'):
            audit(mat,changed,native_fixture())


if __name__=='__main__':
    unittest.main()
