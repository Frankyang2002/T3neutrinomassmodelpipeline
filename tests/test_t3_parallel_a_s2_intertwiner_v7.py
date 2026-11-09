"""V7 invariance checks and fail-closed symbolic/native provenance tests.

Synthetic unit tests execute independently; the all-16 fixture test runs when
installed into the user's full repository containing the v5/v6 tests.
"""
from __future__ import annotations

from hashlib import sha256
from copy import deepcopy
import json
from unittest.mock import patch
from pathlib import Path
import tempfile
import unittest

import numpy as np

from Numerical.parallel_t3.wilson.s2_intertwiner_v7 import (
    S2EvidenceError,TOL,STATUS,canonical_checks,frame_ambiguity,
    generator_checks,induced_s2_generators,intertwiner_from_generators,rel,audit,
    audit_from_paths,
)


class TestS2Math(unittest.TestCase):
    def test_canonical_induced_s2(self):
        c,g,info=canonical_checks()
        self.assertLess(info['maximum'],1e-12)
        np.testing.assert_allclose(g[2],np.diag([1,0,-1]),atol=1e-13)
        np.testing.assert_allclose(np.einsum('aij,bij->ab',c.conj(),c),
                                   np.eye(3),atol=1e-13)

    def test_unitary_frame_change_and_ambiguity(self):
        c,g,_=canonical_checks()
        rng=np.random.default_rng(76023)
        for _ in range(12):
            raw=rng.normal(size=(3,3))+1j*rng.normal(size=(3,3))
            q,r=np.linalg.qr(raw)
            q=q@np.diag(np.exp(-1j*np.angle(np.diag(r))))
            transformed=np.einsum('ux,axy,vy->auv',q,g,q.conj())
            u,checks=intertwiner_from_generators(transformed,g)
            self.assertLess(checks['intertwining'],1e-12)
            self.assertLess(checks['unitarity'],1e-12)
            self.assertLess(rel(u@u.conj().T,np.eye(3)),1e-12)
            witness=frame_ambiguity(g,c,u,transformed)
            self.assertGreater(witness['distinct_generator_frame_distance'],.01)
            self.assertLess(witness['rotated_intertwining'],1e-12)
            self.assertLess(witness['transported_adjoint_tensor_invariance'],1e-12)

    def test_negative_lie_sign_and_nonrepresentation_rejected(self):
        _,g,_=canonical_checks()
        self.assertGreater(generator_checks(-g)['lie_bracket'],.2)
        broken=g.copy();broken[1,0,0]+=0.15
        with self.assertRaises(S2EvidenceError):
            intertwiner_from_generators(broken,g)

    def test_transpose_of_triplet_generator_is_not_equivalent_fundamental_action(self):
        c,g,_=canonical_checks()
        # S2 triplet is complex in the chosen weight basis (T2 antisymmetric).
        self.assertGreater(rel(g,g.transpose(0,2,1)),.1)
        self.assertLess(generator_checks(g)['maximum'],1e-12)


class TestV7Fixture(unittest.TestCase):
    @staticmethod
    def prepared_v6():
        try:
            from tests.test_t3_parallel_a_portal_orientation_v6 import prepared
            from Numerical.parallel_t3.wilson.portal_orientation_v6 import (
                audit as v6_audit,SOURCES_FIERZ)
        except ImportError as exc:
            raise unittest.SkipTest('Full repository v5/v6 test fixtures unavailable') from exc
        mat,rgb,v5,witness=prepared()
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)
            for name in SOURCES_FIERZ:
                p=root/name;p.parent.mkdir(parents=True,exist_ok=True)
                p.write_text(name,encoding='utf-8')
            witness['source_sha256']={p:sha256((root/p).read_bytes()).hexdigest()
                                      for p in SOURCES_FIERZ}
            v6=v6_audit(mat,rgb,v5,witness,root,input_sha256={
                'matchete':'a'*64,'rgbeta':'b'*64,
                'v5_audit':'c'*64,'native_fierz':'d'*64})
        return mat,v6

    def test_all16_historical_registered_cg_fixture(self):
        mat,v6=self.prepared_v6()
        report=audit(mat,v6)
        self.assertEqual(report['status'],STATUS)
        self.assertEqual(report['triplet_sectors_conditional_comparison_verified'],14)
        self.assertFalse(report['independent_native_RGBeta_S2_to_Matchete_intertwiner_certified'])
        self.assertLess(report['maximum_comparison_residual'],1e-12)
        self.assertEqual(set(report['models']),set(v6['models']))
        for row in (s for a in report['models'].values() for s in a):
            self.assertFalse(row['independent_native_RGBeta_to_Matchete_intertwiner_certified'])
            self.assertTrue(row['source_generator_origin'].startswith('RECONSTRUCTED_'))
            self.assertGreater(row['alternative_frame_checks']['distinct_generator_frame_distance'],.01)

    def test_v6_triplet_schema_not_v5_schema(self):
        mat,v6=self.prepared_v6()
        triplets=[s for model in v6['models'].values() for s in model['sectors']
                  if s['sector']=='triplet']
        self.assertEqual(len(triplets),14)
        self.assertTrue(all('candidate_rgbeta_from_matchete_real_matrix' in s
                            for s in triplets))
        self.assertTrue(all('candidate_rgbeta_from_native_matchete_real_matrix' not in s
                            for s in triplets))
        self.assertEqual(audit(mat,v6)['triplet_sectors_conditional_comparison_verified'],14)

    def test_saved_v6_source_recheck_uses_v6_schema_and_detects_tampering(self):
        # Test v7's saved-v6 recheck separately from the expensive native
        # source revalidation, which is exercised by the existing v6 tests.
        mat,v6=self.prepared_v6()
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)
            inputs={'matchete.json':mat,'rgbeta.json':{},'v5.json':{},
                    'fierz.json':{},'v6.json':v6}
            for filename,obj in inputs.items():
                (root/filename).write_text(json.dumps(obj),encoding='utf-8')
            args=(root/'matchete.json',root/'rgbeta.json',root/'v5.json',
                  root/'fierz.json',root/'v6.json',root)
            with patch('Numerical.parallel_t3.wilson.portal_orientation_v6.audit_from_paths',
                       return_value=v6):
                report=audit_from_paths(*args)
            self.assertEqual(report['triplet_sectors_conditional_comparison_verified'],14)
            edited=deepcopy(v6)
            sector=next(s for model in edited['models'].values()
                        for s in model['sectors'] if s['sector']=='triplet')
            sector['candidate_rgbeta_from_matchete_real_matrix'][1][0]+=0.1
            (root/'v6.json').write_text(json.dumps(edited),encoding='utf-8')
            with patch('Numerical.parallel_t3.wilson.portal_orientation_v6.audit_from_paths',
                       return_value=v6):
                with self.assertRaisesRegex(S2EvidenceError,'saved v6 coupling map differs'):
                    audit_from_paths(*args)

    def test_historical_native_symbolic_status_accepted_but_not_certified(self):
        mat,v6=self.prepared_v6()
        witness={'status':'NATIVE_SYMBOLIC_GENERATOR_RULES_CAPTURED_NUMERIC_AXES_PENDING',
                 'RGBeta_numeric_eps_components_certified':False,
                 'RGBeta_numeric_delS2_components_certified':False,
                 'RGBeta_to_Matchete_basis_identified':False,
                 'physical_coupling_dictionary_certified':False}
        fixture={'status':'candidate_from_repository_test_fixtures_not_live_export',
                 'generator_basis_is_live_export':False}
        result=audit(mat,v6,witness,fixture)
        self.assertFalse(result['matchete_generator_fixture_is_live_native'])
        self.assertTrue(result['native_witness_is_historical_diagnostic_only'])
        self.assertFalse(result['independent_native_RGBeta_S2_to_Matchete_intertwiner_certified'])
        witness['RGBeta_numeric_delS2_components_certified']=True
        with self.assertRaisesRegex(S2EvidenceError,'unexpectedly changed'):
            audit(mat,v6,witness,fixture)
        witness['RGBeta_numeric_delS2_components_certified']=False
        fixture['generator_basis_is_live_export']=True
        with self.assertRaisesRegex(S2EvidenceError,'Unknown generator candidate provenance'):
            audit(mat,v6,witness,fixture)


if __name__=='__main__':
    unittest.main()
