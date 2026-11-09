"""v6 SU2 Fierz and conditional Lie algebra sign tests.

Fixture CGs were exported from historical registered Matchete five-class data.
They are NOT substitutes for a fresh local all-16 native v5 run.
"""
from __future__ import annotations

from copy import deepcopy
from hashlib import sha256
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from tests.test_t3_parallel_a_portal_v5 import fixture
from Numerical.parallel_t3.wilson.portal_tensor_audit import audit as v5_audit
from Numerical.parallel_t3.wilson.portal_orientation_v6 import (
    FIERZ_FORMAT,FIERZ_STATUS,SOURCES_FIERZ,STATUS,_pauli,_fundamental_fierz,
    _triplet_sign_checks,_extract_triplet_generators,_identity,
    audit,audit_from_paths,
)
from Numerical.parallel_t3.wilson.matchete_rgbeta_dictionary import DictionaryEvidenceError
from Numerical.parallel_t3.wilson.portal_tensor_audit import SOURCE_MAT, SOURCE_RGB


def prepared():
    m,r=fixture()
    v=v5_audit(m,r,sha_m='a'*64,sha_r='b'*64,samples=3)
    v['input_sha256']={'matchete_portals':'a'*64,'rgbeta_portals':'b'*64}
    witness={'format':FIERZ_FORMAT,'status':FIERZ_STATUS,
             'native_fundamental_Fierz_exact':True,
             'native_symmetric_S2_rep_dimension_three':True,
             'native_S2_trace_normalization_two':True,
             'native_triplet_generators_symbolic':True,
             'native_barred_generator_minus_transpose':True,
             'native_S2_component_intertwiner_certified':False}
    return m,r,v,witness


class TestPortalV6(unittest.TestCase):
    def test_native_fierz_agrees_with_pauli(self):
        s=_pauli()
        z=np.einsum('aij,akl->ijkl',s,s)
        np.testing.assert_allclose(z,_fundamental_fierz(),atol=1.e-14)
        self.assertAlmostEqual(np.vdot(z,z).real,.75)

    def test_all16_doublet_and_triplet_orientation_fixture(self):
        m,r,v,w=prepared()
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)
            for p in SOURCES_FIERZ:
                path=root/p;path.parent.mkdir(parents=True,exist_ok=True)
                path.write_text(p,encoding='utf-8')
            w['source_sha256']={p:sha256((root/p).read_bytes()).hexdigest()
                                for p in SOURCES_FIERZ}
            out=audit(m,r,v,w,root,input_sha256={'matchete':'a'*64,
                        'rgbeta':'b'*64,'v5_audit':'c'*64,'native_fierz':'d'*64})
        self.assertEqual(out['status'],STATUS)
        self.assertEqual(out['doublet_sectors_Fierz_verified'],12)
        self.assertEqual(out['triplet_sectors_Pauli_algebra_conditional'],14)
        self.assertEqual(out['singlet_sectors'],6)
        self.assertLess(out['max_relative_component_and_potential_residual'],1.e-12)
        for cls in ('B','C'):
            for alpha in (-3,-1,1):
                for sec in out['models'][f'{cls}_{alpha:+d}']['sectors']:
                    np.testing.assert_allclose(sec['rgbeta_from_matchete_real_matrix'],
                         [[.5,np.sqrt(3)/2],[2.,-2/np.sqrt(3)]],atol=1.e-12)
                    self.assertTrue(sec['native_RGBeta_fundamental_fierz_verified'])
        for model in out['models'].values():
            for sec in model['sectors']:
                if sec['sector']=='triplet':
                    self.assertEqual(sec['selected_Pauli_frame_sign_relative_to_v5_adjoint'],-1)
                    self.assertLess(sec['selected_generator_checks']['max_residual'],1.e-12)
                    self.assertGreater(sec['wrong_sign_lie_bracket_residual'],1.e-3)
                    self.assertFalse(sec['native_RGBeta_S2_intertwiner_certified'])

    def test_minus_sign_violates_triplet_algebra(self):
        m,r,v,_=prepared()
        from Numerical.parallel_t3.wilson.portal_orientation_v6 import _projected_portals
        proj=_projected_portals(m['models']['A_-4']['portals'][1:],2,3)
        t=proj[0]
        candidate=3**.5*t/np.linalg.norm(t)
        yes=_triplet_sign_checks(candidate)
        no=_triplet_sign_checks(-candidate)
        self.assertLess(min(yes['max_residual'],no['max_residual']),1.e-12)
        self.assertGreater(max(yes['lie_bracket_residual'],no['lie_bracket_residual']),1.e-2)
        g=_extract_triplet_generators(-candidate)
        np.testing.assert_allclose(g[0]@g[1]-g[1]@g[0],1j*g[2],atol=1.e-12)
        np.testing.assert_allclose(sum(x@x for x in g),2*np.eye(3),atol=1.e-12)

    def test_fierz_witness_missing_identity_rejected(self):
        m,r,v,w=prepared();w['native_fundamental_Fierz_exact']=False
        with tempfile.TemporaryDirectory() as temp:
            with self.assertRaisesRegex(DictionaryEvidenceError,'Incomplete native Fierz'):
                audit(m,r,v,w,Path(temp),input_sha256={'matchete':'a'*64,
                      'rgbeta':'b'*64,'v5_audit':'c'*64,'native_fierz':'d'*64})

    def test_cannot_promote_triplet_to_native_certificate(self):
        m,r,v,w=prepared()
        w['native_S2_component_intertwiner_certified']=True
        with tempfile.TemporaryDirectory() as temp:
            with self.assertRaisesRegex(DictionaryEvidenceError,'Unsupported native S2'):
                audit(m,r,v,w,Path(temp),input_sha256={'matchete':'a'*64,
                      'rgbeta':'b'*64,'v5_audit':'c'*64,'native_fierz':'d'*64})

    def test_v5_input_hash_mismatch_rejected(self):
        m,r,v,w=prepared()
        v['input_sha256']['matchete_portals']='f'*64
        with tempfile.TemporaryDirectory() as temp:
            with self.assertRaisesRegex(DictionaryEvidenceError,'not derived'):
                # witness guard is prior; make it valid first.
                root=Path(temp)
                for p in SOURCES_FIERZ:
                    path=root/p;path.parent.mkdir(parents=True,exist_ok=True)
                    path.write_text(p,encoding='utf-8')
                w['source_sha256']={p:sha256((root/p).read_bytes()).hexdigest()
                                    for p in SOURCES_FIERZ}
                audit(m,r,v,w,root,input_sha256={'matchete':'a'*64,
                      'rgbeta':'b'*64,'v5_audit':'c'*64,'native_fierz':'d'*64})

    def test_edit_to_native_matrix_is_rejected_when_real_data_available(self):
        m,r,v,w=prepared()
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)
            for p in SOURCES_FIERZ:
                path=root/p;path.parent.mkdir(parents=True,exist_ok=True)
                path.write_text(p,encoding='utf-8')
            w['source_sha256']={p:sha256((root/p).read_bytes()).hexdigest()
                                for p in SOURCES_FIERZ}
            v['models']['B_-1']['sectors'][0][
                'candidate_rgbeta_from_native_matchete_real_matrix'][1][0]=100
            with self.assertRaisesRegex(DictionaryEvidenceError,'v5 oriented adjoint'):
                audit(m,r,v,w,root,input_sha256={'matchete':'a'*64,
                        'rgbeta':'b'*64,'v5_audit':'c'*64,'native_fierz':'d'*64})

    def test_source_witness_hash_change_rejected(self):
        m,r,v,w=prepared()
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)
            for p in SOURCES_FIERZ:
                path=root/p;path.parent.mkdir(parents=True,exist_ok=True)
                path.write_text(p,encoding='utf-8')
            w['source_sha256']={p:sha256((root/p).read_bytes()).hexdigest()
                                for p in SOURCES_FIERZ}
            (root/SOURCES_FIERZ[-1]).write_text('stale',encoding='utf-8')
            with self.assertRaisesRegex(DictionaryEvidenceError,'stale'):
                audit(m,r,v,w,root,input_sha256={'matchete':'a'*64,
                      'rgbeta':'b'*64,'v5_audit':'c'*64,'native_fierz':'d'*64})

    def test_actual_path_recheck_blocks_source_drifts(self):
        m,r,v,w=prepared()
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)
            for p in set(SOURCE_MAT+SOURCE_RGB+SOURCES_FIERZ):
                path=root/p;path.parent.mkdir(parents=True,exist_ok=True)
                path.write_text(p,encoding='utf-8')
            m['source_sha256']={p:sha256((root/p).read_bytes()).hexdigest()
                                for p in SOURCE_MAT}
            r['source_sha256']={p:sha256((root/p).read_bytes()).hexdigest()
                                for p in SOURCE_RGB}
            w['source_sha256']={p:sha256((root/p).read_bytes()).hexdigest()
                                for p in SOURCES_FIERZ}
            for n,obj in (('mat.json',m),('rgb.json',r),('fierz.json',w)):
                (root/n).write_text(json.dumps(obj),encoding='utf-8')
            v['input_sha256']={'matchete_portals':sha256((root/'mat.json').read_bytes()).hexdigest(),
                                'rgbeta_portals':sha256((root/'rgb.json').read_bytes()).hexdigest()}
            (root/'v5.json').write_text(json.dumps(v),encoding='utf-8')
            out=audit_from_paths(root/'mat.json',root/'rgb.json',root/'v5.json',
                                 root/'fierz.json',root)
            self.assertEqual(out['doublet_sectors_Fierz_verified'],12)
            (root/SOURCE_MAT[-1]).write_text('BAD',encoding='utf-8')
            with self.assertRaisesRegex(DictionaryEvidenceError,'Stale'):
                audit_from_paths(root/'mat.json',root/'rgb.json',root/'v5.json',
                                 root/'fierz.json',root)


if __name__=='__main__': unittest.main()
