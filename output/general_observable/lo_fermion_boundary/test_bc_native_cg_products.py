"""Offline tests for numerical B/C native CG tensor decomposition.
Project path: output/general_observable/lo_fermion_boundary/test_bc_native_cg_products.py
"""
from __future__ import annotations
import unittest
import math
import json
import tempfile
from pathlib import Path
from copy import deepcopy
import numpy as np
import check_bc_native_cg_products as cg


def record(name,x):
    return {'name':name,'dimensions':list(x.shape),
            'representations_inputform':'{SU2L[fund], ...}',
            'real_flat':x.real.flatten().tolist(),
            'imag_flat':x.imag.flatten().tolist()}


def c_basis():
    q=np.zeros((2,2,2,2),dtype=complex)
    for i,j,a,b,v in [(0,0,0,0,1),(0,1,0,1,.5),(0,1,1,0,.5),
                        (1,0,0,1,.5),(1,0,1,0,.5),(1,1,1,1,1)]:
        q[i,j,a,b]=v
    eps=np.array([[0,1],[-1,0]],dtype=complex)
    bridge=eps.copy()
    e=np.einsum('ij,ba->ijab',eps,bridge)
    return {'T3ExtraCGQ020020I1':q,'T3MEBridgeBB':bridge,'eps[SU2L]':eps},q,e


def triplet_fixture():
    registry,q,e=c_basis()
    k=q/math.sqrt(3)+(math.sqrt(3)/2)*e
    # Takagi factor: symmetric (i,a),(j,b) matrix of rank three.
    a=np.array([[k[i,j,u,v] for j in range(2) for v in range(2)]
                for i in range(2) for u in range(2)])
    w,v=np.linalg.eigh(a)
    cols=[v[:,n]*np.sqrt(abs(w[n]))*(1 if w[n]>=0 else 1j)
          for n in range(4) if abs(w[n])>1e-10]
    assert len(cols)==3
    ybar=np.stack(cols,axis=1).reshape(2,2,3).transpose(0,2,1)
    y=ybar.conj()
    rec={'fermion_su2_dimension':3,
         'native_yukawa1_cg':record('T3MENativeY1CG',y),
         'native_yukawa2_cg':record('T3MENativeY2CG',y),
         'old_y2_epsilon_transport_alpha_inputform':'1',
         'uv_y2_explicit_prefactor_inputform':'-1',
         'native_y1_to_new_basis_alpha_inputform':'1'}
    return rec,registry


class BCProductTests(unittest.TestCase):
    def test_c_exact_two_tensor_decomposition(self):
        rec,reg=triplet_fixture()
        out=cg.analyze('C_-1',rec,reg)
        self.assertEqual(out['status'],'C_YukawaProductsMatchSignedQAndBridgeFactors')
        for v in out['native_yukawa_product_projections'].values():
            self.assertLess(v['relative_residual'],1e-11)
            self.assertAlmostEqual(v['real'][0],1/math.sqrt(3),places=10)
            self.assertAlmostEqual(v['real'][1],math.sqrt(3)/2,places=10)

    def test_c_wrong_bridge_orientation_rejected(self):
        rec,reg=triplet_fixture()
        reg['T3MEBridgeBB']=-reg['T3MEBridgeBB']
        with self.assertRaises(ValueError):cg.analyze('C_-1',rec,reg)

    def test_c_wrong_yukawa_phase_rejected(self):
        rec,reg=triplet_fixture()
        rec['native_yukawa2_cg']['real_flat']=[-x for x in rec['native_yukawa2_cg']['real_flat']]
        rec['native_yukawa2_cg']['imag_flat']=[-x for x in rec['native_yukawa2_cg']['imag_flat']]
        with self.assertRaises(ValueError):cg.analyze('C_-1',rec,reg)

    def test_c_corrupt_triplet_tensor_rejected(self):
        rec,reg=triplet_fixture()
        rec['native_yukawa1_cg']['real_flat'][0]+=0.1
        with self.assertRaises(ValueError):cg.analyze('C_-1',rec,reg)

    def test_b_identity_yukawa_cgs(self):
        x=np.eye(2,dtype=complex)
        rec={'fermion_su2_dimension':1,
             'native_yukawa1_cg':record('T3MENativeY1CG',x),
             'native_yukawa2_cg':record('T3MENativeY2CG',x),
             'old_y2_epsilon_transport_alpha_inputform':'-1',
             'uv_y2_explicit_prefactor_inputform':'1',
             'native_y1_to_new_basis_alpha_inputform':'1'}
        out=cg.analyze('B_-1',rec,{'T3MENativeY1CG':x,'T3MENativeY2CG':x})
        self.assertTrue(out['native_yukawa_cg_equal_to_matched_source'])

    def test_b_wrong_y2_tensor_rejected(self):
        x=np.eye(2,dtype=complex)
        rec={'fermion_su2_dimension':1,
             'native_yukawa1_cg':record('T3MENativeY1CG',x),
             'native_yukawa2_cg':record('T3MENativeY2CG',-x),
             'old_y2_epsilon_transport_alpha_inputform':'-1',
             'uv_y2_explicit_prefactor_inputform':'1',
             'native_y1_to_new_basis_alpha_inputform':'1'}
        with self.assertRaises(ValueError):
            cg.analyze('B_-1',rec,{'T3MENativeY1CG':x,'T3MENativeY2CG':x})

    def test_parse_actual_b_compressed_row_sparse_identity(self):
        # Exact serialization of native B_-1 T3MENativeY1CG/Y2CG registry.
        form='SparseArray[Automatic, {2, 2}, 0, {1, {{0, 1, 2}, {{1}, {2}}}, {1, 1}}]'
        a=cg.parse_matchete_sparse_array(form,tensor_name='B-matched')
        np.testing.assert_array_equal(a,np.eye(2))

    def test_parse_actual_c_four_index_compressed_row_cg(self):
        # Exact C_-1 Matchete registry serialization, 6 entries across 2 rows.
        form=('SparseArray[Automatic, {2, 2, 2, 2}, 0, '
              '{1, {{0, 3, 6}, {{1, 1, 1}, {2, 1, 2}, {2, 2, 1}, '
              '{1, 1, 2}, {1, 2, 1}, {2, 2, 2}}}, '
              '{1, 1/2, 1/2, 1/2, 1/2, 1}}]')
        a=cg.parse_matchete_sparse_array(form,tensor_name='C-Q')
        expected,_,_=c_basis()
        np.testing.assert_array_equal(a,expected['T3ExtraCGQ020020I1'])
        self.assertEqual(np.count_nonzero(a),6)

    def test_parse_actual_c_bridge_and_epsilon(self):
        form='SparseArray[Automatic, {2, 2}, 0, {1, {{0, 1, 2}, {{2}, {1}}}, {1, -1}}]'
        a=cg.parse_matchete_sparse_array(form,tensor_name='C-bridge')
        np.testing.assert_array_equal(a,np.array([[0,1],[-1,0]]))

    def test_parse_wrong_row_pointer_rejected(self):
        form='SparseArray[Automatic, {2, 2}, 0, {1, {{0, 1, 3}, {{1}, {2}}}, {1, 1}}]'
        with self.assertRaises(ValueError):
            cg.parse_matchete_sparse_array(form)
        form='SparseArray[Automatic, {2, 2}, 0, {1, {{0, 2, 1}, {{1}, {2}}}, {1, 1}}]'
        with self.assertRaises(ValueError):
            cg.parse_matchete_sparse_array(form)

    def test_parse_out_of_range_suffix_rejected(self):
        form='SparseArray[Automatic, {2, 2}, 0, {1, {{0, 1, 2}, {{1}, {3}}}, {1, 1}}]'
        with self.assertRaisesRegex(ValueError,'outside dimensions'):
            cg.parse_matchete_sparse_array(form)

    def test_parse_duplicate_sparse_element_rejected(self):
        form='SparseArray[Automatic, {2, 2}, 0, {1, {{0, 2, 2}, {{1}, {1}}}, {1, 1}}]'
        with self.assertRaisesRegex(ValueError,'duplicate'):
            cg.parse_matchete_sparse_array(form)

    def test_parse_unsupported_expression_rejected(self):
        form='SparseArray[Automatic, {2, 2}, 0, {1, {{0, 1, 2}, {{1}, {2}}}, {1, Sqrt[2]}}]'
        with self.assertRaises(ValueError):
            cg.parse_matchete_sparse_array(form)

    def test_reject_nonfinite_tensor(self):
        x=np.eye(2,dtype=complex)
        raw=record('T3MENativeY1CG',x)
        raw['real_flat'][0]=float('nan')
        with self.assertRaises(ValueError):cg.array(raw,'T3MENativeY1CG',(2,2))

    def test_reject_tensor_schema_mismatch(self):
        x=np.eye(2,dtype=complex)
        raw=record('T3MENativeY1CG',x)
        with self.assertRaises(ValueError):cg.array(raw,'T3MENativeY2CG',(2,2))

    def test_fresh_kernel_reports_join_without_phase_claim(self):
        import tempfile
        with tempfile.TemporaryDirectory() as tmp:
            paths={}
            for model in ('B_-1','C_-1'):
                path=Path(tmp)/(model+'.json')
                payload={
                    'format':'t3-bc-independent-yukawa-cg-components-v1',
                    'status':'IndependentNativeYukawaCGComponentsExtracted_SingleModel',
                    'model':model,
                    'record':{'fermion_su2_dimension':1 if model=='B_-1' else 3},
                    'master_weyl_phase_certified':False,
                    'native_BC_full_exchange_amplitude_computed':False,
                    'previous_release_changed':False,
                    'native_cg_product_computed_in_python':False,
                }
                path.write_text(json.dumps(payload),encoding='utf8')
                paths[model]=path
            obj,joined_sha,by_model=cg.load_separate_reports(paths)
            self.assertEqual(set(obj['models']),{'B_-1','C_-1'})
            self.assertEqual(len(joined_sha),64)
            self.assertEqual(set(by_model),set(paths))
            self.assertFalse(obj['master_weyl_phase_certified'])
            broken=json.loads(paths['C_-1'].read_text(encoding='utf8'))
            broken['model']='B_-1'
            paths['C_-1'].write_text(json.dumps(broken),encoding='utf8')
            with self.assertRaises(ValueError):
                cg.load_separate_reports(paths)

    def test_missing_fresh_kernel_model_rejected(self):
        with self.assertRaises(ValueError):
            cg.load_separate_reports({'B_-1':Path('nonexistent')})

    def test_master_phase_remains_open(self):
        rec,reg=triplet_fixture()
        self.assertFalse(cg.analyze('C_-1',rec,reg)['complete_feynman_amplitude_certified'])


if __name__=='__main__':unittest.main(verbosity=2)
