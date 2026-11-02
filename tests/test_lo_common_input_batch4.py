"""Regression tests for the isolated A_-4 dense-fitted common-input batch."""
from __future__ import annotations
from copy import deepcopy
from pathlib import Path
import tempfile
import unittest
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0,str(ROOT))

from Numerical.lo_certification.benchmark import common_input_batch4 as m


def cx(a):
    return {"real": np.asarray(a,float).tolist(),
            "imag": np.zeros_like(np.asarray(a,float)).tolist()}


def frozen_fixture():
    sm = {'gY':.38265, 'g2':.59493, 'g3':.7649, 'lambdaH':.25,
          'yu':cx(np.diag([1e-5,.0022,.6462])),
          'yd':cx(np.diag([1e-5,.0002,.0096])),
          'ye':cx(np.diag([2e-6,.0006,.01]))}
    ordinary={'y1':cx([[1,2,3],[4,5,6],[7,8,9]]),
              'y2':cx([[2,3,4],[5,6,7],[8,9,10]]),
              'MF':cx(np.diag([100000,110000,120000])),
              'mS1Sq':1.306e6,'mS2Sq':1e6,'lambdaS1':.1,'lambdaS2':.12,
              'lambdaH1':.02,'lambdaH2':.03,'lambda12':.01,
              'lambdaT3':{'real':.1,'imag':0.},
              'lambdaH2Adj':0.,'lambdaS2Adj':0.}
    return {'representation':{'d_s1':1,'d_s2':3,'d_f':2,'alpha':-4},
            'base_state':{'mu_gev':1e7,'sm':sm,'ordinary':ordinary},
            'scales':{'mu_fermion_threshold_gev':1e5,'mu_scalar_threshold_gev':1e3}}


class TestCommonInputBatch4(unittest.TestCase):
    def test_01_all_twenty_uv_couplings_pinned(self):
        f=frozen_fixture()
        keys=set(f['base_state']['sm'])|set(f['base_state']['ordinary'])
        self.assertEqual(len(keys),20)
        o=m.overrides_from_frozen(f,keys)
        self.assertEqual(set(o),keys)
        self.assertEqual(o['y1'][0][1],{'real':2.,'imag':0.})

    def test_02_missing_export_source_coupling_fails(self):
        f=frozen_fixture(); keys=set(f['base_state']['sm'])|set(f['base_state']['ordinary'])
        with self.assertRaisesRegex(m.Batch4Error,'export only'):
            m.overrides_from_frozen(f,keys|{'inventedQuartic'})

    def test_03_omitted_frozen_coupling_fails(self):
        f=frozen_fixture(); keys=set(f['base_state']['sm'])|set(f['base_state']['ordinary'])
        with self.assertRaisesRegex(m.Batch4Error,'frozen only'):
            m.overrides_from_frozen(f,keys-{'lambdaH'})

    def test_04_wrong_model_fails(self):
        f=frozen_fixture();f['representation']['alpha']=-2
        keys=set(f['base_state']['sm'])|set(f['base_state']['ordinary'])
        with self.assertRaisesRegex(m.Batch4Error,'representation'):
            m.overrides_from_frozen(f,keys)

    def test_05_wrong_threshold_fails(self):
        f=frozen_fixture();f['scales']['mu_scalar_threshold_gev']=100
        keys=set(f['base_state']['sm'])|set(f['base_state']['ordinary'])
        with self.assertRaisesRegex(m.Batch4Error,'threshold'):
            m.overrides_from_frozen(f,keys)

    def test_06_unsupported_complex_matrix_fails(self):
        with self.assertRaises(m.Batch4Error):
            m.z_matrix({'real':[[1,2]],'imag':[[0,0]]},'bad')

    def test_07_matrix_round_trip(self):
        a=np.eye(3,dtype=complex)*(1-2j);a[0,2]=.2+.5j
        self.assertTrue(np.array_equal(m.z_matrix(m.z_json(a),'matrix'),a))

    def test_08_diagonal_input_makes_diagonal_flavour(self):
        a=np.diag([1e-4,2e-4,3e-4]).astype(complex)
        w=m.invariant_flavour(a,a,np.array([1e5,1.1e5,1.2e5]))
        self.assertEqual(np.count_nonzero(w-np.diag(np.diag(w))),0)

    def test_09_dense_input_gives_offdiagonal_flavour(self):
        a=np.array([[1,2,3],[4,5,6],[7,8,9]],complex)*1e-4
        b=np.array([[2,3,4],[5,6,7],[8,9,10]],complex)*1e-4
        w=m.invariant_flavour(a,b,np.array([1e5,1.1e5,1.2e5]))
        self.assertGreater(abs(w[0,2]),1e-14)

    def test_10_conjugate_orientation_not_plain(self):
        a=np.eye(3,dtype=complex)*(1+2j)*1e-4
        b=np.eye(3,dtype=complex)*(3-1j)*1e-4
        w=m.invariant_flavour(a,b,np.array([1e5,1.1e5,1.2e5]))
        self.assertEqual(np.sign(w[0,0].imag),-np.sign(((1+2j)*(3-1j)).imag))
        self.assertTrue(np.allclose(w,(a.conj()/np.array([1e5,1.1e5,1.2e5])[None,:])@b.conj().T,rtol=0,atol=1e-25))
        self.assertFalse(np.allclose(w,(a/np.array([1e5,1.1e5,1.2e5])[None,:])@b.T,rtol=0,atol=1e-25))

    def test_11_zero_mass_rejected(self):
        with self.assertRaisesRegex(m.Batch4Error,'Nonpositive'):
            m.invariant_flavour(np.eye(3),np.eye(3),np.array([1.,0.,3.]))

    def test_12_bad_running_encoding_rejected(self):
        with self.assertRaises(m.Batch4Error):
            m.decode_state([[{'real':1,'imag':0}]],'bad')

    def test_13_boundary_is_new_uncertified_copy(self):
        old={'summed_independent_symmetric_coordinates':[{'coordinate':[1,1,1,1,5,7]}],
             'sources':[{'independent_symmetric_coordinates':[],
                         'flavour_permutation_lepton_order':[1,0]}],
             'spinor_normalization':{'certified':False,
                                     'gamma_cc_to_master_weyl_factor':None},
             'completeness':{'one_loop_F_hard_matching_included':False}}
        rows=[{'coordinate':[1,1,1,1,5,7]},{'coordinate':[1,1,1,2,5,7]}]
        doc=m.make_boundary(old,rows,{'source_term_inputform':'native'},'r','t','q')
        self.assertIsNone(doc['spinor_normalization']['gamma_cc_to_master_weyl_factor'])
        self.assertFalse(doc['completeness']['one_loop_F_hard_matching_included'])
        self.assertEqual(len(old['summed_independent_symmetric_coordinates']),1)
        self.assertEqual(len(doc['summed_independent_symmetric_coordinates']),2)
        self.assertNotIn('batch4_provenance',old)

    def test_14_no_old_boundary_overwrite(self):
        self.assertNotIn('releases/source_',str(m.OWNED))
        self.assertNotEqual(m.BENCH,m.OWNED)

    def test_15_frozen_dense_reference_numerical_witness(self):
        a=np.array([[-.00002618833696469664,.00005192426489666104,.000042472246941179035],
              [-.0000795598204666997,.00006316598653793336,-.0001337733481710395],
              [-.00004550523936687876,.00006617305918789167,.000037284407112747435]])
        b=np.array([[.00009910122072324159,.0000793926564976573,-.000021200984437018638],
              [.00004159613623782387,2.3139864206313102e-7,.000014121611182189383],
              [.000004319550376385441,-.00010891978349536659,.000044861549511551865]])
        w=m.invariant_flavour(a,b,np.array([1e5,1.1e5,1.2e5]))
        self.assertAlmostEqual(w[1,2],-1.159929776161694e-13,delta=5e-26)
        self.assertEqual(np.count_nonzero(w-np.diag(np.diag(w))),6)

    def test_16_missing_local_checkout_skips_integration(self):
        if not (ROOT/m.BENCH).is_file():
            self.skipTest('Full Windows repository not available in this build container')
        # Read-only preflight: tests canonical source exports and retained operators.
        result=m.preflight(ROOT)
        self.assertEqual(result['status'],'PREFLIGHT_OK_NO_OUTPUT_WRITTEN')
        self.assertEqual(len(result['complete_frozen_UV_parameter_keys']),20)


if __name__=='__main__':unittest.main()
