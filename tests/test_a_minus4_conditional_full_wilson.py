"""Pure conditional boundary tests; production integration runs separately.

A local standalone run requires only numpy and this newly delivered module.
The actual full 45-Weyl ODE is exercised by the batch runner on the user's
current repository after 17-beta/provenance checks.
"""
from __future__ import annotations

import copy
import importlib.util
from pathlib import Path
import unittest

import numpy as np

HERE=Path(__file__).resolve().parents[1]
SOURCE=HERE/'Numerical/running/Aminus4ConditionalFullWilson.py'
spec=importlib.util.spec_from_file_location('a4cond',SOURCE)
lib=importlib.util.module_from_spec(spec)
spec.loader.exec_module(lib)


def synthetic_report():
    pos=[];neg=[]
    for p in (1,2,3):
        for i,j in ((1,1),(1,2),(2,2)):
            local_i,local_j=2*(p-1)+i,2*(p-1)+j
            for a,b in ((5,7),(5,8),(6,7),(6,8)):
                z=complex((p*10+i*3+j+a+b)*1e-14,
                          (i-j+1)*(b-6)*1e-15)
                base={'i':local_i,'j':local_j,'a':a,'b':b,
                      'distinct_orbit_size':lib._orbit_size(local_i,local_j,a,b)}
                pos.append({**base,'value_re_gev_inv':z.real,'value_im_gev_inv':z.imag})
                neg.append({**base,'value_re_gev_inv':-z.real,'value_im_gev_inv':-z.imag})
    return {'format':'t3-a-minus4-declared-weyl-uv-coupling-dictionary-v1',
            'model':'A_-4','independent_coordinate_count':36,
            'native_matchete_to_declared_Weyl_factor':None,
            'physical_C_ijab_at_MF':None,'physical_certification':False,
            'master_full_shape':[45,45,12,12],
            'declared_gamma_cc_choice':{'native_GammaCC_to_C0_measured':None},
            'declared_rho_plus_one_independent_coordinates':pos,
            'rho_minus_one_fixed_native_y1_alternative_coordinates':neg}


def indices():
    return {f'L[{p},{w}]':36+2*(p-1)+(w-1)
            for p in (1,2,3) for w in (1,2)}


class ConditionalBoundaryTests(unittest.TestCase):
    def test_full_45_weyl_embedding_and_pair_symmetries(self):
        c=lib.declared_tensor(synthetic_report(),indices())
        self.assertEqual(c.shape,(45,45,12,12))
        self.assertEqual(c[36,36,4,6],complex(2.6e-13,1e-15))
        self.assertTrue(np.array_equal(c,c.swapaxes(0,1)))
        self.assertTrue(np.array_equal(c,c.swapaxes(2,3)))
        self.assertFalse(np.any(c[:36]))
        self.assertFalse(np.any(c[:,:,0:4]))

    def test_minus_sign_branch_exact(self):
        doc=synthetic_report()
        cp=lib.declared_tensor(doc,indices(),rho=1)
        cm=lib.declared_tensor(doc,indices(),rho=-1)
        self.assertTrue(np.array_equal(cp,-cm))

    def test_nontrivial_cross_flavour_coordinate_maps_by_labels(self):
        doc=synthetic_report()
        # Swap the actual global position of first and last lepton coordinate.
        idx=indices();idx['L[1,1]'],idx['L[3,2]']=idx['L[3,2]'],idx['L[1,1]']
        original=lib.declared_tensor(doc,indices())
        c=lib.declared_tensor(doc,idx)
        self.assertEqual(c[41,41,4,6],original[36,36,4,6])
        self.assertEqual(c[36,36,4,6],original[41,41,4,6])
        self.assertNotEqual(c[36,36,4,6],original[36,36,4,6])

    def test_source_master_candidate_algebra(self):
        doc=synthetic_report();fi=indices()
        cp=lib.declared_tensor(doc,fi)
        rows=[{'i':r['i'],'j':r['j'],'a':r['a'],'b':r['b'],
               'distinct_orbit_size':r['distinct_orbit_size'],
               'coefficient_real_if_rho_1':r['value_re_gev_inv'],
               'coefficient_imag_if_rho_1':r['value_im_gev_inv']}
              for r in doc['declared_rho_plus_one_independent_coordinates']]
        master={'model':'A_-4','physical_certification':False,
                'physical_C_ijab_at_MF':None,
                'conditional_master_weight_basis_independent_entries':rows}
        self.assertEqual(lib.compare_master_source(doc,master,cp,fi),0)
        master['conditional_master_weight_basis_independent_entries'][0][
            'coefficient_real_if_rho_1']*=1.3
        with self.assertRaisesRegex(lib.ConditionalBlocked,'differs'):
            lib.compare_master_source(doc,master,cp,fi)

    def test_duplicate_source_coordinate_rejected(self):
        doc=synthetic_report()
        doc['declared_rho_plus_one_independent_coordinates'][2]=copy.deepcopy(
            doc['declared_rho_plus_one_independent_coordinates'][1])
        with self.assertRaisesRegex(lib.ConditionalBlocked,'Duplicated'):
            lib.declared_tensor(doc,indices())

    def test_inconsistent_rho_alternative_rejected(self):
        doc=synthetic_report()
        doc['rho_minus_one_fixed_native_y1_alternative_coordinates'][2][
            'value_re_gev_inv']+=1e-13
        with self.assertRaisesRegex(lib.ConditionalBlocked,'sign alternatives'):
            lib.declared_tensor(doc,indices())

    def test_refuse_physically_certified_claim_without_native_evidence(self):
        doc=synthetic_report();doc['physical_certification']=True
        with self.assertRaisesRegex(lib.ConditionalBlocked,'normalization status'):
            lib.declared_tensor(doc,indices())

    def test_wrong_orbit_normalisation_rejected(self):
        doc=synthetic_report()
        doc['declared_rho_plus_one_independent_coordinates'][0]['distinct_orbit_size']=4
        with self.assertRaisesRegex(lib.ConditionalBlocked,'orbit normalization'):
            lib.declared_tensor(doc,indices())

    def test_invalid_non_lepton_scalar_rejected(self):
        doc=synthetic_report()
        doc['declared_rho_plus_one_independent_coordinates'][0]['a']=1
        with self.assertRaisesRegex(lib.ConditionalBlocked,'coordinate keys'):
            lib.declared_tensor(doc,indices())


if __name__=='__main__':unittest.main()
