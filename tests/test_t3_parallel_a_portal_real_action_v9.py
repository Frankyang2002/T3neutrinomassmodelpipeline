"""v9 conditional real-scalar action: source-factor and spectral negative controls.

The 16-model integration tests require the user's local repository/exports.
Standalone representation and factor tests do not pretend to run RGBeta Wolfram.
"""
from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
import unittest

import numpy as np

from Numerical.parallel_t3.wilson.portal_real_action_v9 import (
    V9EvidenceError, fundamental_fierz, identity, canonical_lam, complex_action,
    lam_action, source_half_plus_hc, spectral_triplet_adjoint,
    audit_sector, audit, read,
)


def metadata(which:int, dim:int, kind:str) -> dict:
    sym=f'lambdaH{which}'+('Adj' if kind=='adjoint' else '')
    return {'symbol':sym,'which':which,'kind':kind,
            'dimensions':[2,2,dim,dim],'native_self_conjugate':True,
            'native_permutation_count':24,
            'native_coupling_inputform':f'6*Matrix[{sym}][##1] & ',
            'native_invariant_inputform':('tGen[SU2L[fund],A,#1,#2]*tGen[SU2L[S2],A,#3,#4] &'
                 if kind=='adjoint' else 'del[SU2L[fund], #1,#2]*del[SU2L[S2],#3,#4] &'),
            'native_numerical_generator_components_exported':False}


def source(which:int,dim:int,k:int,phase:str='1') -> dict:
    return {'which':which,'invariant':k,'coupling':f'lambdaH{which}Inv{k}',
            'phase_inside_plus_hc':phase,'dimensions':[2,2,dim,dim]}


def spin_one_adjoint() -> np.ndarray:
    # Declared highest-weight basis; NOT an RGBeta numeric S2 frame witness.
    p=np.array([[0.,1.],[1.,0.]],complex)/2
    y=np.array([[0.,-1j],[1j,0.]],complex)/2
    z=np.diag([.5,-.5]).astype(complex)
    rise=np.array([[0.,np.sqrt(2),0.],[0.,0.,np.sqrt(2)],[0.,0.,0.]],complex)
    g=np.array([(rise+rise.conj().T)/2,(rise-rise.conj().T)/(2j),
                np.diag([1.,0.,-1.])])
    return np.einsum('aij,auv->ijuv',np.array([p,y,z]),g)


class TestRealActionMath(unittest.TestCase):
    def test_canonical_real_lam_for_all_three_dimensions(self):
        rng=np.random.default_rng(901)
        for d in (1,2,3):
            t=identity(d)
            for _ in range(4):
                x=rng.normal(size=4+2*d)
                self.assertLess(abs(lam_action(canonical_lam(t),x)-
                                    complex_action(t,x)),1.e-11)

    def test_missing_six_prefactor_is_not_equal_to_complex_action(self):
        t=identity(3)
        x=np.array([.3,.1,-.2,.7,1.,0.,0.,.9,.5,-.1])
        wanted=complex_action(t,x)
        self.assertGreater(abs(wanted),.01)
        bad=lam_action(canonical_lam(t)/6,x)
        self.assertGreater(abs(wanted-bad),.01)

    def test_doublet_fierz_normalization(self):
        k=fundamental_fierz()
        self.assertAlmostEqual(np.vdot(k,k).real,.75,places=12)
        self.assertAlmostEqual(np.vdot(k,identity(2)).real,0,places=12)

    def test_native_like_doublet_sector_real_action(self):
        d=2;which=1
        k=fundamental_fierz();i=identity(d)
        tensors=[.5*i+2*k,np.sqrt(3)/2*i-2/np.sqrt(3)*k]
        r=[source(which,d,j+1) for j in range(2)]
        regs=[metadata(which,d,'singlet'),metadata(which,d,'adjoint')]
        v6={'sector':'doublet','rgbeta_from_matchete_real_matrix':
             [[.5,np.sqrt(3)/2],[2.,-2/np.sqrt(3)]]}
        row=audit_sector(tensors,tensors,r,regs,v6,tag='B_-1',which=which,d=d)
        self.assertLess(row['max_relative_residual'],1.e-12)
        self.assertFalse(row['native_RGBeta_Lam_projector_executed'])

    def test_triplet_phase_and_spectral_sign_without_v6_choice(self):
        d=3;which=2
        k=spin_one_adjoint();i=identity(d)
        phased=[-np.sqrt(2)*k,i]
        raw=[-1j*phased[0],phased[1]]
        g,sg,err=spectral_triplet_adjoint(phased)
        self.assertEqual(sg,-1)
        self.assertLess(err,1.e-12)
        np.testing.assert_allclose(g,k,atol=1.e-13)
        records=[source(which,d,1,'I'),source(which,d,2)]
        regs=[metadata(which,d,'singlet'),metadata(which,d,'adjoint')]
        v6={'sector':'triplet','candidate_rgbeta_from_matchete_real_matrix':
             [[0.,1.],[-np.sqrt(2),0.]]}
        v8={'source_CG_sign_selected_without_using_v6_map':True,
            'sign_relative_to_first_native_traceless_CG':-1,
            'selected_spectral_checks':{'max_residual':1.e-15},
            'RGBeta_numeric_triplet_frame_certified':False}
        result=audit_sector(raw,phased,records,regs,v6,tag='A_-4',which=which,d=d,
                            v8_sector=v8)
        self.assertLess(result['max_relative_residual'],1.e-12)
        self.assertEqual(result['native_Matchete_phase_inside_half_PlusHc'],['I','1'])

    def test_incorrect_native_rgbeta_prefactor_rejected(self):
        i=identity(1)
        rg=metadata(1,1,'singlet');rg['native_coupling_inputform']='12*Matrix[lambdaH1][##1] &'
        with self.assertRaisesRegex(V9EvidenceError,'6\\*Matrix'):
            audit_sector([i],[i],[source(1,1,1)],[rg],
                {'sector':'singlet','rgbeta_from_matchete_real_matrix':[[1.]]},
                tag='A_-4',which=1,d=1)

    def test_incorrect_perm_count_rejected(self):
        i=identity(1);rg=metadata(1,1,'singlet');rg['native_permutation_count']=12
        with self.assertRaisesRegex(V9EvidenceError,'unique arrangements'):
            audit_sector([i],[i],[source(1,1,1)],[rg],
                {'sector':'singlet','rgbeta_from_matchete_real_matrix':[[1.]]},
                tag='A_-4',which=1,d=1)

    def test_missing_triplet_i_phase_rejected(self):
        k=spin_one_adjoint()
        raw=[1j*k,identity(3)]
        phased=[k,identity(3)]
        with self.assertRaisesRegex(V9EvidenceError,'PlusHc projection mismatch'):
            audit_sector(raw,phased,[source(2,3,1),source(2,3,2)],
                [metadata(2,3,'singlet'),metadata(2,3,'adjoint')],
                {'sector':'triplet','candidate_rgbeta_from_matchete_real_matrix':[[0,1],[1,0]]},
                tag='A_-4',which=2,d=3)

    def test_wrong_v8_sign_rejected(self):
        k=spin_one_adjoint();i=identity(3)
        phased=[-np.sqrt(2)*k,i];raw=[-1j*phased[0],i]
        with self.assertRaisesRegex(V9EvidenceError,'saved v8 independent spectral sign'):
            audit_sector(raw,phased,[source(2,3,1,'I'),source(2,3,2)],
                [metadata(2,3,'singlet'),metadata(2,3,'adjoint')],
                {'sector':'triplet','candidate_rgbeta_from_matchete_real_matrix':
                    [[0,1],[-np.sqrt(2),0]]},tag='A_-4',which=2,d=3,
                v8_sector={'source_CG_sign_selected_without_using_v6_map':True,
                           'sign_relative_to_first_native_traceless_CG':1,
                           'selected_spectral_checks':{'max_residual':0},
                           'RGBeta_numeric_triplet_frame_certified':False})

    def test_wrong_v6_matrix_rejected(self):
        i=identity(1)
        with self.assertRaisesRegex(V9EvidenceError,'v6 transformation'):
            audit_sector([i],[i],[source(1,1,1)],[metadata(1,1,'singlet')],
                {'sector':'singlet','rgbeta_from_matchete_real_matrix':[[2.]]},
                tag='A_-4',which=1,d=1)


class TestRepositoryNativeExports(unittest.TestCase):
    def setUp(self):
        self.root=Path(__file__).resolve().parents[1]
        self.directory=self.root/'output/parallel_t3/wilson'
        needed=['matchete_portals_native_v5.json','rgbeta_portals_native_v5.json',
                'portal_orientation_conditional_v6.json','s2_spectrum_conditional_v8.json']
        if not all((self.directory/f).is_file() for f in needed):
            self.skipTest('Current local full-repository v5-v8 JSON not mounted')
        self.mat=read(self.directory/needed[0]);self.rgb=read(self.directory/needed[1])
        self.v6=read(self.directory/needed[2]);self.v8=read(self.directory/needed[3])

    def test_16_models_58_native_portal_actions(self):
        result=audit(self.mat,self.rgb,self.v6,self.v8,samples=2)
        self.assertEqual(result['native_Matchete_portal_coupling_slots'],58)
        self.assertEqual(result['triplet_sectors_spectral_normalization_verified'],14)
        self.assertEqual(result['anti_Hermitian_Matchete_phase_slots_checked'],14)
        self.assertLess(result['max_constructed_real_action_relative_residual'],2.e-10)
        self.assertFalse(result['certification']['native_RGBeta_Lam_projector_numerically_executed'])

    def test_live_registration_tampering_rejected(self):
        altered=deepcopy(self.rgb)
        altered['models']['A_-4']['portals'][0]['native_coupling_inputform'] = (
            '12*Matrix[lambdaH1][##1] &')
        with self.assertRaisesRegex(V9EvidenceError,'6\\*Matrix'):
            audit(self.mat,altered,self.v6,self.v8,samples=1)


if __name__=='__main__':unittest.main()
