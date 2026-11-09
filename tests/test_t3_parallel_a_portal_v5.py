"""Synthetic ABI tests. The actual v5 Wolfram exporters require native Windows runs.

Triplet and doublet fixture values are drawn from previously saved registered
Matchete five-class CG component evidence, not manufactured RGBeta generators.
"""
from __future__ import annotations

from copy import deepcopy
from hashlib import sha256
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from Numerical.parallel_t3.wilson.portal_tensor_audit import (
    ALL16, MAT_FORMAT, RGB_FORMAT, REP, SOURCE_MAT, SOURCE_RGB,
    audit, audit_from_paths,
)
from Numerical.parallel_t3.wilson.matchete_rgbeta_dictionary import DictionaryEvidenceError


# A_+0: H2Inv1 anti-Hermitian 2x2x3x3 registered Matchete CG, and
# H2Inv2 real identity. This is a fixture, NOT an actual native v5 export.
TRIPLET = [
    '0','0','1/Sqrt[2]','0','0','0','-(1/Sqrt[2])','0','0',
    '0','I/Sqrt[2]','0','(-I)/Sqrt[2]','0','-(1/Sqrt[2])','0','1/Sqrt[2]','0',
    '0','(-I)/Sqrt[2]','0','I/Sqrt[2]','0','-(1/Sqrt[2])','0','1/Sqrt[2]','0',
    '0','0','-(1/Sqrt[2])','0','0','0','1/Sqrt[2]','0','0',
]
DOUBLE = [
    ['1','0','0','0','0','0','1','0','0','1','0','0','0','0','0','1'],
    ['1/Sqrt[3]','0','0','2/Sqrt[3]','0','0','-(1/Sqrt[3])','0',
     '0','-(1/Sqrt[3])','0','0','2/Sqrt[3]','0','0','1/Sqrt[3]'],
]


def fixture():
    m={'format':MAT_FORMAT,'models':{}}
    r={'format':RGB_FORMAT,'models':{}}
    count=0
    for tag in sorted(ALL16):
        cls=tag[0];d1,d2,df,alphas=REP[cls]
        alpha=int(tag.split('_')[1]);assert alpha in alphas
        mr=[];rr=[]
        for which,d in ((1,d1),(2,d2)):
            n=1 if d==1 else 2
            if d == 1:
                raw = [['1','0','0','1']]
            elif d == 2:
                raw = DOUBLE
            else:
                raw = [TRIPLET,[str(int(i==j and a==b))
                       for i in range(2) for j in range(2)
                       for a in range(3) for b in range(3)]]
            for k in range(n):
                mat_phase='I' if d==3 and k==0 else '1'
                coupling=f'lambdaH{which}Inv{k+1}'
                mr.append({'which':which,'invariant':k+1,'coupling':coupling,
                           'cg':f'T3HiggsPortal{which}CGInv{k+1}',
                           'dimensions':[2,2,d,d],
                           'native_dimensions':[2,2,d,d] if d>1 else [2,2],
                           'components_inputform':raw[k],
                           'phase_inside_plus_hc':mat_phase,
                           'native_hermiticity':('AntiHermitian' if mat_phase=='I'
                                                 else 'Hermitian'),
                           'native_tensor_accessor':'Matchete`GetCGTensor',
                           'assembled_lagrangian_contains_cg_and_coupling':True})
            for kind in (['singlet'] if d==1 else ['singlet','adjoint']):
                sym=f'lambdaH{which}' + ('Adj' if kind=='adjoint' else '')
                rr.append({'symbol':sym,'which':which,'kind':kind,
                           'dimensions':[2,2,d,d],
                           'fields_inputform':f'{{Bar[H], H, Bar[S{which}], S{which}}}',
                           'native_invariant_inputform':
                           ('tGen[SU2L[fund],A,#1,#2]*tGen[SU2L[S2],A,#3,#4]&'
                            if kind=='adjoint' else
                            'del[SU2L[fund],#1,#2]*del[SU2L[S2],#3,#4]&'),
                           'native_coupling_inputform':f'6*Matrix[{sym}][##1] &',
                           'native_self_conjugate':True,
                           'native_permutation_count':12})
            count+=n
        mixed = tag in {'B_-1','C_-1'}
        m['models'][tag]={'dimensions':[d1,d2,df],'alpha':alpha,
                          'native_build_status': ('RequiresMixedScalarMassDiagonalization'
                                                  if mixed else 'Success'),
                          'full_lagrangian_validation':not mixed,
                          'matching_mass_basis_ready':not mixed,
                          'mixed_scalar_mass_basis_blocked':mixed,
                          'additional_quadratics_heavy_mass_basis_blocked':int(mixed),
                          'portal_cg_operator_inventory_verified':True,
                          'portals':mr}
        r['models'][tag]={'dimensions':[d1,d2,df],'alpha':alpha,'portals':rr}
    m['models_requiring_mixed_scalar_mass_diagonalization']=['B_-1','C_-1']
    m['mixed_scalar_mass_basis_matching_authorized']=False
    assert count==58
    return m,r


class TestT3ParallelAPortalV5(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mat,cls.rgb=fixture()
        cls.result=audit(cls.mat,cls.rgb,sha_m='a'*64,sha_r='b'*64,samples=3)

    def test_all16_portals_fullrank(self):
        out=self.result
        self.assertEqual(out['portal_coupling_slots'],58)
        self.assertEqual(out['model_count'],16)
        self.assertEqual(out['native_anti_Hermitian_CGs_given_I_inside_PlusHc'],14)
        self.assertLess(out['max_relative_tensor_and_potential_residual'],1e-12)
        for model in out['models'].values():
            for sec in model['sectors']:
                self.assertEqual(sec['raw_CG_rank'],sec['corrected_real_Hermitian_rank'])

    def test_Aminus4_anti_Hermitian_correction_keeps_second_direction(self):
        row=self.result['models']['A_-4']['sectors'][1]
        self.assertEqual(row['native_matchete_phases_inside_half_PlusHc'],['I','1'])
        self.assertEqual(row['corrected_real_Hermitian_rank'],2)
        np.testing.assert_allclose(
            row['candidate_rgbeta_from_native_matchete_real_matrix'],
            [[0.,1.],[np.sqrt(2),0.]],atol=1e-12)
        self.assertEqual(row['candidate_orientation_sign_relative_to_native_RGBeta_tGen'],
                         'UNDETERMINED_PLUS_OR_MINUS')

    def test_B_doublet_mixing_is_not_one_to_one(self):
        row=self.result['models']['B_-1']['sectors'][0]
        np.testing.assert_allclose(
            row['candidate_rgbeta_from_native_matchete_real_matrix'],
            [[.5, np.sqrt(3)/2],[2.,-2/np.sqrt(3)]],atol=1e-12)

    def test_singlet_only_has_no_adjoint(self):
        row=self.result['models']['A_-4']['sectors'][0]
        self.assertEqual(row['scalar_dimension'],1)
        self.assertEqual(row['candidate_rgbeta_from_native_matchete_real_matrix'],[[1.]])
        self.assertIsNone(row['candidate_orientation_sign_relative_to_native_RGBeta_tGen'])

    def test_missing_I_is_rejected(self):
        m=deepcopy(self.mat)
        row=m['models']['A_-4']['portals'][1]
        row['phase_inside_plus_hc']='1'
        with self.assertRaisesRegex(DictionaryEvidenceError,'anti-Hermitian|PlusHc|phase'):
            audit(m,self.rgb,sha_m='a'*64,sha_r='b'*64,samples=1)

    def test_missing_portal_is_rejected(self):
        m=deepcopy(self.mat)
        m['models']['D_-2']['portals'].pop()
        with self.assertRaises(DictionaryEvidenceError):
            audit(m,self.rgb,sha_m='a'*64,sha_r='b'*64,samples=1)

    def test_rank_deficiency_is_rejected(self):
        m=deepcopy(self.mat)
        src=m['models']['A_-4']['portals']
        src[2]['components_inputform']=src[1]['components_inputform']
        src[2]['phase_inside_plus_hc']='I'
        src[2]['native_hermiticity']='AntiHermitian'
        with self.assertRaises(DictionaryEvidenceError):
            audit(m,self.rgb,sha_m='a'*64,sha_r='b'*64,samples=1)

    def test_native_rgbeta_changed_factor_is_rejected(self):
        r=deepcopy(self.rgb)
        r['models']['E_+2']['portals'][0]['native_coupling_inputform']='12*Matrix[lambdaH1][##1]&'
        with self.assertRaisesRegex(DictionaryEvidenceError,'prefactor'):
            audit(self.mat,r,sha_m='a'*64,sha_r='b'*64,samples=1)

    def test_native_rgbeta_changed_field_order_is_rejected(self):
        r=deepcopy(self.rgb)
        r['models']['B_-1']['portals'][0]['fields_inputform']='{H,Bar[H],Bar[S1],S1}'
        with self.assertRaisesRegex(DictionaryEvidenceError,'order'):
            audit(self.mat,r,sha_m='a'*64,sha_r='b'*64,samples=1)

    def test_source_sha256_contract_fails_stale(self):
        m=deepcopy(self.mat);r=deepcopy(self.rgb)
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)
            for path in SOURCE_MAT+SOURCE_RGB:
                dest=root/path;dest.parent.mkdir(parents=True,exist_ok=True)
                dest.write_text(f'{path}\n',encoding='utf-8')
            m['source_sha256']={p:sha256((root/p).read_bytes()).hexdigest()
                                for p in SOURCE_MAT}
            r['source_sha256']={p:sha256((root/p).read_bytes()).hexdigest()
                                for p in SOURCE_RGB}
            mat=root/'mat.json';rgb=root/'rgb.json'
            mat.write_text(json.dumps(m),encoding='utf-8')
            rgb.write_text(json.dumps(r),encoding='utf-8')
            result=audit_from_paths(mat,rgb,root,samples=1)
            self.assertEqual(result['portal_coupling_slots'],58)
            (root/SOURCE_MAT[-1]).write_text('CHANGED',encoding='utf-8')
            with self.assertRaisesRegex(DictionaryEvidenceError,'Stale'):
                audit_from_paths(mat,rgb,root,samples=1)

    def test_BC_mixed_scalar_mass_is_tensor_only(self):
        self.assertEqual(
            self.result['models_requiring_mixed_scalar_mass_diagonalization'],
            ['B_-1','C_-1'])
        self.assertFalse(self.result['mixed_scalar_mass_basis_matching_authorized'])
        for tag in ('B_-1','C_-1'):
            row=self.result['models'][tag]['matchete_native_build_validation']
            self.assertEqual(row['native_build_status'],
                             'RequiresMixedScalarMassDiagonalization')
            self.assertTrue(row['portal_CG_extraction_only'])
            self.assertFalse(row['matching_mass_basis_ready'])

    def test_unexpected_mass_basis_failure_is_rejected(self):
        m=deepcopy(self.mat)
        x=m['models']['A_-4']
        x.update(native_build_status='FullValidationFailed',
                 full_lagrangian_validation=False,
                 matching_mass_basis_ready=False)
        with self.assertRaisesRegex(DictionaryEvidenceError,'native Lagrangian status'):
            audit(m,self.rgb,sha_m='a'*64,sha_r='b'*64,samples=1)

    def test_BC_claim_of_mass_basis_readiness_is_rejected(self):
        m=deepcopy(self.mat)
        m['models']['B_-1']['matching_mass_basis_ready']=True
        with self.assertRaisesRegex(DictionaryEvidenceError,'mixed-mass'):
            audit(m,self.rgb,sha_m='a'*64,sha_r='b'*64,samples=1)

    def test_stale_mixed_model_inventory_is_rejected(self):
        m=deepcopy(self.mat)
        m['models_requiring_mixed_scalar_mass_diagonalization']=[]
        with self.assertRaisesRegex(DictionaryEvidenceError,'mixed-mass status inventory'):
            audit(m,self.rgb,sha_m='a'*64,sha_r='b'*64,samples=1)

    def test_cannot_claim_physical_certification(self):
        scope=self.result['scope']
        self.assertFalse(scope['physical_RGBeta_adjoint_orientation_sign_certified'])
        self.assertFalse(scope['native_RGBeta_numeric_generator_components_extracted'])
        self.assertFalse(scope['GammaCC_to_Weyl_normalization_certified'])


if __name__=='__main__': unittest.main()
