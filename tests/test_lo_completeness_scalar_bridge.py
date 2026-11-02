"""Standalone scalar-input conversion tests; no Matchete or project files required."""
from __future__ import annotations
import copy
import unittest
import numpy as np

from Numerical.lo_certification.completeness.scalar_bridge import (
    parse_scalar_background, conditional_scalar_payload, SourceBridgeError)


class ScalarBridgeTests(unittest.TestCase):
    def setUp(self):
        n=8
        self.labels=tuple(['NewScalar1']*2+['NewScalar2']*6)
        self.m={'format':'t3-joint-scalar-conditional-background-v1',
                'model':'A_-4','source_boundary_sha256':'sourcehash',
                'source_scalar_real_labels':list(self.labels),
                'renormalization_scheme':'MSbar',
                'tensor_convention':'T3_REAL_LLSS_HHSS_WEINBERG_v1',
                'no_separate_scalar_threshold':True,
                'same_muS_quartic_and_mass':True,
                'scalar_mass_origin':'unbroken_electroweak',
                'physical_C5_certified':False,'benchmark_not_physical':True,
                'absolute_GammaCC_to_Weyl_factor':None,
                'operator_inventory_complete':False,
                'provenance':{'source':'v4 native seed and muS running',
                              'quartic_potential_normalization':
                              'V=lambda_abcd*phi_a*phi_b*phi_c*phi_d/4!'},
                'mu_f_gev':100000.,'mu_s_gev':1000.,'reference_scale_gev':1000.,
                'quartic_channel_names':['portal-a','portal-b']}
        q=np.zeros((4,4,n,n),float)
        q[0,0,0,0]=.5
        q[1,1,2,2]=.125
        q[0,1,0,2]=.3
        q[1,0,2,0]=.3
        q[1,0,0,2]=.3
        q[0,1,2,0]=.3
        x=np.zeros_like(q);x[0,0,0,0]=.5
        self.z={'source_to_matching_real_U':np.eye(n),
                'mass_squared_gev2_muS':np.diag(np.linspace(1e6,1.2e6,n)),
                'lambda_HHSS_muS':q,
                'lambda_HHSS_channels_muS':np.stack([x,q-x])}

    def parse(self,m=None,z=None):
        return parse_scalar_background(self.m if m is None else m,
                     self.z if z is None else z,model='A_-4',
                     source_sha256='sourcehash',labels=self.labels,
                     mu_f_gev=100000.,meta_hash='metahash',arrays_hash='npzhash')

    def test_valid_independent_quartics_no_missing_channels(self):
        r=self.parse()
        self.assertEqual(r.lambda_channels_source.shape,(2,4,4,8,8))
        np.testing.assert_allclose(r.lambda_source,r.lambda_channels_source.sum(axis=0))
        meta,data=conditional_scalar_payload(r)
        self.assertFalse(meta['physical_C5_certified'])
        self.assertIsNone(meta['provenance']['GammaCC_to_Weyl_signed_factor'])
        self.assertEqual(meta['source_channel_coverage_scope'],
                         'complete_cached_native_HHSS_only_NOT_general_UV_closure')
        np.testing.assert_allclose(data['scalar_mass_squared_gev2'],self.z['mass_squared_gev2_muS'])

    def test_nontrivial_real_basis_pullback_rotates_mass_and_every_channel(self):
        z=copy.deepcopy(self.z)
        U=np.eye(8);U[[0,1],:]=U[[1,0],:]
        U[4,4]=-1
        z['source_to_matching_real_U']=U
        d=self.parse(z=z)
        np.testing.assert_allclose(d.mass2_source,U@z['mass_squared_gev2_muS']@U.T)
        np.testing.assert_allclose(d.lambda_source,np.einsum('ABab,ia,jb->ABij',z['lambda_HHSS_muS'],U,U))
        np.testing.assert_allclose(d.lambda_source,d.lambda_channels_source.sum(axis=0))

    def test_refuses_unmatched_source_or_flavour_field_basis(self):
        m=copy.deepcopy(self.m);m['source_boundary_sha256']='wrong'
        with self.assertRaises(SourceBridgeError):self.parse(m=m)
        m=copy.deepcopy(self.m);m['source_scalar_real_labels'][0]='guess'
        with self.assertRaises(SourceBridgeError):self.parse(m=m)

    def test_refuses_two_scalar_thresholds_or_inconsistent_mu(self):
        m=copy.deepcopy(self.m);m['no_separate_scalar_threshold']=False
        with self.assertRaises(SourceBridgeError):self.parse(m=m)
        m=copy.deepcopy(self.m);m['reference_scale_gev']=1200
        with self.assertRaises(SourceBridgeError):self.parse(m=m)

    def test_rejects_quartic_channel_omission_even_when_small(self):
        z=copy.deepcopy(self.z)
        z['lambda_HHSS_channels_muS'][1]*=0.9999
        with self.assertRaises(SourceBridgeError):self.parse(z=z)

    def test_rejects_unverified_physical_claims(self):
        m=copy.deepcopy(self.m);m['physical_C5_certified']=True
        with self.assertRaises(SourceBridgeError):self.parse(m=m)
        m=copy.deepcopy(self.m);m['absolute_GammaCC_to_Weyl_factor']=1
        with self.assertRaises(SourceBridgeError):self.parse(m=m)

    def test_rejects_bad_mass_matrix_and_rotation(self):
        z=copy.deepcopy(self.z);z['mass_squared_gev2_muS'][0,0]=-1
        with self.assertRaises(SourceBridgeError):self.parse(z=z)
        z=copy.deepcopy(self.z);z['source_to_matching_real_U'][0,1]=0.2
        with self.assertRaises(SourceBridgeError):self.parse(z=z)

    def test_rejects_missing_channels_and_EWSB_resummation(self):
        z=copy.deepcopy(self.z);z.pop('lambda_HHSS_channels_muS')
        with self.assertRaises(SourceBridgeError):self.parse(z=z)
        m=copy.deepcopy(self.m);m['scalar_mass_origin']='EWSB_resummed'
        with self.assertRaises(SourceBridgeError):self.parse(m=m)


    def test_extracted_bundle_to_existing_hard_interface_complex_flavour(self):
        """True numeric end-to-end adapter test with no upstream UV claims."""
        from Numerical.lo_certification.completeness.source_bridge import (
            FormalDecomposition, prepare_conditional_input)
        from Numerical.lo_certification.completeness.hard_matching import weinberg_tensor
        from scripts.RunLOCompletenessAudit import compute_conditional
        from pathlib import Path
        import tempfile, json
        n=8
        m=copy.deepcopy(self.m)
        z=copy.deepcopy(self.z)
        q=np.zeros((4,4,n,n),float)
        c=np.zeros((3,3,3,2,2,n,n),complex)
        t=weinberg_tensor()
        for a in range(4):
            for b in range(4):
                q[a,b,a,b]+=.5
                q[b,a,a,b]+=.5
                for k in range(3):
                    # Complex diagonal Yukawa benchmark; no phase is lost.
                    c[k,k,k,:,:,a,b]=t[:,:,a,b]*((k+1)+1j*(k+.1))*1e-8
        z['lambda_HHSS_muS']=q
        z['lambda_HHSS_channels_muS']=np.stack([.001*q,.999*q])
        z['mass_squared_gev2_muS']=np.eye(n)*1e6
        bundle=self.parse(m=m,z=z)
        md,arrays=conditional_scalar_payload(bundle)
        formal=FormalDecomposition('A_-4',c.sum(axis=0),c,
                np.array([100000.,110000.,120000.]),self.labels,(1,2),
                100000.,'sourcehash','runhash')
        meta,out=prepare_conditional_input(formal,md,arrays,'npzhash')
        self.assertFalse(meta['physical_C5_certified'])
        self.assertEqual(out['C_LLSS_by_fermion'].shape,(3,3,3,2,2,8,8))
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)
            np.savez_compressed(p/'A_m4.npz',**out)
            (p/'A_m4.json').write_text(json.dumps(meta),encoding='utf-8')
            diagnostic=compute_conditional('A_-4',p)
        self.assertEqual(diagnostic['status'],'CONDITIONAL_NUMERICAL_NOT_PHYSICS_CERTIFIED')
        self.assertTrue(diagnostic['projection_passed'])
        self.assertFalse(diagnostic['physical_C5_certified'])
        self.assertEqual(len(diagnostic['quartic_channels']),2)
        hard=diagnostic['matrices']['hard_muF']
        self.assertGreater(abs(hard['imag'][0][0]),0.)
        self.assertEqual(hard['real'][0][1],0.)


if __name__=='__main__':unittest.main()
