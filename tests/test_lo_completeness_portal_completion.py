"""Workstream3 V3/V4 cross-stream portal completion tests (independent fixtures)."""
from __future__ import annotations

import copy
import tempfile
import unittest
from pathlib import Path
import json
import numpy as np

from Numerical.lo_certification.completeness.portal_completion_bridge import (
    read_ws2_completed_portals, original_to_boundary, rotate_hhss,
    merge_completed_baseline, parse_extra_reference, crosscheck_v4_projections,
    extra_unit_response, response)
from Numerical.lo_certification.completeness.source_bridge import (
    SourceBridgeError, FormalDecomposition)
from Numerical.lo_certification.completeness.hard_matching import weinberg_tensor


class PortalCompletionTest(unittest.TestCase):
    def setUp(self):
        self.n=8
        self.names=('p1','p2','missing-s2')
        self.native=('p1','p2')
        self.arr={}
        q=np.zeros((3,4,4,8,8),float)
        q[0,0,0,0,0]=.1
        q[1,1,1,2,2]=.2
        q[2,0,0,3,3]=.3
        for scale in ('muF','muS'):
            self.arr[scale+'_channel_labels']=np.array(self.names)
            self.arr[scale+'_channel_tensors']=q.copy()*(2 if scale=='muF' else 1)
            self.arr[scale+'_portal_sum']=self.arr[scale+'_channel_tensors'].sum(0)
        species=[{'species':'S1','status':'exact_full_portal_tensor_basis_in_declared_original_component_frame',
                  'native_symbols':['p1'],'added_complement_symbols':[],
                  'augmented_symbols':['p1'],'completion_count':0},
                 {'species':'S2','status':'exact_full_portal_tensor_basis_in_declared_original_component_frame',
                  'native_symbols':['p2'],'added_complement_symbols':['missing-s2'],
                  'augmented_symbols':['p2','missing-s2'],'completion_count':1}]
        self.status={'format':'t3-ws2-complete-baseline-portal-and-extra-reference-batch-v3',
           'model':'A_-4','status':'exact_baseline_portal_completed_extra_reference_basis_not_physically_mapped',
           'source_boundary_sha256':'s','source_sha256':{'running':'r'},
           'baseline_portal_exact_tensor_reconstruction_at_both_scales':True,
           'safe_for_existing_joint_scalar_matching_as_full_input':False,
           'full_HHSS_operator_tensor_and_couplings_certified':False,
           'physical_C5_calculated':False,
           'baseline_portal_completion':{
              'format':'t3-ws2-portal-completed-native-plus-rgbeta-v3',
              'higher_HHSS_and_lambdaT3_included':False,
              'source_scalar_basis':'original S1/S2 canonical real components, interleaved Re/Im',
              'quartic_normalization':'V4=lambda_abcd phi_a phi_b phi_c phi_d/4!',
              'potential_sign':'V4=-L4 quartic interaction',
              'baseline_missing_native_directions_filled':1,'species':species,
              'scales':{k:{'status':'exact_portal_numeric_reconstruction_conditional_component_frame',
                           'missing_rgbeta_couplings':[],
                           'channels_included':list(self.names)} for k in ('muF','muS')}}}

    def get(self,status=None,arrays=None):
        return read_ws2_completed_portals(status or self.status,arrays or self.arr,
            model='A_-4',source_hash='s',running_hash='r',n=8,native_names=self.native)

    def test_two_scale_exact_components_and_14_baseline_direction_semantics(self):
        p=self.get()
        self.assertEqual(p.baseline_missing_native_directions,1)
        self.assertEqual(p.names,self.names)
        np.testing.assert_allclose(p.mu_f_channels_original,2*p.mu_s_channels_original)

    def test_fail_missing_additional_portal_channel(self):
        bad=copy.deepcopy(self.status)
        bad['baseline_portal_completion']['species'][1]['added_complement_symbols']=[]
        with self.assertRaises(SourceBridgeError):self.get(status=bad)

    def test_fail_wrong_running_or_source_hash(self):
        bad=copy.deepcopy(self.status);bad['source_sha256']['running']='other'
        with self.assertRaises(SourceBridgeError):self.get(status=bad)
        bad=copy.deepcopy(self.status);bad['source_boundary_sha256']='other'
        with self.assertRaises(SourceBridgeError):self.get(status=bad)

    def test_fail_v3_quartic_channel_sum(self):
        bad=copy.deepcopy(self.arr);bad['muS_portal_sum']=bad['muS_portal_sum'].copy()
        bad['muS_portal_sum'][0,0,3,3]+=.01
        with self.assertRaises(SourceBridgeError):self.get(arrays=bad)

    def test_reject_physical_claims_and_missing_two_scales(self):
        for key in ('physical_C5_calculated','full_HHSS_operator_tensor_and_couplings_certified'):
            bad=copy.deepcopy(self.status);bad[key]=True
            with self.assertRaises(SourceBridgeError):self.get(status=bad)
        bad=copy.deepcopy(self.arr);bad.pop('muF_channel_tensors')
        with self.assertRaises(KeyError):self.get(arrays=bad)

    def test_special_model_basis_requires_express_opt_in(self):
        with self.assertRaises(SourceBridgeError):
            original_to_boundary('B_-1',8,R=np.eye(8),allow_exploratory_special=False)
        with self.assertRaises(SourceBridgeError):
            original_to_boundary('B_-1',8,R=None,allow_exploratory_special=True)
        self.assertEqual(original_to_boundary('A_-4',8,R=None,allow_exploratory_special=False).shape,(8,8))
        bad=np.eye(8);bad[1,2]=.1
        with self.assertRaises(SourceBridgeError):
            original_to_boundary('B_-1',8,R=bad,allow_exploratory_special=True)

    def test_rotate_original_to_boundary_and_preserve_nonportal(self):
        portal=self.get()
        old=np.zeros((3,4,4,8,8),float)
        old[0,0,0,0,0]=.2
        old[1,1,1,2,2]=.1
        old[2,0,0,7,7]=.4
        R=np.eye(8);R[[0,1]]=R[[1,0]]
        names,assembled,change=merge_completed_baseline(('p1','p2','mixing'),old,portal,R)
        self.assertEqual(names,('p1','p2','missing-s2','mixing'))
        self.assertEqual(change['preserved_nonportal_names'],['mixing'])
        self.assertEqual(change['added_baseline_directions'],1)
        np.testing.assert_allclose(assembled[-1],old[2])
        np.testing.assert_allclose(assembled[0],rotate_hhss(portal.mu_s_channels_original[:1],R)[0])
        self.assertNotAlmostEqual(assembled.sum(),(old.sum()+portal.mu_s_channels_original.sum()))

    def test_unknown_portal_native_must_not_be_double_counted(self):
        portal=self.get();old=np.zeros((2,4,4,8,8),float)
        with self.assertRaises(SourceBridgeError):
            merge_completed_baseline(('p1','mixing'),old,portal,np.eye(8))

    def test_extra_reference_exact_sqrt_and_bosonic_symmetry(self):
        e={'format':'t3-ws2-extra-hhss-exact-su2-reference-space-v1',
           'model':'A_-4','physical_C5_certified':False,
           'numerical_extra_couplings_available':False,
           'registry_HHSS_extra_symbols':1,
           'monomial_sectors':[{'counts_H_BarH_S1_BarS1_S2_BarS2':[1,1,0,0,0,2],
            'independent_complex_invariant_count':1,
            'source_registry_symbols':[{'rgbeta_symbol':'lambda_extra'}],
            'exact_reference_real_component_tensor_basis':[{'id':'reference_invariant_1',
             'Re_c_unit_tensor_components':[
                {'H1_H2_S1_S2_0based':[0,1,0,2],'value_exact':'sqrt(2)/2'}],
             'Im_c_unit_tensor_components':[
                {'H1_H2_S1_S2_0based':[1,1,2,2],'value_exact':'-2'}]}]}]}
        refs=parse_extra_reference(e,model='A_-4',n=8)
        self.assertEqual(len(refs),1)
        a=refs[0]['tensors']['Re']
        self.assertAlmostEqual(a[0,1,0,2],2**-.5)
        self.assertAlmostEqual(a[1,0,2,0],2**-.5)
        self.assertEqual(refs[0]['tensors']['Im'][1,1,2,2],-2.)
        v4={'format':'t3-ws2-source-named-rgbeta-HHSS-projection-v4','model':'A_-4',
            'two_thresholds_only':True,'physical_named_Matchete_RGBeta_maps_certified':0,
            'full_physical_eft1_quartics_certified':False,
            'source_projection_count':1,'orientation_ambiguous_symbols':1,
            'entries':[{'rgbeta_symbol':'lambda_extra',
                        'candidate_orientation_projection_maps':[{'in_exact_canonical_SU2_reference_span':True}]}]}
        self.assertEqual(crosscheck_v4_projections(v4,refs,model='A_-4')['extra_named_symbols'],1)
        v4['physical_named_Matchete_RGBeta_maps_certified']=1
        with self.assertRaises(SourceBridgeError):crosscheck_v4_projections(v4,refs,model='A_-4')
        e['numerical_extra_couplings_available']=True
        with self.assertRaises(SourceBridgeError):parse_extra_reference(e,model='A_-4',n=8)

    def test_zero_extra_reference_inventory(self):
        e={'format':'t3-ws2-extra-hhss-exact-su2-reference-space-v1',
           'model':'A_-4','physical_C5_certified':False,
           'numerical_extra_couplings_available':False,
           'registry_HHSS_extra_symbols':0,'monomial_sectors':[]}
        self.assertEqual(parse_extra_reference(e,model='A_-4',n=8),[])

    def test_unit_response_not_inserted_as_named_coefficient(self):
        n=8; c=np.zeros((3,3,3,2,2,n,n),complex)
        t=weinberg_tensor()
        for k in range(3):
            for a in range(4):
                for b in range(4):
                    c[k,k,k,:,:,a,b]=t[:,:,a,b]*1e-8/(k+1)
        formal=FormalDecomposition('A_-4',c.sum(0),c,np.array([1e5,1.1e5,1.2e5]),
           ('a',)*8,(1,),1e5,'s','r')
        q=np.zeros((4,4,n,n),float)
        for a in range(4):
            for b in range(4):q[a,b,a,b]+=.5;q[b,a,a,b]+=.5
        refs=[{'reference_id':'unit1','unassigned_source_symbols':['UNVERIFIED'],
               'sector_counts':None,'tensors':{'Re':q,'Im':-q}}]
        entries=extra_unit_response(formal,np.eye(8)*1e6,1e3,refs,np.eye(8))
        self.assertEqual(entries[0]['unit_parts']['Re']['status'],'CONDITIONAL_UNIT_REFERENCE_RESPONSE')
        self.assertFalse(entries[0]['physical_coefficient_assigned'])
        self.assertFalse(entries[0]['reference_basis_to_named_RGBeta_map_certified'])
        src=response(formal,np.eye(8)*1e6,1e3,('x',),q[None,...])
        self.assertTrue(src['projection_passed'])
        self.assertGreater(src['dim5_norm'],0)
        self.assertFalse(src['physical_C5_certified'])


if __name__=='__main__':unittest.main()


class RunnerSourceInventoryTests(unittest.TestCase):
    """Validate run-gate contracts without relying on GitHub or the local checkout."""
    def test_exact_all16_summary_requirements_and_no_overwrite(self):
        from scripts.RunLOPortalCompletionResponse import source_inventory, allowed_out
        from Numerical.lo_certification.completeness.inventory import ALL_MODELS
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp);v3=p/'v3';v4=p/'v4';v3.mkdir();v4.mkdir()
            def out():
                (v3/'summary.json').write_text(json.dumps({
                  'format':'t3-ws2-complete-baseline-portal-and-extra-reference-batch-v3',
                  'all16_model_audits_passed':True,'models_processed':16,
                  'baseline_portal_directions_completed':14,'extra_HHSS_registry_entries':24,
                  'physical_C5_calculated':False,
                  'models':[{'model':m} for m in ALL_MODELS]}))
                (v4/'summary.json').write_text(json.dumps({
                  'format':'t3-ws2-source-named-rgbeta-HHSS-projection-v4',
                  'all16_source_projections_passed':True,'models_processed':16,
                  'extra_HHSS_source_symbols_examined':24,
                  'physical_named_maps_certified':0,'physical_C5_calculated':False,
                  'models':[{'model':m} for m in ALL_MODELS]}))
            out()
            a=source_inventory(v3,v4,ALL_MODELS)
            self.assertEqual(a['v3_baseline_directions_added'],14)
            self.assertEqual(a['v3_extra_named_registry_symbols'],24)
            j=json.loads((v4/'summary.json').read_text())
            j['physical_named_maps_certified']=1
            (v4/'summary.json').write_text(json.dumps(j))
            with self.assertRaises(ValueError):source_inventory(v3,v4,ALL_MODELS)
        with self.assertRaises(ValueError):allowed_out(Path('scripts'))
        with self.assertRaises(ValueError):allowed_out(Path('output/general_observable/lo_certification_completeness'))
