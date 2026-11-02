"""Independent fixed-input and source-evidence gates, no project files needed."""
from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest
import numpy as np

from Numerical.lo_certification.completeness.consistency_audit import (
    audit_model, assemble, crossworkstream_evidence, ConsistencyError,
    relative_norm_error, decode_matrix,
)
from Numerical.lo_certification.completeness.hard_matching import (
    evaluate_hard_scalar_channel, weinberg_tensor,
)
from Numerical.lo_certification.completeness.inventory import ALL_MODELS
from Numerical.lo_certification.completeness.source_bridge import digest


def example_tensors():
    n = 8
    c = np.zeros((1,3,3,2,2,n,n), complex)
    q = np.zeros((4,4,n,n), float)
    t = weinberg_tensor()
    for a in range(4):
        for b in range(4):
            q[a,b,a,b] += 0.5
            q[b,a,a,b] += 0.5
            c[0,:,:,:,:,a,b] = (1e-8+1j*2e-9)*t[:,:,a,b][None,None,:,:]*np.eye(3)[:,:,None,None]
    return c, np.array([1e5]), q, np.eye(n)*1e6


def enc(a):
    return {'real': a.real.tolist(), 'imag': a.imag.tolist()}


def write_json(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2) + '\n', encoding='utf-8')


class ConsistencyAuditTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.run = Path(self.temp.name)
        m = 'A_m4'
        c, masses, q, ms2 = example_tensors()
        self.arrays = {
            'C_LLSS_by_fermion': c,
            'C_LLSS_total_muF': c.sum(axis=0),
            'fermion_masses_gev': masses,
            'lambda_HHSS': q,
            'lambda_HHSS_channels': np.stack([q*.4, q*.6]),
            'scalar_mass_squared_gev2': ms2,
        }
        self.bgmeta = self.run/'scalar_backgrounds'/f'{m}.json'
        self.bgnpz = self.run/'scalar_backgrounds'/f'{m}.npz'
        self.bgmeta.parent.mkdir(parents=True,exist_ok=True)
        write_json(self.bgmeta, {'model':'A_-4'})
        np.savez_compressed(self.bgnpz,x=np.eye(8))
        self.input_path = self.run/'conditional_inputs'/f'{m}.json'
        self.input_path.parent.mkdir(parents=True,exist_ok=True)
        self.input_npz = self.run/'conditional_inputs'/f'{m}.npz'
        np.savez_compressed(self.input_npz,**self.arrays)
        self.meta = {
            'model':'A_-4','scheme':'MSbar',
            'tensor_convention':'T3_REAL_LLSS_HHSS_WEINBERG_v1',
            'source_normalization_certified':False,
            'physical_C5_certified':False,
            'mu_f_gev':1e5,'mu_s_gev':1e3,
            'tensor_npz':f'{m}.npz',
            'quartic_channel_names':['q1','q2'],
            'source_scalar_background_metadata_sha256':digest(self.bgmeta),
            'source_scalar_background_npz_sha256':digest(self.bgnpz),
            'provenance':{'scalar_bundle_sha256':digest(self.bgnpz),
                          'fermion_source_sha256':'independent_llss_source',
                          'running_sha256':'independent_running_source'},
        }
        write_json(self.input_path,self.meta)
        r = evaluate_hard_scalar_channel(c,masses,q,ms2,1e5,1e3,
                lambda_channels=self.arrays['lambda_HHSS_channels'],
                channel_names=('q1','q2'),C_total_muF=c.sum(axis=0))
        self.saved = {
            'model':'A_-4','projection_passed':r['projection_passed'],
            'physical_C5_certified':False,
            'all_fermions_one_muF':True,
            'both_scalars_one_muS':True,
            'other_six_EFT1_RGE_sources_included':False,
            'complete_hard_one_loop_matching_certified':False,
            'matrices': {k:enc(x) for k,x in r['matrices'].items()},
        }
        self.diagnostic = self.run/'conditional_hard'/f'{m}.json'
        write_json(self.diagnostic,self.saved)

    def test_model_hard_rg_soft_partition_and_power_scaling(self):
        r = audit_model('A_-4',self.run)
        self.assertEqual(r['status'],'PORTAL_NUMERIC_CHECKS_PASS')
        self.assertLess(r['matching_scale_muF_relative_drift_fixed_inputs'],1e-10)
        self.assertLess(r['matching_scale_muS_relative_drift_fixed_inputs'],1e-10)
        self.assertGreater(r['relative_exact_triangle_truncation_error'],0)
        self.assertTrue(.5<r['scalar_mass2_power_exponent_1_to_1quarter']<1.5)
        self.assertFalse(r['physical_C5_certified'])

    def test_changed_scalar_source_sha_is_rejected(self):
        self.bgnpz.write_bytes(self.bgnpz.read_bytes()+b'x')
        with self.assertRaisesRegex(ConsistencyError,'changed'):
            audit_model('A_-4',self.run)

    def test_changed_conditional_matrix_is_rejected(self):
        self.saved['matrices']['hard_muF']['real'][0][0] += 0.00001
        write_json(self.diagnostic,self.saved)
        with self.assertRaisesRegex(ConsistencyError,'disagrees'):
            audit_model('A_-4',self.run)

    def test_wrong_fermion_threshold_flag_is_rejected(self):
        self.saved['both_scalars_one_muS']=False
        write_json(self.diagnostic,self.saved)
        with self.assertRaisesRegex(ConsistencyError,'threshold'):
            audit_model('A_-4',self.run)

    def test_signed_physical_claim_is_rejected(self):
        self.meta['source_normalization_certified']=True
        write_json(self.input_path,self.meta)
        with self.assertRaisesRegex(ConsistencyError,'normalization claim'):
            audit_model('A_-4',self.run)

    def test_nonfinite_saved_matrix_is_rejected(self):
        self.saved['matrices']['hard_muF']['real'][0][0] = float('nan')
        write_json(self.diagnostic,self.saved)
        with self.assertRaisesRegex(ConsistencyError,'finite'):
            audit_model('A_-4',self.run)

    def test_decode_complex_and_error_function(self):
        a=decode_matrix(enc(np.eye(3,dtype=complex)*1j))
        self.assertEqual(a[1,1],1j)
        self.assertEqual(relative_norm_error(a,a),0)


class CrossWorkstreamTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.d = Path(self.temp.name)
        self.portal = self.d/'portal.json'
        self.legacy = self.d/'legacy.json'
        self.gamma = self.d/'gamma.json'
        write_json(self.portal, {
            'format':'t3-ws2-exact-portal-span-audit-v2',
            'all16_audit_passed':True,
            'models':[{'model':m,
                       'RGBeta_portal_directions_missing_from_native_seed':1 if m=='A_-4' else 0,
                       'extra_HHSS_candidate_count':2 if m=='B_-1' else 0,
                       'blockers':['uncertified']}
                      for m in ALL_MODELS],
        })
        write_json(self.legacy, {
            'format':'t3-exploratory-all16-conditional-c5-benchmark-v6',
            'models':[{'model':m,
                       'status':'ERROR_unverified_source_or_physics_consistency' if m=='A_-2' else 'conditional',
                       'error':'full LLhh amplitude violates the Weinberg Ward projection' if m=='A_-2' else ''}
                      for m in ALL_MODELS],
        })
        write_json(self.gamma, {
            'format':'t3-bc-conditional-majorana-mass-normalization-v1',
            'native_GammaCC_to_Warsaw_prefactor_certified':False,
            'installed_source_probe':{'status':'GammaCCSourceLocated_NotYetInterpreted'},
        })

    def test_provenance_and_blockers_do_not_promote_physical_c5(self):
        e=crossworkstream_evidence(self.portal,self.gamma,self.legacy)
        self.assertFalse(e['gammacc']['signed_factor_certified'])
        rows=[{'model':m,'status':'PORTAL_NUMERIC_CHECKS_PASS'} for m in ALL_MODELS]
        s=assemble(rows,e)
        self.assertFalse(s['physical_C5_certified'])
        self.assertEqual(s['models_with_missing_native_portal_directions'],['A_-4'])
        self.assertEqual(s['previous_scalar_projection_failures'],['A_-2'])
        a=next(r for r in s['models'] if r['model']=='A_-2')
        self.assertTrue(any('Earlier independent' in x for x in a['physical_blockers']))

    def test_missing_external_evidence_stays_uncertified(self):
        e=crossworkstream_evidence(None,None,None)
        self.assertEqual(e['portal']['status'],'NOT_PROVIDED')
        s=assemble([{'model':m,'status':'PORTAL_NUMERIC_CHECKS_PASS'} for m in ALL_MODELS],e)
        self.assertEqual(s['portal_numeric_consistency_passed'],16)
        self.assertFalse(s['physical_full_one_loop_certified'])

    def test_reject_outdated_portal_inventory(self):
        p=json.loads(self.portal.read_text())
        p['models'].pop()
        write_json(self.portal,p)
        with self.assertRaisesRegex(ConsistencyError,'all 16'):
            crossworkstream_evidence(self.portal,None,None)

    def test_gammacc_source_locator_cannot_claim_physical_signed_map(self):
        p=json.loads(self.gamma.read_text())
        p['native_GammaCC_to_Warsaw_prefactor_certified']=True
        write_json(self.gamma,p)
        with self.assertRaisesRegex(ConsistencyError,'changed'):
            crossworkstream_evidence(None,self.gamma,None)

    def test_unrecognised_old_scalar_format_rejected(self):
        p=json.loads(self.legacy.read_text())
        p['format']='unknown'
        write_json(self.legacy,p)
        with self.assertRaisesRegex(ConsistencyError,'format changed'):
            crossworkstream_evidence(None,None,self.legacy)


if __name__=='__main__':
    unittest.main()
