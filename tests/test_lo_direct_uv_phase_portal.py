"""Numerical basis invariance and explicitly conditional coupling-map hypotheses.

Synthetic fixtures are NOT native Matchete source certificates. All 16 catalogued
charge assignments are exercised, including Majorana, Dirac and extra charges.
"""
from __future__ import annotations

import json
from pathlib import Path
import tempfile
import hashlib
import unittest
from contextlib import redirect_stdout
import io

import numpy as np

from Numerical.direct_one_loop.core import evaluate, encoded, _complex
from Numerical.lo_certification.direct_uv.certification import CLASSES, all_models, charge_bases
from Numerical.lo_certification.direct_uv.phase_portal import (
    rephase, portal_geometry, portal_candidate, convention_scan,
    inspect_phase_orbit, ConventionError,
)
from Numerical.lo_certification.direct_uv.runner import run
from test_lo_direct_uv_certification import synthetic_all16, toy_B


def restricted(short):
    d=synthetic_all16(short)
    ds1,ds2=CLASSES[short[0]][:2]
    A=_complex(d['scalar_A_gev2'],'A')
    v=246.22
    m1,m2=(900.**2,1020.**2)
    h1,h2=(.03,-.02)
    for s in range(ds1):A[s,s]=m1-h1*v*v/2
    for s in range(ds1,ds1+ds2):A[s,s]=m2-h2*v*v/2
    d['scalar_A_gev2']=encoded(A)
    bench={'numerical':{'vev_gev':v},
           'base_state':{'ordinary':{'mS1Sq':m1,'mS2Sq':m2,'lambdaH1':h1,'lambdaH2':h2}}}
    return d,bench


def portal_stub(short):
    ds1,ds2=CLASSES[short[0]][:2]
    rows=[]
    for index,d in enumerate((ds1,ds2)):
        rows.append({'scalar':f'S{index+1}','dimension':d,
              'full_tensor_identity_in_real_native_span':True,
              'does_not_certify_coupling_sign_or_normalization':True,
              'native_coupling_names':[f'lambdaH{index+1}Inv1',f'lambdaH{index+1}Inv2'],
              'native_coupling_coefficients_for_unit_identity':[0.,1.],
              'native_vev_projector_diagonals':[[0.]*d,[1.]*d],
              'max_abs_four_index_residual':0.})
    return {'model':short,'scalars':rows}


class PhaseAndPortal(unittest.TestCase):
    def test_complex_weyl_and_scalar_rephasing_including_bilinears(self):
        d=toy_B()
        B=np.zeros((4,4),complex)
        B[0,0]=30.+20j;B[3,3]=-25.+13j;B[0,3]=B[3,0]=10.-5j
        d['scalar_B_gev2']=encoded(B)
        pf=np.exp(1j*np.array([.25,.74,-.47]))
        ps=np.exp(1j*np.array([.31,.07,-.37,1.25]))
        r=rephase(d,pf,ps)
        before=evaluate(d,retain_individual_terms=False)
        after=evaluate(r,retain_individual_terms=False)
        a=_complex(before['total']['complex_symmetric_matrix_gev'],'m')
        b=_complex(after['total']['complex_symmetric_matrix_gev'],'m')
        self.assertLess(np.linalg.norm(a-b)/max(np.linalg.norm(a),1e-80),2e-6)
        # Check all three independent tensor transport laws explicitly.
        np.testing.assert_allclose(_complex(r['fermion_mass_gev'],'M'),
                   pf[:,None]*_complex(d['fermion_mass_gev'],'M')*pf[None,:])
        np.testing.assert_allclose(_complex(r['scalar_B_gev2'],'B'),
                   ps[:,None]*B*ps[None,:])
        self.assertFalse(r.get('provenance',{}).get('validated_uv_component_conventions'))

    def test_every_representation_exact_rephasing(self):
        counted=0
        for short,_ in all_models():
            with self.subTest(model=short):
                d=synthetic_all16(short)
                base=evaluate(d,retain_individual_terms=False)
                orbit=inspect_phase_orbit(d,base)
                self.assertEqual(len(orbit['full_loop_sign_rephasing_cases']),2)
                self.assertTrue(all(c['passed'] for c in orbit['full_loop_sign_rephasing_cases']))
                self.assertFalse(orbit['physical_source_vertex_normalisation_certified'])
                counted+=1
        self.assertEqual(counted,16)

    def test_native_portal_identity_reconstructed_all16(self):
        for short,_ in all_models():
            with self.subTest(model=short):
                res=portal_geometry(portal_stub(short),short[0])
                self.assertEqual(res['status'],'native_portal_geometry_verified_physical_sign_pending')
                self.assertEqual(len(res['sector_witnesses']),2)
                self.assertFalse(res['portal_dictionary_certified'])

    def test_portal_sign_and_exact_mixing_rephasing_on_all16(self):
        counted=0
        for short,_ in all_models():
            with self.subTest(model=short):
                data,b=restricted(short)
                baseline=evaluate(data,retain_individual_terms=False)
                scan=convention_scan(data,short[0],int(short[2:]),b,portal_stub(short),baseline)
                self.assertEqual(len(scan['variants']),5)
                self.assertLess(scan['relative_matrix_discrepancy_mix_flip_vs_y2_flip'],4e-5)
                self.assertEqual(scan['status'],'conditional_exact_mass_variants_not_physical_certification')
                self.assertFalse(scan['physical_certified'])
                self.assertLess(scan['variants']['S2_rephasing_with_all_vertices_transported']
                                ['mass_matrix_difference_frobenius_relative_to_baseline'],4e-5)
                self.assertGreaterEqual(scan['relative_shift_from_portal_hypothesis'],0.)
                self.assertEqual(len(scan['conditional_minus_portal_charge_and_suppression']
                                     ['per_fermion_eigenstate']),3*CLASSES[short[0]][2]*(1 if short in ('B_-1','C_-1') else 2))
                counted+=1
        self.assertEqual(counted,16)

    def test_source_like_native_doublet_portal_basis_requires_two_invariants(self):
        # Actual existing source audit has these doublet values (not an
        # invented assumption that each native portal is an identity).
        p=portal_stub('B_-1')
        coeff=[.5,np.sqrt(3.)/2]
        proj=[[0.,1.],[2/np.sqrt(3.),1/np.sqrt(3.)]]
        for entry in p['scalars']:
            entry['native_coupling_coefficients_for_unit_identity']=coeff
            entry['native_vev_projector_diagonals']=proj
        g=portal_geometry(p,'B')
        self.assertLess(max(x['maximum_neutral_projector_error']
                            for x in g['sector_witnesses']),2e-14)
        d,b=restricted('B_-1')
        alt,info=portal_candidate(d,b,p,'B',sign=-1)
        a=_complex(alt['scalar_A_gev2'],'A')
        old=_complex(d['scalar_A_gev2'],'A')
        self.assertGreater(np.linalg.norm(a-old),0)
        self.assertFalse(info['physical_cross_package_dictionary_certified'])

    def test_all16_runner_with_source_tied_native_portal_audits(self):
        # Synthetic test-only reports exercise the complete new runner path,
        # all SHA guards, all four variants, summaries, and no overwrites.
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            model_rows=[]
            for short,key in all_models():
                bpath=root/'Numerical/benchmarks/fitted_reference'/(key+'.json')
                bpath.parent.mkdir(parents=True,exist_ok=True)
                d1,d2,df=CLASSES[short[0]][:3]
                data,benchmark=restricted(short)
                benchmark.update({
                    'representation':{'d_s1':d1,'d_s2':d2,'d_f':df,'alpha':int(short[2:])},
                    'base_state':{'mu_gev':1e7,**benchmark['base_state']},
                    'scales':{'mu_fermion_threshold_gev':1e5,'mu_scalar_threshold_gev':1e3},
                })
                bpath.write_text(json.dumps(benchmark))
                hashed=hashlib.sha256(bpath.read_bytes()).hexdigest()
                data['provenance']={'benchmark_sha256':hashed,'validated_uv_component_conventions':False}
                ip=root/'output/general_observable/lo_direct_one_loop/provisional_inputs'/(short+'.json')
                ip.parent.mkdir(parents=True,exist_ok=True)
                ip.write_text(json.dumps(data))
                entry=portal_stub(short);entry['frozen_sha256']=hashed
                model_rows.append(entry)
            auditpath=root/'output/general_observable/lo_direct_one_loop/source_audits/native_portal_identity_decomposition.json'
            auditpath.parent.mkdir(parents=True,exist_ok=True)
            auditpath.write_text(json.dumps({
               'status':'all_32_geometric_identities_representable_physical_dictionary_pending',
               'models_checked':16,'scalar_portal_sectors_checked':32,
               'full_tensor_identities_representable':32,'models':model_rows,
            }))
            out=root/'output/general_observable/lo_certification_direct_uv'
            with redirect_stdout(io.StringIO()): result=run(root,out,convention_scan_enabled=True)
            self.assertEqual(result['counts']['exact_uv_provisional'],16)
            self.assertEqual(result['counts']['convention_scan_completed'],16)
            self.assertEqual(result['counts']['convention_scan_blocked'],0)
            self.assertEqual(result['counts']['physical_source_certified'],0)
            self.assertEqual(len(list((out/'convention_scan').glob('*.json'))),16)
            self.assertTrue((out/'all_16_convention_scan_summary.csv').is_file())
            with self.assertRaises(FileExistsError):run(root,out,convention_scan_enabled=True)

    def test_bad_projector_identity_fails_closed(self):
        p=portal_stub('C_-1');p['scalars'][1]['native_vev_projector_diagonals'][1][0]=2.
        with self.assertRaisesRegex(ConventionError,'does not reconstruct'):
            portal_geometry(p,'C')

    def test_bad_rephasing_modulus_fails_closed(self):
        d=toy_B()
        with self.assertRaisesRegex(ConventionError,'unit modulus'):
            rephase(d,np.array([1.,2.,1.]),np.ones(4))

    def test_source_candidate_diagonal_mismatch_rejected(self):
        d,b=restricted('B_-1')
        A=_complex(d['scalar_A_gev2'],'A');A[0,0]+=1234.
        d['scalar_A_gev2']=encoded(A)
        with self.assertRaisesRegex(ConventionError,'differs from restricted'):
            portal_candidate(d,b,portal_stub('B_-1'),'B',sign=-1)

    def test_runner_reports_missing_portal_without_blocking_uv(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            for short,key in all_models():
                bpath=root/'Numerical/benchmarks/fitted_reference'/(key+'.json')
                bpath.parent.mkdir(parents=True,exist_ok=True)
                d1,d2,df=CLASSES[short[0]][:3]
                b={'representation':{'d_s1':d1,'d_s2':d2,'d_f':df,'alpha':int(short[2:])},
                   'base_state':{'mu_gev':1e7,'ordinary':{}},
                   'scales':{'mu_fermion_threshold_gev':1e5,'mu_scalar_threshold_gev':1e3}}
                bpath.write_text(json.dumps(b))
                data=synthetic_all16(short)
                data['provenance']={'benchmark_sha256':hashlib.sha256(bpath.read_bytes()).hexdigest()}
                ip=root/'output/general_observable/lo_direct_one_loop/provisional_inputs'/(short+'.json')
                ip.parent.mkdir(parents=True,exist_ok=True);ip.write_text(json.dumps(data))
            out=root/'output/general_observable/lo_certification_direct_uv'
            with redirect_stdout(io.StringIO()): summary=run(root,out,convention_scan_enabled=True)
            self.assertEqual(summary['counts']['exact_uv_provisional'],16)
            self.assertEqual(summary['counts']['blocked'],0)
            self.assertEqual(summary['counts']['convention_scan_blocked'],16)
            self.assertEqual(summary['counts']['convention_scan_completed'],0)
            self.assertEqual(summary['counts']['physical_source_certified'],0)


if __name__=='__main__':unittest.main()
