"""B1 source-spurion algebra, adversarial controls, provenance and synthetic audit."""
from __future__ import annotations
import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
import sympy as sp

from Numerical.lo_certification.physical_normalisation.b1_phase_transport import (
    PhaseAuditError, exact_invariant_witness, exponent_matrix, sign_orbit_witness,
    complex_matrix, complex_scalar, source_flavour_numerator,
    change_basis, deterministic_unitary, test_flavour_covariance, relative,
)
from Numerical.lo_certification.physical_normalisation.b1_matchete_weyl import model_specs
from scripts.RunPhysicalB1PhaseTransport import (
    benchmark_filename, run_benchmarks, source_contract, build,
)


def sh(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


class PhaseTransportTests(unittest.TestCase):
    def setUp(self):
        self.rng=np.random.default_rng(443)
        self.a=self.rng.normal(size=(3,3))+1j*self.rng.normal(size=(3,3))
        self.b=self.rng.normal(size=(3,3))+1j*self.rng.normal(size=(3,3))
        self.m=self.rng.normal(size=(3,3))+1j*self.rng.normal(size=(3,3))
        self.lam=0.7+1.2j
        self.uc=deterministic_unitary(self.rng)
        self.ux=deterministic_unitary(self.rng)
        self.un=deterministic_unitary(self.rng)

    def test_exact_dirac_unique_spurion(self):
        d=exact_invariant_witness(False)
        self.assertTrue(d['internal_invariant_witness_exact'])
        self.assertEqual(d['unique_internal_invariant_exponents'],[1,1,-1,-1])
        self.assertEqual(exponent_matrix(False).shape,(4,5))

    def test_exact_majorana_unique_spurion(self):
        d=exact_invariant_witness(True)
        self.assertTrue(d['external_nu_pair_covariance_exact'])
        self.assertEqual(exponent_matrix(True).shape,(4,4))
        self.assertEqual(exponent_matrix(True)[2,1],2)

    def test_sign_orbit_dirac_majorana(self):
        d=sign_orbit_witness()['cases']
        for item in d:
            self.assertEqual(item['allowed_internal_rephase_sign_patterns'],8)
            self.assertFalse(item['single_yukawa_flip_in_internal_orbit'])
            self.assertTrue(item['both_yukawas_flip_in_internal_orbit'])
            self.assertTrue(item['Yukawa_plus_MF_flip_in_internal_orbit'])
            self.assertTrue(item['Yukawa_plus_lambda_flip_in_internal_orbit'])

    def test_dirac_independent_heavy_u3_and_nu_u3(self):
        t=change_basis(self.a,self.b,self.m,self.lam,self.uc,self.ux,self.un,.33,-.92,False)
        p=source_flavour_numerator(*t)
        expected=self.un.T @ source_flavour_numerator(self.a,self.b,self.m,self.lam) @ self.un
        self.assertLess(relative(p,expected),2e-13)

    def test_majorana_same_u3_rotates_symmetric_mass(self):
        mass=self.m+self.m.T
        t=change_basis(self.a,self.b,mass,self.lam,self.uc,None,self.un,.13,.75,True)
        self.assertLess(relative(t[2],t[2].T),1e-14)
        want=self.un.T @ source_flavour_numerator(self.a,self.b,mass,self.lam) @ self.un
        self.assertLess(relative(source_flavour_numerator(*t),want),2e-13)

    def test_majorana_cannot_use_independent_xi(self):
        with self.assertRaises(PhaseAuditError):
            change_basis(self.a,self.b,self.m+self.m.T,self.lam,self.uc,self.ux,self.un,0,0,True)

    def test_majorana_nonsymmetric_mass_rejected(self):
        with self.assertRaises(PhaseAuditError):
            change_basis(self.a,self.b,self.m,self.lam,self.uc,None,self.un,0,0,True)

    def test_independent_covariance_both_fermions(self):
        for majorana in (False,True):
            m=self.m+self.m.T if majorana else self.m
            out=test_flavour_covariance(self.a,self.b,m,self.lam,majorana,seed=500)
            self.assertLess(out['U3_flavour_covariance_residual'],2e-13)
            self.assertLess(out['symmetric_nu_pair_covariance_residual'],2e-13)

    def test_wrong_nonconjugated_mass_not_covariant(self):
        def wrong(a,b,m,lam):return a@m@b.T*complex(lam).conjugate()
        t=change_basis(self.a,self.b,self.m,self.lam,self.uc,self.ux,self.un,0,0,False)
        self.assertGreater(relative(wrong(*t),self.un.T@wrong(self.a,self.b,self.m,self.lam)@self.un),.1)

    def test_wrong_transpose_only_mass_not_covariant(self):
        def wrong(a,b,m,lam):return a@m.T@b.T*complex(lam).conjugate()
        t=change_basis(self.a,self.b,self.m,self.lam,self.uc,self.ux,self.un,0,0,False)
        self.assertGreater(relative(wrong(*t),self.un.T@wrong(self.a,self.b,self.m,self.lam)@self.un),.1)

    def test_wrong_mf_conjugate_without_transpose_not_covariant(self):
        def wrong(a,b,m,lam):return a@m.conj()@b.T*complex(lam).conjugate()
        t=change_basis(self.a,self.b,self.m,self.lam,self.uc,self.ux,self.un,0,0,False)
        self.assertGreater(relative(wrong(*t),self.un.T@wrong(self.a,self.b,self.m,self.lam)@self.un),.1)

    def test_wrong_lambda_not_conjugated(self):
        def wrong(a,b,m,lam):return a@m.conj().T@b.T*lam
        t=change_basis(self.a,self.b,self.m,self.lam,np.eye(3),np.eye(3),np.eye(3),.7,.2,False)
        self.assertGreater(relative(wrong(*t),wrong(self.a,self.b,self.m,self.lam)),.1)

    def test_only_one_yukawa_sign_flips_invariant(self):
        original=source_flavour_numerator(self.a,self.b,self.m,self.lam)
        changed=source_flavour_numerator(-self.a,self.b,self.m,self.lam)
        self.assertLess(relative(changed,-original),1e-14)
        self.assertGreater(relative(changed,original),1.99)

    def test_paired_yukawa_flip_harmless(self):
        a=source_flavour_numerator(self.a,self.b,self.m,self.lam)
        b=source_flavour_numerator(-self.a,-self.b,self.m,self.lam)
        self.assertLess(relative(a,b),1e-14)

    def test_paired_mass_and_yukawa_flip_harmless(self):
        a=source_flavour_numerator(self.a,self.b,self.m,self.lam)
        b=source_flavour_numerator(-self.a,self.b,-self.m,self.lam)
        self.assertLess(relative(a,b),1e-14)

    def test_paired_scalar_mix_and_yukawa_flip_harmless(self):
        a=source_flavour_numerator(self.a,self.b,self.m,self.lam)
        b=source_flavour_numerator(-self.a,self.b,self.m,-self.lam)
        self.assertLess(relative(a,b),1e-14)

    def test_majorana_rephase_by_i_is_required_for_some_sign_patterns(self):
        # chi_old=i chi_new, S1_old=i S1_new, S2_old=i S2_new.
        # Then (y1,y2,MF,lambda) -> (-y1,y2,-MF,lambda).
        m=self.m+self.m.T
        t=change_basis(self.a,self.b,m,self.lam,1j*np.eye(3),None,
                       np.eye(3),np.pi/2,np.pi/2,True)
        self.assertLess(relative(t[0],-self.a),1e-14)
        self.assertLess(relative(t[1],self.b),1e-14)
        self.assertLess(relative(t[2],-m),1e-14)
        self.assertLess(abs(t[3]-self.lam),1e-14)

    def test_external_neutrino_rephase_covariant_not_invariant(self):
        nu=np.exp(1j*.4)*np.eye(3)
        t=change_basis(self.a,self.b,self.m,self.lam,np.eye(3),np.eye(3),nu,0,0,False)
        old=source_flavour_numerator(self.a,self.b,self.m,self.lam)
        self.assertLess(relative(source_flavour_numerator(*t),np.exp(.8j)*old),1e-13)
        self.assertGreater(relative(source_flavour_numerator(*t),old),.1)

    def test_reject_bad_matrices_nonunitaries_nonfinite(self):
        with self.assertRaises(PhaseAuditError):complex_matrix([[0]],'y1')
        with self.assertRaises(PhaseAuditError):complex_scalar(complex(float('nan')),'lambdaT3')
        with self.assertRaises(PhaseAuditError):
            change_basis(self.a,self.b,self.m,self.lam,2*np.eye(3),np.eye(3),self.un,0,0,False)
        with self.assertRaises(PhaseAuditError):
            change_basis(self.a,self.b,self.m,self.lam,self.uc,np.eye(3),self.un,float('nan'),0,False)


class PhaseTransportPipelineTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)
        self.action=self.root/'action.json';self.prior=self.root/'prior.json'
        self.cc=self.root/'cc.json';self.native=self.root/'native'
        self.native.mkdir()
        (self.root/'Numerical/benchmarks/fitted_reference').mkdir(parents=True)
        (self.root/'output/general_observable/lo_direct_one_loop/provisional_inputs').mkdir(parents=True)
        self.natives={}
        for tag in ('A_p0','B_m1','C_m1'):
            p=self.native/(tag+'.json');p.write_text(json.dumps({'model':tag}),encoding='utf-8')
            self.natives[tag]={'sha256':sh(p),'couplings_and_chiral_projectors':{'mass':tag}}
        rows=[]
        for name,dims in model_specs():
            b=self.root/'Numerical/benchmarks/fitted_reference'/benchmark_filename(dims)
            j=self.root/'output/general_observable/lo_direct_one_loop/provisional_inputs'/(name+'.json')
            matrix=lambda n:{'real':n,'imag':np.zeros((3,3)).tolist()}
            p=np.eye(3).tolist()
            payload={'representation':dict(zip(('d_s1','d_s2','d_f','alpha'),dims)),
                     'base_state':{'ordinary':{'y1':matrix(p),'y2':matrix(p),
                                    'MF':matrix((np.eye(3)*100).tolist()),
                                    'lambdaT3':{'real':.1,'imag':.03}}}}
            b.write_text(json.dumps(payload),encoding='utf-8')
            j.write_text(json.dumps({'model':name}),encoding='utf-8')
            rows.append({'model':name,'benchmark_sha256':sh(b),'frozen_input_sha256':sh(j),
                         'physically_verified':False,'residuals':{'mass':0,'Y_after_global_minus':0,'Z_after_global_minus':0}})
        previous={'status':'All16ConditionalCandidateWeylTensorsReconstructed_PhysicalDictionaryPending',
                  'count_conditional_comparison_passed':16,
                  'physical_source_certified':False,'native_action_observed':True,
                  'models':rows}
        self.prior.write_text(json.dumps(previous),encoding='utf-8')
        self.cc.write_text(json.dumps({'format':'t3-b1-native-cc-abstract-v2',
                  'status':'NativeAbstractCCAlgebraVerified_NumericOrientationPending',
                  'model':'A_p0','checks':{str(i):True for i in range(6)},
                  'chosen_numeric_gamma_matrix_test_performed_here':False,
                  'installed_numeric_gamma_cc_sign_measured':False,
                  'cross_package_spinor_phase_certified':False,
                  'physical_source_certified':False,
                  'installed_symbol_manifest':{'symbol_full_name':'Matchete`GammaCC'}}),encoding='utf-8')
        from scripts.RunPhysicalB1NativeCoefficientAudit import check_all16,check_cc
        c=check_all16(previous); cc=check_cc(json.loads(self.cc.read_text()))
        self.guard={'source.wl':{'git_blob_sha1':'123','github_main_blob_match':True}}
        act={'format':'t3-b1-native-source-factor-audit-v2',
             'status':'ThreeNativeActionsCoefficientAndChiralityVerified_AbsoluteCCAndB2Pending',
             'native_actions_structurally_verified_count':3,
             'physical_source_certified':False,'safe_to_integrate_into_production':False,
             'previous_conditional_audit_sha256':sh(self.prior),
             'previous_conditional_16':c,
             'guarded_source_blobs':{'source.wl':{'git_blob_sha1':'123'}},
             'native_action_witnesses':self.natives,
             'native_cc_abstract_probe':cc,
             'chosen_basis_algebra':{'absolute_native_GammaCC_phase_measured':False}}
        self.action.write_text(json.dumps(act),encoding='utf-8')
        self.patch_guard=patch('scripts.RunPhysicalB1PhaseTransport.guard_source',return_value=self.guard)
        self.patch_parse=patch('scripts.RunPhysicalB1PhaseTransport.verify_record',
                               side_effect=lambda rec,model:{'couplings_and_chiral_projectors':{'mass':model}})
        self.patch_guard.start();self.addCleanup(self.patch_guard.stop)
        self.patch_parse.start();self.addCleanup(self.patch_parse.stop)

    def launch(self):
        return build(self.root,self.action,self.prior,self.cc,self.native)

    def test_end_to_end_synthetic_all16(self):
        r=self.launch()
        self.assertEqual(r['source_flavour_benchmarks_checked_count'],16)
        self.assertLess(r['maximum_numerator_covariance_residual'],1e-12)
        self.assertFalse(r['physical_source_certified'])
        self.assertFalse(r['safe_to_integrate_into_production'])

    def test_changed_action_audit_checksum_rejected(self):
        obj=json.loads(self.action.read_text());obj['previous_conditional_audit_sha256']='0'*64
        self.action.write_text(json.dumps(obj))
        with self.assertRaises(PhaseAuditError):self.launch()

    def test_changed_native_file_sha_rejected(self):
        (self.native/'B_m1.json').write_text(json.dumps({'model':'B_m1','changed':True}))
        with self.assertRaises(PhaseAuditError):self.launch()

    def test_changed_cc_manifest_rejected(self):
        obj=json.loads(self.cc.read_text());obj['installed_symbol_manifest']['symbol_full_name']='Wrong`GammaCC'
        self.cc.write_text(json.dumps(obj))
        with self.assertRaises(PhaseAuditError):self.launch()

    def test_false_physical_certificate_rejected(self):
        obj=json.loads(self.action.read_text());obj['physical_source_certified']=True
        self.action.write_text(json.dumps(obj))
        with self.assertRaises(PhaseAuditError):self.launch()

    def test_source_guard_drift_rejected(self):
        with patch('scripts.RunPhysicalB1PhaseTransport.guard_source',return_value={'source.wl':{'git_blob_sha1':'bad','github_main_blob_match':True}}):
            with self.assertRaises(PhaseAuditError):self.launch()

    def test_benchmark_sha_mismatch_rejected(self):
        p=next((self.root/'Numerical/benchmarks/fitted_reference').glob('*.json'))
        p.write_text(p.read_text()+' ',encoding='utf-8')
        with self.assertRaises(PhaseAuditError):self.launch()

    def test_frozen_sha_mismatch_rejected(self):
        p=next((self.root/'output/general_observable/lo_direct_one_loop/provisional_inputs').glob('*.json'))
        p.write_text(p.read_text()+' ',encoding='utf-8')
        with self.assertRaises(PhaseAuditError):self.launch()

    def test_b2_ward_report_is_conditional_only(self):
        b2=self.root/'b2.json';b2.write_text(json.dumps({'physical_source_certified':False,
                    'physical_crosspackage_complete':False,
                    'status':'LITERAL_RGBETA_DOUBLE_DELTA_WARD_MISMATCH_CONDITIONAL'}))
        r=build(self.root,self.action,self.prior,self.cc,self.native,b2)
        self.assertTrue(r['B2_evidence']['supplied'])
        self.assertFalse(r['physical_crosspackage_SU2_dictionary_verified'])

    def test_false_b2_certification_rejected(self):
        b2=self.root/'b2.json';b2.write_text(json.dumps({'physical_source_certified':True,
                    'physical_crosspackage_complete':True}))
        with self.assertRaises(PhaseAuditError):
            build(self.root,self.action,self.prior,self.cc,self.native,b2)

    def test_majorana_mass_not_symmetric_rejected(self):
        p=self.root/'Numerical/benchmarks/fitted_reference'/benchmark_filename((2,2,1,-1))
        record=json.loads(p.read_text())
        record['base_state']['ordinary']['MF']['real'][0][1]=4
        p.write_text(json.dumps(record))
        previous=json.loads(self.prior.read_text())
        for row in previous['models']:
            if row['model']=='B_-1':row['benchmark_sha256']=sh(p)
        self.prior.write_text(json.dumps(previous))
        action=json.loads(self.action.read_text());action['previous_conditional_audit_sha256']=sh(self.prior)
        from scripts.RunPhysicalB1NativeCoefficientAudit import check_all16
        action['previous_conditional_16']=check_all16(previous)
        self.action.write_text(json.dumps(action))
        with self.assertRaises(PhaseAuditError):self.launch()


if __name__=='__main__':unittest.main()
