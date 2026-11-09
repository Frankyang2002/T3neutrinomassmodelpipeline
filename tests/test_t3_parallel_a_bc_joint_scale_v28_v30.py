"""v28–v30: actual numerical root, saved-state, provenance and EFT-order guards."""
from __future__ import annotations
import json
import math
import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np

from Numerical.parallel_t3.wilson import bc_joint_scalar_scale_v28 as v28
from Numerical.parallel_t3.wilson import bc_joint_scalar_scale_audit_v29 as v29
from Numerical.parallel_t3.wilson import bc_fermion_first_common_handoff_v30 as v30


def pair(z): return {'real':float(complex(z).real),'imag':float(complex(z).imag)}
def mat(z):
    z=np.asarray(z,complex)
    return {'real':z.real.tolist(),'imag':z.imag.tolist()}

class DummyBasis:
    def __init__(self,model,m1,m2,mu):
        self.model=model;self.m_s1_sq=m1;self.m_s2_sq=m2;self.mu_sq=mu
        self.A=np.array([[m1,mu],[mu.conjugate(),m2]],complex)
        self.eigenvalues_sq,self.U=np.linalg.eigh(self.A)
        self.tachyonic=bool(self.eigenvalues_sq[0]<0)
    @classmethod
    def build(cls,*args):return cls(*args)
    def info(self):
        return {'model':self.model,'masses_squared_gev2':self.eigenvalues_sq.tolist(),
                'unitary_flavour_rotation_U':mat(self.U),
                'mass_squared_flavour_matrix_gev2':mat(self.A),
                'matchete_field_redefinition_installed':False,
                'threshold_matching_performed':False}

class DummyState:
    def __init__(self,z):
        self.baseline=types.SimpleNamespace(mS1Sq=float(z[0]),mS2Sq=float(z[1]))
        self.extra={'lambdaExtra':complex(.19,.04)}


def setup_numeric(root:Path,backend:str):
    paths={name:root/f'{name}.json' for name in ('registry','manifest','beta')}
    if backend=='quadratic_only': paths['saved']=root/'saved.json'
    else: paths['mass_beta']=root/'mass_beta.json'
    for k,p in paths.items():p.write_text(json.dumps({'source':k}),encoding='utf-8')
    m1=1400**2;m2=1600**2;mu=complex(5000,1200)
    y0=np.array([m1,m2,mu.real,mu.imag]+([-7800.] if backend=='five_mass' else []),float)
    ctx={'rhs':lambda t,y:np.zeros_like(y),'y0':y0,
      'layout':types.SimpleNamespace(size=2),'registry':{},
      'unpack':lambda y,layout,registry,uv,mu_gev:DummyState(y),
      'Basis':DummyBasis,'gauge_check':lambda b,trials,seed:0.,
      'hessian':lambda m,model,m1,m2,mu:{'mass_squared_eigenvalues_gev2':np.repeat(DummyBasis(model,m1,m2,mu).eigenvalues_sq,4).tolist()},
      's34':types.SimpleNamespace(STOP=1e3),'manifest':{},'sources':paths,
      'qeft':types.SimpleNamespace(), 'mu_start':1e5}
    ctx['saved']=(None if backend=='five_mass' else {
         'final_baseline':{'mS1Sq':m1,'mS2Sq':m2},
         'mu2_final_gev2':pair(mu),
         'final_new_quartics':{'lambdaExtra':pair(.19+.04j)}})
    return ctx

class TestJointRoot(unittest.TestCase):
    def test_constant_masses_geometric_and_upper(self):
        fn=lambda t:np.array([1400**2,1600**2])
        g=v28.joint_scale_from_running(math.log(1000),math.log(1e5),fn)
        u=v28.joint_scale_from_running(math.log(1000),math.log(1e5),fn,scheme='upper')
        self.assertAlmostEqual(g['matching_scale_gev'],math.sqrt(1400*1600),places=9)
        self.assertAlmostEqual(u['matching_scale_gev'],1600,places=9)
        self.assertAlmostEqual(g['log_mu_over_mX'],-g['log_mu_over_mY'],places=11)
        self.assertTrue(g['two_scalar_eigenstates_both_retained_in_EFT1_until_matching'])
        self.assertFalse(g['intermediate_single_scalar_EFT2_used'])

    def test_running_masses_nontrivial_selfconsistent_root(self):
        ta=math.log(1500)
        fn=lambda t:np.sort(np.array([1400**2*math.exp(.10*(t-ta)),1600**2*math.exp(-.13*(t-ta))]))
        g=v28.joint_scale_from_running(math.log(1000),math.log(100000),fn)
        t=g['matching_scale_gev']
        self.assertLess(abs(math.log(t)-.25*sum(np.log(fn(math.log(t))))),1e-10)

    def test_unbracketed_tachyon_and_nonunique_rejected(self):
        for fn in [lambda t:np.array([20.,30.]),
                   lambda t:np.array([-1.,2.]),
                   lambda t:np.array([float('nan'),30.])]:
            with self.assertRaises(v28.JointScaleBlocked):
                v28.joint_scale_from_running(math.log(1000),math.log(1e5),fn)
        with self.assertRaises(v28.JointScaleBlocked):
            v28.joint_scale_from_running(10.,10.,lambda t:np.array([1.,2.]))

    def test_endpoint_replay_rejects_forgery(self):
        self.assertAlmostEqual(v28.compare_endpoint([1,2,3],[1,2,3]),0)
        with self.assertRaises(v28.JointScaleBlocked):v28.compare_endpoint([1,2,3],[1,2,4])
        with self.assertRaises(v28.JointScaleBlocked):v28.compare_endpoint([1,2],[1,2,3])

def make_prior_v27_sources(folder:Path, payload:dict)->dict:
    files={'v22':'bc_quartic_mass_readiness_conditional_v22.json',
           'v25':'bc_native_matchete_mass_eigen_concordance_v25.json',
           'v26':'bc_existing_one_scale_matchete_probe_v26.json'}
    digests={}
    for label,name in files.items():
        f=folder/name
        f.write_text(json.dumps({'label':label,'test_only':True}),encoding='utf-8')
        digests[label]=v28.sha(f)
    payload['source_sha256']=digests
    return payload

class TestEvidence(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)
        self.s33=types.ModuleType('scripts.FermionThresholdHandoff')
        self.s33.json_state=lambda value:pair(value) if isinstance(value,complex) else (
            vars(value) if hasattr(value,'__dict__') else value)
        self.scripts=types.ModuleType('scripts')
        self.scripts.__path__=[]
        self.scripts.FermionThresholdHandoff=self.s33

    def produce(self,backend='quadratic_only'):
        ctx=setup_numeric(self.root,backend)
        for model in v28.MODELS:
            tag=model.replace('-','m')
            req=(f'{tag}_eft1_full_mass_v1.json' if backend=='five_mass' else
                 f'{tag}_eft1_quadratic_running.json')
            (self.root/req).write_text('{}',encoding='utf-8')
        with patch.object(v28,'_stage33_or_saved_qeft',return_value=ctx),\
             patch.dict(sys.modules,{'scripts':self.scripts,'scripts.FermionThresholdHandoff':self.s33}):
            a=v28.run_model(self.root,'B_-1',backend=backend,stage33_dir=self.root,
                            qeft_dir=self.root,mass_dir=self.root)
            b=v28.run_model(self.root,'C_-1',backend=backend,stage33_dir=self.root,
                            qeft_dir=self.root,mass_dir=self.root)
        return a,b,ctx

    def test_quadratic_backend_complete_synthetic_BC(self):
        a,b,_=self.produce()
        self.assertEqual(a['backend'],'quadratic_only')
        self.assertIsNone(a['state_at_common_scale']['mHSq_gev2'])
        self.assertLess(a['saved_quadratic_endpoint_relative_residual'],1e-9)
        self.assertEqual(a['source_mass_parameters_at_start'],
                         {'taken_from_saved_candidate_quadratic_eft1':True})
        result=v29.audit(self.root,a,b)
        self.assertEqual(result['models_checked'],2)
        self.assertGreater(result['models']['B_-1']['common_muS_gev'],1400)
        self.assertFalse(result['certification']['physical_Wilson_evolution_authorized'])

    def test_five_mass_full_higgs_source_still_unmatched(self):
        a,b,_=self.produce('five_mass')
        self.assertEqual(a['backend'],'five_mass')
        self.assertEqual(a['state_at_common_scale']['mHSq_gev2'],-7800)
        self.assertIsNone(a['saved_quadratic_endpoint_relative_residual'])
        self.assertTrue(v29.audit(self.root,a,b)['models_checked']==2)
        self.assertFalse(a['certification']['joint_two_scalar_threshold_Wilson_matching_performed'])

    def test_source_modified_or_external_path_fails(self):
        a,_,ctx=self.produce()
        ctx['sources']['beta'].write_text('changed',encoding='utf-8')
        with self.assertRaisesRegex(v28.JointScaleBlocked,'stale source hash'):
            v29.audit_model(self.root,a)
        a['source_files']['beta']='../secret.json'
        with self.assertRaises(v28.JointScaleBlocked):v29.audit_model(self.root,a)

    def test_physics_overclaim_wrong_root_and_corrupt_matrix_rejected(self):
        a,b,_=self.produce()
        a['certification']['physical_Wilson_evolution_authorized']=True
        with self.assertRaises(v28.JointScaleBlocked):v29.audit(self.root,a,b)
        a['certification']['physical_Wilson_evolution_authorized']=False
        a['common_scalar_matching_scale']['matching_scale_gev']+=7
        with self.assertRaises(v28.JointScaleBlocked):v29.audit_model(self.root,a)
        a,b,_=self.produce()
        a['state_at_common_scale']['mass_basis']['unitary_flavour_rotation_U']['real'][0][0]*=1.1
        with self.assertRaises(v28.JointScaleBlocked):v29.audit_model(self.root,a)

    def test_v30_ledger_does_not_claim_any_eft(self):
        a,b,_=self.produce()
        w=v29.audit(self.root,a,b)
        v29path=self.root/'v29.json';v29path.write_text(json.dumps(w),encoding='utf-8')
        # Fixture follows actual previously implemented v27 format/status.
        v27data={'format':'t3-parallel-a-bc-fermion-first-scalar-threshold-itinerary-v27',
                 'status':'BCMatcheteEigenbasisLinkedToSourceEFT1_ThresholdReevaluationRequired_PhysicsBlocked',
                 'models_checked':2,'certification':{'EFT_order':v28.ORDER},
                 'models':{key:{'source_UV_fermion_first_matching_verified':False,
                                'joint_scalar_matching_completed':False,
                                'direct_one_scale_matchete_probe_not_valid_EFT1_scalar_matching':True}
                           for key in v28.MODELS}}
        make_prior_v27_sources(self.root,v27data)
        p27=self.root/'v27.json';p27.write_text(json.dumps(v27data),encoding='utf-8')
        # v30 import calls v27 source which comes from prior Workstream A v25.
        result=v30.audit(p27,v29path,root=self.root)
        self.assertEqual(result['models_checked'],2)
        self.assertFalse(result['certification']['SMEFT_matching'])
        self.assertFalse(result['certification']['physical_Wilson_evolution_authorized'])

    def test_saved_native_beta_windows_keys_and_byte_mutation(self):
        paths={k:self.root/f'{k}.txt' for k in ('native_beta','registry','manifest')}
        for path in paths.values():path.write_text('source')
        saved={'input_sha256':{str(path).replace('/','\\'):v28.sha(path)
                                for path in paths.values()}}
        v28.verify_saved_beta_sources(saved,paths)
        paths['native_beta'].write_text('tampered')
        with self.assertRaises(v28.JointScaleBlocked):
            v28.verify_saved_beta_sources(saved,paths)

    def test_saved_source_duplicate_path_alias_rejected(self):
        path=self.root/'registry.txt';path.write_text('source')
        key=str(path)
        saved={'input_sha256':{key:v28.sha(path),key.replace('/','\\'):v28.sha(path)}}
        with self.assertRaises(v28.JointScaleBlocked):
            v28.verify_saved_beta_sources(saved,{'registry':path})

    def test_reject_tampered_other_common_scheme(self):
        a,b,_=self.produce()
        a['other_joint_scale_scheme']['matching_scale_gev']*=1.01
        with self.assertRaises(v28.JointScaleBlocked):v29.audit(self.root,a,b)

    def test_reject_eft_two_step_and_wrong_operator_basis(self):
        a,b,_=self.produce()
        a['common_scalar_matching_scale']['intermediate_single_scalar_EFT2_used']=True
        with self.assertRaises(v28.JointScaleBlocked):v29.audit(self.root,a,b)
        a,b,_=self.produce()
        a['EFT_order']='UV -> scalar matching -> fermion matching'
        with self.assertRaises(v28.JointScaleBlocked):v29.audit(self.root,a,b)

    def test_higgs_mass_missing_in_full_backend_rejected(self):
        a,b,_=self.produce('five_mass')
        a['state_at_common_scale']['mHSq_gev2']=None
        with self.assertRaises((v28.JointScaleBlocked,TypeError)):v29.audit(self.root,a,b)

    def test_wrong_source_pair_and_shaped_eigenvectors_rejected(self):
        a,b,_=self.produce()
        with self.assertRaises(v28.JointScaleBlocked):v29.audit(self.root,a,a)
        a['state_at_common_scale']['mass_basis']['masses_squared_gev2']=[1,2,3]
        with self.assertRaises(v28.JointScaleBlocked):v29.audit_model(self.root,a)

    def test_v30_rejects_false_prior_joint_scalar_matching(self):
        a,b,_=self.produce()
        checked=v29.audit(self.root,a,b)
        p29=self.root/'v29.json';p29.write_text(json.dumps(checked))
        past={'format':'t3-parallel-a-bc-fermion-first-scalar-threshold-itinerary-v27',
              'status':'BCMatcheteEigenbasisLinkedToSourceEFT1_ThresholdReevaluationRequired_PhysicsBlocked',
              'models_checked':2,'certification':{'EFT_order':v28.ORDER},
              'models':{name:{'source_UV_fermion_first_matching_verified':False,
                              'joint_scalar_matching_completed':name=='B_-1',
                              'direct_one_scale_matchete_probe_not_valid_EFT1_scalar_matching':True}
                        for name in v28.MODELS}}
        make_prior_v27_sources(self.root,past)
        p27=self.root/'v27.json';p27.write_text(json.dumps(past))
        with self.assertRaises(v28.JointScaleBlocked):v30.audit(p27,p29,root=self.root)

    def test_provenance_relative_path_escape_rejected(self):
        a,_,_=self.produce()
        a['source_files']['beta']='/etc/passwd'
        with self.assertRaises(v28.JointScaleBlocked):v29.audit_model(self.root,a)
        a,_,_=self.produce()
        a['source_files']['beta']='C:/Windows/System32/config'
        with self.assertRaises(v28.JointScaleBlocked):v29.audit_model(self.root,a)

    def test_v30_upstream_v27_hash_tampering_rejected(self):
        a,b,_=self.produce();checked=v29.audit(self.root,a,b)
        p29=self.root/'v29.json';p29.write_text(json.dumps(checked))
        prior=make_prior_v27_sources(self.root,{
           'format':'t3-parallel-a-bc-fermion-first-scalar-threshold-itinerary-v27',
           'status':'BCMatcheteEigenbasisLinkedToSourceEFT1_ThresholdReevaluationRequired_PhysicsBlocked',
           'models_checked':2,'certification':{'EFT_order':v28.ORDER},
           'models':{model:{'source_UV_fermion_first_matching_verified':False,
                             'joint_scalar_matching_completed':False,
                             'direct_one_scale_matchete_probe_not_valid_EFT1_scalar_matching':True}
                    for model in v28.MODELS}})
        p27=self.root/'v27.json';p27.write_text(json.dumps(prior))
        self.assertEqual(v30.audit(p27,p29,root=self.root)['models_checked'],2)
        (self.root/'bc_existing_one_scale_matchete_probe_v26.json').write_text('{}')
        with self.assertRaises(v28.JointScaleBlocked):v30.audit(p27,p29,root=self.root)

    def test_run_cli_immutable_file(self):
        from scripts.RunT3ParallelABCJointScaleV28V30 import emit
        target=self.root/'a.json'
        self.assertEqual(emit(target,{'a':1}),'CREATED')
        self.assertEqual(emit(target,{'a':1}),'REUSED')
        with self.assertRaises(v28.JointScaleBlocked):emit(target,{'a':2})
        self.assertEqual(json.loads(target.read_text()),{'a':1})

if __name__=='__main__':unittest.main()
