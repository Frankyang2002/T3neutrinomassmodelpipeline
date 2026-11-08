"""A(-4) Batch-2 post-export verification; does not require archived outputs."""
from __future__ import annotations

from hashlib import sha256
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch
import json
import unittest

import numpy as np

from Numerical.running import Aminus4EFT1Audit as audit
from Numerical.running import Aminus4EFT1Couplings as eft


def point(scale:float) -> dict:
    k=np.log(100000./scale)
    u=np.diag([1e-5, 1e-3, .6 + .003*k]).astype(complex)
    d=np.diag([2e-5, 5e-4, .03]).astype(complex)
    e=np.diag([3e-6, 6e-4, .01]).astype(complex)
    return {'gY':.33,'g2':.65,'g3':.89,
            'yu':u,'yd':d,'ye':e,'lambdaH':.14,'lambdaS1':.12,
            'lambdaS2':.09,'lambdaS2Adj':.02,'lambdaH1':.03,
            'lambdaH2':.04,'lambdaH2Adj':.01,'lambda12':.02,
            'lambdaT3':.08+.005j*k,'mS1Sq':1.1e6,'mS2Sq':1.4e6}


def encode(p:dict) -> dict:
    d={}
    for k,v in p.items():
        if isinstance(v,np.ndarray):
            d[k]={'real':v.real.tolist(),'imag':v.imag.tolist()}
        elif isinstance(v,complex):d[k]={'real':v.real,'imag':v.imag}
        else:d[k]=v
    return d


def beta_payload() -> dict:
    return {'status':'Success','metadata':{
        'SharedScalar':False,'dS1':1,'dS2':3,'dF':2,'alpha':-4,
        'IntegratedField':'F','ActiveBSMFields':['S1','S2'],
        'ReportBetaConvention':eft.BETA_CONVENTION},
        'report_betas':{n:'0' for n in eft.BETA_NAMES}}


def make_fixture(folder:Path):
    scales=np.array([1e5,1e4,1e3])
    points=[point(float(mu)) for mu in scales]
    tensors=[eft.tensor_point(float(mu),p) for mu,p in zip(scales,points)]
    arr={'mu_gev':scales,'Y':np.stack([p.Y for p in tensors]),
         'lambda4':np.stack([p.lambda4 for p in tensors])}
    for kind,sec in tensors[0].gauge.items():
        for rep in ('fermion','scalar'):
            arr[f'generator_{kind}_{rep}']=np.stack(sec[rep])
    np.savez_compressed(folder/audit.F_OUTPUT,**arr)
    bp=folder/'native_beta.json'
    bp.write_text(json.dumps(beta_payload()),encoding='utf-8')
    basis=tensors[0].to_machine()
    basis['sampled_scales_gev']=scales.tolist()
    basis['running_couplings']=[encode(p) for p in points]
    basis['rge_payload_sha256']=sha256(bp.read_bytes()).hexdigest()
    (folder/audit.B_OUTPUT).write_text(json.dumps(basis),encoding='utf-8')
    return bp,basis,arr


class PostExportTests(unittest.TestCase):
    def test_real_potential_beta_self_reconstruction(self):
        x=point(1e5)
        for key in eft.SCALAR_COUPLINGS:x[key]=0j if key=='lambdaT3' else 0.
        x['lambdaH']=.2
        coeff=audit._scalar_coordinate_columns()
        self.assertEqual(coeff.shape,(12**4,10))
        self.assertEqual(np.linalg.matrix_rank(coeff,tol=1e-10),10)
        lam=eft.quartic_tensor(x)
        result=audit.scalar_beta_tensor(lam).ravel()
        pred=np.linalg.lstsq(coeff,result,rcond=None)[0]
        self.assertAlmostEqual(pred[0],12*.2**2,places=11)
        np.testing.assert_allclose(pred[1:],0,rtol=0,atol=3e-13)
        self.assertLess(np.linalg.norm(coeff@pred-result),1e-11)

    def test_isolated_triplet_adjoint_portal_reference_coefficients(self):
        p=point(1e5)
        for k in eft.SCALAR_COUPLINGS:p[k]=0j if k=='lambdaT3' else 0.
        p['lambdaH2Adj']=1.
        full=audit.scalar_beta_tensor(eft.quartic_tensor(p)).ravel()
        C=audit._scalar_coordinate_columns()
        prediction=np.linalg.lstsq(C,full,rcond=None)[0]
        q=dict(zip(audit.SCALAR_NAMES,prediction))
        self.assertAlmostEqual(q['lambdaH'],1.)
        self.assertAlmostEqual(q['lambdaH2'],2.)
        self.assertAlmostEqual(q['lambdaS2'],0.)
        self.assertAlmostEqual(q['lambdaS2Adj'],1.)

    def test_scalar_beta_comparison_uses_native_expression_and_reports_mismatch(self):
        p=point(1e5)
        for k in eft.SCALAR_COUPLINGS:p[k]=0j if k=='lambdaT3' else 0.
        p['lambdaH']=.3
        src=beta_payload()
        src['report_betas']['lambdaH']='12*lambdaH^2'
        def evaluator(expr,env):
            return 12*env['lambdaH']**2 if expr=='12*lambdaH^2' else 0.
        result=audit.scalar_beta_comparison(p,src,evaluator=evaluator)
        self.assertTrue(result['reference_and_native_scalar_loop_match'])
        self.assertAlmostEqual(result['coupling_beta_comparison']['lambdaH']
                               ['reference_scalar_loop'],12*.3**2)
        src['report_betas']['lambdaH']='3*lambdaH^2'
        def wrong(expr,env):return 3*env['lambdaH']**2 if expr=='3*lambdaH^2' else 0.
        incorrect=audit.scalar_beta_comparison(p,src,evaluator=wrong)
        self.assertIn('lambdaH',incorrect['mismatched_coupling_betas'])
        self.assertFalse(incorrect['reference_and_native_scalar_loop_match'])

    def test_three_scale_saved_tensor_audit(self):
        with TemporaryDirectory() as tmp:
            folder=Path(tmp)
            src,_,_=make_fixture(folder)
            dummy={'reference_and_native_scalar_loop_match':False,
                   'mismatched_coupling_betas':['lambdaH2Adj']}
            with patch.object(audit,'scalar_beta_comparison',return_value=dummy), \
                 patch.object(audit,'scalar_unit_channel_audits',return_value=[]):
                report=audit.audit_export(folder,src,check_sources=False,verify_ode=False)
            self.assertEqual(report['sample_count'],3)
            self.assertTrue(report['tensor_reconstruction_and_gauge_ward_passed'])
            self.assertFalse(report['all_native_reference_scalar_beta_coefficients_agree'])
            self.assertFalse(report['ready_for_physically_certified_master_Wilson_RGE'])
            for row in report['samples']:
                self.assertTrue(row['passed'])
                for group in ('U1','SU2','SU3'):
                    self.assertLess(row['gauge_ward'][group]
                                   ['quartic_max_absolute'],1e-11)

    def test_corrupted_saved_yukawa_rejected(self):
        with TemporaryDirectory() as tmp:
            folder=Path(tmp);src,_,arr=make_fixture(folder)
            arr['Y'][1,0,0,0]+=.2
            np.savez_compressed(folder/audit.F_OUTPUT,**arr)
            with self.assertRaisesRegex(audit.AuditError,'reconstruction or Ward'):
                with patch.object(audit,'scalar_beta_comparison',return_value={}):
                    audit.audit_export(folder,src,check_sources=False,verify_ode=False)

    def test_wrong_beta_hash_rejected(self):
        with TemporaryDirectory() as tmp:
            folder=Path(tmp);src,_,_=make_fixture(folder)
            obj=json.loads((folder/audit.B_OUTPUT).read_text())
            obj['rge_payload_sha256']='0'*64
            (folder/audit.B_OUTPUT).write_text(json.dumps(obj))
            with self.assertRaisesRegex(audit.AuditError,'SHA-256'):
                audit.audit_export(folder,src,check_sources=False,verify_ode=False)


if __name__=='__main__':unittest.main()
