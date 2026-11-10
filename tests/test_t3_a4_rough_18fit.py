"""Offline tests for the fixed-background 18-real-Yukawa A(-4) fit."""
from __future__ import annotations

import sys
import unittest
from types import ModuleType, SimpleNamespace
from unittest.mock import patch

import numpy as np

from Numerical.parallel_t3.background.general_batch1_api import encode_native
from Numerical.parallel_t3.complete.a4_rough_18fit import (
    SCALE, PAIRS, LinearSMEFT, UVToFLinearized, FitWorkspace,
    RoughFitError, _validate_target, fit_18_real,
)
from Numerical.parallel_t3.complete.a4_rough_benchmark import estimate_boundary


def fixture():
    a=np.array([[.2,-.4,.3],[.5,.6,.4],[-.2,.5,.7]])
    b=np.array([[.9,.1,-.3],[.6,-.2,.4],[.3,.7,.5]])
    cfg={'schema':'t3-thesis-benchmark-two-model-pilot-v1',
         'scales_GeV':{'UV':1e7,'F':1e5,'S':1e3,'low':100},
         'vev_GeV':246.22,'models':{'A_-4':{
             'y1_UV_times_1e4':a.tolist(),'y2_UV_times_1e4':b.tolist()}}}
    f={'model':'A_-4','scale_GeV':1e5,
       'all_heavy_fermions_to_integrate_out_together':True,
       'native_uv_rgbeta_couplings_at_muF':{
           'MF':encode_native(np.diag([1e5,1.1e5,1.2e5]).astype(complex)),
           'y1':encode_native(a*1.1*SCALE),
           'y2':encode_native(b*.9*SCALE),
           'lambdaT3':encode_native(.1)}}
    s={'model':'A_-4','scale_GeV':1e3,
       'both_heavy_scalar_multiplets_to_integrate_out_together':True,
       'scalar_quadratic_terms_GeV2':{'mS1Sq':1.6e6,'mS2Sq':1.4e6},
       'native_eft1_rgbeta_couplings_at_muS':{
           'gY':.35,'g2':.60,'g3':.90,'lambdaH':.25,
           'ye':encode_native(np.diag([1e-5,5e-4,1e-2]).astype(complex)),
           'yu':encode_native(np.diag([1e-5,.003,.85]).astype(complex)),
           'yd':encode_native(np.diag([1e-5,.0004,.02]).astype(complex))}}
    target={'ordering':'NO',
            'central_values':dict(zip(('delta_m21_sq_ev2','delta_m3l_sq_ev2',
                                        'sin2_theta12','sin2_theta13','sin2_theta23'),
                                       (7.537e-5,.002511,.3088,.02248,.47))),
            'one_sigma_errors':dict(zip(('delta_m21_sq_ev2','delta_m3l_sq_ev2',
                                         'sin2_theta12','sin2_theta13','sin2_theta23'),
                                        (9.7e-7,2.05e-5,.00665,.00057,.0155)))}
    return f,s,cfg,target


class TestA4Rough18Fit(unittest.TestCase):
    def test_uv_seed_maps_to_saved_f_endpoint(self):
        f,s,c,t=fixture()
        transfer=UVToFLinearized.from_sources(f,c)
        y1,y2=transfer.f_yukawas(transfer.original_uv_scaled)
        np.testing.assert_allclose(y1,transfer.saved_f_y1,atol=0,rtol=1e-14)
        np.testing.assert_allclose(y2,transfer.saved_f_y2,atol=0,rtol=1e-14)

    def test_uniform_scaling_gives_quadratic_C5(self):
        f,s,c,t=fixture()
        uv=UVToFLinearized.from_sources(f,c)
        old=estimate_boundary(uv.modified_f_input(f,uv.original_uv_scaled),s).c5_at_ms
        fac=.74705187846
        new=estimate_boundary(uv.modified_f_input(f,fac*uv.original_uv_scaled),s).c5_at_ms
        np.testing.assert_allclose(new,fac**2*old,rtol=1e-12,atol=1e-28)

    def test_18_real_parameters_only(self):
        f,s,c,t=fixture()
        uv=UVToFLinearized.from_sources(f,c)
        for x in (np.ones(17),np.ones(19),np.full(18,float('nan'))):
            with self.assertRaises(RoughFitError):uv.f_yukawas(x)

    def test_zero_uv_seed_rejected(self):
        f,s,c,t=fixture()
        c['models']['A_-4']['y1_UV_times_1e4'][0][0]=0.
        with self.assertRaises(RoughFitError):UVToFLinearized.from_sources(f,c)

    def test_bad_target_rejected(self):
        f,s,c,t=fixture()
        t['one_sigma_errors']['sin2_theta13']=0
        with self.assertRaises(RoughFitError):_validate_target(t)

    def test_symmetric_SMEFT_basis_reconstruction(self):
        k=np.array([[2.,2j,3.],[2j,4.,5j],[3.,5j,-6.]],complex)*1e-15
        images=[]
        for p,q in PAIRS:
            b=np.zeros((3,3),complex);b[p,q]=b[q,p]=1.
            images.append(.9*b)
        tran=LinearSMEFT(np.asarray(images),np.eye(3,dtype=complex),100.,6)
        np.testing.assert_allclose(tran.low(k),.9*k,rtol=1e-13,atol=1e-29)

    def test_smeft_map_uses_repository_solver_six_times(self):
        f,s,c,t=fixture()
        calls=[]
        def evolve(initial,mu0,mu1,**kwargs):
            calls.append((mu0,mu1,kwargs['atol']))
            return SimpleNamespace(K=.9*initial.K,ye=initial.ye,solver_success=True)
        m=ModuleType('Numerical.running.SMWeinbergEvolution')
        m.SMWeinbergInitialConditions=lambda **kw:SimpleNamespace(**kw)
        m.evolve_sm_weinberg=evolve
        with patch.dict(sys.modules,{'Numerical.running.SMWeinbergEvolution':m}):
            obj=LinearSMEFT.from_scalar_input(s)
        self.assertEqual(len(calls),6)
        self.assertTrue(all(x==(1000.,100.,1e-24) for x in calls))
        self.assertEqual(obj.images.shape,(6,3,3))

    def test_rough_fit_optimizes_18_real_variables(self):
        f,s,c,t=fixture()
        uv=UVToFLinearized.from_sources(f,c)
        base=uv.original_uv_scaled
        center,sigma=_validate_target(t)
        # Synthetic analytic five-observable function, no physics claims.
        # Its atmospheric splitting obeys the exact quartic scaling.
        def predict(x):
            ratio=np.asarray(x)/base
            obs={'mixing_defined':True,
                 'dm21_sq_eV2':center[0]*ratio[1]**4,
                 'dm31_sq_eV2':center[1]*4*ratio[0]**4,
                 'sin2_theta12':center[2]+.01*(x[2]-base[2]),
                 'sin2_theta13':center[3]+.001*(x[3]-base[3]),
                 'sin2_theta23':center[4]+.01*(x[4]-base[4]),
                 'masses_eV':[.005,.01,.05]}
            arr={'C5_low':np.eye(3)*1e-15,'C5_MS':np.eye(3)*1e-15,
                 'mnu_low':np.eye(3)*-1e-11}
            return obs,arr
        ws=FitWorkspace(f,s,c,t,uv,
                        LinearSMEFT(np.zeros((6,3,3)),np.eye(3),100.,6),
                        center,sigma,1.)
        ws.predict=predict
        ws.boundary=lambda x: SimpleNamespace(c5_at_ms=np.eye(3)*1e-15)
        fake_summary={'observables':None}
        def fake_full(*args,**kwargs):
            # Test direct solver consistency check only.
            fake_summary['observables']=predict(ws._last)[0]
            return fake_summary,{'C5_low_approx_GeVinv':np.eye(3)*1e-15,
                                 'mnu_low_approx_GeV':np.eye(3)*-1e-11}
        old=ws.predict
        def track(x):
            ws._last=np.asarray(x)
            return old(x)
        ws.predict=track
        path='Numerical.parallel_t3.complete.a4_rough_18fit.run_smeft_observables'
        with patch(path,side_effect=fake_full):
            result,arr=fit_18_real(ws,max_nfev=35)
        self.assertEqual(result['parameters'],18)
        self.assertLess(result['final_fit_chi2'],2.)
        self.assertEqual(np.asarray(result['fitted_y1_UV_times_1e4']).shape,(3,3))
        self.assertFalse(result['full_physical_C5_verified'])

if __name__=='__main__': unittest.main()
