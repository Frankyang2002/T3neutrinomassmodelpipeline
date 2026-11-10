"""Fast algebraic tests for the approximate A(-4) numeric handoff."""
from __future__ import annotations

import unittest
from unittest.mock import patch
from types import SimpleNamespace
import numpy as np

from Numerical.parallel_t3.background.general_batch1_api import encode_native
from Numerical.parallel_t3.complete.a4_rough_benchmark import (
    RoughBenchmarkError, estimate_boundary, run_smeft_observables,
)


def inputs():
    m = np.diag([1e5,1.1e5,1.2e5]).astype(complex)
    y1 = np.array([[1,2,3],[4,5,6],[2,3,4]],complex)*1e-4
    y2 = np.array([[-2,1,3],[3,5,-4],[4,-5,6]],complex)*1e-4
    f={'model':'A_-4','scale_GeV':1e5,
       'all_heavy_fermions_to_integrate_out_together':True,
       'native_uv_rgbeta_couplings_at_muF':{
           'MF':encode_native(m),'y1':encode_native(y1),'y2':encode_native(y2),
           'lambdaT3':encode_native(.1+0j)}}
    s={'model':'A_-4','scale_GeV':1e3,
       'both_heavy_scalar_multiplets_to_integrate_out_together':True,
       'scalar_quadratic_terms_GeV2':{'mS1Sq':1.6e6,'mS2Sq':1.4e6},
       'native_eft1_rgbeta_couplings_at_muS':{
          'gY':.35,'g2':.6,'g3':.9,'lambdaH':.25,
          'ye':encode_native(np.diag([1e-5,.0005,.01]).astype(complex)),
          'yu':encode_native(np.diag([1e-5,.004,.8]).astype(complex)),
          'yd':encode_native(np.diag([1e-5,.0004,.02]).astype(complex))}}
    return f,s


class TestA4RoughBenchmark(unittest.TestCase):
    def test_finite_symmetric_3x3(self):
        f,s=inputs(); result=estimate_boundary(f,s)
        self.assertEqual(result.c5_at_ms.shape,(3,3))
        np.testing.assert_allclose(result.c5_at_ms,result.c5_at_ms.T,atol=1e-30)
        np.testing.assert_allclose(result.c5_at_ms,result.ordered+result.ordered.T)

    def test_operator_factor_is_linear(self):
        f,s=inputs(); a=estimate_boundary(f,s,operator_factor=1)
        b=estimate_boundary(f,s,operator_factor=-2)
        np.testing.assert_allclose(b.c5_at_ms,-2*a.c5_at_ms,rtol=1e-12,atol=1e-28)

    def test_rejects_zero_or_nonfinite_factor(self):
        f,s=inputs()
        for value in (0, float('nan'),float('inf'),True):
            with self.subTest(value=value), self.assertRaises(RoughBenchmarkError):
                estimate_boundary(f,s,operator_factor=value)

    def test_flavour_exchange_stays_symmetric(self):
        f,s=inputs(); v=f['native_uv_rgbeta_couplings_at_muF']
        v['y1'],v['y2']=v['y2'],v['y1']
        b=estimate_boundary(f,s)
        self.assertTrue(np.allclose(b.c5_at_ms,b.c5_at_ms.T,atol=1e-28))

    def test_singular_mass_rejected(self):
        f,s=inputs(); f['native_uv_rgbeta_couplings_at_muF']['MF']=encode_native(np.zeros((3,3),complex))
        with self.assertRaises(ValueError): estimate_boundary(f,s)

    def test_wrong_threshold_rejected(self):
        f,s=inputs();s['both_heavy_scalar_multiplets_to_integrate_out_together']=False
        with self.assertRaises(ValueError): estimate_boundary(f,s)

    def test_full_pipeline_adapts_repo_solvers_without_fake_certification(self):
        f,s=inputs(); boundary=estimate_boundary(f,s)
        # Test only the handoff wiring here. On the user's workstation the
        # actual repository's SMEFT integration/Takagi implementations run.
        def fake_evolve(initial,mu0,mu1,**kwargs):
            self.assertEqual(mu0,1e3);self.assertEqual(mu1,100.)
            self.assertLess(kwargs["atol"],1e-20)
            return SimpleNamespace(solver_success=True,K=initial.K,
                                   ye=initial.ye,nfev=15)
        def fake_neutrino(m,ye):
            self.assertEqual(m.shape,(3,3))
            return ({'mixing_defined':True,'masses_eV':[.01,.012,.05],
                     'dm21_sq_eV2':4.4e-5,'dm31_sq_eV2':.0024,
                     'sin2_theta12':.31,'sin2_theta13':.022,
                     'sin2_theta23':.46}, {})
        def fake_score(obs,target):
            return {'chi2_gaussian_uncorrelated':123.}
        import sys
        from types import ModuleType
        from unittest.mock import patch
        sm = ModuleType('Numerical.running.SMWeinbergEvolution')
        sm.SMWeinbergInitialConditions = lambda **kwargs: SimpleNamespace(**kwargs)
        sm.evolve_sm_weinberg = fake_evolve
        rough = ModuleType('Numerical.parallel_t3.threshold.rough_numerics')
        rough.neutrino_observables = fake_neutrino
        completion = ModuleType('Numerical.parallel_t3.complete.full_general_completion')
        completion.score_nufit = fake_score
        import Numerical.running
        import Numerical.parallel_t3.complete
        with patch.dict(sys.modules, {
            'Numerical.running.SMWeinbergEvolution':sm,
            'Numerical.parallel_t3.threshold': ModuleType('Numerical.parallel_t3.threshold'),
            'Numerical.parallel_t3.threshold.rough_numerics':rough,
            'Numerical.parallel_t3.complete.full_general_completion':completion,
        }):
            doc,arr=run_smeft_observables(boundary,s,{'ordering':'NO'})
        self.assertTrue(doc['SMEFT_one_loop_solver_success'])
        self.assertFalse(doc['full_physical_C5_verified'])
        self.assertFalse(doc['double_counted_separate_formal_RG_term_added'])
        self.assertEqual(doc['nufit_gaussian_approximation']['chi2_gaussian_uncorrelated'],123.)
        self.assertEqual(arr['mnu_low_approx_GeV'].shape,(3,3))


if __name__=='__main__':unittest.main()
