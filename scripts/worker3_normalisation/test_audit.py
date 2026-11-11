"""Run with python -m unittest discover -s scripts/worker3_normalisation -p test_*.py -v."""
import importlib.util
from pathlib import Path
import unittest
import numpy as np
try:
    from scripts.worker3_normalisation import audit
except ModuleNotFoundError:
    import audit

class ScalingTests(unittest.TestCase):
    def setUp(self):
        self.m=audit.synthetic_matrix(101)
        self.a=audit.synthetic_matrix(102)
        self.b=audit.synthetic_matrix(103)

    def test_majorana_matrix_is_symmetric(self):
        np.testing.assert_allclose(self.m,self.m.T)

    def test_real_global_mass_and_mass_squared(self):
        v=audit.invariant_checks(self.m,-3.1,2+0.8j)
        self.assertLess(v['mass_real_error'],1e-12)
        self.assertLess(v['dm_real_error'],1e-12)

    def test_angles_ratios_pmns_global_real_and_complex(self):
        v=audit.invariant_checks(self.m,2.3,0.7+1j)
        for key in ('mass_ratio_real_error','angle_real_error','pmns_abs_real_error','angle_complex_global_error','mass_complex_global_error'):
            self.assertLess(v[key],1e-12,key)

    def test_zero_and_degenerate_masses_require_caution(self):
        with self.assertRaises(ValueError):audit.obs(np.eye(3,dtype=complex))

    def test_real_linear_beta_exact_obstruction(self):
        nr,_=audit.phase_obstruction(self.m,0.2+0.1j,0.4-0.3j,3.2)
        nc,_=audit.phase_obstruction(self.m,0.2+0.1j,0.4-0.3j,2j)
        self.assertLess(nr,1e-12)
        self.assertGreater(nc,1e-3)

    def test_complex_linear_exception(self):
        nc,_=audit.phase_obstruction(self.m,0.2+0.1j,0,2j)
        self.assertLess(nc,1e-12)

    def test_relative_channel_conversion_changes_shape(self):
        v=audit.complex_channel_checks(self.a,self.b,2.5)
        self.assertGreater(v['angles_shift'],1e-6)
        self.assertGreater(v['mass_ratio_shift'],1e-6)

    def test_collinear_channel_exception(self):
        v=audit.complex_channel_checks(self.a,2*self.a,2.5)
        self.assertLess(v['angles_shift'],1e-12)

    def test_stage_scaling(self):
        v=audit.stage_ratio_check(self.a,self.b,-7)
        self.assertLess(v['error'],1e-12)

    def test_higgs_four_real_components(self):
        # Model's neutral H0=(hR+i hI)/sqrt(2); Weinberg L=1/2 C5 (nu H0)^2
        # Canonical master L=C_ijab psi_i psi_j phi_a phi_b + h.c. => C_nunu,hR,hR=C5/4.
        h=np.array([1/np.sqrt(2),1j/np.sqrt(2)])
        c5=0.2+0.3j
        t=0.5*c5*np.outer(h,h)
        self.assertAlmostEqual(abs(t[0,0]-c5/4),0)
        self.assertAlmostEqual(abs(t[1,1]+c5/4),0)
        self.assertAlmostEqual(abs(t[0,1]-1j*c5/4),0)

    def test_production_projector_when_present(self):
        root=Path.cwd()
        path=root/'RGE/running/intermediate/GeneralEFT1C5Mixing.py'
        if not path.is_file():self.skipTest('Production source not installed in this environment')
        spec=importlib.util.spec_from_file_location('worker3_original_projection',path)
        module=importlib.util.module_from_spec(spec)
        import sys
        sys.modules[spec.name]=module
        spec.loader.exec_module(module)
        from types import SimpleNamespace
        c=SimpleNamespace(nf=6,ns=4,leptons=((0,1),(2,3),(4,5)),higgs=(0,1,2,3))
        mat=audit.synthetic_matrix(200)
        tensor=module._weinberg_template(mat,c)
        recovered,ward,sym=module.project_c5(tensor,c)
        np.testing.assert_allclose(recovered,mat)
        self.assertLess(ward,1e-12)
        self.assertLess(sym,1e-12)
        # Negative control: modifying only imaginary neutral component breaks Weinberg orbit.
        corrupted=tensor.copy();corrupted[0,0,3,3]+=0.01
        _,ward,_=module.project_c5(corrupted,c)
        self.assertGreater(ward,0.009)

if __name__=='__main__':unittest.main()

class CoupledIVPTests(unittest.TestCase):
    def test_two_stage_coupled_equivariance(self):
        try:
            from scripts.worker3_normalisation.ivp_regression import integration
        except ModuleNotFoundError:
            from ivp_regression import integration
        x=(audit.synthetic_matrix(121),audit.synthetic_matrix(122))
        original=integration(x)
        multiplied=integration(x,scale=-3)
        for a,b in zip(original,multiplied):
            np.testing.assert_allclose(b,-3*a,rtol=1e-8,atol=1e-8)
