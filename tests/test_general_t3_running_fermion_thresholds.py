"""Workstream D Stage 2: synthetic heavy-mass crossing tests (not T3 predictions)."""
from __future__ import annotations

import math
import unittest
from types import SimpleNamespace

import numpy as np

from Numerical.experimental.GeneralT3RunningFermionThresholds import (
    mass_matrix_from_vector, singular_masses, locate_fermion_crossings,
)


class TestRunningFermionThresholds(unittest.TestCase):
    def test_singular_values_are_takagi_masses_for_symmetric_complex_mf(self):
        phi=np.exp(1j*np.array([.4,1.2,-.6]))
        m=np.diag(np.array([200.,300.,500.])*phi)
        np.testing.assert_allclose(singular_masses(m,majorana=True),[200,300,500],atol=1e-11)
        u=np.array([[1,0,0],[0,np.sqrt(.5),1j*np.sqrt(.5)],
                    [0,1j*np.sqrt(.5),np.sqrt(.5)]],complex)
        m=u.conj()@m@u.conj().T
        self.assertTrue(np.allclose(m,m.T))
        np.testing.assert_allclose(singular_masses(m,majorana=True),[200,300,500],atol=1e-10)

    def test_majorana_rejects_nonsymmetric_but_vectorlike_accepts(self):
        m=np.diag([100.,200.,300.]).astype(complex)
        m[0,1]=4j
        with self.assertRaisesRegex(ValueError,'symmetric'):
            singular_masses(m,majorana=True)
        self.assertEqual(len(singular_masses(m,majorana=False)),3)

    def test_degenerate_and_nonpositive_fail_closed(self):
        with self.assertRaisesRegex(ValueError,'degenerate'):
            singular_masses(np.diag([100.,100.,200.]),majorana=False)
        with self.assertRaisesRegex(ValueError,'nonpositive'):
            singular_masses(np.diag([0.,100.,200.]),majorana=True)

    def test_all_three_constant_crossings_mass_sorted(self):
        m=np.diag([2e5,4e5,6e5])
        r=locate_fermion_crossings(lambda t:m,1e7,1e5,majorana=True)
        self.assertTrue(r['all_three_crossings_bracketed'])
        self.assertFalse(r['sequential_heavy_field_decoupling_performed'])
        self.assertEqual([x['particle'] for x in r['crossings']],
                         ['fermion:2','fermion:1','fermion:0'])
        np.testing.assert_allclose([x['scale_gev'] for x in r['crossings']],
                                   [6e5,4e5,2e5],rtol=1e-12)
        self.assertLess(max(x['relative_mass_scale_residual'] for x in r['crossings']),1e-10)

    def test_running_masses_evaluated_at_root_not_frozen(self):
        ref=1e6
        power=.15
        base=np.array([2e5,4e5,6e5])
        fn=lambda t: np.diag(base*(math.exp(t)/ref)**power)
        r=locate_fermion_crossings(fn,1e7,1e5,majorana=True)
        expected=(base/ref**power)**(1/(1-power))
        np.testing.assert_allclose([e['scale_gev'] for e in r['crossings']],
                                   expected[::-1],rtol=1e-12)

    def test_below_grouped_endpoint_not_extrapolated(self):
        m=np.diag([99999.99962690244,110000.,120000.])
        r=locate_fermion_crossings(lambda t:m,1e7,1e5,majorana=True)
        self.assertEqual(len(r['crossings']),2)
        self.assertFalse(r['all_three_crossings_bracketed'])
        self.assertEqual(r['unresolved_crossings'][0]['particle'],'fermion:0')
        self.assertEqual(r['unresolved_crossings'][0]['status'],
                         'below_saved_UV_stop_undetermined')

    def test_above_uv_start_not_extrapolated(self):
        m=np.diag([2e5,4e5,2e7])
        r=locate_fermion_crossings(lambda t:m,1e7,1e5,majorana=False)
        self.assertEqual([x['particle'] for x in r['crossings']],['fermion:1','fermion:0'])
        self.assertEqual(r['unresolved_crossings'][0]['status'],
                         'above_saved_UV_start_undetermined')

    def test_inverse_crossing_rejected(self):
        # Mass grows toward high scales faster than mu: descending trajectory
        # would carry this state from below to above mu (inverse decoupling).
        def fn(t):
            mu=math.exp(t)
            m0=(mu/1e6)**2*1e6
            return np.diag([m0,1e10,2e10])
        with self.assertRaisesRegex(ValueError,'inverse crossing'):
            locate_fermion_crossings(fn,1e7,1e5,majorana=False)

    def test_double_crossing_rejected(self):
        def fn(t):
            mu=math.exp(t)
            f=1+0.3*math.sin(4*math.pi*(t-math.log(1e5)) / math.log(100)+0.2)
            return np.diag([mu*f,2e8,4e8])
        with self.assertRaises(ValueError):
            locate_fermion_crossings(fn,1e7,1e5,majorana=False,scan_points=129)

    def test_layout_real_imag_matrix_is_read_correctly(self):
        class Layout:
            size=20
            def block(self,name):
                self_name=name
                if self_name!='MF':raise KeyError(name)
                return SimpleNamespace(name='MF',kind='complex_matrix',shape=(3,3),start=2,stop=20)
        m=np.array([[1+2j,3+4j,5+6j],[7+8j,9+10j,11+12j],[13+14j,15+16j,17+18j]])
        vec=np.r_[[-1.,-2.],m.real.ravel(),m.imag.ravel()]
        np.testing.assert_allclose(mass_matrix_from_vector(vec,Layout()),m)
        with self.assertRaisesRegex(ValueError,'packed'):
            mass_matrix_from_vector(np.zeros(19),Layout())

    def test_full_uv_state_with_self_or_mixed_bilinear(self):
        """A full UV ODE has two extra coordinates in all eight bilinear models."""
        class Layout:
            size = 20
            def block(self, name):
                if name != 'MF':
                    raise KeyError(name)
                return SimpleNamespace(name='MF', kind='complex_matrix',
                                       shape=(3,3), start=2, stop=20)

        mass = np.array([[1+2j, 0, 0], [0, 2+3j, 0], [0, 0, 3+4j]])
        packed = np.r_[[5.0, -7.0], mass.real.ravel(), mass.imag.ravel()]
        mu2 = np.array([1234.5, -98.7])  # mass dimension 2 (GeV^2)
        full_ode = np.r_[packed, mu2]

        # The exact bug reported for A_+0, B_-1, C_-1, etc.
        with self.assertRaisesRegex(ValueError, 'Invalid packed RGBeta MF/layout'):
            mass_matrix_from_vector(full_ode, Layout())
        np.testing.assert_allclose(
            mass_matrix_from_vector(full_ode, Layout(), trailing_components=2), mass)
        # The same interpretation is required at every dense-solution sample.
        for z in (0.1, 0.5, 0.9):
            sample = np.r_[packed, z*mu2]
            np.testing.assert_allclose(
                mass_matrix_from_vector(sample, Layout(), trailing_components=2), mass)

    def test_bilinear_suffix_fails_closed_on_shape_or_nonfinite(self):
        class Layout:
            size=20
            def block(self, name):
                if name!='MF':raise KeyError(name)
                return SimpleNamespace(kind='complex_matrix',shape=(3,3),
                                       start=2,stop=20)
        base=np.r_[[0.0,0.0], np.eye(3).ravel(), np.zeros(9)]
        with self.assertRaisesRegex(ValueError, 'Invalid packed'):
            mass_matrix_from_vector(np.r_[base, 1.0], Layout(), trailing_components=2)
        with self.assertRaisesRegex(ValueError, 'Invalid packed'):
            mass_matrix_from_vector(np.r_[base, 1.0, 2.0, 3.0], Layout(),
                                    trailing_components=2)
        with self.assertRaisesRegex(ValueError, 'Invalid packed'):
            mass_matrix_from_vector(np.r_[base, 1.0, float('nan')], Layout(),
                                    trailing_components=2)
        with self.assertRaisesRegex(ValueError, 'Invalid packed'):
            mass_matrix_from_vector(np.r_[base, 1.0, 2.0], Layout(),
                                    trailing_components=3)

    def test_invalid_ode_span_and_scan_rejected(self):
        f=lambda t:np.diag([2e5,4e5,6e5])
        with self.assertRaisesRegex(ValueError,'bounds'):
            locate_fermion_crossings(f,1e5,1e7,majorana=True)
        with self.assertRaisesRegex(ValueError,'scan_points'):
            locate_fermion_crossings(f,1e7,1e5,majorana=True,scan_points=4)


if __name__=='__main__':
    unittest.main()
