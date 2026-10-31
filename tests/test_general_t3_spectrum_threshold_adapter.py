"""Workstream D candidate spectra: unit tests using existing source inventory.

All numbers in this module are SYNTHETIC algorithm probes, NOT T3 benchmarks.
"""
from pathlib import Path
import math
import unittest

import numpy as np

from Numerical.experimental.GeneralT3ScalarQuadraticMasses import (
    REP, load_manifest, quadratic_kind, quadratic_mass_hessian,
)
from Numerical.experimental.GeneralT3MixedScalarMassBasis import MixedDoubletMassBasis
from Numerical.experimental.GeneralT3SpectrumThresholdAdapter import (
    DEFAULT_MANIFEST, RunningMassSample, scalar_spectrum,
    ordered_mass_events, find_mass_crossings,
)


class TestGeneralT3SpectrumThresholdAdapter(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.manifest = load_manifest(DEFAULT_MANIFEST)

    def test_all_16_match_existing_full_hessians(self):
        self.assertEqual(len(self.manifest['models']), 16)
        expected_kinds = {
            'A_-4':'none','A_-2':'self_s2','A_+0':'self_s1',
            'B_-3':'none','B_-1':'mixed','B_+1':'none',
            'C_-3':'none','C_-1':'mixed','C_+1':'none',
            'D_-2':'self_s2','D_+0':'self_s1','D_+2':'none',
            'E_-4':'none','E_-2':'self_s2','E_+0':'self_s1','E_+2':'none',
        }
        self.assertEqual({m:quadratic_kind(self.manifest,m)
                          for m in self.manifest['models']}, expected_kinds)
        for model in self.manifest['models']:
            with self.subTest(model=model):
                kind = expected_kinds[model]
                mu2 = (17_000 + 3_000j) if kind != 'none' else 0j
                spec = scalar_spectrum(model, 300_000., 520_000., mu2,
                                       manifest=self.manifest).require_positive()
                raw = quadratic_mass_hessian(self.manifest,model,300_000.,520_000.,mu2)
                np.testing.assert_allclose(np.sort(np.linalg.eigvalsh(spec.hessian_gev2)),
                                           raw['mass_squared_eigenvalues_gev2'], rtol=2e-12,atol=2e-9)
                self.assertEqual(sum(m.real_dof for m in spec.modes),
                                 2 * sum(REP[model.split('_')[0]][:2]))
                self.assertEqual(spec.kind, kind)
                self.assertTrue(all(x.real_dof >= 1 for x in spec.modes))

    def test_b_c_mixed_basis_is_exact_legacy_object(self):
        for model in ('B_-1','C_-1'):
            with self.subTest(model=model):
                s = scalar_spectrum(model, 300_000.,520_000.,13_000.+4_000j,
                                    manifest=self.manifest)
                legacy = MixedDoubletMassBasis.build(model,300_000.,520_000.,13_000.+4_000j)
                self.assertIsNotNone(s.mixed_doublet_basis)
                np.testing.assert_allclose(s.mixed_doublet_basis.U,legacy.U,rtol=0,atol=1e-14)
                np.testing.assert_allclose([s.by_label('X').mass_squared_gev2,
                                           s.by_label('Y').mass_squared_gev2],
                                          legacy.eigenvalues_sq,rtol=1e-13,atol=1e-9)
                self.assertEqual([m.real_dof for m in s.modes],[4,4])
                x=np.array([1.+0.1j,-0.3+0.2j]);y=np.array([.8-0.3j,.5+0.4j])
                s1,s2=legacy.original_fields(x,y)
                self.assertAlmostEqual(legacy.old_potential(s1,s2),legacy.new_potential(x,y),places=7)

    def test_self_pair_split_is_two_mu2_not_mu2(self):
        s = scalar_spectrum('A_+0',1000.,2000.,100j,manifest=self.manifest)
        self.assertEqual(s.by_label('S1_minus').mass_squared_gev2,800.)
        self.assertEqual(s.by_label('S1_plus').mass_squared_gev2,1200.)
        self.assertEqual(s.by_label('S2').mass_squared_gev2,2000.)

    def test_tachyon_and_forbidden_bilinear_fail_closed(self):
        s=scalar_spectrum('A_+0',100.,800.,60j,manifest=self.manifest)
        self.assertFalse(s.positive)
        with self.assertRaisesRegex(ValueError,'non-positive'):
            ordered_mass_events(s)
        with self.assertRaises(ValueError):
            scalar_spectrum('A_-4',100.,200.,1j,manifest=self.manifest)

    def test_event_order_is_from_mass_not_model(self):
        s=scalar_spectrum('A_-4',300**2,450**2,manifest=self.manifest)
        events=ordered_mass_events(s, [200.,600.])
        self.assertEqual([e.particle for e in events],
                         ['fermion:1','scalar:S2','scalar:S1','fermion:0'])
        np.testing.assert_allclose([e.mass_gev for e in events],[600,450,300,200])
        self.assertTrue(all(e.status=='spectrum_only_matching_not_applied' for e in events))

    def test_constant_running_crossings_all_types(self):
        samples=[RunningMassSample(float(mu),300.**2,450.**2,0j,(200.,600.))
                 for mu in (1000.,700.,520.,400.,250.,100.)]
        roots=find_mass_crossings('A_-4',samples,manifest=self.manifest)
        self.assertEqual([r.particle for r in roots],
                         ['fermion:1','scalar:S2','scalar:S1','fermion:0'])
        np.testing.assert_allclose([r.scale_gev for r in roots],
                                   [600,450,300,200],rtol=1e-12,atol=1e-9)
        self.assertTrue(all(r.relative_mass_scale_residual < 1e-10 for r in roots))

    def test_mixed_crossings_are_legacy_mass_roots(self):
        m1,m2,mu2=300.**2,430.**2,9_000.+7_000j
        legacy=MixedDoubletMassBasis.build('B_-1',m1,m2,mu2)
        samples=[RunningMassSample(t,m1,m2,mu2) for t in (900.,180.)]
        roots=find_mass_crossings('B_-1',samples,manifest=self.manifest)
        self.assertEqual([r.particle for r in roots],['scalar:Y','scalar:X'])
        np.testing.assert_allclose([r.scale_gev for r in roots],
                                   np.sqrt(legacy.eigenvalues_sq[::-1]),
                                   rtol=1e-12,atol=1e-9)

    def test_running_interpolation_diagonalizes_at_root(self):
        samples=[RunningMassSample(850.,300**2,500**2,5000+4000j),
                 RunningMassSample(170.,350**2,440**2,8000-2000j)]
        roots=find_mass_crossings('C_-1',samples,manifest=self.manifest)
        self.assertEqual([r.particle for r in roots],['scalar:Y','scalar:X'])
        for r in roots:
            self.assertLess(r.relative_mass_scale_residual,1e-10)

    def test_existing_report_audit_contract(self):
        from scripts.AuditGeneralT3MassCrossings import audit_record
        rows=[{'mu_gev':mu,'mS1Sq_gev2':300**2,'mS2Sq_gev2':450**2,
               'mu2_gev2':{'real':0.,'imag':0.}} for mu in (1000.,500.,200.)]
        report={'format':'t3-continuous-general-scalar-candidate-v1',
                'status':'CandidateContinuousUVFermionThresholdEFT1Integrated',
                'model':'A_-4','uv_trajectory':rows,'eft1_trajectory':rows,
                'uv_fermion_threshold_state':{'mS1Sq':300**2,'mS2Sq':450**2},
                'quadratic_mu2_at_fermion_threshold_gev2':{'real':0.,'imag':0.},
                'heavy_fermion_threshold_masses_gev':[600.,200.]}
        result=audit_record(report,self.manifest)
        self.assertFalse(result['finite_matching_applied'])
        self.assertFalse(result['uv_fermion_crossings_available'])
        self.assertEqual([e['particle'] for e in result['at_grouped_F_mass_order']],
                         ['fermion:0','scalar:S2','scalar:S1','fermion:1'])
        self.assertEqual([e['particle'] for e in result['uv_scalar_crossings']],
                         ['scalar:S2','scalar:S1'])

    def test_no_crossing_not_assumed(self):
        samples=[RunningMassSample(1000.,100**2,110**2),
                 RunningMassSample(500.,100**2,110**2)]
        self.assertEqual(find_mass_crossings('A_-4',samples,manifest=self.manifest),())

    def test_reject_missing_or_nonphysical_trajectory(self):
        s=[RunningMassSample(1000.,300**2,400**2),
           RunningMassSample(1300.,300**2,400**2)]
        with self.assertRaisesRegex(ValueError,'descending'):
            find_mass_crossings('B_+1',s,manifest=self.manifest)
        with self.assertRaisesRegex(ValueError,'inventory'):
            find_mass_crossings('B_+1',[
                RunningMassSample(1000.,300**2,400**2,0j,(600.,)),
                RunningMassSample(100.,300**2,400**2)],manifest=self.manifest)


if __name__=='__main__':
    unittest.main()
