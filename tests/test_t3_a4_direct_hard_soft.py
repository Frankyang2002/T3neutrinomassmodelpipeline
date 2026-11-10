"""Offline regression tests for the fresh-source direct F/S Laurent bridge."""
from __future__ import annotations

import math
import unittest
import json
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import sympy as sp

from Numerical.parallel_t3.background.general_batch1_api import encode_native
from Numerical.parallel_t3.complete.a4_direct_hard_soft import (
    DirectBridgeError, _mass_log_divided_difference,
    audit_direct_sources, audit_hard_first_threshold, audit_scalar_soft,
    formal_numeric_direct_remainder, split_matchete_sum,
)


def make_hard() -> str:
    terms=[]
    for chi in ('PR','PL'):
        if chi == 'PL':
            ys = 'Bar[Coupling[y1, {Index[p, Flavor], Index[r, NFlavor]}, 0]]*Bar[Coupling[y2, {Index[q, Flavor], Index[r, NFlavor]}, 0]]'
            lam = 'Coupling[lambdaT3, {}, 0]'
            cg = 'CG[Bar[eps[SU2L]], {Bar[Index[i, SU2L[fund]]], Bar[Index[j, SU2L[fund]]]}]'
            higgs = 'Field[H, Scalar, {Index[i, SU2L[fund]]}, {}]*Field[H, Scalar, {Index[j, SU2L[fund]]}, {}]'
            proj='Proj[-1]'
        else:
            ys = 'Coupling[y1, {Index[p, Flavor], Index[r, NFlavor]}, 0]*Coupling[y2, {Index[q, Flavor], Index[r, NFlavor]}, 0]'
            lam = 'Bar[Coupling[lambdaT3, {}, 0]]'
            cg = 'CG[eps[SU2L], {Index[i, SU2L[fund]], Index[j, SU2L[fund]]}]'
            higgs = 'Bar[Field[H, Scalar, {Index[i, SU2L[fund]]}, {}]]*Bar[Field[H, Scalar, {Index[j, SU2L[fund]]}, {}]]'
            proj='Proj[1]'
        for kind in ('constant','pole','log'):
            inside = f'2*hbar*{lam}*{ys}*{cg}*{cg}*{higgs}'
            if kind == 'log':
                inside += '*Log[μbar2/Coupling[MF, {Index[r, NFlavor]}, 0]^2]'
            inside += (f'*NCM[Field[l, Fermion, {{Index[i, SU2L[fund]]}}, {{}}],'
                       f'Field[l, Fermion, {{Index[j, SU2L[fund]]}}, {{}}],'
                       f'DiracProduct[GammaCC, {proj}]]')
            eps='ϵ*' if kind=='pole' else ''
            terms.append(f'({inside})/(Sqrt[3]*{eps}Coupling[MF, {{Index[r, NFlavor]}}, 0])')
    return ' + '.join(terms)


SOFT = (
    '(-2*((Sqrt[3]*hbar*Bar[Coupling[T3EFT1W12, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*Coupling[lambdaT3, {}, 0]*Coupling[MS1, {}, 0]^2)/(Coupling[MS1, {}, 0]^2 - Coupling[MS2, {}, 0]^2) '
    '+ (Sqrt[3]*hbar*Bar[Coupling[T3EFT1W12, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*Coupling[lambdaT3, {}, 0]*Coupling[MS1, {}, 0]^2)/(ϵ*(Coupling[MS1, {}, 0]^2 - Coupling[MS2, {}, 0]^2)) '
    '- (Sqrt[3]*hbar*Bar[Coupling[T3EFT1W12, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*Coupling[lambdaT3, {}, 0]*Coupling[MS2, {}, 0]^2)/(Coupling[MS1, {}, 0]^2 - Coupling[MS2, {}, 0]^2) '
    '- (Sqrt[3]*hbar*Bar[Coupling[T3EFT1W12, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*Coupling[lambdaT3, {}, 0]*Coupling[MS2, {}, 0]^2)/(ϵ*(Coupling[MS1, {}, 0]^2 - Coupling[MS2, {}, 0]^2)) '
    '+ (Sqrt[3]*hbar*Bar[Coupling[T3EFT1W12, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*Coupling[lambdaT3, {}, 0]*Coupling[MS1, {}, 0]^2*Log[μbar2/Coupling[MS1, {}, 0]^2])/(Coupling[MS1, {}, 0]^2 - Coupling[MS2, {}, 0]^2) '
    '- (Sqrt[3]*hbar*Bar[Coupling[T3EFT1W12, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*Coupling[lambdaT3, {}, 0]*Coupling[MS2, {}, 0]^2*Log[μbar2/Coupling[MS2, {}, 0]^2])/(Coupling[MS1, {}, 0]^2 - Coupling[MS2, {}, 0]^2)))/3'
)


def make_f():
    y1 = np.array([[1+1j,2,3],[.1,.2,1j],[.4,.5,.6]],complex)*1e-4
    y2 = np.array([[1,0,.1],[.4j,1,.2],[.7,.9,.3]],complex)*1e-4
    return {'model':'A_-4', 'scale_GeV':1e5,
            'all_heavy_fermions_to_integrate_out_together':True,
            'native_uv_rgbeta_couplings_at_muF': {
                'MF':encode_native(np.diag([1e5,1.1e5,1.2e5]).astype(complex)),
                'y1':encode_native(y1),'y2':encode_native(y2),
                'lambdaT3':encode_native(complex(.1))}}


def make_s(a=1.6e6,b=1.4e6):
    return {'model':'A_-4', 'scale_GeV':1e3,
            'both_heavy_scalar_multiplets_to_integrate_out_together':True,
            'scalar_quadratic_terms_GeV2':{'mS1Sq':a,'mS2Sq':b}}


class DirectBridgeTests(unittest.TestCase):
    def test_balanced_split_not_internal(self):
        self.assertEqual(split_matchete_sum('Log[x+y] + 3 - (a+b)'),
                         ['Log[x+y]','+ 3','- (a+b)'])
        with self.assertRaises(DirectBridgeError):
            split_matchete_sum('Log[x+y')

    def test_hard_six_complete_source_triplets(self):
        count, groups = audit_hard_first_threshold(make_hard())
        self.assertEqual(count, 6)
        self.assertEqual(set(groups),{'holomorphic_HH','anti_holomorphic_HdagHdag'})
        self.assertEqual(len(groups['holomorphic_HH']),3)

    def test_hard_rejects_missing_pole_and_old_CG(self):
        s = make_hard()
        with self.assertRaises(DirectBridgeError):
            audit_hard_first_threshold(s.replace('ϵ*','',1))
        with self.assertRaises(DirectBridgeError):
            audit_hard_first_threshold(s.replace('CG[Bar[eps[SU2L]]', 'CG[eps[SU2L]'))

    def test_native_soft_exact_laurent_expression(self):
        pole, finite = audit_scalar_soft(SOFT)
        self.assertIn('WB',pole)
        self.assertIn('log',finite)
        self.assertTrue(audit_direct_sources(make_hard(),SOFT).formal_ordered_pole_cancelled)

    def test_native_soft_rejects_changed_residue(self):
        with self.assertRaises(DirectBridgeError):
            audit_scalar_soft(SOFT.replace('(-2*','(-3*',1))

    def test_stable_equal_scalar_mass_limit(self):
        m2=np.array([1e10,1.21e10,1.44e10])
        b=1.5e6
        actual=_mass_log_divided_difference(b,b,m2)
        np.testing.assert_allclose(actual,np.log(b/m2)+1,rtol=1e-14,atol=1e-14)
        near=_mass_log_divided_difference(b*(1+1e-11),b,m2)
        np.testing.assert_allclose(near,actual,rtol=1e-10,atol=1e-10)

    def test_three_flavour_numeric_pole_and_finite(self):
        result=formal_numeric_direct_remainder(make_f(),make_s())
        self.assertFalse(result['physical_C5_computed'])
        self.assertEqual(len(result['all_three_heavy_masses_GeV']),3)
        from Numerical.parallel_t3.background.general_batch1_api import decode_native
        p=np.asarray(decode_native(result['ordered_hard_pole_GeVinv']),complex)
        q=np.asarray(decode_native(result['ordered_soft_pole_GeVinv']),complex)
        np.testing.assert_array_equal(p+q,np.zeros((3,3)))
        mat=np.asarray(decode_native(result['ordered_formal_finite_remainder_GeVinv']),complex)
        self.assertTrue(np.isfinite(mat).all())
        self.assertGreater(np.linalg.norm(mat),0)

    def test_offdiagonal_MF_blocks(self):
        f=make_f()
        mf=np.diag([1e5,1.1e5,1.2e5]).astype(complex)
        mf[0,1]=1e-1
        f['native_uv_rgbeta_couplings_at_muF']['MF']=encode_native(mf)
        with self.assertRaises(DirectBridgeError):
            formal_numeric_direct_remainder(f,make_s())

    def test_end_to_end_source_contract(self):
        from scripts.RunT3A4DirectHardSoft import run
        with TemporaryDirectory() as folder:
            root=Path(folder)
            a4=root/'a4'
            b1=root/'batch1'
            (a4/'data').mkdir(parents=True)
            b1.mkdir()
            (a4/'data/eft1_after_F_one_loop_inputform.txt').write_text(
                make_hard(),encoding='utf-8')
            (a4/'A_scalar_soft_only_unrenormalized.txt').write_text(
                SOFT,encoding='utf-8')
            fixtures={
                a4/'native_two_threshold_status.json': {
                    'status':'A_PhysicalScalarSoftOnly_HardInheritanceAndWilsonRGPending',
                    'model':'A_-4','heavy_fermion_flavours':3,
                    'fermion_groups':[['F']], 'scalar_groups':[['S1','S2']],
                    'first_threshold_reused':False,
                    'A_joint_scalar_match_performed':True,
                    'A_physical_scalar_soft_Weinberg_terms':10,
                    'exploratory_source_theory_used':False},
                a4/'A_joint_scalar_lagrangian_audit.json': {
                    'CheckLagrangian_original':'True',
                    'accepted_canonicalisation':'ExactSavedTree',
                    'projection_used_for_physical_match':False},
                a4/'A_joint_scalar_match_diagnostics.json': {
                    'matching_at_one_loop_succeeded':True,
                    'source_theory':'ExactSavedTree',
                    'exploratory_non_equivalent_source':False},
                a4/'A_EFT1_Wilson_definitions.json': {
                    'W12[p,q]':'Sum[y1[p,r]*y2[q,r]/MF[r],{r,1,3}]',
                    'all_three_fermions_integrated_together':True},
                b1/'fermion_matching_input.json': make_f(),
                b1/'scalar_matching_input.json': make_s(),
            }
            for filename, data in fixtures.items():
                filename.write_text(json.dumps(data),encoding='utf-8')
            report=run(a4,b1)
            self.assertEqual(report['status'],
                'SourceVerifiedFormalDirectChannel_PolePairedNotFullC5')
            self.assertFalse(report['physical_C5_computed'])
            self.assertEqual(report['source_audit']['total_first_threshold_terms'],6)
            self.assertEqual(len(report['source_files']),8)

    def test_no_D2_substitution(self):
        with self.assertRaises(DirectBridgeError):
            formal_numeric_direct_remainder(make_f()|{'model':'D_-2'}, make_s())


if __name__=='__main__':
    unittest.main()
