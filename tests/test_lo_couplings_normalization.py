"""Independent algebraic tests for the Workstream 2 portal normalization audit."""
from __future__ import annotations

from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
import json

import numpy as np
import sympy as sp

from Numerical.lo_certification.couplings.portal_normalization import (
    TensorAuditError, su2_generators, real_bilinear_hessian,
    rgbeta_portal_columns, native_portal_columns, portal_species_audit,
    construct_rgbeta_portal_background, real_running_coupling)
from Numerical.lo_certification.couplings.quartic_audit import run_portal_batch

DIMS={'A':(1,3),'B':(2,2),'C':(2,2),'D':(3,1),'E':(3,3)}


def native_mock(model, *, remove=(), shift=None):
    d1,d2=DIMS[model[0]]
    channels={}
    for i,d in enumerate((d1,d2),1):
        for j,(kind,rows) in enumerate(rgbeta_portal_columns(d).items(),1):
            name=f'lambdaH{i}Inv{j}'
            if name in remove:
                continue
            offset=0 if i==1 else 2*d1
            channels[name]={(A,B,a+offset,b+offset):v*sp.Symbol(name)
                            for (A,B,a,b),v in rows.items()}
    if shift:
        name,k,v=shift
        channels[name][k]=sp.simplify(channels[name].get(k,0)+sp.Symbol(name)*v)
    return SimpleNamespace(model=model,by_name=channels)


class TestExactRGBetaPortalNormalization(unittest.TestCase):
    def test_su2_generators_obey_algebra(self):
        for d in (1,2,3):
            a,b,c=su2_generators(d)
            self.assertEqual(a*b-b*a,sp.I*c)
            self.assertEqual(b*c-c*b,sp.I*a)
            self.assertEqual(c*a-a*c,sp.I*b)

    def test_real_bilinear_hessian_exact(self):
        T=su2_generators(2)[1]
        r=sp.symbols('r0:4', real=True)
        z=sp.Matrix([(r[2*k]+sp.I*r[2*k+1])/sp.sqrt(2) for k in range(2)])
        polynomial=sp.expand((z.conjugate().T*T*z)[0])
        k=real_bilinear_hessian(T)
        reconstructed=(sp.Matrix(r).T*k*sp.Matrix(r))[0]/2
        self.assertEqual(sp.simplify(polynomial-reconstructed),0)
        self.assertEqual(k,k.T)

    def test_v4_exact_fourth_derivative_factor(self):
        h=sp.symbols('h0:4',real=True)
        s=sp.symbols('s0:4',real=True)
        H=sp.Matrix([(h[2*k]+sp.I*h[2*k+1])/sp.sqrt(2) for k in range(2)])
        S=sp.Matrix([(s[2*k]+sp.I*s[2*k+1])/sp.sqrt(2) for k in range(2)])
        potential=-sp.expand((H.conjugate().T*H)[0]*(S.conjugate().T*S)[0])
        self.assertEqual(sp.diff(potential,h[0],h[0],s[0],s[0]),-1)
        self.assertEqual(rgbeta_portal_columns(2)['singlet'][(0,0,0,0)],-1)

    def test_all_representation_dimensions_exact_unit_mapping(self):
        for model in ('A_-4','B_-1','C_-1','D_+0','E_-2'):
            with self.subTest(model=model):
                ch=native_mock(model)
                for species in (1,2):
                    out=portal_species_audit(ch,species)
                    self.assertEqual(out['status'],'exact_full_tensor_span')
                    self.assertEqual(out['rgbeta_portal_directions_missing_from_native_span'],0)
                    self.assertTrue(out['normalization_conversion_physically_usable_for_portal_subsector'])

    def test_exact_nontrivial_normalization_rotation(self):
        model='B_-1'
        ch=native_mock(model)
        a,b=ch.by_name['lambdaH1Inv1'],ch.by_name['lambdaH1Inv2']
        x=sp.Symbol('lambdaH1Inv1')
        y=sp.Symbol('lambdaH1Inv2')
        ch.by_name['lambdaH1Inv1']={
                                     k:sp.simplify(2*a.get(k,0)+3*b.get(k,0).subs(y,x))
                                     for k in set(a)|set(b)}
        ch.by_name['lambdaH1Inv2']={k:sp.simplify(-b.get(k,0)/5)
                                    for k in b}
        result=portal_species_audit(ch,1)
        self.assertEqual(result['status'],'exact_full_tensor_span')
        self.assertEqual(result['native_unit_to_rgbeta_invariant_coefficients'],[['2','0'],['3','-1/5']])
        self.assertEqual(result['rgbeta_couplings_to_native_couplings']['lambdaH1Inv1']['lambdaH1'],'1/2')

    def test_missing_baseline_portal_direction_is_detected(self):
        ch=native_mock('A_-4',remove=('lambdaH2Inv2',))
        out=portal_species_audit(ch,2)
        self.assertEqual(out['status'],'exact_partial_tensor_span')
        self.assertEqual(out['rgbeta_portal_directions_missing_from_native_span'],1)
        self.assertFalse(out['normalization_conversion_physically_usable_for_portal_subsector'])

    def test_missing_one_of_two_portals_is_not_certified(self):
        ch=native_mock('B_+1',remove=('lambdaH1Inv2',))
        out=portal_species_audit(ch,1)
        self.assertEqual(out['rgbeta_portal_directions_missing_from_native_span'],1)
        self.assertNotIn('rgbeta_couplings_to_native_couplings',out)

    def test_out_of_span_tensor_is_not_fitted_away(self):
        ch=native_mock('B_-1')
        symbol=sp.Symbol('lambdaH1Inv1')
        k=(0,2,0,0)
        ch.by_name['lambdaH1Inv1'][k]=ch.by_name['lambdaH1Inv1'].get(k,0)+symbol
        out=portal_species_audit(ch,1)
        self.assertEqual(out['status'],'unresolved_generator_or_complex_scalar_basis')
        self.assertTrue(out['native_tensor_outside_rgbeta_span_components'])

    def test_missing_numeric_coupling_is_not_a_zero(self):
        ch=native_mock('A_-4')
        both={'lambdaH1':.01,'lambdaH2':.02}
        report,arr=construct_rgbeta_portal_background(ch,both,both)
        self.assertEqual(report['muS_missing_rgbeta_portals'],['lambdaH2Adj'])
        self.assertEqual(arr['muS_portal_labels'].tolist(),['lambdaH1','lambdaH2'])
        self.assertFalse(report['safe_for_complete_seven_term_or_joint_matching'])

    def test_coupling_value_types(self):
        with self.assertRaises(TensorAuditError): real_running_coupling(True,'lambdaH1')
        with self.assertRaises(TensorAuditError): real_running_coupling({'real':1.0,'imag':.1},'lambdaH1')

    def test_wrong_output_owner_rejected_without_accessing_sources(self):
        with TemporaryDirectory() as t:
            root=Path(t)
            with self.assertRaisesRegex(ValueError,'ownership'):
                run_portal_batch(root,root/'wrong')

    def test_no_overwrite_output(self):
        with TemporaryDirectory() as t:
            root=Path(t)
            old=root/'output/general_observable/lo_certification_couplings/preexisting'
            old.mkdir(parents=True)
            with self.assertRaises(FileExistsError):run_portal_batch(root,old)


if __name__=='__main__': unittest.main()

# Reuse the same exact-portals testing machinery for the complete numerical
# coordinate conversion rather than trusting a candidate alias assignment.
class TestNumericExactPortalMap(unittest.TestCase):
    def test_numeric_full_rank_map_nontrivial(self):
        from Numerical.lo_certification.couplings.portal_normalization import numerically_map_full_portal_subsector
        ch=native_mock('B_-1')
        orig=ch.by_name['lambdaH1Inv1']
        symbol=sp.Symbol('lambdaH1Inv1')
        ch.by_name['lambdaH1Inv1']={k:2*v for k,v in orig.items()}
        state={'lambdaH1':.4,'lambdaH1Adj':.6}
        out=numerically_map_full_portal_subsector(ch,1,state)
        self.assertEqual(out['status'],'exact_portal_only_numeric_conversion_verified')
        self.assertAlmostEqual(out['native_couplings']['lambdaH1Inv1'],.2)
        self.assertAlmostEqual(out['native_couplings']['lambdaH1Inv2'],.6)
        self.assertLess(out['normalized_real_HHSS_tensor_relative_residual'],1e-12)

    def test_numeric_incomplete_map_produces_no_fake_value(self):
        from Numerical.lo_certification.couplings.portal_normalization import numerically_map_full_portal_subsector
        ch=native_mock('A_-4',remove=('lambdaH2Inv2',))
        out=numerically_map_full_portal_subsector(ch,2,{'lambdaH2':.02,'lambdaH2Adj':.01})
        self.assertIsNone(out['native_couplings'])

    def test_numeric_missing_running_does_not_fallback_to_zero(self):
        from Numerical.lo_certification.couplings.portal_normalization import numerically_map_full_portal_subsector
        ch=native_mock('B_-1')
        out=numerically_map_full_portal_subsector(ch,1,{'lambdaH1':.01})
        self.assertEqual(out['status'],'blocked_missing_running_rgbeta_couplings')
        self.assertEqual(out['missing_numeric_couplings'],['lambdaH1Adj'])
