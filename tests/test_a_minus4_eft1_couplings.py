"""Batch 2: A(-4) RGBeta-frame, complete active-EFT1 tensor tests.

The tests do NOT certify the Matchete CG normalization. Source-level
operator coefficients are checked by equivalent fourth derivatives and
Ward contractions in the explicitly declared RGBeta reference frame.
"""
import json
import math
from pathlib import Path
from itertools import product
import unittest

import numpy as np
import sympy as sp

from Numerical.running import Aminus4EFT1Couplings as b


def fixture():
    r=np.arange(9,dtype=float).reshape(3,3)
    return dict(gY=.37,g2=.61,g3=1.05,
        yu=.01*(1+r)+.002j*(r-2),yd=.004*(2+r)-.001j*r,
        ye=.003*(3+r)+.001j*(r-4),lambdaH=.26,
        lambdaS1=.2,lambdaS2=.13,lambdaH1=.07,lambdaH2=.03,
        lambda12=.04,lambdaT3=.04+.025j,
        lambdaH2Adj=.017,lambdaS2Adj=-.012,
        mS1Sq=1.2e6,mS2Sq=1.7e6)


class AMinus4EFT1Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.x=b.tensor_point(1e5,fixture())
        cls.cols=b.quartic_columns()

    def test_complete_45_plus_12_and_active_fermion_only(self):
        self.assertEqual(len(self.x.fermion_basis),45)
        self.assertEqual(len(self.x.real_scalar_basis),12)
        self.assertTrue('L[3,2]' in self.x.fermion_basis)
        self.assertTrue('Q[3,3,2]' in self.x.fermion_basis)
        self.assertFalse(any(s.startswith('F[') for s in self.x.fermion_basis))
        self.assertEqual(self.x.Y.shape,(45,45,12))
        self.assertEqual(self.x.lambda4.shape,(12,12,12,12))
        np.testing.assert_allclose(self.x.Y,self.x.Y.transpose(1,0,2),atol=0,rtol=0)
        self.assertLess(np.max(np.abs(self.x.Y[:,:,4:])),1e-14)
        np.testing.assert_allclose(self.x.lambda4,self.x.lambda4.transpose(3,2,1,0),atol=2e-13,rtol=0)
        for name in ('lambdaH2Adj','lambdaS2Adj','lambdaT3'):
            self.assertGreater(np.max(np.abs(self.cols[name])),1e-5)

    def test_flavour_and_higgs_conjugation(self):
        p=fixture();Y=self.x.Y
        g,h,c=1,3,2
        i=b.FI[f'Q[{g},{c},1]'];j=b.FI[f'uc[{h},{c}]']
        np.testing.assert_allclose(Y[i,j,2],p['yu'][g-1,h-1]/math.sqrt(2))
        np.testing.assert_allclose(Y[i,j,3],1j*p['yu'][g-1,h-1]/math.sqrt(2))
        jj=b.FI[f'dc[{h},{c}]']
        np.testing.assert_allclose(Y[i,jj,1],-1j*p['yd'][g-1,h-1]/math.sqrt(2))
        ll=b.FI['L[2,2]'];ec=b.FI['ec[1]']
        np.testing.assert_allclose(Y[ll,ec,2],p['ye'][1,0]/math.sqrt(2))
        np.testing.assert_allclose(Y[ll,ec,3],-1j*p['ye'][1,0]/math.sqrt(2))

    def test_4th_derivative_combinatorics(self):
        # All checks are exact rational polynomials in selected real components.
        # Native monomials in RGBeta source: 1/2 lambda self; lambda portal.
        phi,eta=sp.symbols('phi eta',real=True)
        vH=sp.Rational(1,2)*(phi**2/sp.Integer(2))**2
        vH1=(phi**2/sp.Integer(2))*(eta**2/sp.Integer(2))
        self.assertAlmostEqual(float(sp.diff(vH,phi,4)),float(self.cols['lambdaH'][0,0,0,0]))
        self.assertAlmostEqual(float(sp.diff(vH1,phi,2,eta,2)),float(self.cols['lambdaH1'][0,0,4,4]))
        self.assertAlmostEqual(float(sp.diff(vH1,phi,2,eta,2)),float(self.cols['lambda12'][4,4,6,6]))
        # Complex coupling contributes its conjugate: derivative of HH*S1*S2*.
        z=b.quartic_tensor(fixture(),columns=self.cols)
        self.assertTrue(np.isfinite(z).all())
        self.assertGreater(abs(z[0,0,4,6]),1e-6)

    def test_invariant_linear_independence_10_real_directions(self):
        # 8 real invariants + Re/Im of one complex coupling.
        samples=[]
        for name in b.SCALAR_COUPLINGS:
            if name!='lambdaT3':samples.append(self.cols[name].ravel().real)
        samples.extend([2*self.cols['lambdaT3'].real.ravel(),
                        -2*self.cols['lambdaT3'].imag.ravel()])
        self.assertEqual(len(samples),10)
        self.assertEqual(np.linalg.matrix_rank(np.stack(samples),tol=1e-10),10)

    def test_nonabelian_gauge_algebra_and_casimir(self):
        g=self.x.gauge
        for a in g['SU2']['fermion']:
            np.testing.assert_allclose(a,a.conj().T,atol=1e-14)
        for a in g['SU3']['fermion']:
            np.testing.assert_allclose(a,a.conj().T,atol=1e-14)
        for rep in ('fermion','scalar'):
            x,y,z=g['SU2'][rep]
            np.testing.assert_allclose(x@y-y@x,1j*z,atol=1e-12)
        x,y,z=g['SU3']['fermion'][:3]
        np.testing.assert_allclose(x@y-y@x,1j*z,atol=1e-12)
        gen_cas=sum((t@t for t in g['SU3']['fermion']),
                    np.zeros((45,45),complex))
        self.assertAlmostEqual(gen_cas[b.FI['Q[2,2,1]'],b.FI['Q[2,2,1]']].real,4/3)
        self.assertAlmostEqual(gen_cas[b.FI['L[1,1]'],b.FI['L[1,1]']].real,0.)
        self.assertEqual(len(g['SU3']['scalar']),8)
        for a in g['SU3']['scalar']:
            self.assertLess(np.max(np.abs(a)),1e-14)

    def test_yukawa_and_quartic_ward_identities(self):
        yy=self.x.Y;ll=self.x.lambda4; g=self.x.gauge
        for group in ('U1','SU2','SU3'):
            for f,s in zip(g[group]['fermion'],g[group]['scalar']):
                a=np.einsum('ik,kja->ija',f.T,yy,optimize=True)
                a+=np.einsum('ika,kj->ija',yy,f,optimize=True)
                a+=np.einsum('ijb,ba->ija',yy,s,optimize=True)
                self.assertLess(np.max(np.abs(a)),3e-13,
                    msg=f'{group} Yukawa Ward failure')
                v=np.einsum('ea,ebcd->abcd',s,ll,optimize=True)
                v+=np.einsum('eb,aecd->abcd',s,ll,optimize=True)
                v+=np.einsum('ec,abed->abcd',s,ll,optimize=True)
                v+=np.einsum('ed,abce->abcd',s,ll,optimize=True)
                self.assertLess(np.max(np.abs(v)),3e-12,
                    msg=f'{group} quartic Ward failure')

    def test_zero_initial_coupling_is_still_in_inventory(self):
        p=fixture();p['lambdaH2Adj']=0.;p['lambdaS2Adj']=0.
        x=b.tensor_point(7e4,p)
        self.assertEqual(set(x.couplings),b.BETA_NAMES)
        with self.assertRaisesRegex(b.EFT1BasisError,'Missing'):
            p.pop('lambdaS2Adj');b.tensor_point(7e4,p)

    def test_basis_mapping_explicit(self):
        self.assertEqual(b.index_map(['L[1,1]','L[2,1]'],sector='fermion').shape,(2,))
        self.assertEqual(b.index_map(['H[1,R]','S2[3,I]'],sector='scalar').tolist(),[0,11])
        with self.assertRaisesRegex(b.EFT1BasisError,'Unknown'):
            b.index_map(['F[1,1]'],sector='fermion')
        with self.assertRaisesRegex(b.EFT1BasisError,'complete'):
            b.index_map(['L[1,1]'],sector='fermion',require_complete=True)

    def test_beta_contract_fails_closed(self):
        meta={'SharedScalar':False,'dS1':1,'dS2':3,'dF':2,'alpha':-4,
            'IntegratedField':'F','ActiveBSMFields':['S1','S2'],
            'ReportBetaConvention':b.BETA_CONVENTION}
        fixture={'status':'Success','metadata':meta,
                 'report_betas':{x:'1' for x in b.BETA_NAMES}}
        b.validate_beta_payload(fixture)
        for missing in ('ye','lambdaS2Adj','lambdaT3'):
            bad={'status':'Success','metadata':meta,
                 'report_betas':dict(fixture['report_betas'])}
            bad['report_betas'].pop(missing)
            with self.assertRaisesRegex(b.EFT1BasisError,'Incomplete'):
                b.validate_beta_payload(bad)
        bad={'status':'Success','metadata':meta,
             'report_betas':dict(fixture['report_betas'],y1='1')}
        with self.assertRaises(b.EFT1BasisError):b.validate_beta_payload(bad)

    def test_scalar_loop_native_beta_normalization_audit(self):
        expected = {'lambdaH':12.,'lambdaH1':2.,'lambdaH2':6.,
                    'lambdaH2Adj':3.,'lambdaT3':4.}
        self.assertEqual(b.field_content_gauge_beta_coefficients(),
                         {'gY':'55/6','g2':'-5/2','g3':'-7'})
        meta={'SharedScalar':False,'dS1':1,'dS2':3,'dF':2,'alpha':-4,
              'IntegratedField':'F','ActiveBSMFields':['S1','S2'],
              'ReportBetaConvention':b.BETA_CONVENTION}
        # Mirrors the literal expression read independently from native
        # output/full16_symbolic/.../eft1_rgbeta_rge.json; checkout test
        # below checks the actual raw file directly when available.
        expression=' + '.join(f'{int(v)}*{k}^2' for k,v in expected.items()
                             if k != 'lambdaT3')
        expression+=' + 4*lambdaT3*Bar[lambdaT3]'
        payload={'status':'Success','metadata':meta,
                 'report_betas':{k:'1' for k in b.BETA_NAMES}}
        payload['report_betas'].update({
            'gY':'(55*gY^3)/6','g2':'(-5*g2^3)/2',
            'g3':'-7*g3^3','lambdaH':expression})
        audit=b.source_normalization_audit(payload,columns=self.cols)
        self.assertEqual(audit['unresolved_normalization_mismatches'],{
            'lambdaH2Adj':{'reference':1.,'native':3.}})
        self.assertFalse(audit['scalar_basis_physically_certified'])
        self.assertEqual(audit['reference_coefficients']['lambdaT3'],4.)
        root=Path(__file__).resolve().parents[1]
        real=root/'output/full16_symbolic/T3_A_alpha_m4/data/eft1_rgbeta_rge.json'
        if real.exists():
            live=json.loads(real.read_text(encoding='utf-8-sig'))
            self.assertEqual(b.source_normalization_audit(live,columns=self.cols)
                             ['unresolved_normalization_mismatches'],
                             audit['unresolved_normalization_mismatches'])

    def test_entire_potential_equals_fourth_derivative_tensor(self):
        # Independent direct complex-field reconstruction, including triplet
        # adjoint portals, real SSSS and complex HHSS topology coupling.
        rng=np.random.default_rng(7159)
        su2h=b.su2_generators(2)
        su2s=b.su2_generators(3)
        cg=np.zeros((2,2,3),float)
        cg[0,0,0]=cg[1,1,2]=1
        cg[0,1,1]=cg[1,0,1]=1/math.sqrt(2)
        p=fixture()
        for _ in range(5):
            x=rng.normal(size=12)
            h=(x[0:4:2]+1j*x[1:4:2])/math.sqrt(2)
            s1=complex(x[4],x[5])/math.sqrt(2)
            s2=(x[6:12:2]+1j*x[7:12:2])/math.sqrt(2)
            hh=np.vdot(h,h).real
            a1=abs(s1)**2;a2=np.vdot(s2,s2).real
            hs=np.array([np.vdot(h,t@h) for t in su2h])
            ss=np.array([np.vdot(s2,t@s2) for t in su2s])
            top=sum(cg[i,j,m]*h[i]*h[j]*s1*s2[m].conjugate()
                    for i,j,m in product(range(2),range(2),range(3)))
            v=(p['lambdaH']*hh**2/2+p['lambdaS1']*a1*a1/2+
               p['lambdaS2']*a2*a2/2+p['lambdaS2Adj']*np.dot(ss,ss)/2+
               p['lambdaH1']*hh*a1+p['lambdaH2']*hh*a2+
               p['lambdaH2Adj']*np.dot(hs,ss)+p['lambda12']*a1*a2+
               2*np.real(p['lambdaT3']*top))
            from_tensor=np.einsum('abcd,a,b,c,d->',self.x.lambda4,x,x,x,x)/24
            self.assertAlmostEqual(float(v.real),float(from_tensor),places=12)

    def test_three_scale_tensor_reconstructions_without_claiming_running(self):
        # This checks the interpolation-ready tensor interface at three
        # distinct external coupling points. Actual beta evolution is
        # exercised on the installed Windows repository via the CLI.
        for scale in (1.e5,1.e4,1.e3):
            p=fixture()
            p['lambdaH2Adj']*=1.+0.15*math.log(1.e5/scale)
            p['lambdaT3']*=1.+0.04*math.log(1.e5/scale)
            p['ye']*=1.+0.08*math.log(1.e5/scale)
            t=b.tensor_point(scale,p,columns=self.cols)
            self.assertTrue(np.isfinite(t.Y).all())
            self.assertTrue(np.isfinite(t.lambda4).all())
            self.assertLess(np.max(abs(t.Y-t.Y.transpose(1,0,2))),1e-13)
            self.assertLess(np.max(abs(t.lambda4-t.lambda4.transpose(2,3,0,1))),1e-12)

    def test_nested_native_complex_three_generation_matrix(self):
        payload=[[{'real':float(i==j),'imag':float(i-j)/100}
                  for j in range(3)] for i in range(3)]
        expected=np.eye(3)+1j*(np.arange(3)[:,None]-np.arange(3)[None,:])/100
        np.testing.assert_array_equal(b._complex_array(payload,(3,3),'yu'),expected)
        payload[1][0].pop('imag')
        with self.assertRaises(b.EFT1BasisError):
            b._complex_array(payload,(3,3),'yu')

    def test_native_seed_inventory_does_not_assert_basis_alias(self):
        seed={'IntegratedFields':['F'],'ActiveHeavyFields':['S1','S2'],
              'ScalarQuarticTermCount':2,'ScalarQuarticTerms':[
                {'TermInputForm':'Coupling[lambdaH2Inv2, {}, 0]'},
                {'TermInputForm':'Coupling[lambdaS2Inv1, {}, 0]'}]}
        audit=b.audit_native_seed(seed)
        self.assertEqual(audit['native_coupling_symbols'],
                         ['lambdaH2Inv2','lambdaS2Inv1'])
        self.assertFalse(audit['native_source_tensor_map_verified'])

    def test_source_contract_against_checkout_when_present(self):
        root=Path(__file__).resolve().parents[1]
        if not (root/'Lagrangian/model/T3ModelCatalog.wl').exists():
            self.skipTest('Exact local executable WL sources not included in delivery ZIP')
        report=b.source_contract(root)
        self.assertEqual(report['quartic_extra_invariants'],0)
        self.assertEqual(report['additional_active_nonSM_yukawas'],0)


if __name__=='__main__':
    unittest.main()
