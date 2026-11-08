"""Batch 3--4 independent numerical tests, including source regression diagnostics.

Run: python -m unittest tests.test_general_eft1_seven_term_flow -v
No physical-A(-4) claims: fixtures deliberately violate gauge selection rules.
"""
from __future__ import annotations

from itertools import permutations
from math import log
import unittest

import numpy as np

from RGE.running.intermediate.GeneralEFT1SevenTermFlow import (
    Couplings, TensorBasis, SymmetricCoordinates, TERMS, LOOP, c5_projection,
    collinear, evolve, master_beta, symmetry_residual, total_beta,
    gauge_ward_residual,
)
from tests.fixtures.Batch34SyntheticAminus4 import fixtures


def raw_array(x):
    if isinstance(x, dict) and set(x) == {'real','imag'}:
        return np.asarray(x['real'],float)+1j*np.asarray(x['imag'],float)
    return np.asarray(x,complex)


def orbit(tensor, indices, val, kind):
    ij = [q-1 for q in indices]
    if kind == 'quartic':
        for key in set(permutations(ij)):
            tensor[key] = val
    elif kind == 'yukawa':
        i,j,a = ij
        tensor[i,j,a] = tensor[j,i,a] = val
    else:
        i,j,a,b = ij
        for fi,fj in ((i,j),(j,i)):
            for sa,sb in ((a,b),(b,a)):
                tensor[fi,fj,sa,sb] = val


def synthetic():
    b, cp = fixtures()
    nf,ns = b['n_fermions'],b['n_real_scalars']
    basis = TensorBasis(nf,ns,tuple(tuple(i-1 for i in x) for x in b['lepton_components']),
                        tuple(i-1 for i in b['higgs_real_components']))
    c = np.zeros((nf,nf,ns,ns),complex)
    for op in b['operators']:
        part=np.zeros_like(c)
        for row in op['components']:
            orbit(part,row['indices'],complex(row['value']['real'],row['value']['imag']),'wilson')
        c += part
    samples=[]
    for row in cp['samples']:
        y=np.zeros((nf,nf,ns),complex); lam=np.zeros((ns,)*4,complex)
        for line in row['yukawa_components']:
            orbit(y,line['indices'],complex(line['value']['real'],line['value']['imag']),'yukawa')
        for line in row['quartic_components']:
            orbit(lam,line['indices'],line['value'],'quartic')
        gg=[]
        for entry in cp['gauge_generators']:
            pairs=tuple((raw_array(t['fermion']),raw_array(t['scalar'])) for t in entry['generators'])
            gg.append((float(row['gauge_couplings'][entry['name']]),pairs))
        samples.append(Couplings(y,lam,tuple(gg)))
    return basis,c,[x['mu_gev'] for x in cp['samples']],samples


def small_case():
    rng=np.random.default_rng(551)
    nf,ns=2,2
    y=rng.normal(size=(nf,nf,ns))+1j*rng.normal(size=(nf,nf,ns))
    y=(y+y.swapaxes(0,1))/2
    lam=rng.normal(size=(ns,)*4)
    lam=sum(np.transpose(lam,p) for p in set(permutations(range(4))))/24
    c=rng.normal(size=(nf,nf,ns,ns))+1j*rng.normal(size=(nf,nf,ns,ns))
    c=0.25*(c+c.swapaxes(0,1)+c.swapaxes(2,3)+c.swapaxes(0,1).swapaxes(2,3))
    t=np.array([[0.4,0.2j],[-0.2j,-0.1]])
    th=np.array([[0,0.3j],[-0.3j,0]])
    return c,Couplings(y,lam,((0.31,((t,th),)),))


def slow_terms(c,p):
    """Independently indexed Eq.(4.85), with FOUR crossed exchanges."""
    y,lam=p.yukawa,p.quartic
    n,_,s=y.shape
    gs,gf=collinear(p)
    res={k:np.zeros_like(c) for k in TERMS}
    for i in range(n):
      for j in range(n):
       for a in range(s):
        for b in range(s):
         for k in range(n):
          for l in range(n):
           for d in range(s):
            res[TERMS[0]][i,j,a,b] += y[i,j,d]*y[k,l,d].conjugate()*c[k,l,a,b]
            res[TERMS[4]][i,j,a,b] += (2*y[i,l,d]*y[j,k,d]
                                            +y[i,j,d]*y[k,l,d])*c[k,l,a,b].conjugate()
         for aa,bb in ((a,b),(b,a)):
          for d in range(s):
           res[TERMS[5]][i,j,a,b] += gs[bb,d]*c[i,j,aa,d]
          for cc in range(s):
           for dd in range(s):
            gauge=0j
            for g,pairs in p.gauge:
             for _,theta in pairs:
              gauge += 2*g*g*theta[aa,cc]*theta[bb,dd]
            res[TERMS[1]][i,j,a,b] -= 0.5*(gauge-lam[aa,bb,cc,dd])*c[i,j,cc,dd]
         for ii,jj in ((i,j),(j,i)):
          for k in range(n):
           res[TERMS[6]][i,j,a,b] += gf[jj,k]*c[ii,k,a,b]
          for aa,bb in ((a,b),(b,a)):
           for l in range(n):
            for cc in range(s):
             gauge=0j
             for g,pairs in p.gauge:
              for tf,theta in pairs:
               gauge += 4*g*g*tf[l,jj]*theta[bb,cc]
             res[TERMS[2]][i,j,a,b] += gauge*c[ii,l,aa,cc]
             for k in range(n):
              res[TERMS[2]][i,j,a,b] += (
                  2*y[jj,k,cc]*y[k,l,bb].conjugate()
                  +y[jj,k,bb]*y[k,l,cc].conjugate())*c[ii,l,aa,cc]
              res[TERMS[3]][i,j,a,b] += y[jj,l,aa]*y[ii,k,cc]*c[l,k,bb,cc]
    return res


class ContractionTests(unittest.TestCase):
    def test_all_seven_independent_index_loops(self):
        c,p=small_case()
        fast=master_beta(c,p)
        slow=slow_terms(c,p)
        for name in TERMS:
            np.testing.assert_allclose(fast[name],slow[name],rtol=1e-11,atol=2e-12,err_msg=name)

    def test_crossed_source_two_exchange_bug_detected(self):
        c,p=small_case()
        wrong=master_beta(c,p,corrected_crossed=False)[TERMS[3]]
        fixed=master_beta(c,p)[TERMS[3]]
        self.assertGreater(np.linalg.norm(wrong-fixed),1e-6)
        self.assertGreater(symmetry_residual(wrong),1e-5)
        self.assertLess(symmetry_residual(fixed),1e-13)

    def test_individual_terms_pair_symmetry(self):
        c,p=small_case()
        for name,b in master_beta(c,p).items():
            self.assertLess(symmetry_residual(b),2e-13,name)

    def test_real_linearity_and_anti_linearity(self):
        c,p=small_case()
        d=np.roll(c,1,axis=0)
        d=0.5*(d+d.swapaxes(0,1))
        lhs=total_beta(master_beta(0.6*c-0.37*d,p))
        rhs=0.6*total_beta(master_beta(c,p))-0.37*total_beta(master_beta(d,p))
        np.testing.assert_allclose(lhs,rhs,rtol=5e-12,atol=5e-12)
        cb=master_beta(1j*c,p)
        ab=master_beta(c,p)
        for name in TERMS:
            if name != TERMS[4]:
                np.testing.assert_allclose(cb[name],1j*ab[name],atol=4e-11)
        np.testing.assert_allclose(cb[TERMS[4]],-1j*ab[TERMS[4]],atol=4e-11)
        self.assertGreater(np.linalg.norm(cb[TERMS[4]]-1j*ab[TERMS[4]]),1e-5)

    def test_real_orthogonal_basis_covariance(self):
        c,p=small_case(); rng=np.random.default_rng(211)
        uf,_=np.linalg.qr(rng.normal(size=(2,2)))
        rr,_=np.linalg.qr(rng.normal(size=(2,2)))
        def rot_c(z):return np.einsum('ip,jq,ar,bs,pqrs->ijab',uf,uf,rr,rr,z)
        new_y=np.einsum('ip,jq,ar,pqr->ija',uf,uf,rr,p.yukawa)
        new_lam=np.einsum('ap,bq,cr,ds,pqrs->abcd',rr,rr,rr,rr,p.quartic)
        gg=tuple((g,tuple((uf@tf@uf.T,rr@ts@rr.T) for tf,ts in pairs))
                 for g,pairs in p.gauge)
        new_p=Couplings(new_y,new_lam,gg)
        old=master_beta(c,p)
        new=master_beta(rot_c(c),new_p)
        for name in TERMS:
            np.testing.assert_allclose(new[name],rot_c(old[name]),rtol=1e-11,atol=2e-11,
                                       err_msg=name)

    def test_su2_u1_weinberg_gauge_ward_preserved(self):
        sigma=[np.array([[0,1],[1,0]],complex)/2,
               np.array([[0,-1j],[1j,0]],complex)/2,
               np.diag([.5,-.5]).astype(complex)]
        def real_scalar_generator(t):
            theta=np.zeros((4,4),complex)
            for p in range(2):
                for q in range(2):
                    z=t[p,q]
                    theta[2*p,2*q] += 1j*z.imag
                    theta[2*p,2*q+1] += 1j*z.real
                    theta[2*p+1,2*q] -= 1j*z.real
                    theta[2*p+1,2*q+1] += 1j*z.imag
            return theta
        vnu=np.array([0,0,1,1j])/np.sqrt(2)
        ve=-np.array([1,1j,0,0])/np.sqrt(2)
        v=[vnu,ve]
        c=np.empty((2,2,4,4),complex)
        for i in range(2):
            for j in range(2):
                c[i,j]=.25*(1+.13j)*(np.outer(v[i],v[j])+np.outer(v[j],v[i]))
        g2=((.65,tuple((t,real_scalar_generator(t)) for t in sigma)))
        hyper_f=-.5*np.eye(2)
        hyper_s=real_scalar_generator(.5*np.eye(2))
        g1=(.36,((hyper_f,hyper_s),))
        p=Couplings(np.zeros((2,2,4),complex),np.zeros((4,)*4), (g2,g1))
        terms=master_beta(c,p)
        for name,coeff in [('C',c)]+list(terms.items()):
            for _,gens in p.gauge:
                for t,theta in gens:
                    self.assertLess(gauge_ward_residual(coeff,t,theta),5e-13,name)

    def test_original_symbolic_master_if_repo_installed(self):
        try:
            import sympy as sp
            from RGE.general.ScalarBasis import ScalarBasis,ComplexScalar
            from RGE.general.WilsonTensorRGE import WilsonRGEInputs,calculate_wilson_tensor_rge
            from RGE.general.GaugeGenerators import GaugeSector
        except ImportError:
            self.skipTest('Current upstream symbolic source is not mounted in this runtime')
        c,p=small_case()
        model=ScalarBasis((ComplexScalar('toy',1,sp.Integer(0)),))
        def sx(v):return sp.Float(complex(v).real)+sp.I*sp.Float(complex(v).imag)
        class Coeff:
            def __getitem__(self,k):
                i,j,a,b=k
                return sx(c[i-1,j-1,a-1,b-1])
        t,th=p.gauge[0][1][0]
        sectors=(GaugeSector(sp.Float(p.gauge[0][0]),(sp.Matrix(th),),(sp.Matrix(t),)),)
        gs,gf=collinear(p)
        ins=WilsonRGEInputs(2,lambda i,j,a:sx(p.yukawa[i-1,j-1,a-1]),
                            lambda a,b,d,e:sx(p.quartic[a-1,b-1,d-1,e-1]),
                            gauge_sectors=sectors,
                            gamma_scalar=lambda a,b:sx(gs[a-1,b-1]),
                            gamma_fermion=lambda i,j:sx(gf[i-1,j-1]))
        symbolic=calculate_wilson_tensor_rge(model,ins,(1,2,1,2),Coeff())
        upstream=master_beta(c,p,corrected_crossed=False)
        for name in TERMS:
            self.assertAlmostEqual(abs(complex(symbolic[name])-upstream[name][0,1,0,1]),0,places=9)


class FlowTests(unittest.TestCase):
    def test_full_three_generation_basis_includes_all_SM_weyl(self):
        basis=TensorBasis(45,12,((0,1),(2,3),(4,5)),(0,1,2,3))
        coords=SymmetricCoordinates(basis)
        self.assertEqual(coords.size,80730)
        self.assertEqual(coords.labels().shape,(80730,4))
        self.assertTrue(np.any(coords.labels()[:,0]==45))
        self.assertTrue(np.any(coords.labels()[:,2]==12))

    def test_matrix_free_generated_fermions_beyond_lepton_sector(self):
        nf,ns=8,12
        y=np.zeros((nf,nf,ns),complex)
        y[0,0,4]=0.19+0.08j; y[7,7,4]=0.23+0.07j
        lam=np.zeros((ns,)*4)
        c=np.zeros((nf,nf,ns,ns),complex)
        orbit(c,[1,1,5,7],1e-10+2e-11j,'wilson')
        b=master_beta(c,Couplings(y,lam,()))
        self.assertGreater(abs(b[TERMS[0]][7,7,4,6]),0)
        self.assertEqual(c[7,7,4,6],0)

    def test_generated_modes_and_full_trajectory(self):
        basis,c,scales,p=synthetic()
        out=evolve(basis,c,1e5,2e3,scales,p,points=4,rtol=1e-8,atol=1e-20)
        coords=SymmetricCoordinates(basis)
        self.assertEqual(out['trajectory'].shape,(4,coords.size))
        self.assertEqual(coords.size,21*78)
        self.assertGreater(len(out['generated_indices']),0)
        self.assertLess(max(out['symmetry_residuals']),3e-12)
        self.assertEqual(len(out['projections']),4)
        self.assertEqual(len(out['beta_final']),7)
        for index in out['generated_indices']:
            self.assertEqual(coords.pack(c)[index],0)
        for name in TERMS:
            self.assertTrue(np.isfinite(out['beta_final'][name]).all())

    def test_fixed_order_limit_and_threshold_direction(self):
        basis,c,scales,samples=synthetic()
        hi=1e5
        out=evolve(basis,c,hi,hi/1.01,scales,samples,points=2,atol=1e-20)
        self.assertLess(out['scales_gev'][-1],out['scales_gev'][0])
        delta=out['final']-c
        fixed=LOOP*log(1/1.01)*total_beta(out['beta_initial'])
        rel=np.linalg.norm(delta-fixed)/np.linalg.norm(fixed)
        self.assertLess(rel,2e-3)
        self.assertGreater(np.linalg.norm(delta),0)

    def test_source_superposition_at_each_scale(self):
        basis,c,_,samples=synthetic()
        mid=samples[1]
        a=c.copy();a[:]=0
        b=c.copy();b[:]=0
        # partition by scalar support, preserving both pair exchanges
        a[:,:,:, :]=np.where(np.abs(c[:,:,:,6])[:, :, :, None]>0,c,0)  # A LLSS block
        b=c-a
        for name in TERMS:
            out=master_beta(c,mid)[name]
            np.testing.assert_allclose(out,master_beta(a,mid)[name]+master_beta(b,mid)[name],
                                       atol=1e-24,rtol=2e-11)

    def test_sole_scalar_pair_direct_projection_in_gaugeless_synthetic_limit(self):
        basis,c,_,samples=synthetic()
        p=samples[-1]
        p=Couplings(np.zeros_like(p.yukawa),p.quartic,())
        terms=master_beta(c,p)
        for name in TERMS:
            projection=c5_projection(terms[name],basis)
            if name==TERMS[1]:self.assertGreater(np.linalg.norm(projection),0)
            else:self.assertLess(np.linalg.norm(projection),1e-23,name)

    def test_ode_convergence_and_reversal_sign(self):
        basis,c,scales,samples=synthetic()
        coarse=evolve(basis,c,1e5,2e3,scales,samples,points=2,
                      rtol=1e-6,atol=1e-18)
        fine=evolve(basis,c,1e5,2e3,scales,samples,points=2,
                    rtol=1e-9,atol=1e-21)
        self.assertLess(np.linalg.norm(coarse['final']-fine['final'])/
                        np.linalg.norm(fine['final']),1e-6)
        self.assertGreater(np.linalg.norm(fine['final']-c),0)
        self.assertLess(np.max(np.abs(fine['scales_gev'][-1]-2e3)),1e-8)

    def test_installed_original_fixed_log_and_evolved_short_log(self):
        try:
            from RGE.running.intermediate.GeneralEFT1C5Mixing import solve_boundary,load_boundary
            from RGE.running.intermediate.GeneralEFT1SevenTermFlow import run_contracts
        except ImportError:
            self.skipTest('Current upstream GeneralEFT1C5Mixing source is not mounted')
        boundary,cp=fixtures()
        # Force the *same* mu_F running tensors into the old frozen-coupling
        # fixed-order source, so the comparison does not silently use other inputs.
        high=cp['samples'][-1]
        boundary['yukawa_components']=high['yukawa_components']
        boundary['quartic_components']=high['quartic_components']
        boundary['gauge_sectors']=[{'name':row['name'],
             'coupling':high['gauge_couplings'][row['name']],
             'generators':row['generators']} for row in cp['gauge_generators']]
        upstream=solve_boundary(boundary,allow_synthetic=True)
        out=run_contracts(boundary,cp,allow_synthetic=True,
                          points=3,rtol=1e-9,atol=1e-21,source_decomposition=True)
        # Exact fixed-order scalar-pair LLHH coefficient agrees. Other terms
        # are permitted to disagree where the crossed-exchange source defect acts.
        ref=raw_array(upstream['summed_by_rge_term_beta_C5_16pi2'][TERMS[1]])
        actual=c5_projection(out['beta_initial'][TERMS[1]],out['basis'])
        np.testing.assert_allclose(actual,ref,rtol=1e-12,atol=1e-24)
        for name in TERMS:
            refval=raw_array(upstream['summed_by_rge_term_beta_C5_16pi2'][name])
            original=raw_array(out['upstream_reference'][name]['direct_LLHH_projection_upstream'])
            np.testing.assert_allclose(refval,original,rtol=1e-12,atol=1e-24)
        self.assertLess(out['source_recomposition_residual'],1e-7)
        self.assertEqual(len(out['source_runs']),2)
        fixed_delta=LOOP*log(2000./1e5)*c5_projection(total_beta(out['beta_initial']),out['basis'])
        np.testing.assert_allclose(c5_projection(out['fixed_order_full']-out['initial'],out['basis']),
                                   fixed_delta,rtol=2e-12,atol=1e-26)

    def test_missing_scale_range_and_non_symmetric_boundary_fail(self):
        basis,c,scales,p=synthetic()
        with self.assertRaisesRegex(ValueError,'bracket'):
            evolve(basis,c,1e5,2e3,[3000,1e5],p[:2],points=2)
        corrupted=c.copy();corrupted[0,2,4,6]+=1e-10
        with self.assertRaisesRegex(ValueError,'symmetries'):
            evolve(basis,corrupted,1e5,2e3,scales,p,points=2)


if __name__=='__main__':
    unittest.main()
