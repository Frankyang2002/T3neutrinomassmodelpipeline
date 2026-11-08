"""Independent basis-covariance and ODE checks for the *unverified* candidate.

Regression goal is mathematical consistency, not proof of loop coefficient.
All tests use artificial tensors and no physical matching inputs.
"""
from __future__ import annotations

from itertools import permutations
import unittest
import numpy as np

from RGE.running.intermediate.GeneralEFT1SevenTermFlow import (
    Couplings, TERMS, TensorBasis, SymmetricCoordinates, evolve, master_beta,
    symmetry_residual, total_beta, gauge_ward_residual,
)
from RGE.running.intermediate.GeneralEFT1CrossedWardCandidate import (
    candidate_beta, crossed_candidate, fermion_casimir_transpose_delta,
)


def symmetric_instance(seed=11):
    rng=np.random.default_rng(seed)
    n=4; s=3
    y=rng.normal(size=(n,n,s))+1j*rng.normal(size=(n,n,s))
    y=.5*(y+y.swapaxes(0,1))
    lam=rng.normal(size=(s,)*4)
    lam=sum(lam.transpose(p) for p in permutations(range(4)))/24
    c=rng.normal(size=(n,n,s,s))+1j*rng.normal(size=(n,n,s,s))
    c=.25*(c+c.swapaxes(0,1)+c.swapaxes(2,3)+c.swapaxes(0,1).swapaxes(2,3))
    tf=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n))
    tf=.5*(tf+tf.conj().T)
    theta=np.zeros((s,s),complex)
    theta[0,1]=1j; theta[1,0]=-1j
    p=Couplings(y,lam,((.34,((tf,theta),)),))
    u,_=np.linalg.qr(rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)))
    r,_=np.linalg.qr(rng.normal(size=(s,s)))
    def transform_c(z):
        return np.einsum('ip,jq,ar,bs,pqrs->ijab',u.conj(),u.conj(),r,r,z,optimize=True)
    new_y=np.einsum('ip,jq,ar,pqr->ija',u.conj(),u.conj(),r,y,optimize=True)
    new_lam=np.einsum('ap,bq,cr,ds,pqrs->abcd',r,r,r,r,lam,optimize=True)
    new_p=Couplings(new_y,new_lam,((.34,((u@tf@u.conj().T,r@theta@r.T),)),))
    return c,p,transform_c(c),new_p,transform_c


def u1_ode_instance():
    """Fully gauge-invariant U1 input embedded in 6 Weyl x 4 real coordinates."""
    n=6; ns=4
    v=np.array([1,1j],complex)/np.sqrt(2)
    t=np.diag([0,0,-1,0,0,0]).astype(complex)
    theta=np.zeros((ns,ns),complex)
    theta[0,1]=1j;theta[1,0]=-1j
    y=np.zeros((n,n,ns),complex)
    y[0,2,:2]=y[2,0,:2]=.35*v
    c=np.zeros((n,n,ns,ns),complex)
    c[2,2,:2,:2]=(.2+.1j)*np.outer(v,v)
    zero=np.zeros((ns,)*4)
    low=Couplings(.9*y,zero,((.27,((t,theta),)),))
    high=Couplings(y,zero,((.30,((t,theta),)),))
    basis=TensorBasis(n,ns,((0,1),(2,3),(4,5)),(0,1,2,3))
    return basis,c,low,high,t,theta


class CandidateCovarianceTests(unittest.TestCase):
    def test_global_weyl_rephasing_proves_conjugate_index_necessity(self):
        c,p,_,_,_=symmetric_instance()
        phase=np.exp(1j*.31)
        # For psi' = phase psi: C' = phase**(-2) C and y' likewise.
        cp=c*phase**(-2)
        pp=Couplings(p.yukawa*phase**(-2),p.quartic,p.gauge)
        good=crossed_candidate(c,p)
        np.testing.assert_allclose(crossed_candidate(cp,pp),phase**(-2)*good,
                                   atol=1e-11,rtol=1e-11)
        orig=master_beta(c,p)[TERMS[3]]
        orig_new=master_beta(cp,pp)[TERMS[3]]
        self.assertGreater(np.linalg.norm(orig_new-phase**(-2)*orig)/np.linalg.norm(orig),.2)

    def test_arbitrary_complex_unitary_and_real_scalar_rotation(self):
        c,p,cp,pp,rot=symmetric_instance()
        before=candidate_beta(c,p)
        after=candidate_beta(cp,pp)
        for name in TERMS:
            np.testing.assert_allclose(after[name],rot(before[name]),rtol=8e-12,atol=4e-11,
                                       err_msg=name)
        self.assertLess(symmetry_residual(total_beta(after)),1e-12)

    def test_transposed_gauge_casimir_is_not_optional_outside_canonical_basis(self):
        c,p,cp,pp,rot=symmetric_instance()
        old=master_beta(c,p)[TERMS[6]]
        old_rot=master_beta(cp,pp)[TERMS[6]]
        relative=np.linalg.norm(old_rot-rot(old))/np.linalg.norm(rot(old))
        self.assertGreater(relative,1e-2)
        self.assertGreater(np.linalg.norm(fermion_casimir_transpose_delta(cp,pp)),.01)
        fixed=candidate_beta(cp,pp)[TERMS[6]]
        np.testing.assert_allclose(fixed,rot(candidate_beta(c,p)[TERMS[6]]),atol=4e-11)

    def test_candidate_anti_linear_terms_and_real_linearity(self):
        c,p,_,_,_=symmetric_instance()
        terms=candidate_beta(c,p)
        shifted=candidate_beta(1j*c,p)
        for name in TERMS:
            factor=-1j if name in (TERMS[3],TERMS[4]) else 1j
            np.testing.assert_allclose(shifted[name],factor*terms[name],atol=2e-11,err_msg=name)
        d=c.swapaxes(0,1)*.17+c.conj()*.2
        for label,f in [('real_linearity',lambda z:candidate_beta(z,p))]:
            z1=total_beta(f(.37*c-.29*d))
            z2=.37*total_beta(f(c))-.29*total_beta(f(d))
            np.testing.assert_allclose(z1,z2,atol=3e-11,err_msg=label)

    def test_rge_candidate_integrates_generated_zero_modes_with_gauge_ward(self):
        basis,c,low,high,t,theta=u1_ode_instance()
        self.assertLess(gauge_ward_residual(c,t,theta),2e-13)
        for sample in (low,high):
            self.assertLess(gauge_ward_residual(total_beta(candidate_beta(c,sample)),t,theta),
                            2e-12)
        out=evolve(basis,c,1e5,1e3,[1e3,1e5],[low,high],
                   points=7,rtol=1e-10,atol=1e-14,beta_mode='ward_covariant_candidate')
        self.assertEqual(out['anti_linear_term_names'],[TERMS[3],TERMS[4]])
        self.assertFalse(out['beta_physics_certified'])
        self.assertGreater(len(out['generated_indices']),0)
        self.assertGreater(np.linalg.norm(out['final'][0,0,:2,:2]),1e-4)
        coords=SymmetricCoordinates(basis)
        for row in out['trajectory']:
            self.assertLess(gauge_ward_residual(coords.unpack(row),t,theta),2e-10)
        beta0=total_beta(out['beta_initial'])
        delta=(np.log(1e3/1e5)/(16*np.pi**2))*beta0
        # Test loop sign and consistency without incorrectly imposing fixed-order equality.
        self.assertGreater(abs(delta[0,0,0,0]),0)
        self.assertLess(abs(out['final'][0,0,0,0]-delta[0,0,0,0])/abs(delta[0,0,0,0]),.12)

    def test_invalid_beta_mode_fails_closed(self):
        basis,c,low,high,_,_=u1_ode_instance()
        with self.assertRaisesRegex(ValueError,'Unknown beta_mode'):
            evolve(basis,c,1e5,1e3,[1e3,1e5],[low,high],beta_mode='auto_physical')


if __name__=='__main__':
    unittest.main()
