"""Exact matrix test for Majorana Takagi redefinition; not cross-package certification."""
import numpy as np

def check(seed=137):
    rng=np.random.default_rng(seed)
    z=rng.normal(size=(3,3))+1j*rng.normal(size=(3,3))
    u,_=np.linalg.qr(z)
    pL=np.diag([1,1,0,0]);pR=np.eye(4)-pL
    # Flavor acts on the 3-index, chirality on the 4-index.
    transform=np.kron(u,pL)+np.kron(u.conj(),pR)
    inverse=np.kron(u.conj().T,pL)+np.kron(u.T,pR)
    unitary=float(np.linalg.norm(inverse@transform-np.eye(12)))
    assert unitary<1e-12
    # A Weyl-basis charge-conjugation operator interchanges L and R.
    flip=np.array([[0,0,1,0],[0,0,0,1],[1,0,0,0],[0,1,0,0]])
    charge_flip=float(np.linalg.norm(np.kron(np.eye(3),flip)@transform.conj()-transform@np.kron(np.eye(3),flip)))
    assert charge_flip<1e-12
    return {'status':'ExactChiralMatrixIdentityPassed','seed':seed,'inverse_residual':unitary,'charge_conjugation_flip_residual':charge_flip,'native_matchete_dictionary_certified':False}
