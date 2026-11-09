"""No-factorial psi psi phi phi coordinates and Weinberg embedding.

L5 = C_ijab psi_i psi_j phi_a phi_b + h.c. (separate pair symmetry)
L5^W = (1/2) C5_pq (L_p dot H)(L_q dot H) + h.c.
=> C[nu_p,nu_q,H[2,R],H[2,R]] = C5[p,q]/4.
"""
from __future__ import annotations
import numpy as np
from .basis import ModelBasis


def symmetrize(c: np.ndarray) -> np.ndarray:
    c = np.asarray(c,dtype=np.complex128)
    return (c+c.swapaxes(0,1)+c.swapaxes(2,3)+c.transpose(1,0,3,2))/4


def symmetry_residual(c: np.ndarray) -> float:
    return float(max(np.max(np.abs(c-c.swapaxes(0,1))),
                     np.max(np.abs(c-c.swapaxes(2,3)))))


def embed_c5(c5: np.ndarray, basis: ModelBasis) -> np.ndarray:
    """Embed the *whole* gauge invariant Weinberg polynomial, not only nu-nu-H0-H0."""
    c5 = np.asarray(c5,np.complex128)
    if c5.shape != (3,3) or not np.isfinite(c5).all() or not np.allclose(c5,c5.T,atol=1e-23,rtol=1e-12):
        raise ValueError('C5 must be finite complex symmetric 3x3')
    out = np.zeros((basis.nf,basis.nf,basis.ns,basis.ns),complex)
    # H^j = sum_a u[j,a] real_phi[a]; U only for first four H coordinates.
    U = np.array([[1,1j,0,0],[0,0,1,1j]],complex)/np.sqrt(2)
    eps = np.array([[0,1],[-1,0]],int)
    for p in range(3):
        for q in range(3):
            for i in range(2):
                for k in range(2):
                    f,g=basis.leptons[p][i],basis.leptons[q][k]
                    x = np.einsum('j,l,ja,lb->ab',eps[i],eps[k],U,U,optimize=True)
                    # (1/2) prefactor of Weinberg and H=(R+iI)/sqrt(2).
                    out[f,g,:4,:4] += 0.5*c5[p,q]*x
    return symmetrize(out)


def project_c5(c: np.ndarray, basis: ModelBasis) -> tuple[np.ndarray,float]:
    """Extract C5=4*C[nu,nu,h0R,h0R]; residual of full LLHH polynomial."""
    c=np.asarray(c,np.complex128)
    if c.shape != (basis.nf,basis.nf,basis.ns,basis.ns):
        raise ValueError('Incorrect Wilson tensor geometry')
    nu=[pair[0] for pair in basis.leptons]
    h = basis.higgs[2]
    proj=4*c[np.ix_(nu,nu,[h],[h])][:,:,0,0]
    proj=(proj+proj.T)/2
    reconstructed=embed_c5(proj,basis)
    inds = [i for pair in basis.leptons for i in pair]
    a=c[np.ix_(inds,inds,range(4),range(4))]
    b=reconstructed[np.ix_(inds,inds,range(4),range(4))]
    residual=float(np.max(np.abs(a-b)))
    return proj,residual


def embed_formal_llss(raw: np.ndarray, basis: ModelBasis,
                      source_from_canonical: np.ndarray, signed_factor: complex) -> np.ndarray:
    """raw[p,q,weak_i,weak_j,heavy_source_a,heavy_source_b].

    phi_source[a] = sum_A R[a,A] phi_canonical_heavy[A]. The spinor
    factor is an EXTERNAL unproved assertion; never inferred here.
    """
    raw=np.asarray(raw,np.complex128)
    nh=basis.ns-4
    R=np.asarray(source_from_canonical,float)
    if raw.shape != (3,3,2,2,nh,nh):
        raise ValueError(f'Formal LLSS has unexpected axes: {raw.shape}')
    if R.shape != (nh,nh) or not np.isfinite(R).all() or not np.allclose(R@R.T,np.eye(nh),atol=1e-10):
        raise ValueError('Require explicit real orthogonal heavy-source-from-canonical map')
    if not np.isfinite(raw).all() or not np.isfinite(signed_factor):
        raise ValueError('Nonfinite formal tensor/spinor factor')
    if not np.allclose(raw,raw.transpose(1,0,3,2,4,5),atol=1e-22,rtol=1e-10):
        raise ValueError('Formal lepton-field interchange fails')
    if not np.allclose(raw,raw.swapaxes(-1,-2),atol=1e-22,rtol=1e-10):
        raise ValueError('Formal scalar exchange fails')
    transformed=np.einsum('pqijab,aA,bB->pqijAB',raw,R,R,optimize=True)
    out=np.zeros((basis.nf,basis.nf,basis.ns,basis.ns),complex)
    for p in range(3):
        for q in range(3):
            for i in range(2):
                for j in range(2):
                    f,g=basis.leptons[p][i],basis.leptons[q][j]
                    out[f,g,4:,4:] += signed_factor*transformed[p,q,i,j]
    return symmetrize(out)
