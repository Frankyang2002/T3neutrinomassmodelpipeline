"""Independent, broken-phase, p=0 one-loop Weyl-fermion/real-scalar self-energy.

All vertices MUST be supplied in an explicitly normalized two-component UV basis.
No EFT coefficient, heavy-mass expansion, model-dependent CG, or Yukawa fit is used.
The input convention is L_Y = -nu_a psi_i (Y[a,i,s] phi_s +
Z[a,i,s] phi_s^*) + h.c., L_m = -psi^T M psi/2 + h.c.
and V_2 = phi^dag A phi + (phi^T B phi + h.c.)/2.

Warning: a source-specific UV -> component-vertex adapter must validate CG,
charge-conjugation phases, and the sign of scalar-potential couplings first.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping
import math
import numpy as np

LOOP = 1.0/(16.0*math.pi**2)


def _complex(a: Any, name: str, shape: tuple[int, ...] | None = None) -> np.ndarray:
    if isinstance(a, Mapping) and set(a) == {'real', 'imag'}:
        value = np.asarray(a['real'], dtype=np.float64) + 1j*np.asarray(a['imag'], dtype=np.float64)
    else:
        value = np.asarray(a, dtype=np.complex128)
    if shape is not None and value.shape != shape:
        raise ValueError(f'{name}: expected {shape}, got {value.shape}')
    if not np.isfinite(value.real).all() or not np.isfinite(value.imag).all():
        raise ValueError(f'{name}: nonfinite values')
    return value.astype(np.complex128)


def encoded(a: Any) -> dict[str, Any]:
    z = np.asarray(a, dtype=np.complex128)
    return {'real':z.real.tolist(), 'imag':z.imag.tolist()}


def takagi(m: np.ndarray, *, rtol: float = 4e-13) -> tuple[np.ndarray, np.ndarray, float]:
    """Return U, masses>=0, residual for U.T @ m @ U=diag(masses).

    Positive eigenvectors of the real symmetric Takagi embedding avoid the
    degeneracy-sensitive phase-fixing of an ordinary complex SVD. Numerical
    null directions are selected using a complex null-space basis.
    """
    m = _complex(m, 'complex symmetric mass')
    if m.ndim != 2 or m.shape[0] != m.shape[1]:
        raise ValueError('Takagi matrix must be square')
    n = m.shape[0]
    scale = float(np.linalg.norm(m, ord=2)) if n else 0.
    if np.linalg.norm(m-m.T) > 2e-11*max(scale, 1.e-80):
        raise ValueError('Takagi input is not complex symmetric')
    if not n:
        return np.zeros((0,0), complex), np.zeros(0), 0.
    if scale == 0.:
        return np.eye(n,dtype=complex), np.zeros(n), 0.
    z = m/scale
    a, b = z.real, z.imag
    h = np.block([[a,-b],[-b,-a]])
    ev, w = np.linalg.eigh(h)
    # Positive eigenvalues are exactly the nonzero Takagi singular values.
    choices = [i for i in range(2*n) if ev[i] > rtol]
    if len(choices) > n:
        raise RuntimeError('Takagi embedding yielded too many positive directions')
    ucols = [w[:n,i] + 1j*w[n:,i] for i in choices]
    svals = [max(ev[i],0.)*scale for i in choices]
    if len(ucols) < n:
        # The SVD right null vectors satisfy m u=0. All null vectors are
        # orthogonal to the positive Takagi subspace for symmetric m.
        _left, singular, vh = np.linalg.svd(z)
        null = vh.conj().T[:,singular <= rtol]
        if null.shape[1] != n-len(ucols):
            raise RuntimeError('Takagi zero-mode rank is numerically ambiguous')
        ucols.extend(null[:,k] for k in range(null.shape[1]))
        svals.extend([0.] * null.shape[1])
    U = np.stack(ucols,axis=1)
    masses = np.asarray(svals)
    order = np.argsort(masses)
    masses = masses[order]
    U = U[:,order]
    test = U.T @ m @ U
    relative = max(float(np.linalg.norm(test-np.diag(masses))/scale),
                   float(np.linalg.norm(U.conj().T@U-np.eye(n))))
    if relative > 2e-9:
        raise RuntimeError(f'Takagi failed: relative residual {relative:.3g}')
    return U, masses, relative


def scalar_real_basis(A: np.ndarray, B: np.ndarray) -> tuple[np.ndarray,np.ndarray,float]:
    """phi=C R; V2=R^T diag(ms2) R /2, real canonical R."""
    A = _complex(A,'A')
    B = _complex(B,'B',A.shape)
    if A.ndim != 2 or A.shape[0] != A.shape[1]:
        raise ValueError('A and B must be square')
    scale = max(float(np.linalg.norm(A)),float(np.linalg.norm(B)),1.)
    if np.linalg.norm(A-A.conj().T)>1e-11*scale:
        raise ValueError('A must be Hermitian')
    if np.linalg.norm(B-B.T)>1e-11*scale:
        raise ValueError('B must be symmetric')
    # phi=(x+i y)/sqrt(2).  Note the antisymmetric imaginary part of A.
    h = np.block([[A.real+B.real,-A.imag-B.imag],
                  [A.imag-B.imag,A.real-B.real]])
    if np.linalg.norm(h-h.T)>1e-11*scale:
        raise RuntimeError('Scalar real Hessian unexpectedly asymmetric')
    val,O = np.linalg.eigh((h+h.T)*0.5)
    if np.min(val) < -1e-11*scale:
        raise ValueError('Tachyonic scalar eigenstate')
    if np.min(val) <= 0:
        raise ValueError('Zero or negative scalar eigenmass: logarithmic IR case requires separate treatment')
    n = A.shape[0]
    C = (O[:n,:] + 1j*O[n:,:])/math.sqrt(2.)
    relative = float(np.linalg.norm(O.T@h@O-np.diag(val))/scale)
    return C,val,relative


def b0_finite(mf2: float, ms2: float, mu2: float) -> float:
    """MSbar finite B0(0;mf2,ms2), dimensionless, stable at ms2=mf2."""
    if min(mf2,ms2,mu2)<=0 or not all(map(math.isfinite,(mf2,ms2,mu2))):
        raise ValueError('B0 requires positive finite squared masses and reference scale')
    x = (ms2-mf2)/mf2
    if abs(x)<1e-5:
        # F=(1+x)ln(1+x)/x; B0=1-ln(mf2/mu2)-F
        f=1+x*(0.5+x*(-1/6+x*(1/12+x*(-1/20+x/30))))
    else:
        f=(1+x)*math.log1p(x)/x
    return 1-math.log(mf2/mu2)-f


def _matrix_diag(m: np.ndarray) -> np.ndarray:
    return np.diag(np.diag(m))


def _full_report_record(mat: np.ndarray, mass_gev_to_ev: float = 1.e9) -> dict[str,Any]:
    U, masses, res = takagi(mat)
    return {'complex_symmetric_matrix_gev':encoded(mat),
            'neutrino_masses_ev':(masses*mass_gev_to_ev).tolist(),
            'takagi_rotation':encoded(U), 'takagi_residual':res,
            'matrix_relative_asymmetry':float(np.linalg.norm(mat-mat.T)/max(np.linalg.norm(mat),1e-100))}


def evaluate(data: Mapping[str,Any], *, retain_individual_terms: bool=True) -> dict[str,Any]:
    """Exact LO UV mass from ALL declared Weyl fermions, real scalars and vertices.

    Individual contributions are reported per (vertex-pair,fermion,scalar)
    and retained even for cancellations. The full calculation is at p=0.
    """
    mm=_complex(data['fermion_mass_gev'],'fermion_mass_gev')
    aa=_complex(data['scalar_A_gev2'],'scalar_A_gev2')
    bb=_complex(data['scalar_B_gev2'],'scalar_B_gev2')
    nf=mm.shape[0]
    ns=aa.shape[0]
    uf,mf,err_f=takagi(mm)
    cs,ms2,err_s=scalar_real_basis(aa,bb)
    mu=float(data.get('renormalization_scale_gev',max(mf)))
    if mu<=0 or not math.isfinite(mu): raise ValueError('Invalid renormalization scale')
    vertices=data.get('vertices',[])
    if not isinstance(vertices,list) or not vertices:
        raise ValueError('Source-derived component vertices required; never infer Yukawa contractions')
    names=[]; transformed=[]
    for row in vertices:
        name=str(row['name'])
        if name in names:raise ValueError('Duplicate vertex name '+name)
        names.append(name)
        y=_complex(row.get('field',np.zeros((3,nf,ns))),name+'.field',(3,nf,ns))
        z=_complex(row.get('conjugate',np.zeros((3,nf,ns))),name+'.conjugate',(3,nf,ns))
        # fermion psi_i=sum_k U_ik chi_k; scalar phi_s=sum_r C_sr R_r
        g=np.einsum('ais,ik,sr->akr',y,uf,cs,optimize=True)
        g+=np.einsum('ais,ik,sr->akr',z,uf,cs.conj(),optimize=True)
        transformed.append(g)
    # Coefficients of the MSbar UV pole and of ln(mu^2). These MUST cancel
    # for a finite radiatively generated Majorana mass with no counterterm.
    uvpole=np.zeros((3,3),complex)
    total=np.zeros((3,3),complex)
    pairmat={}
    individual=[]
    abscontributions=0.
    maxterm=0.
    for a in range(len(names)):
        for b in range(a,len(names)):
            label=f'{names[a]} x {names[b]}'
            vmat=np.zeros((3,3),complex)
            for k,heavy in enumerate(mf):
                for r,scalar2 in enumerate(ms2):
                    vv=np.outer(transformed[a][:,k,r],transformed[b][:,k,r])
                    if b!=a:
                        vv=vv+vv.T
                    uvpole+=LOOP*heavy*vv
                    term=(LOOP*heavy*b0_finite(float(heavy**2),float(scalar2),mu**2))*vv
                    vmat+=term
                    termsize=float(np.linalg.norm(term))
                    abscontributions+=termsize
                    maxterm=max(maxterm,termsize)
                    if retain_individual_terms:
                        individual.append({'pair':label,'fermion_eigenstate':int(k),
                                           'scalar_eigenstate':int(r),
                                           'fermion_mass_gev':float(heavy),
                                           'scalar_mass_gev':float(math.sqrt(scalar2)),
                                           'matrix_gev':encoded(term)})
            pairmat[label]=vmat
            total+=vmat
    norm=float(np.linalg.norm(total))
    possible_ratios=(np.sqrt(ms2)[:,None]/mf[None,:]).ravel()
    pole_norm=float(np.linalg.norm(uvpole))
    div_relative=pole_norm/max(LOOP*sum(mf)*sum(np.linalg.norm(g)**2 for g in transformed),1e-100)
    return {
        'status':'uv_mass_evaluated_unverified_source_mapping',
        'model':str(data.get('model','unknown')),
        'input_basis':'two_component_Weyl_plus_complex_scalars',
        'loop_order':1,
        'loop_momentum':'p=0',
        'subtraction':'MSbar finite B0; renormalization-scale cancellation checked',
        'renormalization_scale_gev':mu,
        'source_provenance':data.get('provenance'),
        'mass_hierarchy':{'min_scalar_to_fermion_mass_ratio':float(np.min(possible_ratios)),
                          'max_scalar_to_fermion_mass_ratio':float(np.max(possible_ratios)),
                          'max_scalar_squared_to_fermion_squared_ratio':float(np.max(possible_ratios)**2),
                          'eft_expansion_small_parameter_is_ms2_over_mf2':True},
        'fermion_eigenmasses_gev':mf.tolist(),
        'scalar_eigenmasses_gev':np.sqrt(ms2).tolist(),
        'fermion_takagi_residual':err_f,
        'scalar_mass_diagonalization_residual':err_s,
        'vertex_pair_matrices_gev':{k:encoded(v) for k,v in pairmat.items()},
        'individual_terms':individual if retain_individual_terms else None,
        'total':_full_report_record(total),
        'diagnostics':{
            'max_single_term_frobenius_gev':maxterm,
            'sum_individual_frobenius_gev':abscontributions,
            'total_frobenius_gev':norm,
            'cancellation_ratio_total_over_sum':norm/max(abscontributions,1e-100),
            'cancellation_factor_sum_over_total':abscontributions/max(norm,1e-100),
            'relative_uv_pole_coefficient':div_relative,
            'uv_pole_matrix_gev':encoded(uvpole),
            'uv_finite': bool(div_relative<1e-8),
            'neutrino_mass_matrix_symmetric':float(np.linalg.norm(total-total.T))<1e-10*max(norm,1e-90)
        }
    }
