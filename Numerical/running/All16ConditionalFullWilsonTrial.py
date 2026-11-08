"""All-16 *exploratory* full-tensor Wilson-RGE trial and joint scalar matching.

Uses actual all-16 F-first LLSS source release, existing v6 scalar mass and
HHSS tensors, and the actual two EFT1 renormalisable snapshots at muF/muS.
No new scalar couplings are invented: the RGE uses the supplied HHSS tensor
and Higgs self quartic, with all other quartic components set to zero.
The latter is an explicit reduced-potential approximation, NOT the full
model-specific RGBeta quartic tensor used in the A(-4) benchmark.

Convention: C0_master = (1/2) * C0_legacy_joint_matching_formal.  This
exactly matches the legacy joint-scalar (O5, no explicit 1/2) projector to
the master (1/2 C5 O5) projector, for a fixed source and quartic.  It is
an explicit numerical dictionary, not a native Matchete CC measurement.

EFT order: UV -> integrate out ALL heavy fermions -> EFT1(S1,S2)
           -> integrate out BOTH heavy scalars simultaneously -> SMEFT.
"""
from __future__ import annotations

import csv
from hashlib import sha256
import json
from math import log, sqrt
from pathlib import Path
from typing import Any, Mapping
from itertools import permutations

import numpy as np

from Numerical.joint_scalar_matching.adapter import ALL_MODELS, CLASS_DIMS, normalize_model
from Numerical.joint_scalar_matching.fermion_boundary_bridge import (
    load_release_manifest, load_source, slug,
)
from Numerical.joint_scalar_matching.large_batch import load_background
from Numerical.joint_scalar_matching.numerical_benchmark_v5 import load_running_snapshot
from Numerical.joint_scalar_matching.special_mass_v6 import SPECIAL, build_special_basis
from Numerical.running import Aminus4EFT1Couplings as eft
from Numerical.joint_scalar_matching.Aminus4ConditionalSMEFT import (
    project_master, neutrino_mass_proxy, maxabs, relative_error,
)
from Numerical.joint_scalar_matching.core import scalar_bubble_matching, project_weinberg as project_joint
from RGE.running.intermediate.GeneralEFT1SevenTermFlow import Couplings, total_beta, LOOP
from RGE.running.intermediate.GeneralEFT1SixTermReferenceCandidate import six_term_beta

FORMAT = 't3-all16-conditional-full-wilson-reduced-quartic-trial-v2'
ROOT = Path(__file__).resolve().parents[2]
DEFAULT_RELEASE = Path('output/general_observable/lo_fermion_boundary/releases/source_021dcc1e8fe8ba0f')
DEFAULT_V6 = Path('output/general_observable/lo_joint_scalar_matching/exploratory_all16_v6_native_batch1_local')
DEFAULT_RUNNING = Path('output/general_observable/all16_v1')
# Only the existing 45 Weyl SM fermions survive F-first integration.
LEPTON_POSITIONS = np.array([[eft.FI[f'L[{p},{i}]'] for i in (1,2)]
                             for p in (1,2,3)], dtype=int)
SCALAR_CONVENTION_FACTOR = 0.5


class All16TrialBlocked(ValueError):
    pass


def require(ok: Any, why: str) -> None:
    if not bool(ok):
        raise All16TrialBlocked(why)


def checked_path(root: Path, relative: Path) -> Path:
    p = Path(relative)
    return p if p.is_absolute() else root / p


def read_json(path: Path) -> dict:
    require(path.is_file(), f'missing JSON file: {path}')
    result = json.loads(path.read_text(encoding='utf-8-sig'))
    require(isinstance(result,dict), f'not a JSON object: {path}')
    return result


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def enc(z: np.ndarray) -> dict:
    a = np.asarray(z, complex)
    require(np.isfinite(a).all(), 'nonfinite complex matrix')
    return {'real':a.real.tolist(), 'imag':a.imag.tolist()}


def decode_3(obj: Mapping[str,Any]) -> np.ndarray:
    require(isinstance(obj,dict) and 'real' in obj and 'imag' in obj,
            'missing complex 3x3 matrix real/imag parts')
    z=np.asarray(obj['real'],float)+1j*np.asarray(obj['imag'],float)
    require(z.shape==(3,3) and np.isfinite(z).all(), 'not a finite complex 3x3 matrix')
    return z


def yukawa_at_endpoint(state: Mapping[str,Any], nscalar: int) -> np.ndarray:
    """Reuse verified SM 45-Weyl/Yukawa contraction from the A(-4) adapter."""
    require(nscalar in (12,16), 'T3 model requires 12 or 16 real scalars')
    point = {name:0. for name in eft.BETA_NAMES}
    for name in ('gY','g2','g3','yu','yd','ye','lambdaH'):
        require(name in state, 'missing EFT1 SM coupling: '+name)
        point[name] = state[name]
    point['lambdaT3'] = 0j
    y12 = eft.yukawa_tensor(point)
    require(y12.shape==(45,45,12), 'unexpected common SM Yukawa basis')
    y=np.zeros((45,45,nscalar),complex)
    y[:,:,:4] = y12[:,:,:4]
    return y


def native_triplet_su2_generators() -> tuple[np.ndarray, ...]:
    """SU(2) spin-one matrices in the Matchete LLSS source component basis.

    The release's Cartesian triplet components differ from the ordered
    m=+1,0,-1 spherical basis of eft.su2_generators(3). The intertwiner
    was determined from the A(-4) Matchete matched LL-triplet CG, and
    independently checked here against actual A/D/E LLSS Ward identities.

      S_weight = N S_native;  T_native = N^dagger T_weight N.

    N is unitary. Its matrix also explains why the native triplet gauge
    generators are pure-imaginary antisymmetric (a real representation).
    Applying the spherical generators directly to the Matchete source was
    the A/D/E 2.3--2.6% false Weinberg-projection failure in trial v1.
    """
    z = 1.0 / sqrt(2.0)
    N = np.array([[-z, 0, 1j*z],
                  [ 0, 1, 0],
                  [ z, 0, 1j*z]], dtype=complex)
    require(np.max(np.abs(N.conj().T@N-np.eye(3))) < 1e-14,
            'Matchete-to-weight triplet map is not unitary')
    result=tuple(N.conj().T@t@N for t in eft.su2_generators(3))
    for t in result:
        require(np.max(abs(t-t.conj().T)) < 1e-14 and
                np.max(abs(t+t.T)) < 1e-14,
                'native Cartesian triplet generator not Hermitian/pure imaginary')
    return result


def original_scalar_generators(model: str) -> dict[str,tuple[np.ndarray,...]]:
    """H+S1+S2 SU2/U1/SU3 in the current native complex-component order."""
    model=normalize_model(model)
    d1,d2,_=CLASS_DIMS[model[0]]
    alpha=int(model.split('_')[1])
    n=4+2*(d1+d2)
    generators={}
    for group in ('SU2','U1','SU3'):
        if group=='SU2':
            g=np.zeros((3,n,n),complex)
            blocks=((2,.5,0),(d1,alpha/2,4),(d2,(alpha+2)/2,4+2*d1))
            for dim,hyper,off in blocks:
                native_t = (native_triplet_su2_generators() if dim == 3
                            else eft.su2_generators(dim))
                for k,t in enumerate(native_t):
                    dim2=2*dim
                    g[k,off:off+dim2,off:off+dim2]=eft.real_scalar_generator(t)
        elif group=='U1':
            g=np.zeros((1,n,n),complex)
            for dim,hyper,off in ((2,.5,0),(d1,alpha/2,4),(d2,(alpha+2)/2,4+2*d1)):
                dim2=2*dim
                g[0,off:off+dim2,off:off+dim2]=eft.real_scalar_generator(hyper*np.eye(dim))
        else:
            g=np.zeros((8,n,n),complex)
        generators[group]=tuple(g)
    return generators


def rotated_scalar_generators(model:str, rotation: np.ndarray) -> dict[str,tuple[np.ndarray,...]]:
    n=4+rotation.shape[0]
    require(rotation.shape==(n-4,n-4) and np.isrealobj(rotation)
            and np.isfinite(rotation).all()
            and np.allclose(rotation.T@rotation,np.eye(n-4),atol=1e-10,rtol=1e-10),
            'invalid native-to-matching real scalar rotation')
    full=np.eye(n)
    full[4:,4:]=rotation
    return {group:tuple(full.T @ t @ full for t in ts)
            for group,ts in original_scalar_generators(model).items()}


def source_ward_residual(source_matching: np.ndarray,
                         scalar_generators: Mapping[str,tuple[np.ndarray,...]]) -> dict:
    """SU(2)xU(1) Ward test of the actual 3-family LLSS source at muF.

    Uses the 2 weak components (all 3 flavours) instead of constructing a
    45-Weyl dense test tensor. This is algebraically the same gauge action
    on its nonzero support; no numerical tolerance on GeV^-1 coefficients.
    """
    source=np.asarray(source_matching,complex)
    require(source.ndim==6 and source.shape[:4]==(3,3,2,2) and
            source.shape[4]==source.shape[5] and np.isfinite(source).all(),
            'invalid LLSS tensor for Ward audit')
    norm=maxabs(source)
    require(norm>0,'empty source cannot pass source Ward audit')
    result={}
    for group, weak_gens in (('SU2',eft.su2_generators(2)),
                             ('U1',(np.eye(2,dtype=complex)*(-0.5),))):
        trip=scalar_generators[group]
        require(len(weak_gens)==len(trip),'SU2/U1 generator inventory mismatch')
        residual=0.0
        for tf,ts in zip(weak_gens,trip):
            theta=np.asarray(ts)[4:,4:]
            delta=(np.einsum('ki,pqkjab->pqijab',tf,source,optimize=True)
                   +np.einsum('kj,pqikab->pqijab',tf,source,optimize=True)
                   +np.einsum('ca,pqijcb->pqijab',theta,source,optimize=True)
                   +np.einsum('cb,pqijac->pqijab',theta,source,optimize=True))
            residual=max(residual,maxabs(delta)/norm)
        result[group]=float(residual)
    return result


def quartic_ward_residual(q:np.ndarray,
                          scalar_generators:Mapping[str,tuple[np.ndarray,...]]) -> dict:
    """Verify the provided reduced quartic is an invariant in this exact basis."""
    q=np.asarray(q,float)
    require(q.ndim==4 and len(set(q.shape))==1 and np.isfinite(q).all(),
            'invalid reduced scalar quartic for Ward audit')
    norm=maxabs(q)
    require(norm>0,'empty quartic cannot pass Ward audit')
    result={}
    for group in ('SU2','U1'):
        maximum=0.
        for theta in scalar_generators[group]:
            delta=(np.einsum('ea,ebcd->abcd',theta,q,optimize=True)
                   +np.einsum('eb,aecd->abcd',theta,q,optimize=True)
                   +np.einsum('ec,abed->abcd',theta,q,optimize=True)
                   +np.einsum('ed,abce->abcd',theta,q,optimize=True))
            maximum=max(maximum,maxabs(delta)/norm)
        result[group]=float(maximum)
    return result


def higgs_quartic(n: int) -> np.ndarray:
    """Fourth derivative of (1/2)(Hdag H)^2, lambdaH=1."""
    out=np.zeros((n,n,n,n),float)
    for a in range(4):
        for b in range(4):
            for c in range(4):
                for d in range(4):
                    out[a,b,c,d]=int(a==b and c==d)+int(a==c and b==d)+int(a==d and b==c)
    return out


def hhss_full_quartic(hhss: np.ndarray) -> np.ndarray:
    """Fill all permutations of exactly two Higgs and two heavy scalar legs."""
    q=np.asarray(hhss,float)
    require(q.ndim==4 and q.shape[:2]==(4,4) and q.shape[2]==q.shape[3]
            and np.isfinite(q).all(), 'invalid source-derived real HHSS quartic')
    n=4+q.shape[2]
    full=np.zeros((n,n,n,n),float)
    for a in range(4):
        for b in range(4):
            for s in range(n-4):
                for t in range(n-4):
                    v=q[a,b,s,t]
                    inds=(a,b,s+4,t+4)
                    # Every index placement with exactly 2 Higgs and 2 heavy slots.
                    for perm in set(permutations(inds)):
                        full[perm]=v
    require(np.max(abs(full-full.swapaxes(0,1)))<1e-13 and
            np.max(abs(full-full.swapaxes(1,2)))<1e-13 and
            np.max(abs(full-full.swapaxes(2,3)))<1e-13,
            'quartic symmetry reconstruction failed')
    return full


def embedded_boundary(source:np.ndarray, rotation:np.ndarray,
                      *,factor:float=SCALAR_CONVENTION_FACTOR) -> np.ndarray:
    """Embed 3x3x2x2xNheavy xNheavy in the actual 45-Weyl master basis."""
    source=np.asarray(source,complex)
    nheavy=rotation.shape[0]
    require(source.shape==(3,3,2,2,nheavy,nheavy),
            'source boundary not complete 3-generation LLSS')
    require(float(factor)>0 and np.isfinite(factor), 'invalid declared normalization')
    rotated=np.einsum('pqijab,ak,bl->pqijkl',source,rotation,rotation,optimize=True)
    n=4+nheavy
    c=np.zeros((45,45,n,n),complex)
    for p in range(3):
        for q in range(3):
            for i in range(2):
                for j in range(2):
                    c[LEPTON_POSITIONS[p,i],LEPTON_POSITIONS[q,j],4:,4:]=factor*rotated[p,q,i,j]
    require(np.max(abs(c-c.swapaxes(0,1)))<1e-25 and
            np.max(abs(c-c.swapaxes(2,3)))<1e-25,
            'embedded source lost separate pair symmetry')
    return c


def master_projection(amplitude:np.ndarray, *,tol:float=1e-6) -> tuple[np.ndarray,float]:
    return project_master(amplitude, strict=True, tol=tol)


def regression_fixed_order(c0:np.ndarray, bg, known:Mapping[str,Any]) -> dict:
    """Demand exact agreement with earlier all16 matching at zero running.

    The old joint projector is twice the master template.  The chosen 1/2
    source conversion must therefore reproduce BOTH old log and bubble.
    """
    h=generic_leptons(c0,True)
    b=np.einsum('pqijab,ABab->pqijAB',h,bg.quartic,optimize=True)
    beta,rb=master_projection(b)
    oldbeta,_=project_joint(b,strict=False)
    require(relative_error(beta,2*oldbeta,floor=1e-55)<3e-9,
            'old joint and master Weinberg projectors differ by more than factor 2')
    log_c5=LOOP*log(bg.mu_s/bg.mu_f)*beta
    bubble=scalar_bubble_matching(h,bg.quartic,bg.mass2,bg.mu_s)
    bubble_c5,rfinite=master_projection(bubble['amplitude'])
    old_log=decode_3(known['C5_EFT1_induced_leading_log_per_N_gev_inverse'])
    old_bubble=decode_3(known['C5_joint_scalar_finite_bubble_per_N_gev_inverse'])
    lerr=relative_error(log_c5,old_log,floor=1e-55)
    berr=relative_error(bubble_c5,old_bubble,floor=1e-55)
    require(lerr<2e-8 and berr<2e-8,
            f'old v6 C5 convention regression fails: log={lerr:.4g}, bubble={berr:.4g}')
    return {'leading_log_relative_error':lerr,'finite_bubble_relative_error':berr,
            'leading_log_max_abs_gev_inv':maxabs(log_c5),
            'finite_bubble_max_abs_gev_inv':maxabs(bubble_c5),
            'log_weinberg_residual':rb,'finite_weinberg_residual':rfinite}


def generic_leptons(c:np.ndarray, heavy:bool)->np.ndarray:
    """Generic n=12/16 scalar analogue of the tested A(-4) extraction."""
    c=np.asarray(c,complex)
    n=c.shape[-1]
    require(c.shape==(45,45,n,n) and n in (12,16), 'wrong full C shape')
    sl=slice(4,n) if heavy else slice(0,4)
    dim=n-4 if heavy else 4
    out=np.empty((3,3,2,2,dim,dim),complex)
    for p in range(3):
        for q in range(3):
            for i in range(2):
                for j in range(2):
                    out[p,q,i,j]=c[LEPTON_POSITIONS[p,i],LEPTON_POSITIONS[q,j],sl,sl]
    return out


def six_term_at(c:np.ndarray, y:np.ndarray, quartet:np.ndarray,
                gauge_pairs:tuple) -> np.ndarray:
    inp=Couplings(y,quartet,gauge_pairs)
    terms=six_term_beta(c,inp)
    require(terms['crossed_yukawa'].shape==c.shape and
            np.max(abs(terms['crossed_yukawa']))==0,
            'requested six-term beta unexpectedly added crossed term')
    return total_beta(terms)


def midpoint_transport(c0:np.ndarray, y_f:np.ndarray, y_s:np.ndarray,
                       g_f:Mapping[str,float], g_s:Mapping[str,float],
                       hhss:np.ndarray, lambda_h_f:float, lambda_h_s:float,
                       fermionic:Mapping[str,tuple[np.ndarray,...]],
                       scalars:Mapping[str,tuple[np.ndarray,...]],
                       mu_f:float,mu_s:float,steps:int) -> tuple[np.ndarray,dict]:
    """Fixed log(mu) step RK2, full 45x45xnsxns complex tensor.

    Couplings linearly interpolate in log(mu) between two actual snapshots;
    HHSS quartic is held at muS because only its muS tensor is supplied.
    Scalar self/mixed quartics not supplied by v6 are absent in THIS trial.
    """
    require(1<=steps<=64 and 0<mu_s<mu_f,'invalid midpoint integration interval')
    n=c0.shape[-1]
    require(c0.shape==(45,45,n,n),'wrong 45-Weyl tensor')
    qbase=hhss_full_quartic(hhss)
    qh=higgs_quartic(n)
    dt=log(mu_s/mu_f)/steps
    c=c0.copy()
    evals=0
    def rhs(cval:np.ndarray, frac:float) -> np.ndarray:
        nonlocal evals
        yf=(1-frac)*y_f+frac*y_s
        la=(1-frac)*lambda_h_f+frac*lambda_h_s
        q=qbase+la*qh
        g=[]
        for name,k in (('SU3','g3'),('SU2','g2'),('U1','gY')):
            coupling=(1-frac)*float(g_f[k])+frac*float(g_s[k])
            pairs=tuple(zip(fermionic[name],scalars[name]))
            g.append((coupling,pairs))
        evals+=1
        return LOOP*six_term_at(cval,yf,q,tuple(g))
    for i in range(steps):
        # k1= rhs(mu(t_i)), k2=rhs(midpoint, C + dt*k1/2)
        t=(i/steps)
        k1=rhs(c,t)
        k2=rhs(c+.5*dt*k1,(i+.5)/steps)
        c+=dt*k2
        require(np.isfinite(c).all(),f'nonfinite Wilson at RK2 interval {i}')
    residual=max(maxabs(c-c.swapaxes(0,1)),maxabs(c-c.swapaxes(2,3)))
    require(residual<1e-8*max(maxabs(c),1e-100),
            f'full Wilson exchange symmetry failed (abs={residual:g})')
    return c,{'scheme':'explicit_midpoint_RK2_2_endpoint_interpolation',
              'steps':steps,'beta_rhs_evaluations':evals,'mu_F_gev':mu_f,
              'mu_S_gev':mu_s,'delta_t_per_step':dt,
              'quartic_tensor':'HHSS at muS frozen + Higgs self lambdaH interpolated; OTHER quartics OMITTED',
              'complete_scalar_quartic_RGE_included':False,
              'full_complex_C_relative_change':float(np.linalg.norm(c-c0)/max(np.linalg.norm(c0),1e-100)),
              'exchange_symmetry_max_abs_residual':residual}


def c5_at_ms(c_ms:np.ndarray,bg,*,vev:float=246.22)->dict:
    hh=generic_leptons(c_ms,False)
    heavy=generic_leptons(c_ms,True)
    direct,rd=master_projection(hh)
    bubble=scalar_bubble_matching(heavy,bg.quartic,bg.mass2,bg.mu_s)
    cb,rb=master_projection(bubble['amplitude'])
    require(len(bubble['pair_data'])==(bg.mass2.shape[0]*(bg.mass2.shape[0]+1))//2,
            'joint-scalar pair inventory incomplete')
    cbjoint,rjoint=project_joint(bubble['amplitude'],strict=False)
    require(relative_error(cb,2*cbjoint,floor=1e-55)<2e-9,
            'master and joint projector factor 2 fails')
    total=direct+cb
    masses=neutrino_mass_proxy(total,vev)
    return {'c5_direct':enc(direct),'c5_bubble':enc(cb),'c5_total':enc(total),
            'c5_max_abs_gev_inv':maxabs(total),'masses_eV':masses['takagi_singular_values_eV'],
            'mnu_matrix_eV':enc(-0.5*vev**2*1e9*total),
            'scalar_mass_eigenvalues_gev':np.sqrt(np.linalg.eigvalsh(bg.mass2)).tolist(),
            'all_scalar_eigenstate_pairs':len(bubble['pair_data']),
            'direct_weinberg_residual':rd,'bubble_weinberg_residual':rb,
            'physical_neutrino_masses_certified':False}


def _load_old_response(v6_root:Path, model:str)->dict:
    path=v6_root/'C5_per_unknown_NGamma_exploratory'/(slug(model)+'.json')
    obj=read_json(path)
    require(obj.get('model')==model and isinstance(obj.get('response'),dict),
            'legacy all16 response does not match model')
    return obj['response']


def run_one(root:Path, model:str, *,release_dir:Path, v6_dir:Path,
            running_dir:Path,steps:int=4,vev:float=246.22,
            save_tensor_to:Path|None=None)->dict:
    """Read-only inputs; return conditional metrics for this exact model."""
    model=normalize_model(model)
    paths=load_release_manifest(checked_path(root,release_dir))
    src=load_source(paths[model],expected_model=model)
    bgpath=checked_path(root,v6_dir)/'numeric_backgrounds'/(slug(model)+'.json')
    bg=load_background(bgpath,src)
    runpath=checked_path(root,running_dir)/slug(model)/'running.json'
    running=read_json(runpath)
    require(running.get('model')==model, 'running snapshot model does not match source')
    f=running['eft1_at_fermion_threshold']
    s=running['eft1_at_scalar_threshold']
    mf=src.mu_f_gev
    ms=bg.mu_s
    require(abs(float(running['scales_gev']['F'])-mf)<1e-6 and
            abs(float(running['scales_gev']['S'])-ms)<1e-6,
            'renormalisable EFT1 snapshots do not bracket source and matching')
    n=4+len(src.labels)
    require(n in (12,16),'unsupported scalar dimension')
    from Numerical.joint_scalar_matching.numerical_benchmark_v5 import REGULAR
    if model in SPECIAL:
        snapshot=load_running_snapshot(runpath,model)
        special=build_special_basis(model,snapshot,src)
        native_to_source=special['R_original_from_boundary']
    else:
        require(model in REGULAR,'model not classified in existing v6 scalar pipeline')
        native_to_source=np.eye(n-4)
    # Real source field rotation into background matching basis. For v6 it
    # is identity, but preserve the general existing bg.U contract.
    native_to_matching=native_to_source@bg.U
    generator_scalars=rotated_scalar_generators(model,native_to_matching)
    original_gauge=eft.gauge_generators()
    generator_fermions={name:original_gauge[name]['fermion'] for name in ('SU3','SU2','U1')}
    c0=embedded_boundary(src.raw,bg.U)
    src_matching=np.einsum('pqijab,ak,bl->pqijkl',src.raw,bg.U,bg.U,optimize=True)
    source_ward=source_ward_residual(src_matching,generator_scalars)
    # Both the heavy-source Wilson tensor and the supplied quartic must use
    # the SAME native Matchete representation. A/D/E are triplet-bearing.
    qtest=hhss_full_quartic(bg.quartic)+float(s['lambdaH'])*higgs_quartic(n)
    potential_ward=quartic_ward_residual(qtest,generator_scalars)
    require(max(source_ward.values())<2e-9,
            f'{model}: fermion-threshold SU2/U1 source basis mismatch: {source_ward}')
    require(max(potential_ward.values())<2e-9,
            f'{model}: input reduced potential not SU2/U1 invariant: {potential_ward}')
    legacy=_load_old_response(checked_path(root,v6_dir),model)
    controls=regression_fixed_order(c0,bg,legacy)
    yf=yukawa_at_endpoint(f,n)
    ys=yukawa_at_endpoint(s,n)
    cf,transport=midpoint_transport(c0,yf,ys,f,s,bg.quartic,
                                    float(f['lambdaH']),float(s['lambdaH']),
                                    generator_fermions,generator_scalars,mf,ms,steps)
    endpoint=c5_at_ms(cf,bg,vev=vev)
    legacy_mat=decode_3(legacy['C5_LLSS_total_one_loop_per_N_gev_inverse'])
    masslegacy=np.sort(np.linalg.svd(legacy_mat,compute_uv=False))*vev**2*0.5*1e9
    masses=np.asarray(endpoint['masses_eV'])
    report={'format':FORMAT,'model':model,'status':'TRIAL_FULL_WILSON_REDUCED_POTENTIAL',
            'source_boundary_sha256':src.hash,
            'source_running_sha256':digest(runpath),
            'source_scalar_background_sha256':digest(bgpath),
            'mu_F_gev':mf,'mu_S_gev':ms,'n_weyl':45,'n_real_scalars':n,
            'n_complex_symmetric_C_coordinates':(45*46//2)*(n*(n+1)//2),
            'source_normalization_for_master':SCALAR_CONVENTION_FACTOR,
            'source_normalization_reason':'legacy joint projector has twice master Weinberg template; half C0 matches legacy finite+log at zero running',
            'source_and_quartic_ward_audit':{'source':source_ward,'reduced_quartic':potential_ward},
            'fixed_order_legacy_regression':controls,
            'legacy_fixed_order_masses_eV':masslegacy.tolist(),
            'transport':transport,
            'endpoint':endpoint,
            'reduced_potential_quartic_coverage':'Higgs self + source HHSS; no heavy-only quartics, no extra HHHH operators',
            'background_method':'actual F and S snapshots; two-point log interpolation for SM gauge/Yukawa/lambdaH; HHSS at S held constant',
            'matchete_signed_CC_factor_measured':False,
            'full_model_specific_quartic_tensor_included':False,
            'full_one_loop_matched_neutrino_mass_prediction':False}
    if save_tensor_to is not None:
        save_tensor_to.parent.mkdir(parents=True,exist_ok=True)
        require(not save_tensor_to.exists(),'refusing to overwrite saved tensor')
        np.savez_compressed(save_tensor_to,
            full_C_at_MF_gev_inv=c0,full_C_at_MS_gev_inv=cf,
            model=np.array(model),
            c5_combined_gev_inv=decode_3(endpoint['c5_total']),
            neutrino_mass_eV=np.asarray(masses))
        report['full_tensor_npz']=str(save_tensor_to.name)
    return report


def run_all(root:Path, output:Path, *,models:list[str]|None=None,
            release_dir:Path=DEFAULT_RELEASE, v6_dir:Path=DEFAULT_V6,
            running_dir:Path=DEFAULT_RUNNING,steps:int=4,vev:float=246.22,
            save_tensors:bool=False) ->dict:
    root=Path(root).resolve()
    out=checked_path(root,output)
    require(not out.exists(), 'refusing to overwrite output directory: '+str(out))
    models=list(ALL_MODELS) if models is None else [normalize_model(m) for m in models]
    require(len(set(models))==len(models) and models,'duplicate or missing model labels')
    rows=[]
    out.mkdir(parents=True)
    for model in models:
        print('[all16 full Wilson trial] '+model,flush=True)
        try:
            tensor_path=out/'tensors'/(slug(model)+'.npz') if save_tensors else None
            result=run_one(root,model,release_dir=release_dir,v6_dir=v6_dir,
                           running_dir=running_dir,steps=steps,vev=vev,
                           save_tensor_to=tensor_path)
            path=out/'per_model'/(slug(model)+'.json')
            path.parent.mkdir(parents=True,exist_ok=True)
            path.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
            m=result['endpoint']['masses_eV']
            ref=result['legacy_fixed_order_masses_eV']
            rows.append({'model':model,'status':'TRIAL_FULL_WILSON_REDUCED_POTENTIAL',
                         'scalar_type':('mixed_or_split' if model in SPECIAL else 'ordinary'),
                         'm1_eV':m[0],'m2_eV':m[1],'m3_eV':m[2],
                         'legacy_m1_eV':ref[0],'legacy_m2_eV':ref[1],'legacy_m3_eV':ref[2],
                         'C_relative_change':result['transport']['full_complex_C_relative_change'],
                         'source_Ward_max':max(result['source_and_quartic_ward_audit']['source'].values()),
                         'quartic_Ward_max':max(result['source_and_quartic_ward_audit']['reduced_quartic'].values()),
                         'scalar_pairs':result['endpoint']['all_scalar_eigenstate_pairs'],
                         'max_weinberg_residual':max(result['endpoint']['direct_weinberg_residual'],
                                                     result['endpoint']['bubble_weinberg_residual']),
                         'blocker':''})
            print('  masses (eV): '+str([float(v) for v in m]),flush=True)
        except Exception as exc:
            # Preserve independent diagnostic for each model; no fake C5/mass.
            rows.append({'model':model,'status':'BLOCKED_'+type(exc).__name__,
                         'scalar_type':('mixed_or_split' if model in SPECIAL else 'ordinary'),
                         'm1_eV':None,'m2_eV':None,'m3_eV':None,
                         'legacy_m1_eV':None,'legacy_m2_eV':None,'legacy_m3_eV':None,
                         'C_relative_change':None,'source_Ward_max':None,
                         'quartic_Ward_max':None,'scalar_pairs':None,
                         'max_weinberg_residual':None,
                         'blocker':str(exc)})
            print('  BLOCKED: '+str(exc),flush=True)
    columns=list(rows[0])
    with (out/'all16_trial_comparison.csv').open('w',newline='',encoding='utf-8') as stream:
        writer=csv.DictWriter(stream,fieldnames=columns)
        writer.writeheader();writer.writerows(rows)
    success=sum(x['status']=='TRIAL_FULL_WILSON_REDUCED_POTENTIAL' for x in rows)
    summary={'format':FORMAT,'status':'COMPLETE' if success==len(models) else 'PARTIAL',
             'models_checked':len(models),'models_success':success,
             'models_blocked':len(models)-success,
             'input_revision_note':'trial-v2 triplet source basis fix against GitHub f3c12e8 (2026-11-08); local data revalidated per model',
             'integrator':'six-term Fonseca Eq22 gauge-complete candidate, 2-endpoint explicit midpoint RK2',
             'steps':steps,'mu_F_to_mu_S_all_fermions_then_all_scalars':True,
             'scalar_matching':'one simultaneous scalar bubble for each model, no separate scalar EFT2',
             'quartic_coverage_complete':False,
             'conditional_rough_numerics':True,'physical_mass_prediction':False,
             'model_results':rows,
             'file':'all16_trial_comparison.csv'}
    (out/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print('All16 reduced-quartic full-Wilson trial:',success,'/',len(models),'models;',
          'output:',out,flush=True)
    return summary
