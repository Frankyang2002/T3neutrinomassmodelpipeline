"""Deterministic synthetic algebra fixtures. NOT native A/B physics results.

Generates a full, nontrivial real quartic, complex full flavour Wilson,
non-degenerate scalar mass eigenvalues, explicit special-state phases,
and exact Weinberg-template contraction by an engineered synthetic HHSS
channel. This is a test of numerical plumbing, NOT gauge completeness.
"""
from __future__ import annotations
from hashlib import sha256
from pathlib import Path
from itertools import permutations
import json
import numpy as np
from .models import ALL_MODELS, from_parameters, normalize_model, n_real, scalar_labels
from .matching import master_weinberg_template, b0bar
from .pipeline import W_FORMAT,B_FORMAT,W_CONVENTION,Q_CONVENTION


def _sha(path:Path)->str:return sha256(path.read_bytes()).hexdigest()


def _assign_hhss_full(full: np.ndarray, hhss:np.ndarray)->None:
    nh=hhss.shape[-1]
    for A in range(4):
        for B in range(4):
            for a in range(nh):
                for b in range(nh):
                    val=hhss[A,B,a,b]
                    if val==0:continue
                    for idx in set(permutations((A,B,a+4,b+4))):
                        full[idx]=val


def _quartic_sectors(n:int)->np.ndarray:
    """Small extra HHHH, SSSS, S1/S2 terms, all real symmetric."""
    size=n+4
    q=np.zeros((size,)*4)
    eye=np.eye(size)
    # Synthetic O(N) invariant quartic; all internal cross-scalar channels kept.
    for a in range(size):
        for b in range(size):
            for c in range(size):
                for d in range(size):
                    q[a,b,c,d]=.008*(eye[a,b]*eye[c,d]+eye[a,c]*eye[b,d]+eye[a,d]*eye[b,c])
    return q


def _write_npz(path:Path, meta:dict, arrays:dict)->Path:
    path.parent.mkdir(parents=True,exist_ok=True)
    npz=path.parent/meta['tensor_npz']
    np.savez_compressed(npz,**arrays)
    meta['tensor_npz_sha256']=_sha(npz)
    path.write_text(json.dumps(meta,indent=2)+'\n',encoding='utf-8')
    return path


def make_fixture(model:str, directory:Path, *, mu_f:float=100000.,
                 mu_s:float=1000.) -> tuple[Path,Path]:
    """Make A/B-like full inputs that are coherently PROVEN synthetic only."""
    model=normalize_model(model)
    n=n_real(model)
    special=model in ('B_-1','C_-1','A_-2','A_+0','D_-2','D_+0','E_-2','E_+0')
    mu=(1.1e5 + 0.7e5j) if special else None
    mass=from_parameters(model,1.8e6,2.5e6,mu_sq=mu)
    M=mass.mass2_canonical_gev2
    ev,U=np.linalg.eigh(M)
    if len(ev)<4:raise RuntimeError('Need >=4 real heavy fields')
    B=np.array([[b0bar(ev[a],ev[b],mu_s) for b in range(4)] for a in range(4)])
    if np.any(abs(B)<1e-10):raise RuntimeError('Synthetic sample accidental B0 zero')
    T=master_weinberg_template()
    Kheavy=np.array([[.72+.11j,.17-.08j,.05+.02j],
                     [.17-.08j,.55-.09j,.09+.04j],
                     [.05+.02j,.09+.04j,.48+.03j]])*1e-12
    Kevolved=Kheavy*(1.03+.02j)
    Cmass=np.zeros((3,3,2,2,n,n),complex)
    Cmass_tree=np.zeros_like(Cmass)
    Cmass[..., :4,:4]=np.einsum('fg,ijab->fgijab',Kevolved,T)
    Cmass_tree[..., :4,:4]=np.einsum('fg,ijab->fgijab',Kheavy,T)
    # Rotating the two scalar legs to the background real-field basis.
    Cback=np.einsum('fgijab,ka,lb->fgijkl',Cmass,U,U,optimize=True)
    Cback_tree=np.einsum('fgijab,ka,lb->fgijkl',Cmass_tree,U,U,optimize=True)
    # A source for special models uses named native coordinates. B is Cartesian.
    R=mass.R_canonical_from_native.T
    Csrc=np.einsum('fgijab,ka,lb->fgijkl',Cback,R,R,optimize=True)
    Csrc_tree=np.einsum('fgijab,ka,lb->fgijkl',Cback_tree,R,R,optimize=True)
    # Engineer a B0-weighted synthetic quartic to give EXACT Weinberg shape.
    qmass=np.zeros((4,4,n,n))
    strength=.19
    for A in range(4):
        for Bidx in range(4):
            for a in range(4):
                for b in range(4):
                    qmass[A,Bidx,a,b]=strength*((A==a and Bidx==b)+(A==b and Bidx==a))/(2*B[a,b])
    qback=np.einsum('ABab,ka,lb->ABkl',qmass,U,U,optimize=True)
    portal=np.zeros((4+n,)*4)
    _assign_hhss_full(portal,qback)
    baseline=_quartic_sectors(n)
    fullq=portal+baseline
    # Extra baseline q_ABab is flavour-blind and adds to the bubble only
    # via a trace. It cannot violate the Higgs/SU2 tensor projection of
    # this synthetically chosen Cmass (its scalar-trace picks T_ijaa).
    # Remove this baseline HHSS to keep the SU2 residual exactly zero.
    baseline[:4,:4,4:,4:]=0
    baseline[:4,4:,:4,4:]=0
    baseline[:4,4:,4:,:4]=0
    baseline[4:,:4,:4,4:]=0
    baseline[4:,:4,4:,:4]=0
    baseline[4:,4:,:4,:4]=0
    fullq=portal+baseline
    nferm=8
    Cfull=np.zeros((nferm,nferm,4+n,4+n),complex)
    Ctree=np.zeros_like(Cfull)
    K5=np.array([[1.1+.12j,.14-.07j,.065+.025j],
                 [.14-.07j,.96+.05j,.12-.03j],
                 [.065+.025j,.12-.03j,.81-.04j]])*1e-15
    W=np.einsum('fg,ijAB->fgijAB',K5,T)
    for f in range(3):
        for g in range(3):
            for i in range(2):
                for j in range(2):
                    Cfull[2*f+i,2*g+j,:4,:4]=W[f,g,i,j]
                    Ctree[2*f+i,2*g+j,:4,:4]=W[f,g,i,j]*.95
                    Cfull[2*f+i,2*g+j,4:,4:]=Csrc[f,g,i,j]
                    Ctree[2*f+i,2*g+j,4:,4:]=Csrc_tree[f,g,i,j]
    # Other Wilson components still reside in the full tensor (not truncated).
    Cfull[6,6,0,0]=2.3e-16+1j*.4e-16
    Ctree[6,6,0,0]=Cfull[6,6,0,0]
    yc=np.diag([2.9e-6,6.1e-4,1.02e-2]).astype(complex)
    yu=np.diag([1.1e-5,6.4e-3,0.91]).astype(complex)
    yd=np.diag([2.7e-5,5.4e-4,2.0e-2]).astype(complex)
    base=Path(directory)
    slug=model.replace('_-','_m').replace('_+','_p')
    loc=base/slug
    label=f'fixture:synthetic-state:{model}:v1'
    model_basis=f'{model}:native_named' if special else f'{model}:canonical_interleaved'
    wmeta={'format':W_FORMAT,'tensor_npz':'wilson_at_MS.npz','model':model,
           'matching_context_id':label,'matched_seed_id':label,
           'lagrangian_source_id':'SYNTHETIC_NOT_ACTUAL_LAGRANGIAN',
           'mu_f_gev':mu_f,'mu_s_gev':mu_s,'units':'GeV','scheme':'MSbar',
           'canonical_wilson_convention':W_CONVENTION,
           'wilson_evolved_to_ms':True,'includes_eft1_weinberg':True,
           'legacy_direct_llss_to_c5_log_added':False,
           'scalar_basis_id':model_basis,'scalar_source_real_labels':list(mass.native_labels),
           'lepton_indices':[[0,1],[2,3],[4,5]],
           'higgs_source_indices':[0,1,2,3],
           'heavy_source_indices':list(range(4,4+n)),
           'synthetic_fixture':True,'wilson_adm_certified':False,
           'signed_GammaCC_to_Weyl_certified':False,
           'wilson_dim5_closure_certified':False,
           'strict_1loop_endpoint_at_ms':True,
           'strict_1loop_origin':'SYNTHETIC_A_full_one_loop_truncated_EFT1_at_MS'}
    bmeta={'format':B_FORMAT,'tensor_npz':'background.npz','model':model,
           'matching_context_id':label,'matched_seed_id':label,
           'lagrangian_source_id':'SYNTHETIC_NOT_ACTUAL_LAGRANGIAN',
           'mu_f_gev':mu_f,'mu_s_gev':mu_s,'units':'GeV','scheme':'MSbar',
           'quartic_potential_convention':Q_CONVENTION,'quartic_inventory_status':'full_tensor',
           'wilson_scalar_basis_id':model_basis,
           'background_scalar_basis_id':f'{model}:canonical_interleaved',
           'scalar_background_real_labels':scalar_labels(model),
           'single_common_scalar_threshold':True,'synthetic_fixture':True,
           'scalar_mass_origin':'synthetic_algebra_only',
           'scalar_operator_inventory_certified':False,
           'native_matchete_rgbeta_normalization_certified':False,
           'quartic_channel_names':['synthetic_hhss','synthetic_other_quartics'],
           'sm_couplings_at_ms':{'gY':.359,'g2':.646,'g3':1.08,'lambdaH':.256},
           'mass_parameter_info':mass.extra}
    _write_npz(loc/'wilson.json',wmeta,{
        'C_psi2phi2_ms':Cfull,'C_psi2phi2_tree_mf':Ctree,
        'C_psi2phi2_ms_strict_1loop':Cfull.copy(),
        'C5_EFT1_ms':K5})
    _write_npz(loc/'background.json',bmeta,{
        'mass_squared_gev2_ms':M,
        'R_wilson_from_background':R,
        'quartic_real_full_ms':fullq,
        'quartic_real_channels_ms':np.stack([portal,baseline]),
        'ye':yc,'yu':yu,'yd':yd})
    return loc/'wilson.json',loc/'background.json'


def make_all_fixtures(directory:Path) -> dict[str,tuple[Path,Path]]:
    return {model:make_fixture(model,directory) for model in ALL_MODELS}
