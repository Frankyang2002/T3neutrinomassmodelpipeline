"""Fermion-first rough numerical adapters; never evolve Wilsons or match scalars.

Read ONLY B's saved, checksum-protected UV/EFT1 trajectories and F-tree boundary.
Each production output records any non-source-certified coordinate convention.
The intentionally selected reduced potential is a *projection* of B's numerical
trajectory, NOT the full original RGBeta scalar potential or a re-solved RGE.
"""
from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
from itertools import permutations
import json
from pathlib import Path
from typing import Mapping

import numpy as np

from Numerical.parallel_t3.background.api import Background, load_formal_boundary
from Numerical.parallel_t3.background.contracts import MODELS, hash_file, read_json, tag
from Numerical.parallel_t3.wilson.basis import build_basis, spin_generators
from Numerical.parallel_t3.wilson.coordinates import embed_c5, embed_formal_llss, symmetry_residual
from Numerical.parallel_t3.threshold.models import from_parameters, special_kind
from Numerical.parallel_t3.wilson.interchange import BASIS_ID

B_RUN = 'B_all16_qualitative_v1'
RUN_ID = 'rough_fermion_first_v1'
A_FORMAT = 't3-parallel-a-input-v1'
C_FORMAT = 't3-parallel-c-background-input-v1'
W_CONVENTION = 'C_ijab_full_sum_no_factorial_C5_1over2_H_R_interleaved_v1'
Q_CONVENTION = 'V4_lambda_real_symmetric_over_4factorial_v1'
C_SCALAR_ID = 'T3_RGBeta_reference_Weyl45_real_interleaved_v1'
MF = 1e5
MS = 1e3

class InputBlocked(ValueError):
    """No invented input or unrecorded coordinate mapping is permitted."""


def _save_json(path: Path, data: Mapping) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    p = path.with_suffix(path.suffix + '.partial')
    p.write_text(json.dumps(data, indent=2, sort_keys=True, allow_nan=False) + '\n', encoding='utf-8')
    p.replace(path)


def _save_npz(path: Path, arrays: dict) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    p = path.with_suffix('.partial.npz')
    np.savez_compressed(p, **arrays)
    p.replace(path)
    return hash_file(path)


def _nonzero(x: np.ndarray, tolerance=1e-15) -> bool:
    return bool(np.max(np.abs(x)) > tolerance)


def _real(x: np.ndarray | complex, name: str, tol=1e-9) -> float:
    v = complex(np.asarray(x).item())
    if not np.isfinite([v.real,v.imag]).all() or abs(v.imag)>tol*max(1.,abs(v.real)):
        raise InputBlocked(f'{name}: expected finite real RGBeta running coupling, got {v}')
    return v.real


def _complex_json(z: complex) -> dict:
    return {'real': float(np.real(z)), 'imag': float(np.imag(z))}


def _source_ids(bg: dict, boundary: dict, folder: Path) -> dict:
    src = bg.get('input_provenance',{}).get('local_sources',{})
    sources = {name: value['sha256'] for name,value in src.items()
               if isinstance(value,dict) and value.get('required_for_B_reuse') is True
               and isinstance(value.get('sha256'),str)}
    if not sources: raise InputBlocked('B background lacks actual source SHA256 provenance')
    return {'source_sha256':sources,
            'matching_context_id': ('B:'+bg['UV_endpoint_sha256']+':'+boundary['source_signature_sha256']),
            'matched_seed_id':boundary['source_matchete_seed_sha256'],
            'lagrangian_source_id': sources.get('Lagrangian/interactions/ScalarPotential.wl'),
            'B_background_manifest_sha256':hash_file(folder/'background_manifest.json'),
            'B_boundary_manifest_sha256':hash_file(folder/'fermion_boundary_manifest.json'),
            'B_background_npz_sha256':hash_file(folder/'background.npz'),
            'B_boundary_npz_sha256':hash_file(folder/'fermion_boundary_formal.npz')}


def _u(ns: int, offset: int, index: int, barred: bool=False) -> np.ndarray:
    """Complex component as a linear form in interleaved real fields."""
    a = np.zeros(ns, complex)
    a[offset+2*index] = 1/np.sqrt(2)
    a[offset+2*index+1] = (-1j if barred else 1j)/np.sqrt(2)
    return a


def _bilinear(ns: int, offset1: int, offset2: int, M: np.ndarray) -> np.ndarray:
    """z1^dag M z2 = phi.T @ B @ phi (B symmetric, possibly complex)."""
    n1,n2 = M.shape
    B = np.zeros((ns,ns),complex)
    for i in range(n1):
        a=_u(ns,offset1,i,True)
        for j in range(n2):
            if M[i,j] != 0:
                b=_u(ns,offset2,j)
                B += M[i,j]*np.outer(a,b)
    return (B+B.T)/2


def _product(A:np.ndarray, B:np.ndarray) -> np.ndarray:
    """Fourth derivative of (phi.T A phi)(phi.T B phi)."""
    out=np.zeros(A.shape*2,complex)
    for a,b,c,d in ((0,1,2,3),(0,2,1,3),(0,3,1,2),
                    (2,3,0,1),(1,3,0,2),(1,2,0,3)):
        # assign A to two positions and B to the other two
        p=np.einsum('ij,kl->ijkl',A,B,optimize=True)
        out+=4*p.transpose(np.argsort((a,b,c,d)))
    return out


def _monomial(a:np.ndarray,b:np.ndarray,c:np.ndarray,d:np.ndarray)->np.ndarray:
    """Fourth derivative of a product of four scalar linear forms."""
    raw=np.einsum('a,b,c,d->abcd',a,b,c,d,optimize=True)
    out=np.zeros_like(raw)
    for perm in permutations(range(4)):
        out+=raw.transpose(perm)
    return out


def _cg() -> np.ndarray:
    # Explicit CONDITIONAL orthonormal symmetric-spinor SU(2) frame.
    g=np.zeros((3,2,2))
    g[0,0,0]=1
    g[1,0,1]=g[1,1,0]=1/np.sqrt(2)
    g[2,1,1]=1
    return g


def _triple(g:np.ndarray)->np.ndarray:
    # Source T3RGBetaS2TripleInvariant from T3RGBetaModel.wl lines 462-470.
    eps=np.array([[0,1],[-1,0]],float)
    return np.einsum('aij,bkl,cmn,ik,jm,ln->abc',
                     g,g,g,eps,eps,eps,optimize=True)


def quartic_unit_tensors(model: str) -> dict[str,np.ndarray]:
    """Explicit source-basic and adjoint RGBeta invariants + conditional T3 portal.

    These component maps do NOT certify the Matchete or RGBeta native CG frame.
    Unsupported further native invariants must be selected as explicit imposed zeros.
    """
    basis=build_basis(model)
    d1,d2,_=basis.dims
    ns=basis.ns
    offsets={'H':0,'S1':4,'S2':4+2*d1}
    dims={'H':2,'S1':d1,'S2':d2}
    labels={}
    identity={key:_bilinear(ns,offsets[key],offsets[key],np.eye(dim))
              for key,dim in dims.items()}
    for key,kind in [('lambdaH','H'),('lambdaS1','S1'),('lambdaS2','S2')]:
        labels[key]=_product(identity[kind],identity[kind]).real/2
    for key,a,b in [('lambdaH1','H','S1'),('lambdaH2','H','S2'),
                    ('lambda12','S1','S2')]:
        labels[key]=_product(identity[a],identity[b]).real
    generators={key:spin_generators(dim) for key,dim in dims.items()}
    for key,a,b,factor in [('lambdaH1Adj','H','S1',1),
                           ('lambdaH2Adj','H','S2',1),
                           ('lambdaS1Adj','S1','S1',.5),
                           ('lambdaS2Adj','S2','S2',.5),
                           ('lambda12Adj','S1','S2',1)]:
        if key.startswith('lambdaH') and dims[b]==1: continue
        if key.startswith('lambdaS') and dims[a]!=3: continue
        if key=='lambda12Adj' and min(d1,d2)==1: continue
        val=np.zeros((ns,)*4,complex)
        for G1,G2 in zip(generators[a],generators[b]):
            val += _product(_bilinear(ns,offsets[a],offsets[a],G1),
                            _bilinear(ns,offsets[b],offsets[b],G2))
        labels[key]=factor*val.real
    if d1==3 and d2==3:
        val=np.zeros((ns,)*4,complex)
        for i in range(3):
            for j in range(3):
                # del[i,k]del[j,l]+del[i,l]del[j,k]
                val+=_product(_bilinear(ns,offsets['S1'],offsets['S1'],
                                          np.eye(3)[i:i+1].T@np.eye(3)[j:j+1]),
                              _bilinear(ns,offsets['S2'],offsets['S2'],
                                          (np.eye(3)[i:i+1].T@np.eye(3)[j:j+1]
                                           +np.eye(3)[j:j+1].T@np.eye(3)[i:i+1])))
        labels['lambda12Cross']=val.real
    g=_cg()
    if (d1,d2)==(1,3):
        T=np.einsum('aij->ija',g).reshape(2,2,1,3)
    elif (d1,d2)==(3,1):
        T=np.einsum('aij->ija',g).reshape(2,2,3,1)
    elif (d1,d2)==(2,2):
        T=np.einsum('aij,akl->ijkl',g,g,optimize=True)
    elif (d1,d2)==(3,3):
        T=np.einsum('aij,abc->ijbc',g,_triple(g),optimize=True)
    else:
        raise InputBlocked('No source-declared lambdaT3 invariant for '+model)
    lam=np.zeros((ns,)*4,complex)
    for i,j,k,l in np.argwhere(np.abs(T)>1e-14):
        lam += T[i,j,k,l] * _monomial(_u(ns,0,int(i)),_u(ns,0,int(j)),
                  _u(ns,4,int(k)),_u(ns,4+2*d1,int(l),True))
    # lambdaT3 is complex and the source has explicitly + h.c.
    labels['lambdaT3']=lam
    return labels


def assemble_quartic(units: Mapping[str,np.ndarray],
                     state: Mapping[str, np.ndarray], *, reduced: bool) -> tuple[np.ndarray,dict]:
    ns=next(iter(units.values())).shape[0]
    lam=np.zeros((ns,)*4,float)
    omissions={}
    for key,val in state.items():
        if not key.startswith('lambda'): continue
        z=np.asarray(val).item()
        if key in units:
            if key=='lambdaT3':
                lam+=2*np.real(z*units[key])
            else:
                coeff=_real(z,key)
                lam+=coeff*units[key]
        else:
            if _nonzero(np.asarray(z)):
                if not reduced:
                    raise InputBlocked(f'Unmapped nonzero native quartic {key}={z}; select --reduced-potential')
                omissions[key]=_complex_json(z)
            else:
                omissions[key]=_complex_json(z)
    if not np.isfinite(lam).all():raise InputBlocked('Nonfinite scalar quartic tensor')
    # Up to floating-point round-off, all four real indices are symmetric.
    for p in ((1,0,2,3),(0,2,1,3),(0,1,3,2)):
        if not np.allclose(lam,lam.transpose(p),rtol=1e-10,atol=1e-10):
            raise InputBlocked('Constructed quartic tensor is not fully symmetric')
    return lam,omissions


def sm_yukawa_tensor(basis, yu: np.ndarray,yd: np.ndarray,ye: np.ndarray)->np.ndarray:
    """Conditional convention -YQuc(Q.eps.H)-YQdc(Q.H*)-YLec(L.H*) +h.c.

    Returned Y symmetric in fermion indices, L_Y=-Y_ija psi_i psi_j phi_a/2.
    Exact overall UV RGBeta-to-Weyl signed relation is NOT independently fixed.
    """
    ns=basis.ns
    y=np.zeros((basis.nf,basis.nf,ns),complex)
    f={lab:i for i,lab in enumerate(basis.fermion_labels)}
    eps=((0,1),(-1,0))
    for p in range(3):
        for q in range(3):
            for w in range(2):
                for col in range(1,4):
                    Q=f[f'Q[{p+1},{col},{w+1}]']
                    u=f[f'uc[{q+1},{col}]']
                    d=f[f'dc[{q+1},{col}]']
                    # epsilon Q_i H_j, epsilon[1,2]=+1
                    for v in range(2):
                        vcoef=eps[w][v]
                        if vcoef:
                            for ri in range(2):
                                h=(1 if ri==0 else 1j)/np.sqrt(2)
                                y[Q,u,2*v+ri]+=yu[p,q]*vcoef*h
                    for ri in range(2):
                        h=(1 if ri==0 else -1j)/np.sqrt(2)
                        y[Q,d,2*w+ri]+=yd[p,q]*h
                L=f[f'L[{p+1},{w+1}]']
                e=f[f'ec[{q+1}]']
                for ri in range(2):
                    h=(1 if ri==0 else -1j)/np.sqrt(2)
                    y[L,e,2*w+ri]+=ye[p,q]*h
    return y+y.swapaxes(0,1)


@dataclass
class Options:
    root: Path
    B_run_id: str = B_RUN
    out_run_id: str = RUN_ID
    gamma_factor: float = 1.
    assume_formal_identity: bool = False
    assume_sm_weyl_yukawa: bool = False
    reduced_potential: bool = False
    hard_c5_zero_benchmark: bool = False
    hard_c5_root: Path | None = None
    overwrite: bool = False


def _c5_hard(folder:Path, opt:Options, arrays:dict, meta:dict):
    if opt.hard_c5_root is not None:
        if opt.hard_c5_zero_benchmark:
            raise InputBlocked('Do not select both an external and zero hard-C5 input')
        fp=opt.root/opt.hard_c5_root/tag(meta['model'])/'hard_C5.npz'
        fmeta=fp.with_name('hard_C5_manifest.json')
        d=read_json(fmeta)
        if (d.get('model')!=meta['model'] or d.get('mu_F_gev')!=MF or
                d.get('same_run_uv_state_sha256')!=meta['UV_endpoint_sha256'] or
                d.get('hard_c5_file_sha256')!=hash_file(fp)):
            raise InputBlocked('Supplied hard C5 metadata not same-run or hash-invalid')
        with np.load(fp,allow_pickle=False) as z:c5=z['C5_hard'].astype(complex)
        origin='supplied_independent_file'
        evidence={'file':str(fp.relative_to(opt.root)),'sha256':hash_file(fp)}
    elif 'C5_hard_external' in arrays:
        c5=np.asarray(arrays['C5_hard_external'],complex)
        origin='B_supplied_external';evidence={'boundary_manifest':str(folder/'fermion_boundary_manifest.json')}
    elif opt.hard_c5_zero_benchmark:
        c5=np.zeros((3,3),complex)
        origin='explicit_zero_benchmark_not_physical_vanishing';evidence={}
    else:
        raise InputBlocked('Hard C5 is unknown: use --hard-c5-zero-benchmark or --hard-c5-root')
    if (c5.shape!=(3,3) or not np.isfinite(c5).all()
            or not np.allclose(c5,c5.T,atol=1e-23,rtol=1e-10)):
        raise InputBlocked('Hard C5 must be symmetric finite complex 3x3')
    return c5,origin,evidence


def _source_hash_check(root:Path, provenance:dict)->None:
    # B's provenance is already immutable in B file manifests; detect source drift.
    # Only required-for-B files are hard gates (Wilson-only may have changed).
    for rel,info in provenance.get('local_sources',{}).items():
        if info.get('required_for_B_reuse',False):
            path=root/rel
            if not path.is_file() or hash_file(path)!=info['sha256']:
                raise InputBlocked(f'Source changed from saved B input: {rel}')


def _scalar_mass(model:str,state:dict):
    a=_real(state['mS1Sq'],'mS1Sq')
    b=_real(state['mS2Sq'],'mS2Sq')
    kind=special_kind(model)
    if kind=='ordinary': return from_parameters(model,a,b)
    mu_names=sorted(n for n in state if n.startswith('muT3'))
    if len(mu_names)!=1:raise InputBlocked(f'{model}: expected one complex mixing/splitting input, got {mu_names}')
    mu=complex(np.asarray(state[mu_names[0]]).item())
    return from_parameters(model,a,b,mu_sq=mu)


def _a_manifest(model:str,basis,provenance:dict,path:Path,kind:str,options:Options,
                *, npz_sha256:str, is_boundary:bool, hard_origin:str | None=None)->dict:
    data={'schema':A_FORMAT,'format':A_FORMAT,'model':model,'basis_id':BASIS_ID,
          'fermion_labels':list(basis.fermion_labels),
          'scalar_labels':list(basis.scalar_labels),
          'full_SM_three_generation_Weyl_basis':True,
          'MF_gev':MF,'MS_gev':MS,'npz_sha256':npz_sha256,
          'source_sha256':provenance['source_sha256'],
          'operator_basis':'psi2phi2_no_factorial',
          'numeric_input_kind':'saved_B_projected_reduced_potential' if options.reduced_potential else 'saved_B_full_native_mapping',
          'gamma_cc_to_weyl_signed_factor':options.gamma_factor,
          'gamma_cc_to_weyl_factor_independently_certified':False,
          'Matchete_RGBeta_scalar_map_certified':False,
          'sm_weyl_yukawa_assumed':options.assume_sm_weyl_yukawa,
          'run_id':options.out_run_id,'B_run_id':options.B_run_id,
          'source_B_manifest_sha256':provenance['B_background_manifest_sha256'],
          'matching_context_id':provenance['matching_context_id'],
          'matched_seed_id':provenance['matched_seed_id'],
          'lagrangian_source_id':provenance['lagrangian_source_id'],
          'scalar_basis_id':model+':canonical_interleaved',
          'scheme':'MSbar','units':'GeV'}
    if is_boundary:
        data.update({'C5_hard_origin':hard_origin,'C_full_is_canonical':False,
                    'native_source_path':str(path).replace('\\','/'),
                    'corrected_Aminus4_portal_evidence':model=='A_-4'})
    return data


def build_model(model:str, opt:Options)->dict:
    if model not in MODELS:raise InputBlocked('Not an all-16 production model: '+model)
    if not opt.assume_formal_identity:
        raise InputBlocked('No certified formal scalar map: explicit --assume-formal-identity needed')
    if not opt.assume_sm_weyl_yukawa:
        raise InputBlocked('No certified SM native-to-Weyl Yukawa map: explicit --assume-sm-weyl-yukawa needed')
    if not np.isfinite(opt.gamma_factor):raise InputBlocked('Nonfinite GammaCC factor')
    root=Path(opt.root)
    folder=root/'output/parallel_t3/background'/opt.B_run_id/tag(model)
    # B API verifies checksum, UV/EFT1 inventories, MF continuity, finite data.
    bg=Background.load(folder)
    boundary, fmeta=load_formal_boundary(folder)
    bmeta=bg.manifest
    if fmeta.get('model')!=model or bmeta.get('model')!=model:
        raise InputBlocked('Model mismatch across B files')
    if (not fmeta.get('all_three_heavy_fermions_decoupled_at_single_MF',False)
            or fmeta.get('same_run_uv_state_sha256')!=bmeta.get('UV_endpoint_sha256')):
        raise InputBlocked('All-fermion common-MF or same-state provenance mismatch')
    if (bmeta['input_scales_gev']['F']!=MF or bmeta['input_scales_gev']['S']!=MS or
           float(fmeta.get('mu_F_gev',-1))!=MF):
        raise InputBlocked('B save has wrong MF/MS')
    _source_hash_check(root,bmeta.get('input_provenance',{}))
    seed_rel=str(fmeta.get('source_matchete_seed','')).replace('\\','/')
    if not seed_rel or Path(seed_rel).is_absolute() or '..' in Path(seed_rel).parts:
        raise InputBlocked('Invalid or absent Matchete fermion source path')
    seed=root/seed_rel
    if not seed.is_file() or hash_file(seed)!=fmeta.get('source_matchete_seed_sha256'):
        raise InputBlocked('Saved Matchete F-boundary source is missing or SHA-mismatched: '+seed_rel)
    if model=='A_-4' and not seed_rel.startswith(
            'output/general_observable/fermion_first_hermitian_portal_v1/A_m4/'):
        raise InputBlocked('A(-4): source is not the corrected Hermitian portal seed')
    if model in ('B_-1','C_-1') and fmeta.get('mixed_scalar_mass_matrix_scale')!='MF':
        raise InputBlocked('B/C(-1): formal X/Y source was not constructed at common MF')
    basis=build_basis(model)
    nh=basis.ns-4
    raw=boundary['C_LLSS_formal']
    if raw.shape!=(3,3,2,2,nh,nh):raise InputBlocked('Wrong F-boundary real scalar geometry')
    c5, hard_origin,hard_evidence=_c5_hard(folder,opt,boundary,bmeta)
    # Formal-source scalar slots mapped by explicitly selected identity convention.
    c=embed_formal_llss(raw,basis,np.eye(nh),opt.gamma_factor) + embed_c5(c5,basis)
    if symmetry_residual(c)>max(1e-24,1e-10*np.max(np.abs(c))):
        raise InputBlocked('Full canonical Wilson symmetry mismatch')
    t=np.asarray(bg.arrays['mu_eft1_gev'],float)
    nstep=len(t)
    units=quartic_unit_tensors(model)
    y=np.empty((nstep,basis.nf,basis.nf,basis.ns),complex)
    lam=np.empty((nstep,)+(basis.ns,)*4,float)
    omitted={}
    for k in range(nstep):
        state={name:bg.arrays['eft1__'+name][k] for name in bg.manifest['beta_symbol_order_eft1']}
        y[k]=sm_yukawa_tensor(basis,state['yu'],state['yd'],state['ye'])
        lam[k],lost=assemble_quartic(units,state,reduced=opt.reduced_potential)
        for name,z in lost.items():
            item=omitted.setdefault(name,{'native_at_mf':None,'native_at_ms':None,
                                          'max_abs_over_saved_trajectory':0.,
                                          'imposed_value_in_tensor':_complex_json(0j)})
            if k==0:item['native_at_mf']=z
            if k==nstep-1:item['native_at_ms']=z
            item['max_abs_over_saved_trajectory']=max(item['max_abs_over_saved_trajectory'],
                                                        abs(complex(state[name].item())))
    msstate={name:bg.arrays['eft1__'+name][-1] for name in bg.manifest['beta_symbol_order_eft1']}
    adapter=_scalar_mass(model,msstate)
    # In our selected source=canonical=background real frame, C requires identity.
    R=np.eye(nh)
    mass=np.asarray(adapter.mass2_canonical_gev2,float)
    if mass.shape!=(nh,nh) or not np.allclose(mass,mass.T,atol=1e-9):
        raise InputBlocked('Scalar Hessian mismatch')
    mu=bg.arrays['mu_eft1_gev']; src=_source_ids(bmeta,fmeta,folder)
    out=root/'output/parallel_t3/rough_inputs'/opt.out_run_id
    a=out/'a'/tag(model);cb=out/'c'/tag(model)
    if not opt.overwrite and (a/'boundary.npz').exists():
        raise InputBlocked('Output already exists; use --overwrite to replace this rough run')
    # Store exact B native trajectories and history in separate names, no re-solving.
    native={name:bg.arrays['eft1__'+name] for name in bg.manifest['beta_symbol_order_eft1']}
    source={'mu_uv_gev':bg.arrays['mu_uv_gev'],
            'mu_eft1_gev':bg.arrays['mu_eft1_gev'],
            't_uv':bg.arrays['t_uv'],'t_eft1':bg.arrays['t_eft1']}
    for name in bg.manifest['beta_symbol_order_uv']:
        source['uv__'+name]=bg.arrays['uv__'+name]
    for name in bg.manifest['beta_symbol_order_eft1']:
        source['eft1__'+name]=bg.arrays['eft1__'+name]
    ssha=_save_npz(a/'saved_native_trajectory.npz',source)
    a_sha=_save_npz(a/'background.npz',{
        'mu_gev':t,'Y':y,'lambda4':lam,'gY':np.real(native['gY']),
        'g2':np.real(native['g2']),'g3':np.real(native['g3'])})
    b_sha=_save_npz(a/'boundary.npz',{
        'C_full':c,'C5_hard':c5,'C_LLSS_formal_muF':raw,
        'heavy_source_from_canonical':np.eye(nh)})
    # A interchange.load_boundary can read C_full ONLY when canonical certification
    # is claimed. Export the *formal* path for that loader; keep C_full as extra
    # for direct consumers, and mark it provisional rather than certified.
    # Remove C_full from the loader-bound file; retain full tensor in C_full_provisional.
    b_sha=_save_npz(a/'boundary.npz',{
        'C_full_provisional':c, 'C5_hard':c5,
        'C_LLSS_formal_muF':raw, 'heavy_source_from_canonical':np.eye(nh)})
    a_bg=_a_manifest(model,basis,src,folder/'background_manifest.json','background',opt,npz_sha256=a_sha,is_boundary=False)
    a_bd=_a_manifest(model,basis,src,folder/'fermion_boundary_manifest.json','boundary',opt,npz_sha256=b_sha,is_boundary=True,hard_origin=hard_origin)
    a_bd['gamma_cc_to_weyl_signed_factor']=float(opt.gamma_factor)
    a_bd['operator_basis']='psi2phi2_no_factorial'
    a_bd['canonical_provisional_tensor_npz_key']='C_full_provisional'
    a_bd['native_source_path']=str(fmeta['source_matchete_seed']).replace('\\','/')
    a_bd['source_sha256'][str(fmeta['source_matchete_seed']).replace('\\','/')]=fmeta['source_matchete_seed_sha256']
    a_bd['corrected_Aminus4_portal_evidence']=(model=='A_-4')
    _save_json(a/'background.json',a_bg)
    _save_json(a/'boundary.json',a_bd)
    _save_json(a/'assumptions.json',{
        'model':model,'run_id':opt.out_run_id,'original_saved_B_trajectories_untouched':True,
        'no_fermions_below_MF':True,'both_scalars_retained_to_MS':True,
        'full_45_Weyl_tensor_created':True,
        'GammaCC_to_Weyl_signed_factor':opt.gamma_factor,'factor_verified':False,
        'formal_scalar_index_map':'identity_selected_not_certified',
        'SM_Yukawa_Weyl_embedding':'explicit_provisional_convention',
        'scalar_quartic_frame':'RGBeta declared SU2 conditional CG',
        'reduced_potential':opt.reduced_potential,
        'omitted_native_nonzero_quartic_coordinates':omitted,
        'reduced_potential_trajectory_status':'projection_of_full_RGBeta_running_not_reintegrated_reduced_RGEs',
        'hard_C5':{'origin':hard_origin,'evidence':hard_evidence,'physical_zero_claim':False},
        'mass_basis_phase_and_native_scalar_map_certified':False,
        'model_mass_kind':adapter.kind,
        'source_B_sha256':src})
    # Workstream C: heavy-only Hessian; H included only in quartic tensor.
    c_sha=_save_npz(cb/'background.npz',{
        'mass_squared_gev2_ms':mass,'R_wilson_from_background':R,
        'quartic_real_full_ms':lam[-1],
        'ye':np.asarray(msstate['ye'],complex),'yu':np.asarray(msstate['yu'],complex),
        'yd':np.asarray(msstate['yd'],complex)})
    sm={key:_real(msstate[key],key) for key in ('gY','g2','g3','lambdaH')}
    cmeta={'format':C_FORMAT,'tensor_npz':'background.npz','tensor_npz_sha256':c_sha,
        'model':model,'mu_f_gev':MF,'mu_s_gev':MS,'scheme':'MSbar','units':'GeV',
        'matching_context_id':src['matching_context_id'],
        'matched_seed_id':src['matched_seed_id'],
        'lagrangian_source_id':src['lagrangian_source_id'],
        'quartic_potential_convention':Q_CONVENTION,
        'quartic_inventory_status':'full_tensor',
        'original_RGBeta_potential_fully_mapped':not bool(omitted),
        'reduced_potential_explicit_zeros':bool(omitted),
        'wilson_scalar_basis_id':model+':canonical_interleaved',
        'background_scalar_basis_id':model+':canonical_interleaved_provisional',
        'scalar_background_real_labels':list(basis.scalar_labels[4:]),
        'single_common_scalar_threshold':True,
        'scalar_mass_origin':'unbroken_electroweak','synthetic_fixture':False,
        'sm_couplings_at_ms':sm, 'quartic_channel_names':[],
        'scalar_operator_inventory_certified':False,
        'native_matchete_rgbeta_normalization_certified':False,
        'gauge_ward_certified':False,
        'omitted_native_quartic_coordinates':omitted,
        'scalar_source_coordinate_assumption':'provisional Matchete/RGBeta identity',
        'mass_hessian_source':'threshold.models.from_parameters plus actual B MS quadratic coordinates',
        'scalar_eigenvalues_squared_gev2':np.linalg.eigvalsh(mass).tolist(),
        'log_scalar_mass2_over_MS2':np.log(np.linalg.eigvalsh(mass)/MS**2).tolist(),
        'common_MS_not_individual_eigenmass':True,
        'retain_mass_dependent_matching_logs':True,
        'source_B_sha256':src,
        'A_boundary_manifest':'a/'+tag(model)+'/boundary.json'}
    _save_json(cb/'background.json',cmeta)
    return {'model':model,'slug':tag(model),
            'status':'written_provisional_A_C_inputs',
            'a_background':str((a/'background.json').relative_to(root)),
            'a_boundary':str((a/'boundary.json').relative_to(root)),
            'c_background':str((cb/'background.json').relative_to(root)),
            'native_uv_samples':int(len(bg.arrays['mu_uv_gev'])),
            'native_eft1_samples':nstep,
            'canonical_wilson_shape':list(c.shape),
            'quartic_shape':list(lam.shape),
            'omitted_nonzero_running_quartic_names':sorted(n for n,item in omitted.items()
                                if item['max_abs_over_saved_trajectory']>1e-15),
            'scalar_mass_eigenvalues_squared_gev2':np.linalg.eigvalsh(mass).tolist(),
            'source_ids':src}


def run_batch(models:list[str],opt:Options)->dict:
    root=Path(opt.root)
    if not (root/'Numerical/parallel_t3/background/api.py').exists():
        raise InputBlocked('Run from the actual repository root, with Workstream B installed')
    output=root/'output/parallel_t3/rough_inputs'/opt.out_run_id
    results={}
    for model in models:
        try:
            results[model]=build_model(model,opt)
            print(f'{model}: A/C rough inputs written')
        except (ValueError,KeyError,TypeError,FileNotFoundError,OSError,RuntimeError) as e:
            results[model]={'model':model,'status':'blocked','blocker':type(e).__name__+': '+str(e)}
            print(f'{model}: BLOCKED: {e}')
    summary={'schema':'t3-parallel-rough-stream1-summary-v1',
        'run_id':opt.out_run_id,'B_run_id':opt.B_run_id,
        'MF_gev':MF,'MS_gev':MS,
        'requested_models':models,
        'count_written':sum(v['status']=='written_provisional_A_C_inputs' for v in results.values()),
        'count_blocked':sum(v['status']=='blocked' for v in results.values()),
        'models':results,
        'physical_normalization_certified':False,
        'extra_fermion_matching_performed':False,
        'Wilson_ODE_performed':False,'scalar_matching_performed':False,
        'EFT_order':'UV -> integrate out ALL fermions at MF -> retain BOTH scalars -> integrate out BOTH at MS -> SMEFT',
        'explicit_reduced_potential':opt.reduced_potential,
        'GammaCC_factor_benchmark':opt.gamma_factor}
    _save_json(output/'summary.json',summary)
    return summary
