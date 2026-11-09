"""Stream 2: canonical fermion-first EFT1 -> coupled full Wilson -> C endpoint.

Additive adapter over the existing 45-Weyl Wilson solver. This module neither
integrates out scalars nor inserts a separate LLSS->C5 logarithm. It does not
convert undocumented native Matchete/RGBeta coordinates by assumption.

Expected directory, per model:
  output/parallel_t3/rough_inputs/rough_fermion_first_v1/a/<slug>/
Either (background.npz, boundary.npz + matching JSON manifests) or a single
canonical_input.npz and canonical_input.json. See README_ROUGH_STREAM2.md.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any
import json
import math
import subprocess
import shutil
import tempfile

import numpy as np

from .adm import SYMMETRIC_TERMS
from .basis import ALL_MODELS, ModelBasis, build_basis
from .coordinates import (embed_c5, embed_formal_llss, project_c5,
                          symmetry_residual)
from .interchange import Background, BlockedInput, BASIS_ID, sha256
from .workflow import run_wilson

INPUT_ROOT = Path('output/parallel_t3/rough_inputs/rough_fermion_first_v1/a')
OUTPUT_ROOT = Path('output/parallel_t3/rough_wilson/rough_fermion_first_v1')
CONVENTION = 'C_ijab_full_sum_no_factorial_C5_1over2_H_R_interleaved_v1'
Q_CONVENTION = 'V4_lambda_real_symmetric_over_4factorial_v1'
W_FORMAT = 't3-parallel-c-wilson-input-v1'
BASIS_FORMAT = 't3-parallel-a-input-v1'
REQUIRED_IDS = ('matching_context_id', 'matched_seed_id', 'lagrangian_source_id',
                'scalar_basis_id')
ALLOWED_ADMS = ('fonseca_eq22_six_term', 'ward_covariant_candidate')
LOOP = 1/(16*math.pi**2)


def slug(model: str) -> str:
    if model not in ALL_MODELS:
        raise BlockedInput('Unsupported model ' + model)
    fam, alpha = model.split('_')
    return f'{fam}_{"m" if int(alpha)<0 else "p"}{abs(int(alpha))}'


def _json(path: Path) -> dict:
    obj = json.loads(path.read_text(encoding='utf-8-sig'))
    if not isinstance(obj, dict):
        raise BlockedInput(f'Expected JSON object: {path}')
    return obj


def _find(folder: Path, names: tuple[str, ...], *, required=True) -> Path | None:
    matches = [folder/name for name in names if (folder/name).is_file()]
    if len(matches)>1:
        raise BlockedInput(f'Ambiguous input manifests in {folder}: {[m.name for m in matches]}')
    if not matches and required:
        raise BlockedInput(f'Missing {names} in {folder}')
    return matches[0] if matches else None


def _npz(path: Path) -> dict[str, np.ndarray]:
    with np.load(path, allow_pickle=False) as data:
        return {k:np.asarray(data[k]).copy() for k in data.files}


def _hash_checked(path: Path, meta: dict, *, project: Path) -> None:
    digest = sha256(path)
    for key in ('tensor_npz_sha256', 'npz_sha256',
                'boundary_npz_sha256', 'background_npz_sha256'):
        value = meta.get(key)
        if value is not None and value != digest:
            raise BlockedInput(f'{key} mismatch: {path}')
    declared = meta.get('tensor_npz')
    if declared is not None and (not isinstance(declared,str) or declared!=path.name):
        raise BlockedInput('Manifest tensor_npz does not name adjacent file')
    refs = meta.get('source_sha256', {})
    if refs and not isinstance(refs, dict):
        raise BlockedInput('source_sha256 must be a dictionary')
    for name, expect in refs.items():
        rel = Path(name)
        if rel.is_absolute() or '..' in rel.parts:
            raise BlockedInput('Unsafe source_sha256 relative path')
        source = project / rel
        if not source.is_file() or sha256(source) != expect:
            raise BlockedInput('Original source changed or unavailable: '+name)


def _same_value(a: dict, b: dict, key: str) -> Any:
    av,bv = a.get(key), b.get(key)
    if av is not None and bv is not None and av != bv:
        raise BlockedInput(f'Incompatible Stream1 {key}: boundary vs background')
    value = av if av is not None else bv
    if value is None:
        raise BlockedInput(f'Missing Stream1 {key}; cannot invent Workstream C context ID')
    return value


def _float_meta(a: dict, b: dict, low: str, upper: str) -> float:
    values = [m[k] for m in (a,b) for k in (low,upper) if k in m]
    if not values:
        raise BlockedInput('Missing threshold '+low)
    x=float(values[0])
    if not math.isfinite(x) or any(abs(float(v)-x)>1e-10*max(x,1) for v in values):
        raise BlockedInput('Disagreeing/nonfinite scale '+low)
    return x


def _validate_meta(basis: ModelBasis, bm: dict, mm: dict, *, mf: float, ms: float) -> dict:
    if not 0 < ms < mf:
        raise BlockedInput('Threshold order must be MF > MS > 0')
    if abs(mf-1e5)>1e-6 or abs(ms-1e3)>1e-6:
        raise BlockedInput('Rough Stream2 requires MF=1e5 GeV and MS=1e3 GeV')
    for key,value in [('model',basis.model),('basis_id',BASIS_ID)]:
        for meta in (bm,mm):
            if key in meta and meta[key]!=value:
                raise BlockedInput(f'Canonical {key} mismatch')
    # An explicit, canonical common field frame must be shared by both inputs.
    for label,expected in (('fermion_labels',basis.fermion_labels),
                           ('scalar_labels',basis.scalar_labels)):
        if not any(label in m for m in (mm,bm)):
            raise BlockedInput('Missing explicit canonical '+label)
        for m in (mm,bm):
            if label in m and tuple(m[label])!=expected:
                raise BlockedInput('Noncanonical '+label+'; requires explicit transform before Stream2')
    if not any(m.get('full_SM_three_generation_Weyl_basis') is True for m in (mm,bm)):
        raise BlockedInput('Full 45-Weyl canonical metadata not declared')
    for m in (mm,bm):
        if 'scheme' in m and m['scheme']!='MSbar':
            raise BlockedInput('Only MSbar background permitted')
        if 'synthetic_fixture' in m and m['synthetic_fixture'] is True:
            raise BlockedInput('Synthetic fixture is not a genuine rough model boundary')
    if mm.get('operator_basis','psi2phi2_no_factorial')!='psi2phi2_no_factorial':
        raise BlockedInput('Wilson normalisation not psi2phi2_no_factorial')
    if bm.get('quartic_potential_convention',Q_CONVENTION)!=Q_CONVENTION:
        raise BlockedInput('Quartic normalisation not lambda/4!')
    ids = {key:_same_value(mm,bm,key) for key in REQUIRED_IDS}
    for key,x in ids.items():
        if not isinstance(x,str) or not x.strip():
            raise BlockedInput('Invalid Stream1 '+key)
    return ids


def _background(basis: ModelBasis, data: dict, bm:dict, source_info:dict,
                *, mf:float,ms:float) -> Background:
    required = ('mu_gev','Y','lambda4','g2','gY','g3')
    missing = [k for k in required if k not in data]
    if missing:
        # Workstream B native-only arrays MUST NOT be converted by position.
        if any(k.startswith('eft1__') for k in data):
            raise BlockedInput('Native Workstream B RGBeta coefficients require the '\
                'Stream1 canonical tensor dictionary; cannot infer Y/lambda4')
        raise BlockedInput('Canonical background missing '+repr(missing))
    mu=np.asarray(data['mu_gev'],float)
    y=np.asarray(data['Y'],complex)
    q=np.asarray(data['lambda4'],float)
    gs=np.column_stack([np.asarray(data[k],float) for k in ('g2','gY','g3')])
    n=len(mu)
    if mu.ndim!=1 or n<2 or y.shape!=(n,basis.nf,basis.nf,basis.ns) or \
       q.shape!=(n,)+(basis.ns,)*4 or gs.shape!=(n,3):
        raise BlockedInput('Canonical background dimensions do not cover 45 Weyl and all scalars')
    if not all(np.isfinite(v).all() for v in (mu,y,q,gs)) or np.min(mu)<=0:
        raise BlockedInput('Nonfinite/invalid background')
    ix=np.argsort(mu)
    mu,y,q,gs=mu[ix],y[ix],q[ix],gs[ix]
    if np.any(np.diff(mu)<=0) or not np.isclose(mu[0],ms,rtol=1e-10,atol=1e-8) \
       or not np.isclose(mu[-1],mf,rtol=1e-10,atol=1e-8):
        raise BlockedInput('Background samples must span exact MF and MS endpoints')
    if np.max(np.abs(y-y.swapaxes(1,2))) > max(1e-12,1e-9*np.max(np.abs(y))):
        raise BlockedInput('Yukawa tensor not symmetric in Weyl indices')
    for p in ((0,2,1,3,4),(0,1,3,2,4),(0,1,2,4,3)):
        if np.max(np.abs(q-q.transpose(p)))>max(1e-12,1e-9*np.max(np.abs(q))):
            raise BlockedInput('Quartic not fully symmetric')
    generators=basis.generators()
    for group,count in [('SU2',3),('U1',1),('SU3',8)]:
        fk=f'generator_{group}_fermion'; sk=f'generator_{group}_scalar'
        if (fk in data)!=(sk in data):
            raise BlockedInput('Incomplete supplied gauge generator pair: '+group)
        if fk in data:
            ff=np.stack([a for a,_ in generators[group]])
            ss=np.stack([b for _,b in generators[group]])
            if np.shape(data[fk])!=(count,basis.nf,basis.nf) or \
               np.shape(data[sk])!=(count,basis.ns,basis.ns) or \
               not np.allclose(data[fk],ff,rtol=1e-10,atol=1e-10) or \
               not np.allclose(data[sk],ss,rtol=1e-10,atol=1e-10):
                raise BlockedInput('Generator frames differ: '+group)
    return Background(basis,np.log(mu),y,q,gs,generators,
        {'source':source_info, 'native_rgbeta_mapping_asserted_by_Stream1':True,
         'interpolation':'linear_in_ln_mu', 'samples':int(n),
         'background_source_sha256':bm.get('source_sha256',{})})


def _load_c5(data:dict, key:str) -> np.ndarray:
    a=np.asarray(data[key],complex)
    if a.shape!=(3,3) or not np.isfinite(a).all() or \
        not np.allclose(a,a.T,atol=1e-22,rtol=1e-10):
        raise BlockedInput('Invalid symmetric complex 3x3 '+key)
    return (a+a.T)/2


def _boundary(basis: ModelBasis,data:dict,mm:dict,*,c5_boundary:str,
              assume_formal_identity:bool=False) -> tuple[np.ndarray,dict]:
    keys=[k for k in ('C_psi2phi2_tree_mf','C_full','C_LLSS_formal_muF') if k in data]
    if len(keys)!=1:
        if 'C_LLSS_formal' in data:
            raise BlockedInput('Legacy B formal GammaCC boundary needs Stream1 canonical map')
        raise BlockedInput('Require exactly one complete tree tensor or explicitly mapped formal LLSS')
    k=keys[0]
    if k=='C_LLSS_formal_muF':
        if 'heavy_source_from_canonical' in data:
            r=np.asarray(data['heavy_source_from_canonical'],float)
            map_source='explicit_Stream1_heavy_real_rotation'
        elif assume_formal_identity:
            r=np.eye(basis.ns-4)
            map_source='EXPLICITLY_ASSUMED_IDENTITY_QUALITATIVE_ONLY'
        else:
            raise BlockedInput('Formal LLSS requires heavy_source_from_canonical')
        fac=mm.get('gamma_cc_to_weyl_signed_factor')
        if fac is None and assume_formal_identity:fac=1.
        if fac is None:raise BlockedInput('Formal LLSS requires explicit GammaCC->Weyl factor')
        if isinstance(fac,dict) and set(fac)=={'real','imag'}:
            fac=complex(fac['real'],fac['imag'])
        try:factor=complex(fac)
        except (ValueError,TypeError):raise BlockedInput('Invalid GammaCC->Weyl factor')
        c=embed_formal_llss(data[k],basis,r,factor)
    else:
        if mm.get('C_full_is_canonical') is not True and mm.get('C_psi2phi2_tree_mf_is_canonical') is not True:
            raise BlockedInput('Full tensor lacks canonical Weyl/real-scalar declaration')
        c=np.asarray(data[k],complex).copy()
        map_source='canonical_full_tensor_supplied_by_Stream1'
    if c.shape!=(basis.nf,basis.nf,basis.ns,basis.ns) or not np.isfinite(c).all():
        raise BlockedInput('Full Wilson tensor geometry/nonfinite error')
    if symmetry_residual(c)>max(1e-24,1e-9*np.max(np.abs(c))):
        raise BlockedInput('Full Wilson violates i<->j or a<->b symmetry')
    c5_before,projection_residual=project_c5(c,basis)
    supplied=[k for k in ('C5_EFT1_mf','C5_hard') if k in data]
    if len(supplied)>1:
        raise BlockedInput('Ambiguous duplicate C5 boundary keys')
    if supplied:
        seed=_load_c5(data,supplied[0])
        if np.linalg.norm(c5_before)>1e-100:
            raise BlockedInput('C5 already in full tensor; do not add a second C5 seed')
        if c5_boundary=='zero' and np.linalg.norm(seed)>0:
            raise BlockedInput('Explicit zero C5 conflicts with nonzero supplied hard C5')
        c += embed_c5(seed,basis)
        c5_origin='supplied_'+supplied[0]
    elif c5_boundary=='zero':
        if np.linalg.norm(c5_before)>1e-100:
            raise BlockedInput('Explicit zero C5 conflicts with nonzero full-tensor C5')
        c5_origin='explicit_zero_boundary_requested'
    elif np.linalg.norm(c5_before)>1e-100:
        c5_origin='already_embedded_in_full_tensor'
    else:
        raise BlockedInput('No C5 boundary supplied; use --c5-boundary zero deliberately')
    c5_mf,c5_res=project_c5(c,basis)
    if supplied and not np.allclose(c5_mf,seed,atol=1e-24,rtol=1e-10):
        raise BlockedInput('C5 full Weinberg injection/projection mismatch')
    return c,{'tensor_source_key':k,'basis_map':map_source,
              'initial_C5_source':c5_origin,
              'MF_C5_HH_projection_residual':float(c5_res),
              'input_LLHH_residual_before_C5_injection':float(projection_residual),
              'C5_at_MF_frobenius':float(np.linalg.norm(c5_mf)),
              'LLSS_at_MF_frobenius':_llss_norm(c,basis)}


def _llss_norm(c:np.ndarray,basis:ModelBasis)->float:
    f=[i for pair in basis.leptons for i in pair]
    return float(np.linalg.norm(c[np.ix_(f,f,range(4,basis.ns),range(4,basis.ns))]))


def _other_than_llss_norm(c:np.ndarray,basis:ModelBasis)->float:
    f=[i for pair in basis.leptons for i in pair]
    outside=np.asarray(c,complex).copy()
    outside[np.ix_(f,f,range(4,basis.ns),range(4,basis.ns))]=0
    return float(np.linalg.norm(outside))


def _runtime_source_state()->dict:
    """Hash the actual LOCAL code used, not a claimed public revision."""
    root=Path(__file__).resolve().parents[3]
    rels=(
        'Numerical/parallel_t3/wilson/rough_stream2.py',
        'Numerical/parallel_t3/wilson/workflow.py',
        'Numerical/parallel_t3/wilson/adm.py',
        'Numerical/parallel_t3/wilson/coordinates.py',
        'Numerical/parallel_t3/wilson/interchange.py',
        'Numerical/parallel_t3/wilson/basis.py',
        'RGE/running/intermediate/GeneralEFT1SixTermReferenceCandidate.py',
        'RGE/running/intermediate/GeneralEFT1SevenTermFlow.py',
        'RGE/running/intermediate/GeneralEFT1CrossedWardCandidate.py')
    state={name:(sha256(root/name) if (root/name).is_file() else None)
           for name in rels}
    try:
        git=subprocess.run(['git','rev-parse','HEAD'],cwd=root,
                           capture_output=True,text=True,check=False,timeout=4)
        head=git.stdout.strip() if git.returncode==0 else None
    except (OSError,subprocess.TimeoutExpired):
        head=None
    return {'local_source_sha256':state,'local_git_head':head,
            'inspected_public_main_at_development':'25eb5e13d09dd2f138176d4e36c4fc65edfa4bb6',
            'local_code_identical_to_public_not_assumed':True}


@dataclass
class CanonicalInput:
    basis: ModelBasis
    initial: np.ndarray
    background: Background
    context:dict
    source:dict
    mf:float
    ms:float


def load_canonical(project:Path,folder:Path,model:str,*,c5_boundary:str='require',
                   assume_formal_identity:bool=False) -> CanonicalInput:
    """Read canonical tensors, not native invariant scalar parameter vectors."""
    project,folder=Path(project),Path(folder)
    basis=build_basis(model)
    if c5_boundary not in ('require','zero'):
        raise ValueError('c5_boundary must be require or zero')
    combined=_find(folder,('canonical_input.npz','rough_input.npz'),required=False)
    if combined:
        p=_find(folder,('canonical_input.json','rough_input.json','manifest.json'))
        bm=mm=_json(p)
        _hash_checked(combined,bm,project=project)
        arr=_npz(combined)
        backdata=seeddata=arr
        refs={'canonical_input_npz':{'path':str(combined),'sha256':sha256(combined)},
              'canonical_input_json':{'path':str(p),'sha256':sha256(p)}}
    else:
        b_npz=folder/'background.npz'
        s_npz=_find(folder,('boundary.npz','fermion_boundary_canonical.npz',
                             'fermion_boundary_formal.npz'))
        if not b_npz.is_file():
            raise BlockedInput('Canonical background.npz not found in '+str(folder))
        bp=_find(folder,('background.json','background_manifest.json'))
        sp=_find(folder,('boundary.json','boundary_manifest.json',
                          'fermion_boundary_manifest.json'))
        bm,mm=_json(bp),_json(sp)
        _hash_checked(b_npz,bm,project=project)
        _hash_checked(s_npz,mm,project=project)
        backdata=_npz(b_npz)
        seeddata=_npz(s_npz)
        refs={name:{'path':str(file),'sha256':sha256(file)} for name,file in
              [('background_npz',b_npz),('background_json',bp),
               ('boundary_npz',s_npz),('boundary_json',sp)]}
    mf=_float_meta(bm,mm,'mu_f_gev','MF_gev')
    ms=_float_meta(bm,mm,'mu_s_gev','MS_gev')
    context=_validate_meta(basis,bm,mm,mf=mf,ms=ms)
    back=_background(basis,backdata,bm,refs,mf=mf,ms=ms)
    c,origin=_boundary(basis,seeddata,mm,c5_boundary=c5_boundary,
                       assume_formal_identity=assume_formal_identity)
    return CanonicalInput(basis,c,back,context,
                          {'file_sha256':refs,'boundary':origin,
                           'background_metadata':{'schema':bm.get('schema'),
                             'format':bm.get('format'),'model':bm.get('model')},
                           'boundary_metadata':{'schema':mm.get('schema'),
                             'format':mm.get('format'),'model':mm.get('model')}},mf,ms)


def _matrix(a:np.ndarray)->dict:
    a=np.asarray(a)
    return {'real':a.real.tolist(),'imag':a.imag.tolist()}


def _channels(result:dict,basis:ModelBasis)->tuple[dict,dict]:
    summary={}
    arrays={}
    for label,key in [('MF','beta_channels_at_MF'),('MS','beta_channels_at_MS')]:
        summary[label]={}
        for name in SYMMETRIC_TERMS:
            ch=result[key][name]
            c5,_=project_c5(ch,basis)
            k=f'beta_C5_{name}_at_{label}_16pi2'
            arrays[k]=c5
            summary[label][name]={
                'C5_projected_16pi2':_matrix(c5),
                'C5_beta_norm_16pi2':float(np.linalg.norm(c5)),
                'LLSS_beta_norm_16pi2':_llss_norm(ch,basis),
                'full_beta_norm_16pi2':float(np.linalg.norm(ch))}
    return summary,arrays


def write_c_endpoint(case:CanonicalInput,result:dict,destination:Path,
                     *,adm_variant:str,rtol:float,atol:float,steps:int,
                     overwrite:bool=False)->dict:
    """Atomic pair export. Never label an unrun/zero fixture as genuine data."""
    destination=Path(destination)
    if destination.exists() and any(destination.iterdir()) and not overwrite:
        raise BlockedInput('Output exists; pass --overwrite to replace this model only: '+str(destination))
    basis=case.basis
    a=np.asarray(result['C_initial'],complex)
    b=np.asarray(result['C_final'],complex)
    if a.shape != b.shape or a.shape!=(basis.nf,basis.nf,basis.ns,basis.ns):
        raise BlockedInput('Solver did not return complete tensor')
    if not np.isfinite(b).all() or symmetry_residual(b)>max(1e-24,1e-8*np.max(np.abs(b))):
        raise BlockedInput('Nonfinite/non-symmetric MS tensor')
    c5_mf,rr_mf=project_c5(a,basis)
    c5_ms,rr_ms=project_c5(b,basis)
    if np.linalg.norm(c5_ms-result['C5_final'])>max(1e-24,1e-9*np.linalg.norm(c5_ms)):
        raise BlockedInput('Evolved C5 must be projected from full tensor')
    chs,beta_arrays=_channels(result,basis)
    f=[i for pair in basis.leptons for i in pair]
    projection=b[np.ix_(f,f,range(4,basis.ns),range(4,basis.ns))]
    xi=a[np.ix_(f,f,range(4,basis.ns),range(4,basis.ns))]
    other_mf=_other_than_llss_norm(a,basis)
    initial_C5_is_zero=bool(np.linalg.norm(c5_mf)<=1e-100)
    generated_C5_is_significant=bool(
        np.linalg.norm(c5_ms)>max(10*atol,1e-14*np.linalg.norm(a)))
    trajectory=result['metadata']['projections']
    if not trajectory or not np.isclose(trajectory[-1]['mu_gev'],case.ms):
        raise BlockedInput('Missing DOP853 scale diagnostics')
    arrays={
        'C_psi2phi2_ms':b,
        'C_psi2phi2_tree_mf':a,
        'C5_EFT1_ms':c5_ms,
        'C5_EFT1_mf':c5_mf,
        'C_LLSS_ms':projection,
        'C_LLSS_mf':xi,
        'beta_C_at_MF_16pi2':result['beta_at_MF'],
        'beta_C_at_MS_16pi2':result['beta_at_MS'],
        'trajectory_mu_gev':np.asarray([q['mu_gev'] for q in trajectory],float),
        'trajectory_C5_norm_gev_inverse':np.asarray(
            [q['C5_frobenius_gev_inverse'] for q in trajectory],float),
        'trajectory_full_norm_gev_inverse':np.asarray(
            [q['C_full_frobenius_gev_inverse'] for q in trajectory],float),
        **beta_arrays,
    }
    # C's strict one-loop optional array is intentionally NOT invented.
    # Matching gets both exact RGE-improved and order-zero boundaries only.
    manifest={
        'format':W_FORMAT,'tensor_npz':'wilson_at_MS.npz',
        'model':basis.model,**case.context,
        'mu_f_gev':case.mf,'mu_s_gev':case.ms,
        'units':'GeV','scheme':'MSbar',
        'canonical_wilson_convention':CONVENTION,
        'wilson_evolved_to_ms':True,'includes_eft1_weinberg':True,
        'legacy_direct_llss_to_c5_log_added':False,
        'lepton_indices':[list(x) for x in basis.leptons],
        'higgs_source_indices':list(basis.higgs),
        'heavy_source_indices':list(range(4,basis.ns)),
        'synthetic_fixture':False,
        'wilson_adm_certified':False,
        'wilson_dim5_closure_certified':False,
        'signed_GammaCC_to_Weyl_certified':False,
        'strict_1loop_endpoint_at_ms':False,
        'output_status':'NumericalCandidateNotCertified',
        'adm_variant':adm_variant,
        'adm_source':('RGE.running.intermediate.GeneralEFT1SixTermReferenceCandidate.six_term_beta'
            if adm_variant=='fonseca_eq22_six_term' else
            'RGE.running.intermediate.GeneralEFT1CrossedWardCandidate.candidate_beta'),
        'adm_reference':('Fonseca-Olgoso-Santiago JHEP07(2026)094 Eq.(22) candidate; '
                         'full-gauge completion not independently certified' if
             adm_variant=='fonseca_eq22_six_term' else
             'unitary-Ward covariant crossed-Cstar candidate; loop factor not certified'),
        'physics_equation':'16*pi^2*dC_ijab/dln(mu)=beta_C(C,C*,Y,lambda,g)',
        'beta_input_convention':'16pi2_times_derivative',
        'solver':{'method':'DOP853','rtol':rtol,'atol':atol,
                  'diagnostic_intervals':steps,'rhs_evaluations':int(result['metadata']['rhs_evaluations']),
                  'backends':result['metadata']['execution_backends']},
        'basis':basis.to_json(),
        'input':case.source,
        'local_runtime':_runtime_source_state(),
        'C5_mf':_matrix(c5_mf),'C5_ms':_matrix(c5_ms),
        'C5_EFT1_mf_norm_gev_inverse':float(np.linalg.norm(c5_mf)),
        'C5_EFT1_ms_norm_gev_inverse':float(np.linalg.norm(c5_ms)),
        'LLSS_ms_norm_gev_inverse':float(np.linalg.norm(projection)),
        'LLSS_mf_norm_gev_inverse':float(np.linalg.norm(xi)),
        'non_LLSS_tree_mf_norm_gev_inverse':other_mf,
        'C_full_mf_norm_gev_inverse':float(np.linalg.norm(a)),
        'C_full_ms_norm_gev_inverse':float(np.linalg.norm(b)),
        'C5_Weinberg_projection_residual_mf':float(rr_mf),
        'C5_Weinberg_projection_residual_ms':float(rr_ms),
        'scale_trajectory':trajectory,
        'beta_channels_16pi2':chs,
        'C5_initial_source':case.source['boundary']['initial_C5_source'],
        'C5_generated_from_zero_boundary':bool(
            initial_C5_is_zero and generated_C5_is_significant),
        'C5_generated_from_LLSS_only_boundary':bool(
            initial_C5_is_zero and generated_C5_is_significant
            and np.linalg.norm(xi)>0
            and other_mf<=1e-12*max(np.linalg.norm(xi),1e-100)),
        'no_scalar_matching_performed':True,
        'no_additional_LLSS_to_C5_leading_log':True,
        'known_limitations':[
            'candidate ADM, not independent one-loop proof',
            'native Matchete/RGBeta and GammaCC maps are upstream Stream1 assertions',
            'other one-loop fermion-hard operators not independently computed here',
            'one-loop scalar threshold and SMEFT running owned by Stream3',
        ],
    }
    # Stage files before moving; another run must not see a half-written pair.
    destination.parent.mkdir(parents=True,exist_ok=True)
    tmp=Path(tempfile.mkdtemp(prefix='.rough_stream2_',dir=destination.parent))
    try:
        f=tmp/'wilson_at_MS.npz'
        np.savez_compressed(f,**arrays)
        manifest['tensor_npz_sha256']=sha256(f)
        j=tmp/'wilson.json'
        j.write_text(json.dumps(manifest,indent=2,allow_nan=False)+'\n',encoding='utf-8')
        if destination.exists():
            if not overwrite:raise BlockedInput('Output already exists')
            if not destination.is_dir():raise BlockedInput('Output path not a directory')
            shutil.rmtree(destination)
        tmp.rename(destination)
    finally:
        if tmp.exists():shutil.rmtree(tmp)
    return manifest


def run_one(project:Path,model:str,*,input_root:Path=INPUT_ROOT,
            output_root:Path=OUTPUT_ROOT,c5_boundary:str='require',
            adm_variant:str='fonseca_eq22_six_term',rtol:float=1e-6,
            atol:float=1e-21,steps:int=6,overwrite:bool=False,
            assume_formal_identity:bool=False)->dict:
    if adm_variant not in ALLOWED_ADMS:
        raise BlockedInput('Only selected covariant candidate ADMs are enabled')
    folder=Path(project)/input_root/slug(model)
    case=load_canonical(Path(project),folder,model,
                        c5_boundary=c5_boundary,
                        assume_formal_identity=assume_formal_identity)
    result=run_wilson(case.basis,case.initial,case.background,case.mf,case.ms,
                      adm_variant=adm_variant,solver='DOP853',steps=steps,
                      rtol=rtol,atol=atol)
    out=Path(project)/output_root/slug(model)
    return write_c_endpoint(case,result,out,adm_variant=adm_variant,rtol=rtol,
                            atol=atol,steps=steps,overwrite=overwrite)


def probe_saved_b(project:Path,b_root:Path,model:str)->dict:
    """Read actual saved Workstream B format, without pretending it is canonical.

    Diagnostic-only; native RGBeta invariant parameters and Matchete GammaCC
    cannot be fed to the 45-Weyl ADM without a legitimate Stream1 map.
    """
    folder=Path(project)/b_root/slug(model)
    bm_path=folder/'background_manifest.json'
    fm_path=folder/'fermion_boundary_manifest.json'
    if not bm_path.is_file() or not fm_path.is_file():
        return {'status':'B_saved_input_missing','folder':str(folder)}
    bm,fm=_json(bm_path),_json(fm_path)
    if bm.get('model')!=model or fm.get('model')!=model:
        raise BlockedInput('Saved B model mismatch')
    if bm.get('UV_endpoint_sha256')!=fm.get('same_run_uv_state_sha256'):
        raise BlockedInput('Saved B fermion and RG background from different UV states')
    obj={'status':'saved_B_native_requires_Stream1_canonical_mapping',
         'folder':str(folder),
         'same_run_uv_state_sha256':bm['UV_endpoint_sha256'],
         'fermion_boundary_available':bool(fm.get('formal_LLsS_tree_available')),
         'hard_C5_status':fm.get('hard_C5_status'),
         'native_scalar_frame':bm.get('scalar_basis'),
         'native_Wilson_frame':fm.get('fermion_matching_source_basis'),
         'canonical_Y_lambda4_available':False}
    p=folder/'fermion_boundary_formal.npz'
    if p.is_file():
        if 'fermion_boundary_npz_sha256' in fm and \
            sha256(p)!=fm['fermion_boundary_npz_sha256']:
            raise BlockedInput('Saved B fermion NPZ hash mismatch')
        with np.load(p,allow_pickle=False) as x:
            if 'C_LLSS_formal' in x:
                c=np.asarray(x['C_LLSS_formal'])
                obj['formal_LLSS_norm_gev_inverse']=float(np.linalg.norm(c))
                obj['formal_LLSS_shape']=list(c.shape)
                obj['nonzero_saved_LLSS']=bool(np.linalg.norm(c)>0)
    return obj
