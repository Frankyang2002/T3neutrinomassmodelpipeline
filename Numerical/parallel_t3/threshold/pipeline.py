"""Workstream C A/B -> simultaneous MS scalar threshold -> SMEFT low scale.

A and B must provide explicit compatible input manifests. No attempt is made
here to recompute/guess their F-threshold or EFT1 Wilson trajectories.
"""
from __future__ import annotations
from dataclasses import dataclass
from hashlib import sha256
from pathlib import Path
from typing import Any
import json
import math
import numpy as np

from .matching import (check_rotation, match_scalar_threshold,
                       mass_eigenbasis, symmetric3, scalar_gauge_ward_residuals)
from .models import normalize_model,n_real,scalar_labels,special_kind

W_FORMAT='t3-parallel-c-wilson-input-v1'
B_FORMAT='t3-parallel-c-background-input-v1'
RESULT_FORMAT='t3-parallel-c-scalar-threshold-result-v1'
W_CONVENTION='C_ijab_full_sum_no_factorial_C5_1over2_H_R_interleaved_v1'
Q_CONVENTION='V4_lambda_real_symmetric_over_4factorial_v1'


def read_json(p: Path) -> dict:
    data=json.loads(Path(p).read_text(encoding='utf-8-sig'))
    if not isinstance(data,dict):raise ValueError('Manifest root must be a JSON object')
    return data


def digest(path: Path) -> str:
    return sha256(Path(path).read_bytes()).hexdigest()


def array_dump(a: np.ndarray) -> dict:
    a=np.asarray(a)
    return {'real':a.real.tolist(),'imag':a.imag.tolist()}


def _npz(path: Path, key: str) -> np.ndarray:
    raise RuntimeError('Use the NPZ context loader instead of isolated open')


def _load_pair(manifest: Path, expected: str) -> tuple[dict,dict,str,str]:
    m=read_json(manifest)
    if m.get('format')!=expected:
        raise ValueError(f'{manifest}: expected format {expected}')
    filename=m.get('tensor_npz')
    if (not isinstance(filename,str) or not filename.endswith('.npz') or
            filename in ('.','..') or Path(filename).name!=filename):
        raise ValueError('tensor_npz must be a sibling .npz basename')
    p=manifest.parent/filename
    if not p.is_file():raise FileNotFoundError(p)
    sha=digest(p)
    if 'tensor_npz_sha256' in m and sha!=m['tensor_npz_sha256']:
        raise ValueError(f'NPZ checksum mismatch: {p}')
    with np.load(p,allow_pickle=False) as data:
        content={name:np.asarray(data[name]).copy() for name in data.files}
    return m,content,digest(manifest),sha


def _scalar_float(x: Any, label: str) -> float:
    x=float(x)
    if not math.isfinite(x):raise ValueError(f'Nonfinite {label}')
    return x


def _contract(w: dict,b: dict) -> tuple[str,float,float]:
    """Fail closed on A/B stale or mixed benchmark states and scales."""
    wm=normalize_model(w.get('model',''))
    bm=normalize_model(b.get('model',''))
    if wm!=bm:raise ValueError('A/B model mismatch')
    required_same=('matching_context_id','matched_seed_id','lagrangian_source_id')
    for key in required_same:
        if not isinstance(w.get(key),str) or not w[key] or w[key]!=b.get(key):
            raise ValueError(f'A/B {key} absent or mismatch: cannot combine different F seeds')
    mf=_scalar_float(w.get('mu_f_gev'),'mu_f_gev')
    ms=_scalar_float(w.get('mu_s_gev'),'mu_s_gev')
    if not (0<ms<mf):raise ValueError('Required EFT sequence: MF > MS > 0')
    for key,val in (('mu_f_gev',mf),('mu_s_gev',ms)):
        if abs(_scalar_float(b.get(key),key)-val)>1e-10*val:
            raise ValueError(f'A/B {key} mismatch')
    if w.get('canonical_wilson_convention')!=W_CONVENTION:
        raise ValueError('A canonical spinor/master Wilson convention not explicitly supplied')
    if b.get('quartic_potential_convention')!=Q_CONVENTION:
        raise ValueError('B full quartic normalization unrecognized')
    if w.get('scheme')!='MSbar' or b.get('scheme')!='MSbar':
        raise ValueError('C implements MSbar only, no implicit scheme conversion')
    if w.get('wilson_evolved_to_ms') is not True or w.get('includes_eft1_weinberg') is not True:
        raise ValueError('A must supply fully evolved EFT1 Weinberg and LLSS at MS')
    if w.get('legacy_direct_llss_to_c5_log_added') is not False:
        raise ValueError('A must explicitly confirm no duplicate fixed-order LLSS->C5 log')
    if b.get('single_common_scalar_threshold') is not True:
        raise ValueError('B must request BOTH scalars integrated simultaneously at MS')
    if w.get('scalar_basis_id')!=b.get('wilson_scalar_basis_id'):
        raise ValueError('B rotation must refer exactly to A Wilson source scalar basis')
    if b.get('scalar_mass_origin') not in ('unbroken_electroweak','EWSB_resummed',
                                          'synthetic_algebra_only'):
        raise ValueError('Unknown scalar mass origin')
    if w.get('synthetic_fixture',False)!=b.get('synthetic_fixture',False):
        raise ValueError('Cannot mix synthetic and real inputs')
    if w.get('synthetic_fixture',False) is not True and any(
            'synthetic' in str(w.get(key,'')).lower()
            for key in ('matching_context_id','matched_seed_id','lagrangian_source_id')):
        raise ValueError('Real A/B input uses a synthetic provenance identifier')
    if b['scalar_mass_origin']=='synthetic_algebra_only' and not w.get('synthetic_fixture',False):
        raise ValueError('Synthetic scalar masses cannot be labelled real')
    if b.get('quartic_inventory_status')!='full_tensor':
        raise ValueError('B must provide the full scalar quartic (no reduced portal)')
    if w.get('units')!='GeV' or b.get('units')!='GeV':
        raise ValueError('A/B units must explicitly be GeV, GeV^-1 for Wilson and GeV^2 for mass')
    return wm,mf,ms


def _load_sm(background: dict) -> tuple[dict,dict]:
    needed=('gY','g2','g3','lambdaH','ye','yu','yd')
    sm=background.get('sm_couplings_at_ms')
    if not isinstance(sm,dict) or any(k not in sm for k in needed[:4]):
        raise ValueError('B must supply four SM couplings at MS')
    scalars={k:_scalar_float(sm[k],k) for k in needed[:4]}
    y={}
    for k in needed[4:]:
        value=background.get(k)
        if value is None:raise ValueError(f'Missing full SM Yukawa matrix {k} in B NPZ')
        a=np.asarray(value,complex)
        if a.shape!=(3,3) or not np.isfinite(a).all():
            raise ValueError(f'{k} must be full complex finite 3x3')
        y[k]=a
    return scalars,y


def evolve_existing_smeft(c5: np.ndarray, b_npz: dict, bmeta: dict,
                          *, mu_s: float, mu_low: float, rtol: float=1e-8,
                          atol: float=1e-11) -> dict:
    """CALL the repository's established one-loop full-flavour SMEFT runner."""
    from Numerical.running.SMWeinbergEvolution import (SMWeinbergInitialConditions,
                                                       evolve_sm_weinberg)
    sm,y=_load_sm({**b_npz, **bmeta})
    if not (0<mu_low<=mu_s):
        raise ValueError('SMEFT low scale must satisfy 0 < mu_low <= MS')
    start=SMWeinbergInitialConditions(gY=sm['gY'],g2=sm['g2'],g3=sm['g3'],
              lambdaH=sm['lambdaH'],ye=y['ye'],yu=y['yu'],yd=y['yd'],K=c5)
    result=evolve_sm_weinberg(start,mu_s,mu_low,rtol=rtol,atol=atol)
    return {'c5_low':symmetric3(result.K,'low-scale C5'),
            'nfev':result.nfev,'solver_success':result.solver_success,
            'running_sm_parameters':{'gY':result.gY,'g2':result.g2,
                 'g3':result.g3,'lambdaH':result.lambdaH}}


def run_case(*, wilson_manifest: Path, background_manifest: Path,
             output_dir: Path, mu_low: float=173., mode: str='rge_improved',
             vev: float=246.22, allow_projection_residual: bool=False,
             match_only: bool=False, overwrite: bool=False) -> dict:
    """Run one A/B pair. All outputs are contained in a new C-owned directory."""
    wp,bp,output_dir=map(Path,(wilson_manifest,background_manifest,output_dir))
    w,wa,w_json_hash,w_npz_hash=_load_pair(wp,W_FORMAT)
    b,ba,b_json_hash,b_npz_hash=_load_pair(bp,B_FORMAT)
    model,mf,ms=_contract(w,b)
    n=n_real(model)
    for name,arr,shape in (('mass_squared_gev2_ms',ba.get('mass_squared_gev2_ms'),(n,n)),
                           ('R_wilson_from_background',ba.get('R_wilson_from_background'),(n,n))):
        if arr is None or np.shape(arr)!=shape:raise ValueError(f'{name}: expected shape {shape}')
    check_rotation(ba['R_wilson_from_background'],n)
    if list(b.get('scalar_background_real_labels',())) and len(b['scalar_background_real_labels'])!=n:
        raise ValueError('B scalar background component labels incomplete')
    if list(w.get('scalar_source_real_labels',())) and len(w['scalar_source_real_labels'])!=n:
        raise ValueError('A source scalar labels incomplete')
    required=('C_psi2phi2_ms',)
    for key in required:
        if key not in wa:raise ValueError(f'A full Wilson endpoint missing {key}')
    for key in ('quartic_real_full_ms','ye','yu','yd'):
        if key not in ba:raise ValueError(f'B full background endpoint missing {key}')
    wc=wa['C_psi2phi2_ms']
    if wc.shape[-2:]!=(4+n,4+n):
        raise ValueError('A scalar count differs from B scalar background')
    ward={'status':'not_supplied_not_certified'}
    if 'scalar_gauge_generators_full' in ba:
        ward=scalar_gauge_ward_residuals(ba['quartic_real_full_ms'],
                  ba['mass_squared_gev2_ms'],ba['scalar_gauge_generators_full'])
        ward['status']='checked_source_supplied_generators'
        if not ward['passed']:
            raise ValueError('Scalar quadratic/quartic Ward check failed on B-supplied generators')
    elif b.get('gauge_ward_certified') is True:
        raise ValueError('B claims gauge Ward certification without the generators')
    tree=wa.get('C_psi2phi2_tree_mf')
    strict_ms=wa.get('C_psi2phi2_ms_strict_1loop')
    if mode=='strict_one_loop':
        if (strict_ms is None or w.get('strict_1loop_endpoint_at_ms') is not True
                or not isinstance(w.get('strict_1loop_origin'),str)
                or not w['strict_1loop_origin'].strip()):
            raise ValueError('Strict mode needs A one-loop-truncated FULL endpoint at MS with explicit strict_1loop_endpoint_at_ms and strict_1loop_origin')
    separately=wa.get('C5_EFT1_ms')
    hard=wa.get('C5_independent_finite_at_MS')
    if hard is not None and w.get('finite_hard_ms_kind')!='independent_finite_only':
        raise ValueError('Hard MS terms must be independently computed finite terms, not old logs')
    if w.get('finite_hard_ms_kind')=='legacy_direct_leading_log':
        raise ValueError('Duplicate fixed-order log is forbidden')
    channel_names=tuple(b.get('quartic_channel_names',()))
    threshold=match_scalar_threshold(full_wilson_ms=wc,
             full_wilson_tree_mf=tree,
             full_wilson_strict_ms=strict_ms,
             lepton_indices=np.asarray(w['lepton_indices'],dtype=int),
             higgs_indices=np.asarray(w['higgs_source_indices'],dtype=int),
             heavy_indices=np.asarray(w['heavy_source_indices'],dtype=int),
             scalar_mass2_background=ba['mass_squared_gev2_ms'],
             full_quartic_ms=ba['quartic_real_full_ms'],
             R_wilson_from_background=ba['R_wilson_from_background'],
             mu_s=ms,mode=mode,separately_supplied_c5_ms=separately,
             independent_hard_c5=hard,
             full_quartic_channels=ba.get('quartic_real_channels_ms'),
             channel_names=channel_names,
             enforce_pure_weinberg=not allow_projection_residual)
    ml=_scalar_float(mu_low,'mu_low')
    v=_scalar_float(vev,'vev')
    if not (0<ml<=ms and 0<v<1000):raise ValueError('Invalid low scale or vev')
    final=threshold['c5_smeft_at_ms']
    if match_only:
        low={'c5_low':final.copy(),'nfev':0,'solver_success':None,
             'running_sm_parameters':{},'smeft_running_performed':False}
    else:
        low=evolve_existing_smeft(final,ba,b,mu_s=ms,mu_low=ml)
        low['smeft_running_performed']=True
    c5low=low['c5_low']
    mass=-.5*v*v*c5low
    singular=np.sort(np.linalg.svd(mass,compute_uv=False))*1e9
    m2,U=mass_eigenbasis(ba['mass_squared_gev2_ms'])
    quality={'synthetic_fixture':bool(w.get('synthetic_fixture',False)),
             'wilson_adm_certified':bool(w.get('wilson_adm_certified',False)),
             'wilson_dim5_closure_certified':bool(w.get('wilson_dim5_closure_certified',False)),
             'scalar_operator_inventory_certified':bool(b.get('scalar_operator_inventory_certified',False)),
             'native_matchete_rgbeta_normalization_certified':bool(b.get('native_matchete_rgbeta_normalization_certified',False)),
             'signed_GammaCC_to_Weyl_certified':bool(w.get('signed_GammaCC_to_Weyl_certified',False)),
             'complete_finite_scalar_matching_certified':False,
             'physical_neutrino_mass_certified':False,
             'low_scale_running_performed':low['smeft_running_performed'],
             'weinberg_projection_pure':max(threshold['rge_input_C5_residual'],
                          threshold['scalar_bubble_projection_residual'])<=1e-8}
    if output_dir.exists() and any(output_dir.iterdir()) and not overwrite:
        raise FileExistsError(f'Output directory not empty: {output_dir}; use --overwrite')
    output_dir.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(output_dir/'scalar_mass_basis.npz',
                mass_squared_background_gev2=ba['mass_squared_gev2_ms'],
                R_wilson_from_background=ba['R_wilson_from_background'],
                mass_squared_eigenvalues_gev2=m2,U_background_from_mass_eigen=U)
    np.savez_compressed(output_dir/'C5_smeft_at_MS.npz',
                C5_EFT1_MS=threshold['c5_eft1'],
                C5_EFT1_MS_rge_improved_reference=threshold['c5_eft1_rge_improved'],
                C5_scalar_bubble_MS=threshold['c5_scalar_bubble'],
                C5_independent_finite_MS=threshold['c5_independent_hard'],
                C5_SMEFT_MS=final)
    np.savez_compressed(output_dir/'C5_smeft_low.npz',C5_SMEFT_low=c5low,
                         m_nu_GeV=mass,mass_singular_values_eV=singular)
    if threshold['channel_c5']:
        np.savez_compressed(output_dir/'C5_per_quartic_channel.npz',
                            **{f'channel_{i:03d}':a for i,a in enumerate(threshold['channel_c5'].values())})
    pair_meta=[]
    for pair in threshold['bubble']['pair_data']:
        pair_meta.append({**{key:val for key,val in pair.items() if key!='projected_c5'},
                          'C5_pair':array_dump(pair['projected_c5'])})
    manifest={'format':RESULT_FORMAT,'model':model,'status':'NumericalCandidateNotCertified',
              'matching_context_id':w['matching_context_id'],'matched_seed_id':w['matched_seed_id'],
              'lagrangian_source_id':w['lagrangian_source_id'],
              'source_commit_readonly':'a944e0030b8ef24b23342d99e0d8f33cef3cffd5',
              'scalar_gauge_ward_diagnostic':ward,
              'threshold_order':['UV','all_fermions_out_at_MF','EFT1_full_dim5_running',
                                 'all_scalar_eigenstates_out_together_at_MS','SMEFT_one_loop_weinberg'],
              'mu_f_gev':mf,'mu_s_gev':ms,'mu_low_gev':ml,'vev_gev':v,
              'scheme':'MSbar','counting_mode':mode,'quartic_convention':Q_CONVENTION,
              'wilson_convention':W_CONVENTION,
              'scalar_mass_origin':b['scalar_mass_origin'],
              'scalar_mass_kind':special_kind(model),
              'scalar_basis_ids':{'A':w['scalar_basis_id'],
                                  'B':b.get('background_scalar_basis_id')},
              'scalar_source_real_labels':w.get('scalar_source_real_labels'),
              'scalar_background_real_labels':b.get('scalar_background_real_labels'),
              'n_real_heavy_scalars':n,'eigenmasses_gev':np.sqrt(m2).tolist(),
              'single_scalar_threshold':True,'scalar_eigenstates_integrated_out':n,
              'all_unordered_eigenstate_pairs_included':threshold['bubble']['pair_count'],
              'c5_contributions':{'already_evolved_EFT1_C5':array_dump(threshold['c5_eft1']),
                   'one_loop_scalar_bubble':array_dump(threshold['c5_scalar_bubble']),
                   'other_independent_finite_hard':array_dump(threshold['c5_independent_hard'])},
              'C5_SMEFT_at_MS':array_dump(final),
              'C5_SMEFT_at_low':array_dump(c5low),
              'm_nu_GeV':array_dump(mass),'mass_singular_values_eV':singular.tolist(),
              'matching_residuals':{'already_evolved_C5':threshold['rge_input_C5_residual'],
                       'strict_C5_one_loop_truncated':threshold['strict_c5_projection_residual'],
                       'finite_scalar_bubble':threshold['scalar_bubble_projection_residual']},
              'normalization_and_status':quality,
              'finite_contributions':{'heavy_scalar_pair_bubble':'implemented_conditional',
                         'independent_hard_MS':'supplied' if hard is not None else 'absent',
                         'other_one_loop_scalar_diagrams':'uncomputed',
                         'scalar_wavefunction_EOM_schemewise':'uncomputed',
                         'other_dim5_classes_at_MS':'not_included_in_bubble'},
              'rge_bookkeeping':{'C5_EFT1_already_evolved':True,
                        'C5_EFT1_order_source':threshold['c5_eft1_power_counting_source'],
                        'strict_C5_one_loop_truncated_endpoint_used':mode=='strict_one_loop',
                        'unexpanded_B_running_quartic_mass_inputs_at_MS':True,
                        'old_LLSS_to_C5_direct_leading_log_added_again':False,
                        'bubble_LLSS_source':threshold['bubble_Wilson_order_source'],
                        'strict_one_loop_requires_MF_orderzero':True},
              'SMEFT_running':{'performed':low['smeft_running_performed'],
                              'backend':'Numerical.running.SMWeinbergEvolution',
                              'nfev':low['nfev'],'success':low['solver_success'],
                              'sm_at_low':low['running_sm_parameters']},
              'input_sha256':{'A_manifest':w_json_hash,'A_npz':w_npz_hash,
                              'B_manifest':b_json_hash,'B_npz':b_npz_hash},
              'output_npz_sha256':{name:digest(output_dir/name) for name in
                       ('scalar_mass_basis.npz','C5_smeft_at_MS.npz','C5_smeft_low.npz')},
              'quartic_channel_names':list(threshold['channel_c5']),
              'pair_contributions':pair_meta,
              'unresolved':[
                   'Canonical signed Matchete GammaCC -> full Weyl Wilson factor (A)',
                   'Matchete/RGBeta complete quartic/mass tensor numerical dictionary (A/B)',
                   'Complete EFT1 dim5 ADM and other operator classes (A)',
                   'Independent complete MS one-loop diagrams/wavefunction/EOM terms',
                   'Realistic EWSB-resummed split spectrum projects to exact SU2 Weinberg only up to higher-dimension terms',
                   'Only computed C5 at fixed v, not pole neutrino masses']}
    (output_dir/'manifest.json').write_text(json.dumps(manifest,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    return manifest
