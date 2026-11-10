"""Batch 16: source-CG -> reference-quartic span analysis for A(-4), D(-2).

This is a CONDITIONAL coordinate audit and re-evolution, not a demonstration
that RGBeta generator conventions, complex phase orientations, the full EFT1
ADM or finite matching have been independently established.

The archived Matchete quartic inventory is incomplete for the RG-closed
potential. Only reference directions lying in the archived native span are
reconstructed; all missing directions remain explicitly reference-generated.
The F threshold is not recomputed and the original Batch 13 files are read-only.
"""
from __future__ import annotations

from collections import OrderedDict
from itertools import permutations
from pathlib import Path
from typing import Any, Mapping
import json
import math
from hashlib import sha256

import numpy as np

MODELS = ('A_-4', 'D_-2')
FORMAT = 't3-native-rgbeta-scalar-coordinate-span-batch16-v1'
ROOT_B13 = Path('output/parallel_t3/thesis_master_wilson_batch13_cgfix_v3')
LOOP = 1. / (16 * math.pi**2)
GROUP_REFERENCES = OrderedDict([
    ('H4', ('lambdaH',)),
    ('S1_4', ('lambdaS1', 'lambdaS1Adj')),
    ('S2_4', ('lambdaS2', 'lambdaS2Adj')),
    ('H2S1_2', ('lambdaH1', 'lambdaH1Adj')),
    ('H2S2_2', ('lambdaH2', 'lambdaH2Adj')),
    ('S1_2S2_2', ('lambda12', 'lambda12Adj')),
])


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError('Batch16: ' + message)


def relative(a: Any, b: Any) -> float:
    x = np.asarray(a, complex)
    y = np.asarray(b, complex)
    return float(np.linalg.norm(x-y)/max(np.linalg.norm(x), np.linalg.norm(y), 1e-100))


def _real(tensor: Any, name: str) -> np.ndarray:
    z = np.asarray(tensor, complex)
    require(z.size > 0 and np.isfinite(z).all(), f'{name}: nonfinite or empty tensor')
    require(np.max(np.abs(z.imag)) < 1e-10*max(np.max(np.abs(z.real)), 1e-100),
            f'{name}: expected a real, Hermitian potential tensor')
    return z.real.astype(float)


def _rank(a: np.ndarray, *, rtol: float=1e-10) -> int:
    s = np.linalg.svd(a, compute_uv=False)
    return int(np.count_nonzero(s > (max(s) * rtol if s.size else 0)))


def fit_sector(reference: Mapping[str, Any], native: Mapping[str, Any],
               *, tolerance: float=1e-8) -> dict:
    """Linear coordinate map R A=N, plus inverse representability check.

    R columns are *potential* fourth-derivative tensors from the RGBeta
    reference definition; N columns are from native Matchete CGs. A converts
    native monomial coefficients to reference monomial coefficients. An
    inverse exists uniquely only for complete, equal-rank square bases.
    """
    require(bool(reference), 'sector has no reference directions')
    ref_names = list(reference)
    native_names = list(native)
    arrays = [_real(reference[k], k) for k in ref_names]
    shape = arrays[0].shape
    require(all(x.shape == shape for x in arrays), 'incompatible reference shapes')
    natives = [_real(native[k], k) for k in native_names]
    require(all(x.shape == shape for x in natives), 'incompatible native shapes')
    R = np.column_stack([x.ravel() for x in arrays])
    N = (np.column_stack([x.ravel() for x in natives]) if natives
         else np.zeros((R.shape[0], 0), float))
    rr = _rank(R)
    nr = _rank(N)
    # A is defined even for non-square or non-unique maps, but that does not
    # imply the reference vectors can all be expressed using the archive.
    A = np.linalg.lstsq(R, N, rcond=1e-12)[0]
    native_fit = (float(np.linalg.norm(R@A-N) / max(np.linalg.norm(N),1e-100))
                  if native_names else None)
    if native_names:
        B = np.linalg.lstsq(N, R, rcond=1e-12)[0]
        reference_fit = float(np.linalg.norm(N@B-R)/max(np.linalg.norm(R),1e-100))
    else:
        B = np.zeros((0,len(ref_names)),float)
        reference_fit = 1.
    native_individual = {
        name:float(np.linalg.norm(R@A[:,i]-N[:,i]) /
                   max(np.linalg.norm(N[:,i]),1e-100))
        for i,name in enumerate(native_names)}
    reference_individual = {
        name:float(np.linalg.norm(N@B[:,i]-R[:,i]) /
                   max(np.linalg.norm(R[:,i]),1e-100))
        for i,name in enumerate(ref_names)}
    native_in_ref = bool(native_names and native_fit <= tolerance)
    ref_in_native = bool(native_names and reference_fit <= tolerance)
    unique = (native_in_ref and ref_in_native and
              len(ref_names)==len(native_names)==rr==nr)
    return {
        'reference_names':ref_names,'native_names':native_names,
        'reference_rank':rr,'native_rank':nr,
        'native_in_reference_span':native_in_ref,
        'reference_in_native_span':ref_in_native,
        'unique_bijective_coordinate_map':unique,
        'native_to_reference_coefficients':A.tolist(),
        'reference_to_native_coefficients':B.tolist() if ref_in_native else None,
        'native_direction_fit_relative':native_fit,
        'native_direction_individual_fit_relative':native_individual,
        'reference_direction_fit_relative':reference_fit,
        'reference_direction_individual_fit_relative':reference_individual,
        'archived_native_columns_missing':max(rr-nr,0),
        'reference_coordinate_provenance':'Aminus4EFT1Couplings RGBeta candidate V4 basis',
        'absolute_matchete_to_rgbeta_generator_normalization_certified':False,
    }


def _group(name: str) -> str | None:
    if name == 'λ': return 'H4'
    if name.startswith('lambdaS1Inv'): return 'S1_4'
    if name.startswith('lambdaS2Inv'): return 'S2_4'
    if name.startswith('lambdaH1Inv'): return 'H2S1_2'
    if name.startswith('lambdaH2Inv'): return 'H2S2_2'
    if name.startswith('lambda12Inv'): return 'S1_2S2_2'
    return None


def _native_phase_tensor(seed: dict, terms: list[dict], model: str,
                         name: str, value: complex) -> np.ndarray:
    """Expand an EXACT Matchete unit topology coupling, including its Bar term.

    Called for real and imaginary orientations; no guessed h.c. sign.
    """
    import sympy as sp
    from RGE.running.intermediate.ScalarOnlyTensorAdapters import build_eft1_quartic_tensor
    d1,d2={'A_-4':(1,3),'D_-2':(3,1)}[model]
    raw=build_eft1_quartic_tensor(
        {'CGRegistry':seed['CGRegistry'],'ScalarQuarticTerms':terms},
        d_s1=d1,d_s2=d2,shared_scalar=False)
    out=np.zeros((12,)*4, float)
    for key,expression in raw.nonzero_items():
        symbols=[s for s in expression.free_symbols if str(s)==name]
        require(len(symbols)<=1, 'ambiguous topology coupling symbol')
        z=expression.subs({symbols[0]:sp.sympify(value)} if symbols else {})
        z=complex(sp.N(z,16))
        require(math.isfinite(z.real) and math.isfinite(z.imag) and
                abs(z.imag)<1e-9*max(abs(z.real),1.),
                'topology potential has non-real component after conjugation')
        for index in set(permutations(tuple(int(i)-1 for i in key))):
            out[index]=z.real
    return out


def transform_tensor(tensor: np.ndarray, U: np.ndarray) -> np.ndarray:
    t=np.asarray(tensor,float)
    u=np.asarray(U,float)
    require(t.shape==(12,)*4 and u.shape==(12,12), 'quartic transform dimension')
    require(np.allclose(u.T@u,np.eye(12),atol=1e-10),'non-orthogonal quartic transform')
    return np.einsum('abcd,aA,bB,cC,dD->ABCD',t,u,u,u,u,optimize=True)


def _phase_audit(seed: dict, grouped: Mapping[str,list[dict]], model:str,
                 native_R: np.ndarray)->dict:
    """Internal real-linearity of the source topology under complex phases.

    Not an *independent* cross-code phase-normalization certificate.
    """
    if 'lambdaT3' not in grouped:
        return {'status':'Blocked_MissingNativeTopology'}
    full_R=np.eye(12)
    full_R[4:,4:]=native_R
    terms=grouped['lambdaT3']
    t_re=transform_tensor(_native_phase_tensor(seed,terms,model,'lambdaT3',1.),full_R)
    t_im=transform_tensor(_native_phase_tensor(seed,terms,model,'lambdaT3',1j),full_R)
    a,b=0.37,-0.61
    mixed=transform_tensor(_native_phase_tensor(seed,terms,model,'lambdaT3',a+1j*b),full_R)
    linearity=relative(mixed,a*t_re+b*t_im)
    norm_re=float(np.linalg.norm(t_re)); norm_im=float(np.linalg.norm(t_im))
    require(norm_re>0 and norm_im>0 and linearity<1e-8,
            'native topology real/imaginary source check fails')
    return {
        'status':'ExactNativeTopologyRealImaginaryLinearityVerified_Only',
        'source_real_norm':norm_re,'source_imaginary_norm':norm_im,
        'linear_combination_relative_residual':linearity,
        'real_imaginary_directions_independent':_rank(np.column_stack((
            t_re.ravel(),t_im.ravel())))==2,
        'independent_RGBeta_topology_phase_normalization_certified':False,
    }


def _input_coupling(couplings: Mapping[str,Any], key: str) -> float:
    require(key in couplings, 'missing native RGBeta value '+key)
    value=couplings[key]
    if isinstance(value,Mapping) and {'real','imag'} <= set(value):
        val=np.asarray(value['real'],float)+1j*np.asarray(value['imag'],float)
    else:
        val=np.asarray(value,complex)
    require(val.shape==() and abs(val.imag)<1e-11*max(1.,abs(val.real)),
            'non-real baseline quartic '+key)
    return float(val.real)


def build_map(seed:dict,model:str,*,tolerance:float=1e-8)->tuple[dict,dict,dict]:
    """Build source-exact archived columns and all candidate reference groups.

    Returns (audit, sector objects [in RAM only], physical basis metadata).
    """
    from Numerical.parallel_t3.complete.two_model_physical_closure_batch14 import (
        _baseline_reference_columns, numeric_quartic_source_channel)
    from Numerical.parallel_t3.complete.native_triplet_basis_batch15 import (
        heavy_source_from_reference)
    import re
    require(model in MODELS,'unsupported model')
    require(seed.get('ScalarQuarticTermCount')==len(seed.get('ScalarQuarticTerms',[])),
            'native quartic seed count inconsistent')
    native_R,basis_meta=heavy_source_from_reference(model,seed)
    full_R=np.eye(12);full_R[4:,4:]=native_R
    grouped:dict[str,list[dict]]={}
    for row in seed['ScalarQuarticTerms']:
        names=re.findall(r'Coupling\[([^,\]]+),',row['TermInputForm'])
        require(len(names)==1,'native scalar monomial has ambiguous coupling')
        grouped.setdefault(names[0],[]).append(row)
    refs=_baseline_reference_columns(model)
    sectors={}
    diagnostics={}
    handled=set()
    for group,candidates in GROUP_REFERENCES.items():
        r={k:np.asarray(refs[k],float) for k in candidates if k in refs}
        if not r:continue
        n={}
        for name,rows in grouped.items():
            if _group(name)!=group:continue
            raw=numeric_quartic_source_channel(seed,rows,model,name)
            n[name]=transform_tensor(raw,full_R)
            handled.add(name)
        info=fit_sector(r,n,tolerance=tolerance)
        sectors[group]={'reference':r,'native':n,'fit':info}
        diagnostics[group]=info
    unhandled=sorted(set(grouped)-handled-{'lambdaT3'})
    diagnostics['source_names_without_reference_sector']=unhandled
    require(not unhandled, 'unhandled archived native quartic channels: '+str(unhandled))
    diagnostics['topology_real_imaginary']=_phase_audit(seed,grouped,model,native_R)
    diagnostics['source_CG_triplet_basis']=basis_meta
    diagnostics['matchete_topology_source_present']='lambdaT3' in grouped
    diagnostics['reference_directions_with_no_archived_span']=[
        name for info in diagnostics.values() if isinstance(info,dict)
        and 'reference_direction_individual_fit_relative' in info
        for name,residual in info['reference_direction_individual_fit_relative'].items()
        if residual>tolerance]
    diagnostics['all_archived_native_channels_fit_reference_span']=all(
        x['native_in_reference_span'] for x in diagnostics.values()
        if isinstance(x,dict) and 'native_in_reference_span' in x
        and x['native_names'])
    diagnostics['all_reference_directions_in_archived_native_span']=all(
        x['reference_in_native_span'] for x in diagnostics.values()
        if isinstance(x,dict) and 'reference_in_native_span' in x)
    diagnostics['complete_independent_native_RGBeta_basis_certified']=False
    return diagnostics,sectors,{'native_R':native_R}


def rebuilt_quartic_delta(sectors: Mapping[str,dict],couplings:Mapping[str,Any]
                          )->tuple[np.ndarray,dict]:
    """Use native source columns only where they SPAN the reference sector.

    If an archive omits a reference direction, retain that reference column;
    never fabricate a Matchete source operator or silently drop a coupling.
    """
    correction=np.zeros((12,)*4,float)
    audit={}
    for group,sector in sectors.items():
        info=sector['fit']
        r=sector['reference'];n=sector['native']
        if not info['reference_in_native_span'] or not info['native_in_reference_span']:
            audit[group]={'mode':'RetainRGBetaReference_BasisIncomplete',
                          'reference_span_in_archive':False}
            continue
        native_names=list(n);ref_names=list(r)
        g=np.array([_input_coupling(couplings,k) for k in ref_names])
        B=np.asarray(info['reference_to_native_coefficients'],float)
        native_coeff=B@g
        native_tensor=sum((native_coeff[i]*n[name] for i,name in enumerate(native_names)),
                          np.zeros((12,)*4,float))
        ref_tensor=sum((g[i]*r[name] for i,name in enumerate(ref_names)),
                       np.zeros((12,)*4,float))
        residual=relative(native_tensor,ref_tensor)
        require(residual<1e-8,
                f'{group}: source-native reconstruction differs from reference')
        correction+=native_tensor-ref_tensor
        audit[group]={'mode':'SourceNativeReconstructionOfRGBetaReference',
                      'reference_span_in_archive':True,
                      'native_coefficients_in_archive_order':native_coeff.tolist(),
                      'reconstruction_relative_residual':residual,
                      'native_unique':info['unique_bijective_coordinate_map']}
    return correction,audit


def _enc(value:Any)->dict:
    from Numerical.parallel_t3.complete.two_model_master_wilson_batch13 import encode_complex
    return encode_complex(value)


def run_one(root:Path,model:str,*,output:Path|None=None,
            include_smeft:bool=True)->dict:
    """Reconstruct resolved quartics; recompute fixed-order seven-channel beta.

    Recompute is a REPRESENTATION-equivalence check, not a new model fit.
    """
    from Numerical.parallel_t3.complete.two_model_joint_scalar_batch9 import load_snapshot
    from Numerical.parallel_t3.complete.two_model_master_wilson_batch13 import (
        B11, decode_complex, gauge_sectors, sm_yukawa_tensor,
        scalar_quartic_from_matchete, add_d2_extra_hhss,
        add_d2_extra_scalar_only, projection, replace_eft1_log)
    from Numerical.running.Aminus4EFT1Couplings import F_LABELS
    from Numerical.parallel_t3.wilson.adm import beta, SYMMETRIC_TERMS
    from Numerical.parallel_t3.complete.native_triplet_basis_batch15 import heavy_source_from_reference
    root=Path(root).resolve()
    require(model in MODELS,'unsupported model')
    b13_file=root/ROOT_B13/model/'result.json'
    npz_file=root/ROOT_B13/model/'full_tensor_candidate.npz'
    b11_file=root/B11/model/'result.json'
    require(b13_file.is_file() and npz_file.is_file() and b11_file.is_file(),
            'missing successful Batch13 or Batch11 source results')
    b13=json.loads(b13_file.read_text(encoding='utf-8-sig'))
    b11=json.loads(b11_file.read_text(encoding='utf-8-sig'))
    require(b13.get('model')==model and
            b13.get('status')=='FullTensorSevenChannelStrictOneLoopConditional_NotPhysicalC5',
            'Batch13 tensor result not successful')
    snap=load_snapshot(root,model)
    seed=json.loads(snap.quartic_seed.read_text(encoding='utf-8-sig'))
    map_info,sectors_fit,basis=build_map(seed,model)
    corr,representation=rebuilt_quartic_delta(
        sectors_fit,snap.trajectory.at_fermion_threshold())
    native_R,basis_info=heavy_source_from_reference(model,seed)
    mass_R=np.eye(8)
    if model=='D_-2':
        from Numerical.running.GeneralT3DynamicFermionWilson import (
            uv_from_json,nonmixed_boundary)
        from Numerical.joint_scalar_matching.special_mass_v6 import build_special_basis
        uv=uv_from_json(model,snap.uv_values,mu_f=snap.running.mu_f)
        source,_=nonmixed_boundary(model,uv,repo_root=root,mu_f=snap.running.mu_f)
        mass_R=np.asarray(build_special_basis(model,snap.running,source)[
            'R_original_from_boundary'],float)
    fstate=snap.trajectory.at_fermion_threshold()
    base_q,base_info=scalar_quartic_from_matchete(
        seed,fstate,model,mass_R,native_source_from_reference=native_R)
    if model=='D_-2':
        add_d2_extra_hhss(base_q,root,fstate,mass_R)
        add_d2_extra_scalar_only(base_q,fstate,mass_R)
    mass_U=np.eye(12); mass_U[4:,4:]=mass_R
    corrected_q=base_q+transform_tensor(corr,mass_U)
    with np.load(npz_file,allow_pickle=False) as data:
        c0=np.asarray(data['C_tree_muF'],complex)
        b13_delta=np.asarray(data['delta_candidate_total'],complex)
    require(c0.shape==b13_delta.shape==(45,45,12,12),
            'Batch13 full tensor size changed')
    yukawa=sm_yukawa_tensor(fstate,tuple(F_LABELS))
    gauges=gauge_sectors(fstate,model,mass_R)
    mode=b13.get('candidate_ADM')
    require(mode in ('fonseca_eq22_six_term','ward_covariant_candidate'),
            'unrecognized Batch13 seven-channel ADM variant')
    # A previously saved source SHA is stricter than assuming the new source
    # has the same contents. Do not do a second expensive ADM calculation.
    source_sha=sha256(snap.quartic_seed.read_bytes()).hexdigest()
    require(source_sha==b13['input_sha256']['source_scalar_quartic'],
            'scalar source SHA changed after Batch13; rerun Batch13 first')
    dt=math.log(snap.running.mu_s/snap.running.mu_f)
    channels,backend=beta(c0,yukawa,corrected_q,gauges,variant=mode)
    require(set(channels)==set(SYMMETRIC_TERMS),'lost seven-channel ADM term')
    delta_by_channel={k:LOOP*dt*v for k,v in channels.items()}
    delta_full=sum(delta_by_channel.values(),np.zeros_like(c0))
    original_residual=relative(delta_full,b13_delta)
    require(original_residual<1e-8,
            'reconstructed quartic changes Batch13 unexpectedly: '
            f'{original_residual:.6g}; inspect map before use')
    c5_log=projection(delta_full,tuple(F_LABELS))
    b13_log=decode_complex(b13['master_candidate_EFT1_log_GeVinv'])
    native_log=decode_complex(b11['source_EFT1_scale_partition_log_GeVinv'])
    native_finite=decode_complex(b11['source_finite_sector_candidate_GeVinv'])
    c5_muS=replace_eft1_log(native_finite,native_log,c5_log,
                            float(b13['master_to_native_factor_assumed']))
    result={
        'format':FORMAT,'model':model,
        'status':'ScalarCoordinateSpanReconstructed_ConditionalNotPhysicalC5',
        'archived_native_vs_RGBeta_coordinate_sectors':map_info,
        'reconstructed_reference_sector_injection':representation,
        'unresolved_reference_directions':map_info['reference_directions_with_no_archived_span'],
        'native_topology_phase_test':map_info['topology_real_imaginary'],
        'scalar_quartic_reference_candidate_info':base_info,
        'scalar_tensor_reconstruction_relative':relative(corrected_q,base_q),
        'strict_fixed_order':'dt=ln(muS/muF), one beta at muF, no double-counting of B11 log',
        'batch13_reproduction_relative':original_residual,
        'full_tensor_beta_change_relative':relative(delta_full,b13_delta),
        'master_EFT1_log_change_relative':relative(c5_log,b13_log),
        'master_EFT1_log_recomputed_GeVinv':_enc(c5_log),
        'conditional_C5_at_muS_GeVinv':_enc(c5_muS),
        'master_native_log_discrepancy_relative':relative(
            float(b13['master_to_native_factor_assumed'])*c5_log,native_log),
        'recomputed_candidate_channel_logs_GeVinv':{
            k:_enc(projection(v,tuple(F_LABELS))) for k,v in delta_by_channel.items()},
        'seven_channel_backend':backend,'candidate_ADM':mode,
        'native_to_RGBeta_complete_coordinate_map_certified':False,
        'full_EFT1_ADM_certified':False,
        'absolute_GammaCC_to_Weyl_certified':False,
        'complete_finite_one_loop_matching_certified':False,
        'physical_C5_certified':False,
        'new_fitting_performed':False,
        'source_CG_triplet_info':basis_info,
    }
    if include_smeft:
        from Numerical.parallel_t3.complete.two_model_conditional_uvfit_batch12 import (
            candidate_c5_as_smeft,BENCHMARK)
        from Numerical.parallel_t3.complete.full_general_completion import (
            finish_smeft,existing_evidence)
        b1=root/Path('output/parallel_t3/full_general_batch1/thesis_full_general_batch1_v1')/model
        b2=root/Path('output/parallel_t3/full_general_batch2/thesis_full_general_batch2_raw_v1')/model
        laurent=root/Path('output/parallel_t3/full_general_batch2/thesis_full_general_batch2_laurent_v2')/model
        evidence=existing_evidence(model,b1,b2,laurent)
        cfg=json.loads((root/BENCHMARK).read_text(encoding='utf-8-sig'))
        target=json.loads((root/cfg['fitting']['target']).read_text(encoding='utf-8-sig'))
        result['conditional_100GeV_neutrinos']=finish_smeft(
            evidence,candidate_c5_as_smeft(c5_muS,
                float(b13.get('signed_factor_assumed',1.))),
            vev=float(cfg['vev_GeV']),target=target)
    if output is not None:
        output=Path(output)
        output.mkdir(parents=True,exist_ok=True)
        (output/'result.json').write_text(json.dumps(result,indent=2,sort_keys=True,
            allow_nan=False,default=str)+'\n',encoding='utf-8')
        np.savez_compressed(output/'recomputed_EFT1_delta.npz',
            C_tree_muF=c0,delta_candidate_total=delta_full,
            delta_EFT1_log_C5=c5_log,quartic_reconstructed=corrected_q,
            quartic_original_batch13=base_q,
            **{'delta_'+k:v for k,v in delta_by_channel.items()})
    return result
