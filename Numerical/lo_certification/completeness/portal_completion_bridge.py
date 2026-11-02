"""Read-only Workstream2 V3/V4 HHSS handoff to Workstream3 conditional portal.

This is NOT a physical full UV -> SMEFT matching certificate.  The exact V3
baseline tensor directions are numerical in a declared original Cartesian frame.
The extra V3 reference invariants lack a certified named coupling dictionary.
The V5/V6 scalar backgrounds are explicitly exploratory.  Never assign a
numerical coefficient to an extra invariant or claim a GammaCC sign.

s_original = R @ s_boundary, hence Q_boundary[ABab] = R[ia] R[jb] Q_original[ABij].
All heavy fermions are integrated together at muF and both scalar multiplets
are retained to one joint muS.  No separate scalar thresholds.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping
import math
import numpy as np
import sympy as sp

from .source_bridge import require, digest
from .hard_matching import evaluate_hard_scalar_channel

V3_FORMAT='t3-ws2-complete-baseline-portal-and-extra-reference-batch-v3'
EXTRA_FORMAT='t3-ws2-extra-hhss-exact-su2-reference-space-v1'
V4_FORMAT='t3-ws2-source-named-rgbeta-HHSS-projection-v4'


@dataclass(frozen=True)
class CompletedPortal:
    model: str
    names: tuple[str,...]
    mu_f_channels_original: np.ndarray
    mu_s_channels_original: np.ndarray
    baseline_missing_native_directions: int
    native_names: tuple[str,...]
    source_hashes: dict[str,str]


def _finite_real_array(obj: Any, shape: tuple[int,...], description: str) -> np.ndarray:
    a=np.asarray(obj)
    require(a.shape==shape and np.isrealobj(a) and np.isfinite(a).all(),
            f'{description}: real finite shape {shape} required')
    return np.asarray(a,dtype=float)


def read_ws2_completed_portals(status: Mapping[str,Any], arrays: Mapping[str,Any], *,
        model: str, source_hash: str, running_hash: str,
        n: int, native_names: tuple[str,...]) -> CompletedPortal:
    """Check exact V3 baseline portal component arrays and both EFT1 scales.

    Existing native portal directions must be listed exactly once; no extras,
    missing quartics or silent zeros. Other native mixing channels live outside
    the baseline portal array and must be preserved separately.
    """
    require(status.get('format')==V3_FORMAT and status.get('model')==model,
            'Workstream2 completion source format/model mismatch')
    require(status.get('status')=='exact_baseline_portal_completed_extra_reference_basis_not_physically_mapped',
            'V3 baseline completion not established')
    require(status.get('source_boundary_sha256')==source_hash and
            status.get('source_sha256',{}).get('running')==running_hash,
            'V3 source boundary or running hash differs from Workstream3')
    require(status.get('baseline_portal_exact_tensor_reconstruction_at_both_scales') is True,
            'V3 exact portal reconstruction has not been verified at both scales')
    require(status.get('safe_for_existing_joint_scalar_matching_as_full_input') is False and
            status.get('full_HHSS_operator_tensor_and_couplings_certified') is False and
            status.get('physical_C5_calculated') is False,
            'Expected CONDITIONAL V3 source; refusing changed physical claims')
    comp=status.get('baseline_portal_completion',{})
    require(comp.get('format')=='t3-ws2-portal-completed-native-plus-rgbeta-v3' and
            comp.get('higher_HHSS_and_lambdaT3_included') is False,
            'Baseline V3 must exclude nonportal lambdaT3/extra HHSS')
    require(comp.get('source_scalar_basis')==
            'original S1/S2 canonical real components, interleaved Re/Im',
            'Unexpected V3 original real scalar basis')
    require(comp.get('quartic_normalization')==
            'V4=lambda_abcd phi_a phi_b phi_c phi_d/4!' and
            comp.get('potential_sign')=='V4=-L4 quartic interaction',
            'V3 scalar potential convention changed')
    species=comp.get('species')
    require(isinstance(species,list) and len(species)==2 and
            {x.get('species') for x in species}=={'S1','S2'},
            'Expected both heavy scalar multiplets in completed portal')
    required_native=[]; all_names=[]; ncomplement=0
    for x in species:
        require(x.get('status')=='exact_full_portal_tensor_basis_in_declared_original_component_frame',
                'Species is not an exact tensor completion')
        native=x.get('native_symbols',[]); complemented=x.get('added_complement_symbols',[])
        augmented=x.get('augmented_symbols',[])
        require(isinstance(native,list) and isinstance(complemented,list) and
                augmented==native+complemented and
                x.get('completion_count')==len(complemented),
                'Native/complemented V3 direction list inconsistent')
        required_native.extend(native);all_names.extend(augmented);ncomplement+=len(complemented)
    require(len(required_native)==len(set(required_native)) and
            set(required_native)==set(native_names),
            'V3 baseline native symbol inventory differs from old V5 native channel set')
    require(ncomplement==comp.get('baseline_missing_native_directions_filled') and
            len(all_names)==len(set(all_names)),
            'V3 complemented direction count/labels inconsistent')
    scales=comp.get('scales',{})
    result={}
    for key in ('muF','muS'):
        metadata=scales.get(key,{})
        require(metadata.get('status')=='exact_portal_numeric_reconstruction_conditional_component_frame' and
                metadata.get('missing_rgbeta_couplings')==[] and
                metadata.get('channels_included')==all_names,
                f'V3 {key} missing or reordered tensor directions')
        lab=np.asarray(arrays[f'{key}_channel_labels'])
        require(lab.ndim==1 and lab.dtype.kind in 'US' and lab.tolist()==all_names,
                f'V3 {key} numeric NPZ channel labels mismatch')
        channels=_finite_real_array(arrays[f'{key}_channel_tensors'],
                                    (len(all_names),4,4,n,n),f'V3 {key} channels')
        total=_finite_real_array(arrays[f'{key}_portal_sum'],(4,4,n,n),f'V3 {key} sum')
        require(np.allclose(channels.sum(axis=0),total,rtol=2e-12,atol=1e-12),
                f'V3 {key} channels do not reconstruct baseline portals')
        require(np.allclose(channels,channels.swapaxes(1,2),atol=1e-12) and
                np.allclose(channels,channels.swapaxes(3,4),atol=1e-12),
                f'V3 {key} HHSS bosonic symmetry failed')
        result[key]=channels
    return CompletedPortal(model,tuple(all_names),result['muF'],result['muS'],
                           ncomplement,tuple(required_native),dict(status.get('source_sha256',{})))


def original_to_boundary(model:str, n:int, *, R:np.ndarray|None,
                         allow_exploratory_special:bool) -> np.ndarray:
    """Never silently assume original and split/mixed named scalar bases agree."""
    special=model in ('A_-2','A_+0','B_-1','C_-1','D_-2','D_+0','E_-2','E_+0')
    if special:
        require(allow_exploratory_special, 'Special model requires explicit exploratory basis opt-in')
        require(R is not None,'Missing actual V6 original-from-boundary scalar rotation')
    else:
        require(R is None or np.allclose(R,np.eye(n),rtol=0,atol=1e-12),
                'Ordinary model cannot silently rotate original component basis')
        R=np.eye(n)
    r=_finite_real_array(R,(n,n),'original-from-boundary R')
    require(np.allclose(r.T@r,np.eye(n),atol=1e-10,rtol=1e-10),
            'Original-to-boundary scalar change must be orthogonal')
    return r


def rotate_hhss(channels:np.ndarray,R:np.ndarray)->np.ndarray:
    """Convert source-original HHSS tensor channels to boundary source fields."""
    return np.einsum('cABij,ia,jb->cABab',channels,R,R,optimize=True)


def merge_completed_baseline(old_names:tuple[str,...],old_channels:np.ndarray,
                             ws2:CompletedPortal,R:np.ndarray)->tuple[tuple[str,...],np.ndarray,dict]:
    """Replace only old baseline portal channels, retain every nonportal channel.

    This is NOT old_total + completed_baseline: that would DOUBLE COUNT the
    original baseline. The replacement deliberately discards provisional V5
    native portal coefficients in favour of V3 numerical baseline tensors.
    """
    require(len(old_names)==len(set(old_names)) and
            np.asarray(old_channels).shape==(len(old_names),4,4,R.shape[0],R.shape[0]),
            'Old V5/V6 channel inventory/shape invalid')
    require(set(ws2.native_names).issubset(set(old_names)),
            'Cannot identify every existing native baseline portal for replacement')
    retained=[(name,np.asarray(old_channels)[i]) for i,name in enumerate(old_names)
              if name not in set(ws2.native_names)]
    require(set(ws2.names).isdisjoint({x for x,_ in retained}),
            'WS2 completed baseline overlaps preserved native mixing channel')
    replacement=rotate_hhss(ws2.mu_s_channels_original,R)
    names=ws2.names+tuple(name for name,_ in retained)
    tensors=np.concatenate([replacement,
                 np.stack([x for _,x in retained]) if retained else
                 np.empty((0,)+replacement.shape[1:],dtype=float)],axis=0)
    q=np.asarray(tensors,dtype=float)
    require(np.isfinite(q).all() and np.allclose(q,q.swapaxes(1,2),atol=1e-12) and
            np.allclose(q,q.swapaxes(3,4),atol=1e-12),
            'Combined completed baseline/nonportal channels have wrong symmetry')
    return names,q,{'replaced_native_baseline_names':list(ws2.native_names),
                    'preserved_nonportal_names':[name for name,_ in retained],
                    'added_baseline_directions':ws2.baseline_missing_native_directions,
                    'channel_count_before':len(old_names),'channel_count_after':len(names),
                    'old_vs_new_tensor_max_absolute':float(np.max(abs(q.sum(axis=0)-old_channels.sum(axis=0))))}


def _exact_value(v:str) -> float:
    require(isinstance(v,str) and len(v)<96, 'Reference tensor entry must be an exact short string')
    import re
    require(re.fullmatch(r'[\d\s+*/().\-sqrt]+',v) is not None,
            'Unsupported reference exact expression')
    x=sp.sympify(v,locals={'sqrt':sp.sqrt},evaluate=True)
    require(not x.free_symbols and x.is_real is True,'Reference expression must be real')
    n=float(x.evalf(22))
    require(math.isfinite(n),'Nonfinite reference entry')
    return n


def parse_extra_reference(doc:Mapping[str,Any],*,model:str,n:int)->list[dict]:
    """Read V3 exact unit Re/Im reference tensors without assigning couplings."""
    require(doc.get('format')==EXTRA_FORMAT and doc.get('model')==model and
            doc.get('physical_C5_certified') is False and
            doc.get('numerical_extra_couplings_available') is False,
            'V3 reference tensor status/model invalid; no physical coupling allowed')
    refs=[]
    sectors=doc.get('monomial_sectors',[])
    require(isinstance(sectors,list),'Missing V3 extra monomial sectors')
    for sector in sectors:
        symbols=sector.get('source_registry_symbols',[])
        rank=sector.get('independent_complex_invariant_count')
        basis=sector.get('exact_reference_real_component_tensor_basis',[])
        require(isinstance(rank,int) and rank>0 and rank==len(basis),
                'V3 source reference basis rank mismatch')
        require(isinstance(symbols,list) and symbols and
                all(isinstance(x,dict) and isinstance(x.get('rgbeta_symbol'),str) for x in symbols),
                'Missing V3 extra registry symbol names')
        for basis_element in basis:
            basis_id=basis_element.get('id')
            require(isinstance(basis_id,str) and basis_id,'Unnamed extra reference direction')
            record={'sector_counts':sector.get('counts_H_BarH_S1_BarS1_S2_BarS2'),
                    'reference_id':basis_id,'unassigned_source_symbols':
                    [x['rgbeta_symbol'] for x in symbols], 'tensors':{}}
            for part,key in (('Re','Re_c_unit_tensor_components'),('Im','Im_c_unit_tensor_components')):
                rows=basis_element.get(key,[])
                require(isinstance(rows,list),'Missing V3 unit component list')
                tensor=np.zeros((4,4,n,n),float)
                for row in rows:
                    idx=row.get('H1_H2_S1_S2_0based')
                    require(isinstance(idx,list) and len(idx)==4 and
                            all(type(x)==int and 0<=x<limit for x,limit in zip(idx,(4,4,n,n))),
                            'V3 reference component out of range')
                    a,b,i,j=idx;v=_exact_value(row.get('value_exact'))
                    for z in {(a,b,i,j),(b,a,i,j),(a,b,j,i),(b,a,j,i)}:
                        if tensor[z]!=0 and not np.isclose(tensor[z],v,rtol=1e-12,atol=1e-12):
                            raise ValueError('Conflicting conjugate/bosonic reference orbit')
                        tensor[z]=v
                record['tensors'][part]=tensor
            refs.append(record)
    require(len(set(x['rgbeta_symbol'] for s in sectors for x in s['source_registry_symbols']))==
            doc.get('registry_HHSS_extra_symbols'),
            'V3 extra HHSS registry symbol count differs from reference sectors')
    return refs


def crosscheck_v4_projections(doc:Mapping[str,Any],refs:list[dict],*,model:str)->dict:
    """Inventory V4 named source projections without identifying physical couplings."""
    require(doc.get('format')==V4_FORMAT and doc.get('model')==model and
            doc.get('two_thresholds_only') is True and
            doc.get('physical_named_Matchete_RGBeta_maps_certified')==0 and
            doc.get('full_physical_eft1_quartics_certified') is False,
            'V4 source projection status must be conditional')
    named={name for ref in refs for name in ref['unassigned_source_symbols']}
    entries=doc.get('entries')
    require(isinstance(entries,list) and {x.get('rgbeta_symbol') for x in entries}==named,
            'V4 and V3 do not describe identical extra HHSS symbols')
    valid=sum(x.get('in_exact_canonical_SU2_reference_span') is True
              for e in entries for x in e.get('candidate_orientation_projection_maps',[]))
    return {'extra_named_symbols':len(named),
            'source_named_symbol_evidence':[{'rgbeta_symbol':e['rgbeta_symbol'],
              'muF_source_coefficient_present':e.get('numerical_running_coefficients',{}).get('muF',{}).get('symbol_present') is True,
              'muS_source_coefficient_present':e.get('numerical_running_coefficients',{}).get('muS',{}).get('symbol_present') is True,
              'physical_coefficient_accepted':False,
              'number_of_candidate_orientation_maps':len(e.get('candidate_orientation_projection_maps',[]))}
              for e in entries],
            'source_projection_count':doc.get('source_projection_count'),
            'valid_orientation_candidate_count':valid,
            'physical_named_maps_certified':0,
            'orientation_ambiguous_symbols':doc.get('orientation_ambiguous_symbols'),
            'source_coefficients_inserted_into_eft1':False}


def response(formal:Any,mass2:np.ndarray,mu_s:float,names:tuple[str,...],channels:np.ndarray,
             *,projection_tolerance:float=1e-8)->dict:
    """Conditional source-portal response, not complete EFT1 seven-term running."""
    ans=evaluate_hard_scalar_channel(formal.c_by_fermion,formal.fermion_masses_gev,
        channels.sum(axis=0),mass2,formal.mu_f_gev,mu_s,
        C_total_muF=formal.c_total,lambda_channels=channels,channel_names=names,
        projection_tolerance=projection_tolerance)
    return {'projection_passed':ans['projection_passed'],
            'residuals':{k:float(v) for k,v in ans['residuals'].items()},
            'mass_hierarchy':ans['mass_hierarchy_max_scalar_over_min_fermion'],
            'matrices':{name:{'real':m.real.tolist(),'imag':m.imag.tolist()}
                        for name,m in ans['matrices'].items()},
            'hard_norm':float(np.linalg.norm(ans['matrices']['hard_muF'])),
            'dim5_norm':float(np.linalg.norm(ans['matrices']['dim5_sum'])),
            'full_triangle_norm':float(np.linalg.norm(ans['matrices']['full_triangle'])),
            'higher_dimension_remainder_norm':float(np.linalg.norm(ans['matrices']['higher_dimension_remainder'])),
            'per_named_channel':{
               label:{'hard_frobenius':float(np.linalg.norm(d['matrices']['hard_muF'])),
                      'dim5_frobenius':float(np.linalg.norm(d['matrices']['dim5_sum'])),
                      'projection_residuals':d['projection_residuals']}
               for label,d in ans['quartic_channels'].items()},
            'physical_C5_certified':False}


def extra_unit_response(formal:Any,mass2:np.ndarray,mu_s:float,references:list[dict],R:np.ndarray)->list[dict]:
    """Unit-strength Re/Im reference sensitivities, not named physical couplings.

    Never sum across references or multiply by a source named coupling: the
    named-to-reference mapping and phase are not source-certified.
    """
    entries=[]
    for ref in references:
        item={'reference_id':ref['reference_id'],
              'unassigned_source_symbols':ref['unassigned_source_symbols'],
              'sector_counts':ref['sector_counts'],'physical_coefficient_assigned':False,
              'reference_basis_to_named_RGBeta_map_certified':False,'unit_parts':{}}
        for part,q in ref['tensors'].items():
            rotated=rotate_hhss(q[None,...],R)[0]
            if not np.any(rotated):
                item['unit_parts'][part]={'status':'IDENTICALLY_ZERO_REFERENCE_TENSOR',
                                         'physical_C5_certified':False}
                continue
            try:
                r=evaluate_hard_scalar_channel(formal.c_by_fermion,formal.fermion_masses_gev,
                          rotated,mass2,formal.mu_f_gev,mu_s,C_total_muF=formal.c_total)
                item['unit_parts'][part]={
                    'status':'CONDITIONAL_UNIT_REFERENCE_RESPONSE' if r['projection_passed'] else
                             'WARD_REJECTED_DO_NOT_USE_AS_C5',
                    'projection_residuals':{k:float(v) for k,v in r['residuals'].items()},
                    'hard_matrix_projected_diagnostic_only':{
                        'real':r['matrices']['hard_muF'].real.tolist(),
                        'imag':r['matrices']['hard_muF'].imag.tolist()},
                    'dim5_matrix_projected_diagnostic_only':{
                        'real':r['matrices']['dim5_sum'].real.tolist(),
                        'imag':r['matrices']['dim5_sum'].imag.tolist()},
                    'physical_C5_certified':False}
            except (ValueError,TypeError,ArithmeticError,AssertionError) as exc:
                item['unit_parts'][part]={'status':'REFERENCE_NUMERICAL_BLOCKED',
                                          'error':f'{type(exc).__name__}: {exc}',
                                          'physical_C5_certified':False}
        entries.append(item)
    return entries
