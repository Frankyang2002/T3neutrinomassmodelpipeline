"""Workstream A v10: verify genuine native RGBeta Lam/projector response evidence.

This verifies installed-package symbolic contractions and independently rechecks
v9's native Matchete CG / *constructed* real-action audit. It does not evaluate
numerical native RGBeta scalar components, or certify the full potential.

The responses R_ij = P_i[Lam_j] form a source-conditioned observer matrix.
Do not assume projectors are mutually orthogonal; offdiagonal entries matter.
"""
from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
import re
from typing import Any

import numpy as np

from .portal_real_action_v9 import (
    EXPECTED_MODELS, DIMENSIONS, STATUS as V9_STATUS,
    audit_from_paths as recheck_v9,
)

FORMAT = 't3-parallel-a-native-rgbeta-quartic-projectors-v10'
NATIVE_STATUS = 'NativeRGBetaRealQuarticLamProjectorResponsesObserved_Conditional'
STATUS = ('NativeRGBetaRealQuarticProjector150DiagonalAnd58PortalResponsesVerified_'
          'NonportalRankAndComplexActionStillConditional')
SOURCE_FILES = (
    'RGE/running/rgbeta/T3RGBetaModel.wl',
    'scripts/RGBetaNativeAccess.wl',
    'scripts/RunT3ParallelANativeQuarticProjectorsV10.wl',
)
SPECIAL = (
    'lambdaHHdagS2S2', 'lambdaHHdagS1barS1bar',
    'lambdaS1bar2S2bar2', 'lambdaS1barS2S2bar2',
    'lambdaS1S1bar2S2bar', 'lambdaHHdagS1barS2barCross',
)
NORM = ('lambdaH','lambdaS1','lambdaS2','lambdaH1','lambdaH2','lambda12')
TOL = 2e-8


class V10Blocked(ValueError):
    """Native source, contraction, rank, or provenance not certified."""


def require(condition: bool, why: str) -> None:
    if not condition:
        raise V10Blocked('Workstream A v10: '+why)


def read(path: Path) -> dict:
    raw=json.loads(path.read_text(encoding='utf-8-sig'))
    require(isinstance(raw,dict),f'{path}: not a JSON object')
    return raw


def file_hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def expected_symbols(tag: str) -> dict[str,bool]:
    """Source-grounded T3RGBetaModel.wl EFT1 AddBasic/RGClosed/Mixing census."""
    require(tag in EXPECTED_MODELS,'unrecognized model '+str(tag))
    d1,d2,_=DIMENSIONS[tag[0]]
    alpha=int(tag.split('_')[1])
    result={name:True for name in NORM}
    result['lambdaT3']=False
    for k,d in ((1,d1),(2,d2)):
        if d>1:result[f'lambdaH{k}Adj']=True
        if d==3:result[f'lambdaS{k}Adj']=True
    if d1>1 and d2>1:result['lambda12Adj']=True
    if d1==3 and d2==3:result['lambda12Cross']=True
    if (d1,d2,alpha)==(2,2,-1):result.update({name:False for name in SPECIAL})
    return result


def _native_value(obj: Any, context: str) -> float:
    require(isinstance(obj,dict) and obj.get('ok') is True,
            context+': native Lam/projector contraction failed: '+str(obj))
    v=obj.get('value')
    require(type(v) in (float,int) and np.isfinite(v),context+': nonnumeric native contraction')
    return float(v)


def validate_native_evidence(native:dict, root:Path|None=None) -> dict:
    """Strictly validate all-16 local Wolfram evidence and registered spectrum."""
    require(native.get('format')==FORMAT,'unrecognized native format')
    require(native.get('status')==NATIVE_STATUS,'native Wolfram run did not certify its contractions')
    for flag in ('native_Lam_evaluated','native_quartic_projectors_evaluated',
                 'native_coupling_response_is_symbolic_not_component_Lam'):
        require(native.get(flag) is True,'native witness missing '+flag)
    for flag in ('native_numeric_S2_frame_measured','native_Matchete_scalar_CG_evaluated',
                 'full_scalar_potential_dictionary_certified',
                 'physical_Wilson_evolution_authorized'):
        require(native.get(flag) is False,'unsupported native certification '+flag)
    hashes=native.get('source_sha256')
    require(isinstance(hashes,dict) and set(hashes)==set(SOURCE_FILES),
            'native source hash inventory inconsistent')
    require(all(isinstance(v,str) and re.fullmatch('[0-9a-f]{64}',v)
                for v in hashes.values()),'invalid native SHA256')
    if root is not None:
        for name in SOURCE_FILES:
            p=root/name
            require(p.is_file() and hashes[name]==file_hash(p),
                    'native source has changed since export: '+name)
        package_path=native.get('installed_RGBeta_main_file')
        package_sha=native.get('installed_RGBeta_main_file_sha256')
        require(isinstance(package_path,str) and Path(package_path).is_file() and
                package_sha==file_hash(Path(package_path)),
                'installed RGBeta entrypoint differs from native evidence')
    models=native.get('models')
    require(native.get('registered_model_count')==16 and isinstance(models,dict)
            and set(models)==EXPECTED_MODELS,'16 EFT1 native models required')
    out={}; registered=real=complex_n=diags=blocks=entries=offdiags=0
    for tag in sorted(EXPECTED_MODELS):
        row=models[tag];expected=expected_symbols(tag)
        d1,d2,df=DIMENSIONS[tag[0]]
        require(row.get('ok') is True and row.get('registry_exact_source_inventory') is True,
                tag+': native model status or quartic inventory failed')
        require(row.get('dimensions')==[d1,d2,df] and
                row.get('alpha')==int(tag.split('_')[1]),tag+': EFT1 dimensions/alpha changed')
        inventory=row.get('inventory')
        require(isinstance(inventory,list) and len(inventory)==len(expected) and
                row.get('registered_quartics_count')==len(expected),
                tag+': native quartic inventory count changed')
        lookup={}
        for rec in inventory:
            symbol=rec.get('symbol')
            require(isinstance(symbol,str) and symbol in expected and symbol not in lookup,
                    tag+': unknown/duplicate quartic '+str(symbol))
            require(rec.get('native_self_conjugate') is expected[symbol],
                    tag+'/'+symbol+': native self-conjugacy differs from source')
            require(rec.get('has_installed_RGBeta_projector') is True,
                    tag+'/'+symbol+': native coupling projector absent')
            count=rec.get('native_unique_arrangement_count')
            require(type(count) is int and 1<=count<=24,
                    tag+'/'+symbol+': invalid native permutation cardinality')
            coup=re.sub(r'\s+','',rec.get('native_coupling_inputform',''))
            require(coup==f'6*Matrix[{symbol}][##1]&',
                    tag+'/'+symbol+': actual RGBeta 6*Matrix convention changed: '+coup)
            require(isinstance(rec.get('fields_inputform'),str) and
                    isinstance(rec.get('native_invariant_inputform'),str),
                    tag+'/'+symbol+': native fields/invariant absent')
            lookup[symbol]=rec
        require(set(lookup)==set(expected),tag+': native registered symbol census changed')
        drows=row.get('native_real_self_projectors')
        require(isinstance(drows,list) and len(drows)==sum(expected.values()),
                tag+': real quartic diagonal response count changed')
        seen=set()
        for dr in drows:
            sym=dr.get('symbol')
            require(sym in expected and expected[sym] and sym not in seen,
                    tag+': unsupported/duplicate real native diagonal '+str(sym))
            seen.add(sym)
            own=_native_value(dr.get('native_response'),tag+'/'+sym)
            require(abs(own-1)<TOL,tag+'/'+sym+': native diagonal projector response != 1')
            diags+=1
        require(seen=={k for k,v in expected.items() if v},
                tag+': native self-projector coverage incomplete')
        parts=row.get('native_portal_response_blocks')
        require(isinstance(parts,list) and len(parts)==2 and
                row.get('native_real_self_projectors_all_unit') is True and
                row.get('native_portal_response_blocks_evaluated') is True,
                tag+': native projector acceptance flags absent')
        portal_data=[]
        for which,dim in ((1,d1),(2,d2)):
            b=parts[which-1]
            symbols=[f'lambdaH{which}']+([f'lambdaH{which}Adj'] if dim>1 else [])
            n=len(symbols)
            require(b.get('which')==which and b.get('dimension')==dim and
                    b.get('symbols')==symbols,
                    tag+f'/S{which}: native portal registration order changed')
            raw=b.get('native_response_matrix')
            require(isinstance(raw,list) and len(raw)==n and
                    all(isinstance(x,list) and len(x)==n for x in raw),
                    tag+f'/S{which}: invalid native response matrix')
            R=np.array([[_native_value(x,tag+f'/S{which}') for x in y] for y in raw])
            require(np.max(abs(np.diag(R)-1))<TOL,
                    tag+f'/S{which}: native portal self-normalization not one')
            sv=np.linalg.svd(R,compute_uv=False)
            require(sv[-1]>1.e-8 and sv[0]/sv[-1]<1.e8,
                    tag+f'/S{which}: native portal response matrix singular')
            portal_data.append({'which':which,'dimension':dim,'symbols':symbols,
                                'native_projector_response_matrix':R.tolist(),
                                'native_projector_response_condition':float(sv[0]/sv[-1])})
            entries+=n*n;offdiags+=n*n-n;blocks+=1
        real+=sum(expected.values());complex_n+=sum(not v for v in expected.values())
        registered+=len(expected)
        out[tag]={'dimensions':[d1,d2,df],
                  'native_quartics_registered':len(expected),
                  'native_real_projector_diagonal_count':len(drows),
                  'native_complex_quartics_inventory_only':sorted(k for k,v in expected.items() if not v),
                  'portal_sectors':portal_data}
    require((registered,real,complex_n,diags,blocks,entries,offdiags)==
            (178,150,28,150,32,110,52),
            'native expected census/response counts changed: '+str(
             (registered,real,complex_n,diags,blocks,entries,offdiags)))
    return {'models':out,'registered_quartic_slots':registered,
            'native_real_projector_diagonals_verified':diags,
            'native_portal_response_blocks':blocks,
            'native_portal_response_matrix_elements':entries,
            'native_portal_offdiagonal_responses':offdiags,
            'complex_quartic_slots_inventory_only':complex_n}


def compare_with_v9(native_rows:dict,v9:dict)->dict:
    require(v9.get('status')==V9_STATUS and
            v9.get('native_Matchete_portal_coupling_slots')==58 and
            v9.get('triplet_sectors_spectral_normalization_verified')==14 and
            v9.get('anti_Hermitian_Matchete_phase_slots_checked')==14,
            'v9 source/conditional portal evidence not certified')
    require(v9.get('certification',{}).get('native_RGBeta_Lam_projector_numerically_executed') is False,
            'v9 incorrectly claims a native RGBeta projector')
    out={};worst=0.0
    for tag in sorted(EXPECTED_MODELS):
        rows=v9.get('models',{}).get(tag,[])
        require(isinstance(rows,list) and len(rows)==2,tag+': v9 sector inventory changed')
        data=[]
        for which in (1,2):
            v9r=rows[which-1];nrow=native_rows['models'][tag]['portal_sectors'][which-1]
            require(v9r.get('which')==which and
                    v9r.get('rgbeta_registered_couplings')==nrow['symbols'],
                    tag+f'/S{which}: v9 vs native registry differs')
            M=np.asarray(v9r.get('conditional_rgbeta_from_matchete_real_matrix'),float)
            R=np.asarray(nrow['native_projector_response_matrix'],float)
            n=len(nrow['symbols'])
            require(M.shape==(n,n) and np.all(np.isfinite(M)) and
                    np.linalg.matrix_rank(M,tol=1.e-8)==n,
                    tag+f'/S{which}: v9 Matchete-to-RGBeta conditional map singular')
            obs=R@M
            error=float(v9r.get('max_relative_residual',1.))
            require(np.isfinite(error) and error<2e-10,
                    tag+f'/S{which}: v9 action residual not certified')
            worst=max(worst,error)
            data.append({'which':which,'source_Matchete_couplings':v9r['native_matchete_couplings'],
              'native_RGBeta_couplings':nrow['symbols'],
              'conditional_RGBeta_from_Matchete':M.tolist(),
              'native_projector_observations_from_Matchete':obs.tolist(),
              'native_projector_deconvolution_matrix':np.linalg.inv(R).tolist(),
              'conditional_not_physical':True})
        out[tag]=data
    return {'models':out,'max_v9_source_action_residual':worst}


def audit(native:dict,v9:dict,root:Path|None=None)->dict:
    coverage=validate_native_evidence(native,root)
    bridge=compare_with_v9(coverage,v9)
    return {
        'format':'t3-parallel-a-native-quartic-projectors-conditional-audit-v10',
        'status':STATUS,
        'registered_quartic_slots':coverage['registered_quartic_slots'],
        'native_real_projector_diagonals_verified':coverage['native_real_projector_diagonals_verified'],
        'native_portal_response_blocks':coverage['native_portal_response_blocks'],
        'native_portal_response_matrix_elements':coverage['native_portal_response_matrix_elements'],
        'native_portal_offdiagonal_responses':coverage['native_portal_offdiagonal_responses'],
        'complex_quartic_slots_inventory_only':coverage['complex_quartic_slots_inventory_only'],
        'max_v9_source_action_residual':bridge['max_v9_source_action_residual'],
        'models':bridge['models'],
        'native_response_metadata':coverage['models'],
        'certification':{
            'installed_RGBeta_symbolic_Lam_and_Projector_executed':True,
            'native_RGBeta_real_quartic_projector_self_normalization_verified':True,
            'native_RGBeta_portal_projector_response_invertible':True,
            'native_RGBeta_numerical_component_Lam_evaluated':False,
            'native_Matchete_to_RGBeta_component_basis_fixed':False,
            'nonportal_Matchete_coupling_dictionary_complete':False,
            'complex_quartic_projector_conjugate_channel_verified':False,
            'BCminus1_mixed_scalar_mass_basis_certified':False,
            'signed_GammaCC_to_Weyl_certified':False,
            'physical_Wilson_evolution_authorized':False,
            'EFT_order':'UV -> integrate out ALL fermions -> integrate out ALL scalars -> SMEFT',
        },
    }


def audit_from_paths(root:Path,native_path:Path,v9_path:Path)->dict:
    """Re-run entire local v9 source chain and verify saved v9 numerical rows."""
    root=root.resolve()
    basedir=root/'output/parallel_t3/wilson'
    v9paths=[basedir/x for x in (
        'matchete_portals_native_v5.json','rgbeta_portals_native_v5.json',
        'portal_dictionary_conditional_v5.json','rgbeta_fundamental_fierz_native_v6.json',
        'portal_orientation_conditional_v6.json','s2_intertwiner_conditional_v7.json',
        'rgbeta_s2_reduction_native_v8.json','s2_spectrum_conditional_v8.json')]
    v9_actual=recheck_v9(*v9paths,root,samples=3)
    saved=read(v9_path)
    require(saved.get('format')==v9_actual['format'] and
            saved.get('status')==v9_actual['status'] and
            saved.get('input_sha256')==v9_actual['input_sha256'],
            'saved v9 provenance/status changed')
    require(saved.get('native_Matchete_portal_coupling_slots')==58 and
            saved.get('max_constructed_real_action_relative_residual',1.)<2e-10,
            'saved v9 has incomplete real-action evidence')
    for tag in EXPECTED_MODELS:
        old=saved['models'][tag];now=v9_actual['models'][tag]
        for a,b in zip(old,now):
            require(a.get('which')==b['which'] and
                    np.allclose(a.get('conditional_rgbeta_from_matchete_real_matrix'),
                                b['conditional_rgbeta_from_matchete_real_matrix'],atol=2e-10),
                    tag+': saved v9 coupling map has been altered')
    out=audit(read(native_path),v9_actual,root)
    out['input_sha256']={
        'native_quartic_projectors_v10':file_hash(native_path),
        'portal_real_action_v9':file_hash(v9_path),
        **{k:file_hash(v) for k,v in zip((
            'matchete_v5','rgbeta_v5','portal_v5','fierz_v6','orientation_v6',
            'intertwiner_v7','reduction_v8','spectrum_v8'),v9paths)},
    }
    return out
