"""Workstream A v18/v19: live Matchete complex CGs + native RGBeta v15 channels.

The baseline EFT1 contains 16 lambdaT3 and 12 B/C(-1) special complex
quartics. In B/C, the OPTIONAL Matchete seven-dimensional extra-quartic
space contains one direction absent from the baseline six-coupling RGBeta
registry. A fitted conditional raw-monomial map is never a physical map.
"""
from __future__ import annotations
from hashlib import sha256
from pathlib import Path
from typing import Any
import json
import numpy as np

from .matchete_rgbeta_dictionary import (
    ALL16, matchete_number, solve_from_paths, DictionaryEvidenceError,
    OUT_FORMAT as V3_FORMAT,
)
from .scalar_topology_action_v12 import (
    topology_constructed_real_lam, topology_source_complex_action,
)
from .native_complex_channel_v15 import audit_from_native as native_v15_audit, STATUS as V15_STATUS
from .native_context_quartic_v14 import read as read_v14
from .native_quartic_projectors_v10 import DIMENSIONS

FORMAT = 't3-parallel-a-complex-scalar-28-source-bridge-v19'
STATUS = 'Matchete16TopologyCGAndBC12BaselineComplexSubspaceChecked_NativeRGBeta28ConjugateChannelsObserved_Conditional'
MAT_FORMAT = 't3-parallel-a-matchete-complex-cg-v18'
MAT_STATUS = 'NativeMatcheteComplexCGSourceOnly_MassAndBasisConditional'
EXPORT = 'scripts/RunT3ParallelAExportComplexMatcheteV18.wl'
SOURCES = (
    'Lagrangian/model/T3ModelCatalog.wl',
    'Lagrangian/model/LagrangianBuilder.wl',
    'Lagrangian/model/SU2Invariants.wl',
    'Lagrangian/model/T3Fields.wl',
    'Lagrangian/interactions/T3Topology.wl',
    'Lagrangian/interactions/ScalarPotential.wl',
    'Lagrangian/interactions/GeneralScalarOperators.wl',
    'Lagrangian/interactions/GeneralScalarVerifiedQuartics.wl',
    EXPORT,
)
BC = ('B_-1', 'C_-1')
SPECIAL = (
    'lambdaHHdagS2S2', 'lambdaHHdagS1barS1bar',
    'lambdaS1bar2S2bar2', 'lambdaS1barS2S2bar2',
    'lambdaS1S1bar2S2bar', 'lambdaHHdagS1barS2barCross',
)

class V19Blocked(ValueError):
    pass

def require(ok: bool, reason: str) -> None:
    if not ok:
        raise V19Blocked('Workstream A v19: '+reason)

def sha(path: Path) -> str:
    return sha256(Path(path).read_bytes()).hexdigest()

def load(path: Path) -> dict:
    obj=json.loads(Path(path).read_text(encoding='utf-8-sig'),parse_constant=lambda s: (_ for _ in ()).throw(V19Blocked('Nonfinite JSON '+s)))
    require(isinstance(obj,dict),'not a JSON object: '+str(path))
    return obj

def native_document(doc: dict, root: Path, mode: str) -> dict:
    require(doc.get('format')==MAT_FORMAT and doc.get('status')==MAT_STATUS and
            doc.get('mode')==mode and doc.get('models_checked')==(16 if mode=='TOPOLOGY' else 2) and
            doc.get('native_complex_CG_slots')==(16 if mode=='TOPOLOGY' else 14),
            'incorrect native Matchete v18 mode/census')
    require(doc.get('independent_RGBeta_native_component_frame_certified') is False and
            doc.get('BC_minus1_mixed_scalar_mass_diagonalized') is False and
            doc.get('physical_Wilson_evolution_authorized') is False,
            'unsafe native Matchete certification claim')
    hashes=doc.get('source_sha256')
    require(isinstance(hashes,dict) and set(hashes)==set(SOURCES),'incomplete native source hash list')
    for p in SOURCES:
        f=root/p
        require(f.is_file() and sha(f)==hashes[p], 'stale native Matchete source: '+p)
    expected=ALL16 if mode=='TOPOLOGY' else set(BC)
    models=doc.get('models')
    require(isinstance(models,dict) and set(models)==set(expected),'wrong native Matchete models')
    observed={}
    for tag,row in models.items():
        d1,d2,df=DIMENSIONS[tag[0]]
        require(row.get('dimensions')==[d1,d2,df] and row.get('alpha')==int(tag.split('_')[1]) and
                row.get('native_builder_status')=='Success' and
                row.get('matching_mass_basis_physically_certified') is False and
                row.get('symbolic_quartic_only') is True,
                tag+': source CG build metadata invalid')
        rr=row.get('records')
        require(isinstance(rr,list) and len(rr)==(1 if mode=='TOPOLOGY' else 7),
                tag+': missing native CG records')
        out={}
        for rec in rr:
            dims=rec.get('dimensions');entries=rec.get('components_inputform')
            require(isinstance(dims,list) and len(dims)==4 and
                    all(type(d)==int and d>0 for d in dims) and
                    rec.get('native_dimensions')==[d for d in dims if d>1] and
                    isinstance(entries,list) and len(entries)==int(np.prod(dims)) and
                    rec.get('assembled_lagrangian_contains_operator') is True,
                    tag+': missing native tensor geometry/operator proof')
            vals=np.array([matchete_number(v) for v in entries],complex).reshape(dims)
            require(np.all(np.isfinite(vals)) and np.linalg.norm(vals)>1e-12,
                    tag+': invalid native CG components')
            symbol=rec.get('symbol')
            require(isinstance(symbol,str) and symbol not in out,
                    tag+': duplicate or invalid CG symbol')
            if mode=='TOPOLOGY':
                require(symbol=='lambdaT3' and rec.get('cg')=='T3MixCG' and
                        rec.get('fields')==['H','H','S1','Bar[S2]'] and
                        dims==[2,2,d1,d2] and
                        np.linalg.norm(vals-vals.swapaxes(0,1))<1e-10,
                        tag+': topology native Higgs pair not symmetric')
            else:
                key=rec.get('counts');species=rec.get('fields')
                require(isinstance(key,list) and len(key)==6 and
                        all(type(x)==int and 0<=x<=4 for x in key) and sum(key)==4 and
                        isinstance(species,list) and len(species)==4 and
                        species==sum(([name]*num for name,num in zip(
                            ('H','Hdag','S1','S1dag','S2','S2dag'),key)),[]) and
                        symbol.startswith('lambdaT3ExtraQ') and
                        rec.get('cg')=='T3ExtraCGQ'+''.join(str(x) for x in key)+'I'+symbol.rsplit('I',1)[-1] and
                        rec.get('operator')=='T3ExtraQ_'+'_'.join(str(x) for x in key)+'_'+symbol.rsplit('I',1)[-1] and
                        dims==[2,2,2,2],
                        tag+': unexpected B/C native extra CG record')
            out[symbol]=(rec,vals)
        observed[tag]=out
    return observed

def topology_action(t: np.ndarray, *, seed: int = 18001, samples: int = 6) -> float:
    """Independent action-level check using the native, not declared, CG."""
    d1,d2=t.shape[2:];rng=np.random.default_rng(seed);worst=0.
    for _ in range(samples):
        x=rng.normal(size=4+2*d1+2*d2)
        c=complex(rng.normal(),rng.normal())
        exact=topology_source_complex_action(t,c,x)
        lam=topology_constructed_real_lam(t,c)
        reconstructed=float(np.einsum('abcd,a,b,c,d->',lam,x,x,x,x,optimize=True)/24.)
        err=abs(exact-reconstructed)/max(1.,abs(exact),abs(reconstructed))
        worst=max(worst,err)
        require(err<1e-10,'native topology PlusHc/real-Lam normalization mismatch')
    return worst

def historical_extra_map(mat:dict,bc_native:dict,old_mat:dict,old_rgb:dict,
                         v3_saved:dict) -> dict:
    """Re-run the source-conditioned 50-slot map and isolate baseline six of seven.

    Native B/C CG components must match the historical 50-slot tensor evidence
    before that old raw-invariant map can be used for current source.
    """
    result={}
    for tag in BC:
        current=bc_native[tag]
        originals=old_mat['models'][tag]['tensors']
        require(len(originals)==len(current)==7,tag+': seven-native-extra-CG inventory failed')
        for old in originals:
            symbol=old['symbol']
            require(symbol in current,tag+': historical extra CG absent from live Matchete')
            freshrec, fresh=current[symbol]
            oldshape=old['dimensions']
            require(list(fresh.shape)==oldshape and
                    freshrec['fields']==old['species'] and
                    freshrec['dimensions']==oldshape,
                    tag+'/'+symbol+': old/live CG field order differs')
            historic=np.array([matchete_number(x) for x in old['components']],complex).reshape(oldshape)
            require(np.linalg.norm(fresh-historic)<2e-9*max(1.,np.linalg.norm(fresh)),
                    tag+'/'+symbol+': current native Matchete CG differs from historical 50-slot evidence')
        row=mat['models'][tag]
        require(row['matchete_symbols']==[x['symbol'] for x in originals],
                tag+': historical Matchete CG ordering differs')
        rgb=row['rgbeta_symbols']
        require(len(rgb)==7 and set(SPECIAL).issubset(rgb) and
                len(set(rgb)-set(SPECIAL))==1 and
                next(iter(set(rgb)-set(SPECIAL))).startswith('lambdaT3GeneralQ'),
                tag+': baseline RGBeta six-special inventory differs')
        x=row['map_rgbeta_from_matchete']
        a=np.asarray(x['real'],dtype=float)+1j*np.asarray(x['imag'],dtype=float)
        require(a.shape==(7,7) and np.all(np.isfinite(a)) and
                np.max(np.abs(a.imag))<1e-9 and np.linalg.matrix_rank(a,tol=1e-10)==7,
                tag+': full optional extra-quartic map singular or non-real')
        ids=[rgb.index(s) for s in SPECIAL]
        b=a[ids]
        require(np.linalg.matrix_rank(b,tol=1e-10)==6,
                tag+': six baseline special quartics are not independent')
        _,_,vh=np.linalg.svd(b,full_matrices=True)
        null=vh[-1].conj()
        null=null/np.linalg.norm(null)
        require(np.linalg.norm(b@null)<1e-9 and
                abs(a[rgb.index(next(iter(set(rgb)-set(SPECIAL))))]@null)>1e-6,
                tag+': no independently missing optional seventh direction')
        result[tag]={'matchete_seven_source_symbols':row['matchete_symbols'],
              'baseline_rgbeta_six_symbols':list(SPECIAL),
              'optional_extra_rgbeta_symbol':next(iter(set(rgb)-set(SPECIAL))),
              'baseline_rgbeta_from_Matchete_seven_complex_matrix':{
                  'real':b.real.tolist(),'imag':b.imag.tolist()},
              'baseline_map_rank':6,'optional_map_rank':7,
              'baseline_kernel_dimension_in_Matchete_seven_space':1,
              'kernel_direction_normalized_complex':{'real':null.real.tolist(),
                                                      'imag':null.imag.tolist()},
              'source_monomial_not_native_projector_normalization':True}
    return result

def audit(root:Path,topology_path:Path,bc_path:Path,native14_path:Path,
          native15_saved:Path,mat50:Path,rgb50:Path,v3_saved:Path,
          *,samples:int=6)->dict:
    root=root.resolve()
    inputs=[topology_path,bc_path,native14_path,native15_saved,mat50,rgb50,v3_saved]
    for p in inputs:require(p.is_file(),'missing required witness '+str(p))
    topo=native_document(load(topology_path),root,'TOPOLOGY')
    bc=native_document(load(bc_path),root,'BCEXTRA')
    native15=native_v15_audit(root,native14_path)
    saved15=load(native15_saved)
    require(native15['status']==V15_STATUS and native15['smoke_only'] is False and
            native15['counts'].get('attempted')==28 and
            native15['counts'].get('unit_direct_zero_conjugate')==28 and
            saved15.get('status')==native15['status'] and
            saved15.get('counts')==native15['counts'] and
            saved15.get('native_v14_input_sha256')==sha(native14_path) and
            saved15.get('native_complex_action_physical_normalization_certified') is False and
            saved15.get('physical_Wilson_evolution_authorized') is False,
            'v15 native RGBeta complex conjugate channels not independently verified')
    for tag,entries in native15['models'].items():
        expected={'lambdaT3'} | (set(SPECIAL) if tag in BC else set())
        require({r['symbol'] for r in entries}==expected and
                all(r['classification']=='native_source_unit_conjugate_zero' for r in entries),
                tag+': baseline EFT1 native complex channel inventory != expected')
    rawmat=load(mat50); rawrgb=load(rgb50); historic=load(v3_saved)
    fresh=solve_from_paths(mat50,rgb50)
    require(fresh.get('format')==historic.get('format')==V3_FORMAT and
            fresh.get('status')==historic.get('status')=='All50RawExtraQuarticInvariantsMatched_ConditionalFrameOnly' and
            fresh.get('source_evidence_sha256')==historic.get('source_evidence_sha256') and
            fresh.get('additional_quartic_slots')==50 and
            historic.get('n_failed_models')==0,
            'saved v3 extra-quartic frame map stale or failed independently')
    for tag in BC:
        x=fresh['models'][tag]['map_rgbeta_from_matchete']
        y=historic['models'][tag]['map_rgbeta_from_matchete']
        require(np.allclose(x['real'],y['real'],atol=2e-9) and
                np.allclose(x['imag'],y['imag'],atol=2e-9),
                tag+': saved conditional 7x7 coupling map altered')
    bc_map=historical_extra_map(fresh,bc,rawmat,rawrgb,historic)
    models={};worst=0.
    for i,tag in enumerate(sorted(ALL16)):
        t=topo[tag]['lambdaT3'][1]
        err=topology_action(t,seed=18001+i,samples=samples)
        worst=max(worst,err)
        models[tag]={'native_Matchete_topology_CG_norm':float(np.linalg.norm(t)),
                     'native_Matchete_topology_source_plus_Hc_action_residual':err,
                     'native_RGBeta_complex_symbols':['lambdaT3']+(list(SPECIAL) if tag in BC else []),
                     'baseline_complex_quartic_slots':7 if tag in BC else 1,
                     'physical_dictionary_certified':False}
    require(sum(x['baseline_complex_quartic_slots'] for x in models.values())==28,
            '28 complex EFT1 slots not accounted for')
    return {'format':FORMAT,'status':STATUS,'models_checked':16,
        'native_Matchete_topology_CGs':16,'native_Matchete_BC_extra_CGs_rechecked':14,
        'BC_baseline_special_complex_quartic_slots':12,
        'native_RGBeta_complex_direct_plus_conjugate_channels':28,
        'baseline_EFT1_complex_slots':28,
        'max_native_topology_6_over_24_real_action_relative_residual':worst,
        'BC_baseline_six_of_seven_conditional_source_maps':bc_map,
        'models':models,
        'input_sha256':{str(p.relative_to(root)) if p.is_relative_to(root) else p.name:sha(p)
                        for p in inputs},
        'certification':{
           'native_Matchete_topology_CG_sources_verified':True,
           'current_BC_Matchete_extra_CGs_identical_to_prior_source_evidence':True,
           'BC_six_baseline_special_subspace_rank_verified':True,
           'RGBeta_native_complex_direct_conjugate_projectors_observed':True,
           'native_RGBeta_complex_component_CG_frame_certified':False,
           'physical_Matchete_RGBeta_complex_coupling_map_certified':False,
           'BC_mixed_scalar_mass_basis_certified':False,
           'signed_GammaCCWeyl_certified':False,
           'physical_Wilson_evolution_authorized':False,
           'EFT_order':'UV -> integrate out ALL fermions -> integrate out ALL scalars -> SMEFT'}}
