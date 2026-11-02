"""Source-backed, explicitly UNVERIFIED 1-TeV scalar benchmark handoff.

Uses existing all16_v1 renormalisable EFT1 numerical snapshots, and the
Matchete-source HHSS quartic tensors. Does not assert that RGBeta coupling
coordinates equal native Matchete invariant coordinates: alias assumptions
below are deliberately exposed in every output, with no physical certification.

UV -> all heavy fermions at muF -> EFT1(S1,S2) -> BOTH scalars at ONE muS
-> SMEFT. The result is a conditional LLSS sector per unknown N_GammaCC.
"""
from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from .fermion_boundary_bridge import SourceBoundary, require, slug
from .source_quartic_v4 import QuarticChannels, evaluate_channels, POTENTIAL_CONVENTION
from .large_batch import load_background, evaluate_conditional

BENCHMARK_FORMAT = 't3-joint-scalar-source-backed-exploratory-benchmark-v5'
RESULT_FORMAT = 't3-joint-scalar-exploratory-c5-per-unknown-gammacc-v5'
SOURCE_SNAPSHOT = 'output/general_observable/all16_v1/{model}/running.json'
REGULAR = frozenset(('A_-4','B_-3','B_+1','C_-3','C_+1','D_+2','E_-4','E_+2'))
SPECIAL = frozenset(('A_-2','A_+0','B_-1','C_-1','D_-2','D_+0','E_-2','E_+0'))

# These aliases encode a *candidate*, NOT a certified Matchete <-> RGBeta map.
# Source native 'Inv1/Inv2' numbering and normalization are not independently
# checked against RGBeta here. This is why all outputs are exploratory.
ALIASES = {
    'lambdaT3': 'lambdaT3',
    'lambdaH1Inv1': 'lambdaH1',
    'lambdaH1Inv2': 'lambdaH1Adj',
    'lambdaH2Inv1': 'lambdaH2',
    'lambdaH2Inv2': 'lambdaH2Adj',
}

@dataclass(frozen=True)
class RunningSnapshot:
    model: str
    source: Path
    sha256: str
    mu_f: float
    mu_s: float
    scalar_mass1_sq: float
    scalar_mass2_sq: float
    couplings: dict[str, Any]
    evidence: dict[str, Any]


def _real_positive(x: Any, label: str) -> float:
    require(type(x) in (int,float), f'{label} must be a real number')
    v=float(x)
    require(math.isfinite(v) and v>0, f'{label} must be finite and positive')
    return v


def load_running_snapshot(path: Path, model: str) -> RunningSnapshot:
    """Require the actual EFT1@S snapshot, not UV / F values or interpolates."""
    path=Path(path)
    doc=json.loads(path.read_text(encoding='utf-8-sig'))
    require(isinstance(doc,dict) and doc.get('model')==model,
            'benchmark running model mismatch')
    require(doc.get('status')=='NumericalRenormalisableRunningOnly',
            'unexpected running pipeline status')
    require(doc.get('eft1_wilson_running_performed') is False,
            'do not reuse an altered Wilson-running provenance')
    require(doc.get('illustrative_benchmark') is True and
            doc.get('physical_mass_parameters_available') is False,
            'this v5 adapter only handles non-certified illustrative benchmark snapshots')
    scales=doc.get('scales_gev')
    require(isinstance(scales,dict) and 'F' in scales and 'S' in scales,
            'missing F/S EFT scales')
    mf=_real_positive(scales['F'],'muF')
    ms=_real_positive(scales['S'],'muS')
    require(ms<mf and math.isclose(mf,1e5,rel_tol=1e-10),
            'incompatible EFT1 threshold scales')
    state=doc.get('eft1_at_scalar_threshold')
    require(isinstance(state,dict),'actual EFT1 scalar-threshold state missing')
    m1=_real_positive(state.get('mS1Sq'),'mS1Sq at muS')
    m2=_real_positive(state.get('mS2Sq'),'mS2Sq at muS')
    require(isinstance(doc.get('symbolic_provenance'),dict),
            'symbolic/source provenance missing')
    return RunningSnapshot(model,path,sha256(path.read_bytes()).hexdigest(),mf,ms,
                           m1,m2,state,doc['symbolic_provenance'])


def select_candidate_couplings(ch: QuarticChannels, snapshot: RunningSnapshot
                               ) -> tuple[dict[str,Any],dict[str,str]]:
    """Never silently fill missing or unmapped source couplings with zero."""
    require(ch.model==snapshot.model,'quartic and running model mismatch')
    values={}
    mapping={}
    for native in sorted(ch.by_name):
        require(native in ALIASES,
                f'UNMAPPED Matchete quartic {native}: no numerical value inferred')
        rhs=ALIASES[native]
        require(rhs in snapshot.couplings,
                f'{native}: missing candidate RGBeta parameter {rhs} at muS')
        raw=snapshot.couplings[rhs]
        require(type(raw) in (float,int) or (isinstance(raw,dict) and
                set(raw)=={'real','imag'}),f'unsupported muS coupling representation: {rhs}')
        values[native]=raw
        mapping[native]=rhs
    return values,mapping


def inspect_model(ch: QuarticChannels,boundary: SourceBoundary,
                  snapshot: RunningSnapshot) -> dict[str,Any]:
    """Report diagnostics for all 16 even where numerical matching is blocked."""
    require(ch.model==boundary.model==snapshot.model,'model mismatch')
    mf_ratio=abs(snapshot.mu_f/boundary.mu_f_gev-1)
    require(mf_ratio<1e-10,'F matching scale differs from source boundary')
    identical=ch.source_basis_labels==boundary.labels
    special=ch.model in SPECIAL
    require(ch.model in SPECIAL|REGULAR,'unknown T3 model')
    require(identical==(not special),
            'special/ordinary basis classification differs from release')
    candidate={}
    try:
        _,candidate=select_candidate_couplings(ch,snapshot)
        mapping_status='candidate_all_native_HHSS_channels_mapped'
    except ValueError as exc:
        mapping_status=f'candidate_coupling_mapping_blocked: {exc}'
    blockers=[
        'absolute Matchete GammaCC to Weyl factor/sign remains unknown',
        'Matchete invariant coupling normalization vs RGBeta parameter map is NOT certified',
        'full scalar UV operator completeness and full hard one-loop matching NOT certified',
    ]
    if special:
        blockers.append('split-real or B/C mixed heavy-scalar mass and original-to-boundary map at muS missing')
    return {'model':ch.model,'mu_f_gev':snapshot.mu_f,'mu_s_gev':snapshot.mu_s,
            'mass_source':'renormalisable EFT1@S illustrative endpoint',
            'source_scalar_basis_matches_fermion_boundary':identical,
            'candidate_coupling_map':candidate,'coupling_mapping_status':mapping_status,
            'scalar_mass_parameters_gev2':{'mS1Sq':snapshot.scalar_mass1_sq,
                                           'mS2Sq':snapshot.scalar_mass2_sq},
            'running_snapshot_sha256':snapshot.sha256,
            'source_quartic_seed_sha256':ch.seed_sha256,
            'fermion_boundary_sha256':boundary.hash,
            'special_mixing_or_split_required':special,
            'physical_predictions_certified':False,'remaining_blockers':blockers}


def build_ordinary_background(ch: QuarticChannels,boundary: SourceBoundary,
                              snapshot: RunningSnapshot,output_dir: Path
                              ) -> dict[str,Any]:
    """Build conditional V3-compatible tensor background for eight regular bases.

    The diagonal scalar mass matrix is the unbroken EFT1 *benchmark* quadratic
    mass.  No EWSB corrections, special mixing, or full-mass claims are made.
    """
    require(ch.model in REGULAR and ch.source_basis_labels==boundary.labels,
            'ordinary ORIGINAL=BOUNDARY basis required; special models are blocked')
    require(abs(snapshot.mu_f/boundary.mu_f_gev-1)<1e-10,'muF mismatch')
    values,mapping=select_candidate_couplings(ch,snapshot)
    channel_tensors,names=evaluate_channels(ch,values)
    n=ch.n
    d1=n//2-{'A':3,'B':2,'C':2,'D':1,'E':3}[ch.model[0]]
    # The dimensions above are fixed to the same scalar blocks as v4.
    require(d1=={'A':1,'B':2,'C':2,'D':3,'E':3}[ch.model[0]],
            'scalar layout incompatible with v4')
    mass=np.diag([snapshot.scalar_mass1_sq]*(2*d1)+
                 [snapshot.scalar_mass2_sq]*(n-2*d1))
    require(np.isfinite(channel_tensors).all() and channel_tensors.shape[1:]==(4,4,n,n),
            'incomplete source-derived numerical quartic array')
    full=channel_tensors.sum(axis=0)
    target=Path(output_dir)
    target.mkdir(parents=True,exist_ok=True)
    basename=slug(ch.model)
    meta_path=target/f'{basename}.json'
    npz_path=target/f'{basename}.npz'
    require(not meta_path.exists() and not npz_path.exists(),
            'refusing to overwrite existing v5 background')
    np.savez_compressed(npz_path,
                        source_to_matching_real_U=np.eye(n),
                        mass_squared_gev2_muS=mass,
                        lambda_HHSS_muS=full,
                        lambda_HHSS_channels_muS=channel_tensors)
    meta={
        'format':'t3-joint-scalar-conditional-background-v1',
        'model':ch.model,'source_boundary_sha256':boundary.hash,
        'source_scalar_real_labels':list(boundary.labels),
        'renormalization_scheme':'MSbar',
        'tensor_convention':'T3_REAL_LLSS_HHSS_WEINBERG_v1',
        'mu_f_gev':snapshot.mu_f,'mu_s_gev':snapshot.mu_s,
        'reference_scale_gev':snapshot.mu_s,
        'scalar_mass_origin':'unbroken_electroweak',
        'same_muS_quartic_and_mass':True,'no_separate_scalar_threshold':True,
        'tensor_npz':npz_path.name,'quartic_channel_names':list(names),
        'provenance':{
            'source':f'Illustrative renormalisable EFT1 running {snapshot.source}; '
                     f'Matchete-source HHSS seed {ch.seed_path}; '
                     'candidate Matchete-to-RGBeta parameter aliases (UNVERIFIED)',
            'quartic_potential_normalization':POTENTIAL_CONVENTION,
            'running_snapshot_sha256':snapshot.sha256,
            'source_quartic_seed_sha256':ch.seed_sha256,
            'source_boundary_sha256':boundary.hash,
            'candidate_native_to_rgbeta_aliases':mapping,
            'coupling_alias_normalization_certified':False,
        },
        'benchmark_not_physical':True,
        'physical_C5_certified':False,
        'absolute_GammaCC_to_Weyl_factor':None,
        'operator_inventory_complete':False,
    }
    meta_path.write_text(json.dumps(meta,indent=2,ensure_ascii=False,allow_nan=False)+'\n',
                         encoding='utf-8')
    # Strictly re-open the exact v3 background format to catch convention errors.
    background=load_background(meta_path,boundary)
    return {'background_json':str(meta_path),'background_npz':str(npz_path),
            'mu_s_gev':background.mu_s,'mass_eigenvalues_gev2':
             [snapshot.scalar_mass1_sq,snapshot.scalar_mass2_sq],
            'source_HHSS_channels':list(names),
            'coupling_aliases_NOT_CERTIFIED':mapping,'physical_C5_certified':False}


def evaluate_ordinary_response(boundary:SourceBoundary,background_metadata:Path
                               ) -> dict[str,Any]:
    bg=load_background(background_metadata,boundary)
    response=evaluate_conditional(boundary,bg)
    return {'format':RESULT_FORMAT,'model':boundary.model,
            'interpretation':'Illustrative LLSS fixed-order leading-log + simultaneous finite scalar bubble, per unit unknown GammaCC/Weyl factor; NOT physical C5',
            'response':response,'physical_C5_certified':False,
            'source_to_master_spinor_factor':None,
            'matchete_rgbeta_quartic_map_certified':False,
            'hard_one_loop_fermion_matching_included':False,
            'complete_dim5_sector_certified':False}


def read_v4_cached_channels(path:Path, seed_path:Path,
                            boundary:SourceBoundary) -> QuarticChannels:
    """Reuse existing v4 source-derived symbolic tensors after exact hash checks."""
    import sympy as sp
    from .source_quartic_v4 import FORMAT, _model_dimensions
    path=Path(path)
    doc=json.loads(path.read_text(encoding='utf-8-sig'))
    seed_path=Path(seed_path)
    sha=sha256(seed_path.read_bytes()).hexdigest()
    require(doc.get('format')==FORMAT and doc.get('model')==boundary.model,
            'cached v4 quartic format/model mismatch')
    require(doc.get('source_boundary_sha256')==boundary.hash and
            doc.get('seed_quartic_sha256')==sha,
            'cached v4 tensor hashes do not match the CURRENT source files')
    d1,d2=_model_dimensions(boundary.model)
    labels=('NewScalar1',)*(2*d1)+('NewScalar2',)*(2*d2)
    require(doc.get('source_real_scalar_labels')==list(labels) and
            doc.get('fermion_boundary_real_scalar_labels')==list(boundary.labels),
            'cached quartic real basis differs from F-boundary')
    require(doc.get('potential_tensor_normalization')==POTENTIAL_CONVENTION,
            'wrong potential factorial convention')
    raw=doc.get('channels')
    require(isinstance(raw,dict) and set(raw)==set(doc.get('couplings_required_at_muS',[])),
            'malformed cached v4 quartic inventory')
    out={}
    n=2*(d1+d2)
    for name,entries in raw.items():
        require(isinstance(name,str) and name and isinstance(entries,list) and entries,
                'empty source-derived HHSS coupling channel')
        found={}
        for item in entries:
            key=item.get('component_zero_based_HA_HB_Sa_Sb')
            expression=item.get('exact_expr')
            require(isinstance(key,list) and len(key)==4 and
                    all(type(k)==int for k in key) and
                    all(0<=k<bound for k,bound in zip(key,(4,4,n,n))) and
                    isinstance(expression,str),
                    f'bad cached HHSS coordinate for {name}')
            key=tuple(key)
            require(key not in found,'duplicate cached quartic component')
            # Expressions are generated by the repository's own exact v4 SymPy
            # emitter, hash-tied to the unchanged Matchete seed. Never accept
            # arbitrary Python execution in symbol names.
            import re
            require(re.fullmatch(r'[A-Za-z_0-9*+/().,\- ]+',expression) is not None,
                    f'invalid cached SymPy expression for {name}')
            allowed={name,'sqrt','conjugate','I','re','im'}
            tokens=set(re.findall(r'[A-Za-z_][A-Za-z_0-9]*',expression))
            require(tokens<=allowed,
                    f'cached source expression uses unrecognised symbols: {sorted(tokens-allowed)}')
            expr=sp.sympify(expression,locals={'sqrt':sp.sqrt,'conjugate':sp.conjugate,
                                                'I':sp.I,'re':sp.re,'im':sp.im})
            require(expr.free_symbols <= {sp.Symbol(name)},
                    f'cached HHSS channel {name} has cross/unknown coupling symbols')
            found[key]=expr
        out[name]=found
    return QuarticChannels(boundary.model,seed_path,sha,labels,out,
                           int(doc['hhss_native_terms']),int(doc['all_quartic_native_terms']),
                           'hash-verified v4 Matchete CG quartic cache')
