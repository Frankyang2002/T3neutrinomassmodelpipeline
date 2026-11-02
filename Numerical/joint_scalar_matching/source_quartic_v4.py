"""Source-backed HHSS quartic export and conditional numerical endpoint assembly.

All source Matchete scalar-quartic terms are processed with the repository's
existing exact CG/real-field parser. Missing couplings, mass matrices, source
basis maps or EFT1 scales are NOT inferred from the UV frozen benchmark.

This module never assigns a physical value to the unknown GammaCC/Weyl factor.
"""
from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
from pathlib import Path
from typing import Any, Callable
import json
import math
import re

import numpy as np
import sympy as sp

from .adapter import CONVENTION, scalar_groups
from .fermion_boundary_bridge import SourceBoundary, require, read_json, digest, slug
from .core import EigenBasis

FORMAT = 't3-source-hhss-quartic-channel-export-v4'
ENDPOINT_FORMAT = 't3-source-quartic-scalar-endpoint-v4'
BACKGROUND_FORMAT = 't3-joint-scalar-conditional-background-v1'
POTENTIAL_CONVENTION = 'V=lambda_abcd*phi_a*phi_b*phi_c*phi_d/4!'


def _model_dimensions(model: str) -> tuple[int, int]:
    dims={'A':(1,3),'B':(2,2),'C':(2,2),'D':(3,1),'E':(3,3)}
    return dims[model[0]]


def source_seed_path(root: Path, model: str) -> Path:
    """The 16 original-field Matchete quartic seeds, not split/mixed EFT1 seeds."""
    t=slug(model).split('_')[1]
    return Path(root)/f'T3_{model[0]}_alpha_{t}'/'data'/'eft1_after_F_scalar_quartic_seed.json'


def _coupling_name(inputform: str) -> str:
    # The source backend also enforces precisely one coupling per term.
    names=re.findall(r'Coupling\[([^,\]]+),',inputform)
    require(len(names)==1 and names[0].strip()!='',
            'every HHSS Matchete source term must contain precisely one coupling')
    return names[0].strip()


@dataclass(frozen=True)
class QuarticChannels:
    model: str
    seed_path: Path
    seed_sha256: str
    source_basis_labels: tuple[str, ...]
    by_name: dict[str, dict[tuple[int,int,int,int], sp.Expr]]
    seed_hhss_term_count: int
    seed_total_term_count: int
    parser: str

    @property
    def n(self) -> int:
        return len(self.source_basis_labels)


def build_source_channels(seed_path: Path, model: str,
                          *, builder: Callable[...,Any] | None = None) -> QuarticChannels:
    """Calculate each exact quartic channel in the ORIGINAL scalar basis.

    `builder` is injection solely for tests. Production uses
    RGE.running.intermediate.ScalarOnlyTensorAdapters.build_eft1_quartic_tensor.
    Expressions retain ALL source CGs, signs and factorials, unevaluated.
    """
    seed_path=Path(seed_path)
    seed=read_json(seed_path)
    terms=seed.get('ScalarQuarticTerms')
    require(isinstance(terms,list) and isinstance(seed.get('CGRegistry'),list),
            f'no Matchete scalar quartic terms/CG registry: {seed_path}')
    require(seed.get('ScalarQuarticTermCount')==len(terms) and
            seed.get('IntegratedFields')==['F'] and
            seed.get('ActiveHeavyFields')==['S1','S2'],
            'not an all-F-first scalar-only EFT1 quartic source')
    if builder is None:
        from RGE.running.intermediate.ScalarOnlyTensorAdapters import (
            build_eft1_quartic_tensor as builder)
    d1,d2=_model_dimensions(model)
    n=2*(d1+d2)
    groups: dict[str,list[dict]]={}
    hhss=0
    for item in terms:
        require(isinstance(item,dict) and isinstance(item.get('TermInputForm'),str),
                'malformed native Matchete quartic term')
        h=item.get('HiggsFieldCount')
        s1=item.get('Scalar1FieldCount')
        s2=item.get('Scalar2FieldCount')
        require(all(type(i)==int and i>=0 for i in (h,s1,s2)) and h+s1+s2==4,
                'quartic field multiplicities invalid')
        if h != 2 or s1+s2 != 2:
            continue
        hhss+=1
        name=_coupling_name(item['TermInputForm'])
        groups.setdefault(name,[]).append(item)
    require(hhss>0, 'there are no HHSS quartic source terms; do not assume zero')
    result: dict[str,dict[tuple[int,int,int,int],sp.Expr]]={}
    for name,source_terms in sorted(groups.items()):
        only={'CGRegistry':seed['CGRegistry'],'ScalarQuarticTerms':source_terms}
        tensor=builder(only,d_s1=d1,d_s2=d2,shared_scalar=False)
        entries={}
        for key,expr in tensor.nonzero_items():
            key=tuple(int(i) for i in key)
            require(len(key)==4 and all(1<=i<=n+4 for i in key),
                    f'{name}: native quartic component outside original-field layout')
            hlegs=[i for i in key if i<=4]
            slegs=[i for i in key if i>=5]
            if len(hlegs)!=2 or len(slegs)!=2:
                continue
            local=(hlegs[0]-1,hlegs[1]-1,slegs[0]-5,slegs[1]-5)
            require(local not in entries, 'duplicate source tensor index')
            entries[local]=sp.simplify(expr)
        require(entries, f'{name}: source parser produced no HHSS quartic coordinates')
        result[name]=entries
    labels=('NewScalar1',)*(2*d1)+('NewScalar2',)*(2*d2)
    return QuarticChannels(model,seed_path,digest(seed_path),labels,result,
                           hhss,len(terms),
                           'RGE.running.intermediate.ScalarOnlyTensorAdapters.build_eft1_quartic_tensor')


def symbolic_export(ch: QuarticChannels, boundary: SourceBoundary) -> dict:
    require(boundary.model==ch.model and len(boundary.labels)==ch.n,
            'quartic release differs in model or scalar dimension')
    original_basis_matches=tuple(boundary.labels)==ch.source_basis_labels
    return {'format':FORMAT,'model':ch.model,'seed_quartic_sha256':ch.seed_sha256,
            'source_boundary_sha256':boundary.hash,
            'source_quartic_seed':str(ch.seed_path),
            'quartic_tensor_basis':'original complex (S1,S2) decomposed into canonical real fields',
            'source_real_scalar_labels':list(ch.source_basis_labels),
            'fermion_boundary_real_scalar_labels':list(boundary.labels),
            'source_basis_matches_fermion_boundary':original_basis_matches,
            'source_to_fermion_basis_rotation_required':not original_basis_matches,
            'potential_tensor_normalization':POTENTIAL_CONVENTION,
            'native_parser':ch.parser,
            'hhss_native_terms':ch.seed_hhss_term_count,
            'all_quartic_native_terms':ch.seed_total_term_count,
            'coupling_count':len(ch.by_name),
            'couplings_required_at_muS':sorted(ch.by_name),
            'channels':{name:[{'component_zero_based_HA_HB_Sa_Sb':list(key),
                               'exact_expr':sp.sstr(expr)}
                              for key,expr in sorted(entries.items())]
                        for name,entries in sorted(ch.by_name.items())},
            'native_scalar_quadratic_mass_at_muS_supplied':False,
            'quartic_couplings_at_muS_supplied':False,
            'normalization_GammaCC_to_Weyl_certified':False,
            'complete_general_scalar_potential_certified':False,
            'physical_C5_certified':False}


def _complex_number(v: Any, name: str) -> complex:
    if isinstance(v,bool): raise ValueError(f'{name}: boolean is not a coupling')
    if isinstance(v,(float,int)):
        z=complex(float(v))
    elif isinstance(v,dict) and set(v)=={'real','imag'}:
        z=complex(float(v['real']),float(v['imag']))
    else:
        raise ValueError(f'{name}: expected number or real/imag object')
    require(math.isfinite(z.real) and math.isfinite(z.imag),f'{name}: nonfinite')
    return z


def evaluate_channels(ch: QuarticChannels, couplings: dict[str,Any]) -> tuple[np.ndarray, tuple[str,...]]:
    """Source CG tensors evaluated at supplied SAME-scale numerical couplings.

    No missing coupling is defaulted to zero. A zero coupling must be explicit.
    """
    require(isinstance(couplings,dict), 'muS quartic coupling dictionary missing')
    names=tuple(sorted(ch.by_name))
    missing=set(names)-set(couplings)
    require(not missing, f'missing muS quartic values for {sorted(missing)}')
    vals={name:_complex_number(couplings[name],name) for name in names}
    # One real scalar potential must produce a real fully symmetric tensor.
    numeric=np.zeros((len(names),4,4,ch.n,ch.n),dtype=float)
    for index,name in enumerate(names):
        substitution={sp.Symbol(k):v for k,v in vals.items()}
        for (a,b,s,t),expr in ch.by_name[name].items():
            value=complex(sp.N(expr.subs(substitution),15))
            require(math.isfinite(value.real) and math.isfinite(value.imag),
                    f'{name}: nonfinite scalar component')
            require(abs(value.imag)<=1e-9*max(abs(value.real),abs(value.imag),1e-15),
                    f'{name}: potential component is complex; coupling reality/HC convention unresolved')
            for A,B in {(a,b),(b,a)}:
                for S,T in {(s,t),(t,s)}:
                    numeric[index,A,B,S,T]=value.real
    return numeric,names


def _validate_orthogonal(U: np.ndarray,n:int) -> None:
    require(np.isrealobj(U) and U.shape==(n,n) and np.isfinite(U).all() and
            np.allclose(U.T@U,np.eye(n),rtol=1e-10,atol=1e-10),
            'source original->boundary scalar map must be an orthogonal real n×n matrix')


def _mass2(value: Any,n:int) -> np.ndarray:
    a=np.asarray(value)
    require(np.isrealobj(a) and a.shape==(n,n) and np.isfinite(a).all() and
            np.allclose(a,a.T,rtol=1e-10,atol=1e-8),
            'muS source-real-basis scalar mass matrix must be finite, real, symmetric')
    a=np.asarray(a,float)
    EigenBasis.from_mass2(a)
    return a


def build_background(ch: QuarticChannels,boundary: SourceBoundary,
                     endpoint_path: Path,output_dir: Path) -> dict:
    """Construct v2/v3 compatible conditional background; no physical C5 claimed.

    Endpoint mass matrix and couplings must ALL be at same muS. For special
    scalar bases, the endpoint supplies and independently sources an explicit
    map: original_real_fields = R @ fermion_boundary_real_fields.
    """
    e=read_json(endpoint_path)
    require(e.get('format')==ENDPOINT_FORMAT and e.get('model')==ch.model,
            'muS endpoint format/model incorrect')
    require(e.get('seed_quartic_sha256')==ch.seed_sha256 and
            e.get('source_boundary_sha256')==boundary.hash,
            'muS scalar endpoint is not hashed to Matchete source AND fermion boundary')
    require(e.get('renormalization_scheme')=='MSbar' and
            e.get('scalar_mass_origin') in ('unbroken_electroweak','EWSB_resummed') and
            e.get('both_scalar_multiplets_active_until_single_muS') is True and
            e.get('masses_and_quartics_at_same_muS') is True,
            'endpoint must declare common muS, scheme, unbroken/broken mass origin')
    mf=float(e.get('mu_f_gev',float('nan')))
    ms=float(e.get('mu_s_gev',float('nan')))
    ref=float(e.get('reference_scale_gev',ms))
    require(all(math.isfinite(x) for x in (mf,ms,ref)) and
            0<ms<=mf and ref>0 and abs(mf/boundary.mu_f_gev-1)<1e-10,
            'muS and muF scales invalid')
    provenance=e.get('provenance')
    require(isinstance(provenance,dict) and
            isinstance(provenance.get('scalar_coupling_trajectory'),str) and
            provenance['scalar_coupling_trajectory'].strip() and
            isinstance(provenance.get('scalar_mass_trajectory'),str) and
            provenance['scalar_mass_trajectory'].strip(),
            'endpoint needs actual numerical scalar mass AND coupling trajectory provenance')
    require(e.get('hhss_source_operator_set_reviewed') is True and
            isinstance(provenance.get('hhss_source_review'),str) and
            provenance['hhss_source_review'].strip(),
            'complete HHSS source operator inventory must be explicitly reviewed')
    n=ch.n
    mass_original=_mass2(e.get('mass_squared_gev2_original_real_basis_muS'),n)
    numeric,names=evaluate_channels(ch,e.get('quartic_couplings_muS'))
    if tuple(boundary.labels)==ch.source_basis_labels:
        require('original_real_equals_boundary_real_R' not in e,
                'ordinary same-basis model must not introduce an implicit field rotation')
        R=np.eye(n)
        transform='source and fermion boundary have identical canonical real component labels'
    else:
        provenance_map=provenance.get('original_to_boundary_basis_map')
        require(isinstance(provenance_map,str) and provenance_map.strip(),
                'special scalar case requires independently sourced original-to-boundary map')
        require(e.get('special_scalar_basis_transform_externally_validated') is True,
                'split/mixed scalar source basis transformation not validated')
        R=np.asarray(e.get('original_real_equals_boundary_real_R'))
        _validate_orthogonal(R,n)
        transform=provenance_map
    # original field vector = R @ boundary field vector
    mass_boundary=R.T@mass_original@R
    channels_boundary=np.einsum('cABij,ia,jb->cABab',numeric,R,R,optimize=True)
    summed=channels_boundary.sum(axis=0)
    require(np.isfinite(summed).all(),'nonfinite generated quartic tensor')
    EigenBasis.from_mass2(mass_boundary)
    output_dir=Path(output_dir)
    output_dir.mkdir(parents=True,exist_ok=True)
    name=slug(ch.model)
    target=output_dir/f'{name}.json'
    npz=output_dir/f'{name}.npz'
    require(not target.exists() and not npz.exists(),
            f'refusing to overwrite existing numerical background for {ch.model}')
    # All arrays now in boundary real component coordinates, so existing
    # conditional matching takes U = identity (no spurious second rotation).
    np.savez_compressed(npz,source_to_matching_real_U=np.eye(n),
                        mass_squared_gev2_muS=mass_boundary,
                        lambda_HHSS_muS=summed,
                        lambda_HHSS_channels_muS=channels_boundary)
    obj={'format':BACKGROUND_FORMAT,'model':ch.model,
         'source_boundary_sha256':boundary.hash,
         'source_scalar_real_labels':list(boundary.labels),
         'renormalization_scheme':'MSbar','tensor_convention':CONVENTION,
         'mu_f_gev':mf,'mu_s_gev':ms,'reference_scale_gev':ref,
         'scalar_mass_origin':e['scalar_mass_origin'],
         'same_muS_quartic_and_mass':True,'no_separate_scalar_threshold':True,
         'tensor_npz':npz.name,'quartic_channel_names':list(names),
         'provenance':{'source':f'Matchete seed {ch.seed_path}; muS endpoint {endpoint_path}',
                       'quartic_potential_normalization':POTENTIAL_CONVENTION,
                       'quartic_source_sha256':ch.seed_sha256,
                       'endpoint_sha256':digest(endpoint_path),
                       'original_to_boundary_basis_map':transform,
                       'scalar_coupling_trajectory':provenance['scalar_coupling_trajectory'],
                       'scalar_mass_trajectory':provenance['scalar_mass_trajectory'],
                       'hhss_source_review':provenance['hhss_source_review']},
         'physical_C5_certified':False,
         'GammaCC_to_Weyl_normalization':None,
         'full_general_one_loop_operator_completeness_certified':False,
         'source_quartic_is_original_basis':True,
         'mass_quartic_rotated_into_fermion_boundary_basis':True}
    target.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    return {'background_json':str(target),'background_npz':str(npz),
            'mu_s_gev':ms,'channels':list(names),
            'max_quartic_abs':float(np.max(abs(summed))),
            'basis_map':transform,'physical_C5_certified':False}
