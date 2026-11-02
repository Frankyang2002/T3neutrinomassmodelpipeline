"""V6: fail-closed provenance and conditional *live Matchete* extra-quartic source audit.

Do not mistake UV builder `PlusHc` operators for an EFT1 matching seed, or an
RGBeta `AddQuartic` declaration for a physical running coefficient. The latter
requires a separate threshold, beta-function and phase/normalisation proof.
"""
from __future__ import annotations

from hashlib import sha256
from pathlib import Path
import json
import re
import subprocess
from typing import Any

import sympy as sp

from .source_cg_v5 import (_exact_number, parse_all_cgs, inspect_all_terms,
                           ExactCG, MatcheteCGError)
from .matchete_rgbeta_dictionary_v5 import candidate_cg_reference_maps

FORMAT = 't3-ws2-live-extra-quartic-source-audit-v6'
OUTPUT_OWNER = Path('output/general_observable/lo_certification_couplings')
MODEL_FAMILIES = {'A':(1,3,2),'B':(2,2,1),'C':(2,2,3),
                  'D':(3,1,2),'E':(3,3,2)}
ALL_MODELS = ('A_-4','A_-2','A_+0','B_-3','B_-1','B_+1','C_-3','C_-1','C_+1',
              'D_-2','D_+0','D_+2','E_-4','E_-2','E_+0','E_+2')
CANONICAL_TYPES = ('H','Bar @ H','S1','Bar @ S1','S2','Bar @ S2')


class V6SourceError(ValueError):
    pass


def file_sha256(path:Path)->str:
    return sha256(path.read_bytes()).hexdigest()


def slug(model:str)->str:
    if model not in ALL_MODELS:raise V6SourceError('Unrecognised ordinary T3 model')
    return model.replace('_-','_m').replace('_+','_p')


def alpha_tag(model:str)->str:
    return slug(model).split('_',1)[1]


def canonical_counts(counts:tuple[int,...])->tuple[int,...]:
    if len(counts)!=6 or any(type(x) is not int or x<0 for x in counts) or sum(counts)!=4:
        raise V6SourceError('Invalid quartic field multiset')
    conjugate=tuple(counts[i] for i in (1,0,3,2,5,4))
    return min(counts,conjugate)


def extra_name(counts:tuple[int,...],index:int,*,cg:bool=False)->str:
    c=canonical_counts(counts)
    if not 1<=index<=10:raise V6SourceError('Invalid invariant index')
    prefix='T3ExtraCGQ' if cg else 'lambdaT3ExtraQ'
    return prefix+''.join(map(str,c))+'I'+str(index)


def model_registry_quartics(model:str,registry:dict)->list[dict]:
    """The 50 generated RGBeta entries are an *inventory*, not a physical map."""
    row=registry['models'][model]
    entries=[]
    by_sector={}
    for group,field in (('new_rgbeta_symbols','name'),('preexisting_rgbeta_symbols','existing_name')):
        for entry in row[group]:
            counts=tuple(entry['counts']); cc=canonical_counts(counts)
            by_sector.setdefault(cc,[]).append((group,entry[field],counts))
    for cc,rows in sorted(by_sector.items()):
        # The Matchete generator enumerates independent invariants per
        # monomial. A registry count of 2 in one sector has 2 native indices,
        # but their orientations/indices must not be identified by name alone.
        for i,(group,name,counts) in enumerate(rows,1):
            entries.append({'rgbeta_symbol':name, 'registry_category':group,
                'RGBeta_field_counts':list(counts),'canonical_field_counts':list(cc),
                'Matchete_generator_candidate_symbol':extra_name(cc,i),
                'Matchete_generator_candidate_CG':extra_name(cc,i,cg=True),
                'Matchete_candidate_index_unverified':True,
                'sector_HHSS':sum(counts[:2])==2,
                'has_physical_conversion':False})
    if len(entries)!=row['total_additional_quartics']:
        raise V6SourceError(f'{model}: registry total differs from candidate count')
    return entries


def generator_contract(root:Path)->dict:
    """Inspect the literal *local* Matchete operator generator, not a guess."""
    p=root/'Lagrangian/interactions/GeneralScalarOperators.wl'
    q=root/'Lagrangian/interactions/GeneralScalarVerifiedQuartics.wl'
    b=root/'Lagrangian/model/LagrangianBuilder.wl'
    t=root/'Lagrangian/interactions/T3Topology.wl'
    for path in (p,q,b,t):
        if not path.is_file():raise V6SourceError(f'Missing local Matchete generator: {path}')
    src=p.read_text(encoding='utf8');builder=b.read_text(encoding='utf8')
    topology=t.read_text(encoding='utf8')
    requirements={
        'defines_extra_complex_couplings':r'DefineCoupling\[coupling,\s*SelfConjugate\s*->\s*False\]',
        'includes_hermitian_conjugate':r'contribution\s*=\s*coupling\[\]\s*product\s*If\[',
        'plus_HC':r'"Expression"\s*->\s*PlusHc\[contribution\]',
        'full_bosonic_symmetrisation':r'T3ExtendedSymmetrizeGroup\[tensor_,\s*positions_List\]',
        'uses_verified_structure_counts':r'T3ExtendedVerifiedStructureCount\[model\]',
        'names_native_extra_CG':r'T3ExtraCG',
        'names_native_extra_couplings':r'lambdaT3Extra',
    }
    missing=[key for key,rx in requirements.items() if not re.search(rx,src)]
    if ('T3ExtendedBuildExtraQuarticCandidates[model]' not in topology or
            'BuildT3InteractionCandidates[model]' not in builder):
        missing.append('generator_used_by_actual_builder')
    if missing:raise V6SourceError('Unreviewed Matchete generator change: '+', '.join(missing))
    return {'type':'experimental_UV_Lagrangian_builder_not_an_EFT1_seed',
        'Matchete_plus_HC_source_pattern_verified':True,
        'complex_coupling_not_a_real_assumption':True,
        'extra_quartics_opt_in':'T3_EXPERIMENTAL_GENERAL_QUARTICS=1',
        'live_native_CG_requires_Matchete_GetCGTensor':True,
        'source_sha256':{str(path.relative_to(root)):file_sha256(path) for path in (p,q,b,t)}}


def _model_tags(model:str)->set[str]:
    s=slug(model)
    return {s,'T3_'+model[0]+'_alpha_'+alpha_tag(model)}


def seed_belongs_to_model(path:Path,model:str)->bool:
    """Accept paths with an explicit model directory, not a substring guess."""
    return bool(_model_tags(model).intersection(path.parts))


def find_seed_files(root:Path,model:str,paths:list[Path])->list[Path]:
    return sorted(p for p in paths if seed_belongs_to_model(p,model))


def indexed_seed_evidence(path:Path,model:str,expected:list[dict])->dict:
    """Read actual EFT1 seeds; exact matching needs matching source terms.

    Name and field sector are checked *independently*. Matching a field sector
    cannot certify a coupling normalisation or CG phase.
    """
    source=json.loads(path.read_text(encoding='utf8'))
    if source.get('IntegratedFields')!=['F'] or source.get('ActiveHeavyFields')!=['S1','S2']:
        raise V6SourceError('Seed does not describe the required EFT1 threshold')
    cgs=parse_all_cgs(source)
    terms=inspect_all_terms(source,cgs)
    by_canonical={}
    for t in terms:
        key=canonical_counts(tuple(t['counts']))
        by_canonical.setdefault(key,[]).append(t)
    rows=[]
    for e in expected:
        cc=tuple(e['canonical_field_counts'])
        sector=by_canonical.get(cc,[])
        exact_name=e['Matchete_generator_candidate_symbol']
        named=[t for t in sector if exact_name in t['native_coupling_symbols']]
        rows.append({'rgbeta_symbol':e['rgbeta_symbol'],
            'matching_field_sector_term_indices':[t['term_index'] for t in sector],
            'matching_named_extra_term_indices':[t['term_index'] for t in named],
            'identifies_physical_RGBeta_map':False})
    return {'path':str(path),'sha256':file_sha256(path),'native_CG_count':len(cgs),
            'native_quartic_term_count':len(terms),'entries':rows,
            'named_extra_HHSS_matches':sum(bool(r['matching_named_extra_term_indices']) for r in rows if next(e for e in expected if e['rgbeta_symbol']==r['rgbeta_symbol'])['sector_HHSS']),
            'matching_operator_does_not_certify_RGBeta_normalization':True}


def _strict_positive_int(x:Any)->int:
    if type(x)!=int or not 1<=x<=3:raise V6SourceError('Invalid CG tensor dimension')
    return x


def parse_live_cg(entry:dict,*,model:str)->ExactCG|None:
    """Exact Matchete live export tensor: no float, roundoff or inferred I."""
    counts=tuple(entry['canonical_field_counts'])
    idx=entry['invariant_index']
    if entry.get('coupling_symbol')!=extra_name(counts,idx):
        raise V6SourceError('Live extra-coupling naming disagrees with generator')
    symbol=extra_name(counts,idx,cg=True)
    if entry.get('CG_symbol')!=symbol:raise V6SourceError('Live CG naming disagrees with generator')
    if not entry.get('present_in_validated_Lagrangian'):
        raise V6SourceError('Extra operator not present in validated Lagrangian')
    if entry.get('implicit_singlet_CG') is True:
        if entry.get('CG_dimensions')!=[] or entry.get('CG_components_inputform')!=['1']:
            raise V6SourceError('Invalid implicit singlet tensor')
        return None
    dims=tuple(_strict_positive_int(x) for x in entry['CG_dimensions'])
    if not dims:raise V6SourceError('Non-singlet CG dimensions missing')
    # No shape guess. Matchete's native GetCGTensor gives the rank and values.
    values=entry['CG_components_inputform']
    from math import prod
    if not isinstance(values,list) or len(values)!=prod(dims):
        raise V6SourceError('Live native CG length does not match shape')
    exact=[_exact_number(v) for v in values]
    if any(v.has(sp.Float) for v in exact):
        raise V6SourceError('Live Matchete CG must contain exact numbers, not machine precision')
    if not any(v!=0 for v in exact):raise V6SourceError('Live native CG is zero')
    from itertools import product
    entries=tuple((k,v) for k,v in zip(product(*(range(d) for d in dims)),exact) if v!=0)
    # Only dimensions are verified here, not the physical bar/dual convention.
    expected_dims=[]
    d1,d2,_=MODEL_FAMILIES[model[0]]
    for i,c in enumerate(counts):expected_dims.extend(([2,2,d1,d1,d2,d2][i],)*c)
    expected_dims=sorted(x for x in expected_dims if x>1)
    if sorted(dims)!=expected_dims:
        raise V6SourceError('Live CG rank/dimensions inconsistent with scalar field multiset')
    digest=sha256(json.dumps({'dims':dims,'exact':[sp.sstr(x) for x in exact]},sort_keys=True).encode()).hexdigest()
    return ExactCG(symbol,tuple(f'dimension-{d}' for d in dims),dims,entries,digest)


def parse_live_operator_export(doc:dict,model:str,expected:list[dict])->dict:
    """Enforce exact operator multiplicities and actual validated Matchete source."""
    if doc.get('format')!='t3-ws2-live-Matchete-extra-quartic-CGs-v6':
        raise V6SourceError('Wrong live Wolfram export format')
    if doc.get('model')!=model or doc.get('status')!='LiveExtraQuarticCGExportVerified':
        raise V6SourceError('Live Matchete model/status mismatch')
    if not doc.get('full_matchete_lagrangian_validation_passed'):
        raise V6SourceError('Matchete did not validate full Lagrangian')
    observed=doc.get('operators')
    if not isinstance(observed,list) or len(observed)!=len(expected):
        raise V6SourceError('Matchete exported incomplete operator inventory')
    expected_by_key={}
    for row in expected:
        k=tuple(row['canonical_field_counts'])
        expected_by_key[k]=expected_by_key.get(k,0)+1
    observed_by_key={}; CGs={}
    for row in observed:
        cc=canonical_counts(tuple(row['canonical_field_counts']))
        idx=row.get('invariant_index')
        if type(idx)!=int or idx<=0:raise V6SourceError('Invalid invariant index')
        observed_by_key.setdefault(cc,set())
        if idx in observed_by_key[cc]:raise V6SourceError('Duplicate exported invariant')
        observed_by_key[cc].add(idx)
        cg=parse_live_cg(row,model=model)
        CGs[(cc,idx)]=cg
    for cc,count in expected_by_key.items():
        if observed_by_key.get(cc)!=set(range(1,count+1)):
            raise V6SourceError(f'Incomplete live Matchete invariant list for {cc}')
    if set(observed_by_key)!=set(expected_by_key):
        raise V6SourceError('Unexpected extra operator sector')
    return {'live_source_verified':True,'live_CG_count':sum(v is not None for v in CGs.values()),
            'live_operator_count':len(observed),'live_operators':[
                {'coupling_symbol':row['coupling_symbol'],
                 'canonical_field_counts':row['canonical_field_counts'],
                 'invariant_index':row['invariant_index'],
                 'CG_symbol':row['CG_symbol'],
                 'CG_dimensions':row['CG_dimensions'],
                 'CG_nonzero_components':(CGs[(tuple(row['canonical_field_counts']),row['invariant_index'])].invariant_summary() if CGs[(tuple(row['canonical_field_counts']),row['invariant_index'])] else None),
                 'physical_RGBeta_conversion_certified':False} for row in observed],
            '_CG_objects':CGs}


def _reference_exact_number(text:str)->sp.Expr:
    """Parse only SymPy's exact rational/sqrt/I numeric printing grammar."""
    if not isinstance(text,str) or not re.fullmatch(r'[0-9\s()+*/.,\-Isqrt]+',text):
        raise V6SourceError('Unsupported reference coefficient grammar')
    names=re.findall(r'[A-Za-z_]+',text)
    if any(n not in ('sqrt','I') for n in names):
        raise V6SourceError('Unknown reference coefficient identifier')
    x=sp.sympify(text,locals={'sqrt':sp.sqrt,'I':sp.I},evaluate=True)
    if x.free_symbols or not x.is_number or x.has(sp.nan,sp.zoo,sp.oo,-sp.oo):
        raise V6SourceError('Invalid reference coefficient')
    return sp.simplify(x)


def proportional_candidate(a:list[str],b:list[str])->str|None:
    if len(a)!=len(b) or not a or not b:raise V6SourceError('Incompatible exact basis vectors')
    av=sp.Matrix([_reference_exact_number(x) for x in a]);bv=sp.Matrix([_reference_exact_number(x) for x in b])
    if not any(x!=0 for x in av) or not any(x!=0 for x in bv):return None
    pos=next(i for i,x in enumerate(bv) if x!=0)
    factor=sp.simplify(av[pos]/bv[pos])
    if all(sp.simplify(x-factor*y)==0 for x,y in zip(av,bv)):
        return sp.sstr(factor)
    return None


def orientation_candidate_ratios(model:str,registry_entries:list[dict],live:dict,v4:dict,*,cap:int=128)->list[dict]:
    """Prove candidate *tensor proportionality*, never physical coupling factors."""
    by_symbol={x['rgbeta_symbol']:x for x in v4['entries']}
    out=[]
    for entry in registry_entries:
        if not entry['sector_HHSS']:continue
        target=by_symbol.get(entry['rgbeta_symbol'])
        if target is None:raise V6SourceError('V4 source projection missing symbol')
        cc=tuple(entry['canonical_field_counts'])
        if tuple(entry['RGBeta_field_counts'])!=cc:
            out.append({'rgbeta_symbol':entry['rgbeta_symbol'],
                'candidate_ratios':[],'candidate_count':0,'physical_map_certified':False,
                'reason':'Conjugated field multiset needs an independently certified dual mapping'})
            continue
        cg_choices=[(k,cg) for k,cg in live['_CG_objects'].items() if k[0]==cc]
        comparisons=[]
        for (key,index),cg in cg_choices:
            if cg is None:continue
            maps=candidate_cg_reference_maps(model,tuple(entry['RGBeta_field_counts']),cg)
            for m in maps:
                for rgb in target['candidate_orientation_projection_maps']:
                    ratio=proportional_candidate(m['exact_reference_coefficients'],
                                                rgb['source_to_reference_polynomial_coefficients'])
                    if ratio is not None:
                        comparisons.append({'Matchete_candidate_invariant_index':index,
                            'Matchete_index_placement':m['CG_axis_to_canonical_field_slots'],
                            'Matchete_dual_map_unverified':m['unverified_dual_map_field_types'],
                            'RGBeta_dual_map_unverified':rgb['dual_map_field_types'],
                            'exact_CG_polynomial_proportionality_ratio_not_a_coupling_factor':ratio})
            if len(comparisons)>cap:raise V6SourceError('Too many orientation candidates')
        out.append({'rgbeta_symbol':entry['rgbeta_symbol'],'candidate_ratios':comparisons,
                    'candidate_count':len(comparisons),'physical_map_certified':False,
                    'RGBeta_AddQuartic_coefficient_and_HC_still_unverified':True})
    return out


def experimental_inputform_evidence(path:Path,model:str,expected:list[dict])->dict:
    """Existing RunFullGeneralScalarBuild.wl UV artifacts are *not* EFT1 seeds.

    Recognises a validated general-scalar UV output and reports whether its
    named extra complex couplings occur. Never manufactures the missing CGs.
    """
    if not path.name.endswith('.interaction.inputform.txt') or not seed_belongs_to_model(path,model):
        raise V6SourceError('Wrong Matchete experimental interaction file/model')
    if path.stat().st_size>20_000_000:
        raise V6SourceError('Experimental input-form export unexpectedly large')
    text=path.read_text(encoding='utf8')
    found=set(re.findall(r'\blambdaT3ExtraQ\d{6}I\d+\b',text))
    names={row['Matchete_generator_candidate_symbol'] for row in expected}
    meta_path=Path(str(path)[:-len('.interaction.inputform.txt')])
    metadata=None
    if meta_path.is_file():
        metadata=json.loads(meta_path.read_text(encoding='utf8'))
    valid=(metadata is not None and metadata.get('full_matchete_check_passed') is True
           and metadata.get('quartic_included')==len(expected)
           and metadata.get('quartic_expected')==len(expected))
    return {'path':str(path),'sha256':file_sha256(path),
            'source_level':'experimental_UV_Lagrangian_not_EFT1_matching_seed',
            'experimental_full_Lagrangian_validation_recorded':valid,
            'metadata_sha256':file_sha256(meta_path) if metadata is not None else None,
            'extra_named_couplings_observed':sorted(found & names),
            'expected_named_extra_coupling_count':len(names),
            'all_expected_extra_names_present':found.issuperset(names),
            'cannot_certify_source_CG_or_EFT1_boundary_from_inputform_alone':True}
