"""V5 exact Matchete CG↔RGBeta HHSS comparison with evidence gates.

This deliberately separates (i) exact source CG tensor parsing, (ii) candidate
SU2 coordinate projections and (iii) a **physical** coupling conversion. (i)
and (ii) can be source-verified without asserting (iii). No running/matching
solvers are reimplemented.
"""
from __future__ import annotations

from itertools import product
from pathlib import Path
import hashlib
import json
import re
import sympy as sp

from .source_cg_v5 import ExactCG, MatcheteCGError, candidate_raw_tensor
from .named_hhss_projection import _dual_variants, _symmetrize, projection_report
from .extra_hhss_geometry import DIMS, reference_invariants, _check_counts

FORMAT='t3-ws2-matchete-rgbeta-named-tensor-dictionary-v5'


def exact_basis_map(source_columns:list[sp.Matrix],target_columns:list[sp.Matrix])->dict:
    """Exact (not least-squares) map: T_target = T_source * M.

    No coupling conversion is certified unless the caller separately proves
    the source columns have exactly the declared *interaction normalization*.
    """
    if not source_columns or not target_columns:
        raise MatcheteCGError('Nonempty tensor column sets required')
    nr=source_columns[0].rows
    if any(x.cols!=1 or x.rows!=nr for x in source_columns+target_columns):
        raise MatcheteCGError('Incompatible flattened tensor shapes')
    A=sp.Matrix.hstack(*source_columns)
    B=sp.Matrix.hstack(*target_columns)
    rA=A.rank();rB=B.rank()
    if rA<len(source_columns):
        return {'status':'dependent_source_tensors', 'source_rank':rA, 'target_rank':rB,
                'coefficients':None,'full_residual_is_zero':False}
    columns=[]
    for j in range(B.cols):
        try:
            solution,parameters=A.gauss_jordan_solve(B[:,j])
            if parameters.rows:raise MatcheteCGError('Nonunique solution')
        except ValueError:
            return {'status':'target_outside_source_span','source_rank':rA,
                    'target_rank':rB,'coefficients':None,'full_residual_is_zero':False}
        if any(sp.simplify(x)!=0 for x in A*solution-B[:,j]):
            raise MatcheteCGError('Inconsistent exact projection residual')
        columns.append(solution)
    M=sp.Matrix.hstack(*columns)
    return {'status':'exact_full_component_identity', 'source_rank':rA, 'target_rank':rB,
            'coefficients':[[sp.sstr(sp.simplify(M[i,j])) for j in range(M.cols)]
                            for i in range(M.rows)],
            'full_residual_is_zero':True,'reversible':M.rows==M.cols and M.det()!=0 if M.rows==M.cols else False}


def candidate_cg_reference_maps(model:str,counts:tuple[int,...],cg:ExactCG,
                                *,limit:int=48)->list[dict]:
    """Find exact **candidate** CG axis/dual maps in the canonical SU2 span.

    Axis assignments based only on dimensions are not a certification of the
    `CG[...]` slot references in Matchete's original `TermInputForm`.
    """
    counts=_check_counts(counts)
    types,dims,base=reference_invariants(*DIMS[model[0]],counts)
    ref=sp.Matrix.hstack(*base)
    dim_by_type=(2,2,DIMS[model[0]][0],DIMS[model[0]][0],
                 DIMS[model[0]][1],DIMS[model[0]][1])
    placements=candidate_raw_tensor(cg,types,dim_by_type)
    rows=[]
    seen=set()
    for placement in placements:
        vec=_symmetrize(placement['raw_tensor'],types,dims)
        for dual_tags,w in _dual_variants(vec,types,dims):
            if len(rows)>=limit:raise MatcheteCGError('Too many CG orientation candidates')
            try:
                sol,pars=ref.gauss_jordan_solve(w)
                matched=pars.rows==0 and all(sp.simplify(x)==0 for x in ref*sol-w)
            except ValueError:
                matched=False
                sol=None
            if not matched:continue
            norm=[sp.sstr(sp.simplify(v)) for v in sol]
            label=(tuple(placement['axis_to_external_slot']),tuple(dual_tags),tuple(norm))
            if label in seen:continue
            seen.add(label)
            rows.append({'CG_axis_to_canonical_field_slots':placement['axis_to_external_slot'],
                         'unverified_dual_map_field_types':[
                            ('H','BarH','S1','BarS1','S2','BarS2')[k] for k in dual_tags],
                         'exact_reference_coefficients':norm,
                         'component_residual_zero':True})
    return rows


def native_seed_sector_matches(native_terms:list[dict],counts:tuple[int,...],
                               expected_symbol:str|None=None)->list[dict]:
    matches=[]
    for row in native_terms:
        if tuple(row['counts'])!=tuple(counts):continue
        names=row['native_coupling_symbols']
        matches.append({'term_index':row['term_index'],
                        'native_coupling_symbols':names,'CGNames':row['CGNames'],
                        'term_source_sha256':row['term_source_sha256'],
                        'is_exact_named_coupling_match':expected_symbol in names if expected_symbol else False,
                        'term_may_be_different_operator':expected_symbol not in names if expected_symbol else True})
    return matches


def per_named_HHSS_evidence(model:str,source_entry:dict,native_terms:list[dict],
                            cgs:dict[str,ExactCG],prior_v4:dict)->dict:
    """Compare source sector and CGs with prior V4 exact RGBeta projections.

    The Matchete source needs a *matching named operator* or an independently
    established operator-basis change. Same field counts do not suffice.
    """
    symbol=source_entry['declaration']['symbol']
    counts=tuple(source_entry['counts'])
    v4=next((x for x in prior_v4['entries'] if x['rgbeta_symbol']==symbol),None)
    if v4 is None or tuple(v4['counts'])!=counts:
        raise MatcheteCGError(f'{model}/{symbol}: V4 source sector mismatch')
    # `new` registry names are explicit candidate Matchete extra couplings.
    match_re=re.match(r'lambdaT3GeneralQ(\d{6})I(\d+)$',symbol)
    expected=(f'lambdaT3ExtraQ{match_re.group(1)}I{match_re.group(2)}'
              if match_re else symbol)
    matches=native_seed_sector_matches(native_terms,counts,expected)
    named=[m for m in matches if m['is_exact_named_coupling_match']]
    candidates=[]
    for row in named:
        for name in row['CGNames']:
            cg=cgs[name]
            projections=candidate_cg_reference_maps(model,counts,cg)
            candidates.append({'term_index':row['term_index'],'CG':name,
                               'candidate_CG_coordinate_maps':projections,
                               'all_candidate_coordinate_maps_unverified':True})
    # A native singlet with no CG still requires explicitly matching term
    # sign, 1/2/HC and Matchete scalar/phase map before physical use.
    statuses=[]
    if not named:statuses.append('matching_named_Matchete_extra_quartic_absent_in_native_EFT1_seed')
    if named and any(not t['CGNames'] for t in named):statuses.append('unnamed_unit_CG_or_singlet_index_requires_source_term_normalization')
    if named and any(not x['candidate_CG_coordinate_maps'] for x in candidates):
        statuses.append('no_candidate_Matchete_CG_orientation_in_reference_span')
    if not v4.get('candidate_orientation_projection_maps'):
        statuses.append('RGBeta_tensor_missing_canonical_Ward_valid_orientation')
    statuses.extend(('RGBeta_AddQuartic_internal_normalization_not_certified',
                     'Wolfram_Matchete_CG_index_and_phase_coordinate_map_not_certified',
                     'complex_coupling_signed_HC_and_two_threshold_running_not_certified'))
    return {'rgbeta_symbol':symbol,'counts_H_BarH_S1_BarS1_S2_BarS2':list(counts),
            'source_category':source_entry['category'],
            'expected_native_symbol_for_direct_name_match':expected,
            'same_field_sector_native_terms':matches,
            'matching_named_native_terms':len(named),
            'candidate_native_CG_orientations':candidates,
            'prior_RGBeta_valid_Ward_orientation_count':v4['valid_canonical_invariant_orientation_count'],
            'source_based_candidate_tensor_reconstruction':bool(candidates),
            'certified_Matchete_RGBeta_coupling_matrix':None,
            'physical_map_certified':False,'physical_numerical_coupling_at_muF':None,
            'physical_numerical_coupling_at_muS':None,
            'blocking_dependencies':statuses}


def summarize_all_extra_registry(model:str,registry:dict,terms:list[dict])->dict:
    """Complete 50-invariant inventory, not just HHSS.

    Structural count deliberately does not imply independent physical CG maps.
    """
    rows=[]
    for category,namekey in (('new','name'),('preexisting','existing_name')):
        for item in registry.get(category+'_rgbeta_symbols',[]):
            counts=tuple(item['counts']);name=item[namekey]
            if sum(counts)!=4:raise MatcheteCGError(f'{model}/{name}: invalid quartic counts')
            refs=native_seed_sector_matches(terms,counts)
            rows.append({'rgbeta_symbol':name,'category':category,
                         'counts_H_BarH_S1_BarS1_S2_BarS2':list(counts),
                         'sector_is_HHSS':sum(counts[:2])==2 and sum(counts[2:])==2,
                         'native_term_same_sector_count':len(refs),
                         'has_explicit_source_normalization':False})
    return {'model':model,'additional_quartic_registry_entries':len(rows),
            'HHSS_entries':sum(x['sector_is_HHSS'] for x in rows),'entries':rows,
            'physical_conversion_status':'source_sector_inventory_only'}
