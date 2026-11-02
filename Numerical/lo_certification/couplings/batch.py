"""All-16 T3 EFT1 coupling inventory, numerical snapshot and conditional basis handoff.

Reads only existing source artifacts, never fills missing quartics with zero and
never promotes the provisional Matchete GammaCC/Weyl factor to physical C5.
Output is exclusive-write and restricted to Workstream 2's output directory.

The existing seven-term and joint matching modules remain sole solvers.
"""
from __future__ import annotations

from dataclasses import replace
from hashlib import sha1, sha256
from pathlib import Path
from typing import Any
import json
import math

import numpy as np
import sympy as sp

from .basis import (checked_eigenbasis, mass2_matrix, relative_basis_rotation,
                    transport_hhss, transport_mass2, transport_llss)

FORMAT = 't3-lo-coupling-dictionary-batch-v1'
MODEL_DIMS = {'A': (1,3,2), 'B': (2,2,1), 'C': (2,2,3),
              'D': (3,1,2), 'E': (3,3,2)}
# Verified source-file identities at GitHub main tree 2c29e5ba2c0350b7e9be292e02ceb8bc7c52c0d3,
# 2026-11-02.  These are *Git blob* SHA1s, not physical certificates.
SOURCE_BLOBS = {
    'Numerical/joint_scalar_matching/source_quartic_v4.py': '5ca01592a1c0e0780e345331b86e684cfa00bbda',
    'Numerical/joint_scalar_matching/numerical_benchmark_v5.py': '7b4fdbb58b2a1c5466264c4c0b9632a19b6508ce',
    'Numerical/joint_scalar_matching/special_mass_v6.py': '2fb97f9eda33b1780336124f5380d085b51ab240',
    'Numerical/joint_scalar_matching/large_batch.py': 'b0c01bd910009cc6a06bf97dec11acb9a081b045',
    'RGE/running/intermediate/GeneralEFT1C5Mixing.py': '2ceafa9afd09375a96b2a54a3eec30468ebe93a5',
    'RGE/running/intermediate/IntegratedEFT1JointScalarBridge.py': 'b0f280aa6bf96dbc4e0a8a3d5994f1d6a5300636',
    'RGE/running/intermediate/ScalarOnlyTensorAdapters.py': 'dc4e6039c1475960be6d6c694fcca8509ec13394',
    'RGE/running/rgbeta/T3RGBetaModel.wl': '7af7a4ddf99a7e38210543a21b3b1cda203eaf02',
    'RGE/general/generated/general_rgbeta_quartic_registry.json': 'bf674314d60b6b4f8fdb7b9a5d200038611cbdf7',
}
# Source-v5 calls these explicitly CANDIDATES; they are NOT physical equalities.
CANDIDATE_ALIAS = {'lambdaT3':'lambdaT3',
                   'lambdaH1Inv1':'lambdaH1', 'lambdaH1Inv2':'lambdaH1Adj',
                   'lambdaH2Inv1':'lambdaH2', 'lambdaH2Inv2':'lambdaH2Adj'}


# HHSS candidates recorded in the inspected *provisional* RGBeta registry.
# This is a regression floor, NOT a certified complete scalar potential.
BASELINE_HHSS_REGISTRY = {
    'A_-4': (),
    'A_-2': ('new:lambdaT3GeneralQ020101I1', 'new:lambdaT3GeneralQ110002I1'),
    'A_+0': ('new:lambdaT3GeneralQ021010I1', 'new:lambdaT3GeneralQ110200I1'),
    'B_-3': ('new:lambdaT3GeneralQ020002I1',),
    'B_-1': ('new:lambdaT3GeneralQ110101I1', 'old:lambdaHHdagS1barS1bar',
             'old:lambdaHHdagS1barS2barCross', 'old:lambdaHHdagS2S2'),
    'B_+1': ('new:lambdaT3GeneralQ022000I1',),
    'C_-3': ('new:lambdaT3GeneralQ020002I1',),
    'C_-1': ('new:lambdaT3GeneralQ110101I1', 'old:lambdaHHdagS1barS1bar',
             'old:lambdaHHdagS1barS2barCross', 'old:lambdaHHdagS2S2'),
    'C_+1': ('new:lambdaT3GeneralQ022000I1',),
    'D_-2': ('new:lambdaT3GeneralQ020101I1', 'new:lambdaT3GeneralQ110002I1'),
    'D_+0': ('new:lambdaT3GeneralQ021010I1', 'new:lambdaT3GeneralQ110200I1'),
    'D_+2': (),
    'E_-4': (),
    'E_-2': ('new:lambdaT3GeneralQ020101I1', 'new:lambdaT3GeneralQ110002I1'),
    'E_+0': ('new:lambdaT3GeneralQ021010I1', 'new:lambdaT3GeneralQ110200I1'),
    'E_+2': (),
}
PARSER_PATH = 'RGE/running/intermediate/ScalarOnlyTensorAdapters.py'
REGISTRY_PATH = 'RGE/general/generated/general_rgbeta_quartic_registry.json'
# Changes to these two files can be rechecked using their actual local content.
# Any other changed source is still blocked unless the user explicitly requests
# the pre-existing diagnostic --allow-local-version-drift mode.
RECHECKABLE_PATHS = frozenset((PARSER_PATH, REGISTRY_PATH))


def _registry_hhss_names(row: dict) -> set[str]:
    found: list[str] = []
    for origin, key, field in (
        ('new', 'new_rgbeta_symbols', 'name'),
        ('old', 'preexisting_rgbeta_symbols', 'existing_name'),
    ):
        entries = row.get(key)
        if not isinstance(entries, list):
            raise ValueError(f'RGBeta registry missing {key}')
        for entry in entries:
            if not isinstance(entry, dict):
                raise ValueError('Invalid RGBeta invariant entry')
            counts = entry.get('counts')
            if not (isinstance(counts, list) and len(counts) == 6 and
                    all(type(x) is int and x >= 0 for x in counts) and
                    sum(counts) == 4):
                raise ValueError('Invalid RGBeta invariant field counts')
            if sum(counts[:2]) != 2 or sum(counts[2:]) != 2:
                continue
            name = entry.get(field)
            if not isinstance(name, str) or not name.strip():
                raise ValueError('Unnamed RGBeta HHSS candidate')
            found.append(f'{origin}:{name}')
    if len(found) != len(set(found)):
        raise ValueError('Duplicate RGBeta HHSS invariant symbol in registry')
    return set(found)


def local_source_recheck(root: Path, models: tuple[str, ...], versions: dict,
                         release: Path, cache: Path, registry: dict) -> dict:
    """Re-audit only the two supported locally changed source files.

    A parser drift is accepted ONLY when recomputing the Matchete CG tensors
    directly from all selected native seeds reproduces each hash-checked
    exported cache component exactly. The candidate registry may add new
    entries but must not silently remove a previously inventoried HHSS term.
    Neither check certifies an independent quartic normalization or C5.
    """
    drift = sorted(rel for rel, v in versions.items()
                   if not v['matches_inspected_main'])
    disallowed = sorted(set(drift) - RECHECKABLE_PATHS)
    if disallowed:
        raise ValueError('Unreviewed source changes outside local re-audit scope: ' +
                         ', '.join(disallowed))
    report = {'status': 'inspected_main_exact' if not drift else 'local_structural_recheck',
              'drift': drift, 'parser_reparsed_models': [],
              'registry_models_checked': [], 'registry_new_HHSS_candidates': {},
              'physical_source_certified': False}
    if REGISTRY_PATH in drift:
        if registry.get('status') != 'GeneratedCandidateBasisNotWolframVerified':
            raise ValueError('Local RGBeta registry status changed; requires a new audit')
        for model in models:
            row = registry.get('models', {}).get(model)
            if not isinstance(row, dict):
                raise ValueError(f'{model}: local RGBeta registry entry missing')
            observed = _registry_hhss_names(row)
            missing = set(BASELINE_HHSS_REGISTRY[model]) - observed
            if missing:
                raise ValueError(f'{model}: prior HHSS registry invariants disappeared: '+
                                 ', '.join(sorted(missing)))
            report['registry_models_checked'].append(model)
            extras = sorted(observed - set(BASELINE_HHSS_REGISTRY[model]))
            if extras:
                report['registry_new_HHSS_candidates'][model] = extras
    if PARSER_PATH in drift:
        from Numerical.joint_scalar_matching.fermion_boundary_bridge import load_source
        from Numerical.joint_scalar_matching.numerical_benchmark_v5 import read_v4_cached_channels
        from Numerical.joint_scalar_matching.source_quartic_v4 import (
            source_seed_path, build_source_channels)
        from Numerical.joint_scalar_matching.fermion_boundary_bridge import slug
        for model in models:
            source = load_source(release/slug(model)/'boundary.json', expected_model=model)
            seed = source_seed_path(root/'output/full16_symbolic', model)
            cached = read_v4_cached_channels(cache/f'{slug(model)}.json', seed, source)
            built = build_source_channels(seed, model)
            if set(built.by_name) != set(cached.by_name):
                raise ValueError(f'{model}: changed parser has different HHSS channel inventory')
            for name in cached.by_name:
                expected, observed = cached.by_name[name], built.by_name[name]
                if set(expected) != set(observed):
                    raise ValueError(f'{model}/{name}: changed parser lost or added scalar components')
                for key in expected:
                    if sp.simplify(expected[key] - observed[key]) != 0:
                        raise ValueError(f'{model}/{name}/{key}: changed parser alters V4 derivative')
            report['parser_reparsed_models'].append(model)
    return report



def digest(path: Path) -> str:
    return sha256(Path(path).read_bytes()).hexdigest()


def git_blob_id(path: Path) -> str:
    data = Path(path).read_bytes()
    return sha1(f'blob {len(data)}\0'.encode() + data).hexdigest()


def read_json(path: Path) -> dict:
    obj = json.loads(Path(path).read_text(encoding='utf-8-sig'))
    if not isinstance(obj, dict):
        raise ValueError(f'Expected JSON object: {path}')
    return obj


def write_json(path: Path, data: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x', encoding='utf-8') as f:
        f.write(json.dumps(data, indent=2, allow_nan=False, sort_keys=True) + '\n')


def scalar_native_catalog(ch) -> dict:
    return {name: [{'indices_HHSS_zero_based':list(key),
                    'V4_derivative_exact':sp.sstr(expr)}
                   for key,expr in sorted(entries.items())]
            for name,entries in sorted(ch.by_name.items())}


def nonzero_sources(model_registry: dict, native_names: set[str]) -> list[dict]:
    """Structural registry is candidate/not-normalized: preserve each HHSS invariant."""
    ret = []
    for item in model_registry.get('new_rgbeta_symbols', []):
        if sum(item['counts'][:2]) == 2 and sum(item['counts'][2:]) == 2:
            ret.append({'rgbeta':item['name'], 'matchete_basis_name':item['matchete_basis_name'],
                        'mapping_status':item.get('mapping_status','NotNormalised'),
                        'represented_as_certified_native_HHSS':False,
                        'required_proof':'independent full real-component quartic tensor and scale/phase normalization'})
    for item in model_registry.get('preexisting_rgbeta_symbols', []):
        if sum(item['counts'][:2]) == 2 and sum(item['counts'][2:]) == 2:
            ret.append({'rgbeta':item['existing_name'],
                        'represented_as_certified_native_HHSS':False,
                        'mapping_status':'NotNormalizedToMatchete',
                        'required_proof':'independent tensor equality, including complex conjugates'})
    return ret


def observed_running_states(path: Path, model: str) -> dict:
    doc = read_json(path)
    if doc.get('model') != model or doc.get('status') != 'NumericalRenormalisableRunningOnly':
        raise ValueError('Unexpected all16 running model/status')
    scales=doc.get('scales_gev', {})
    if not (0 < float(scales['S']) < float(scales['F']) < float(scales['UV'])):
        raise ValueError('Must have exactly one grouped muF and one grouped muS')
    if doc.get('illustrative_benchmark') is not True or doc.get('physical_mass_parameters_available') is not False:
        raise ValueError('Numerical snapshot physical status unexpectedly changed; re-audit')
    if doc.get('eft1_wilson_running_performed') is not False:
        raise ValueError('Unexpected EFT1 Wilson evolution in renormalisable snapshot')
    f=doc['eft1_at_fermion_threshold']; s=doc['eft1_at_scalar_threshold']
    for k in ('mS1Sq','mS2Sq','lambdaT3','gY','g2','g3','ye','yu','yd'):
        if k not in f or k not in s:
            raise ValueError(f'missing EFT1 running field {k} at muF or muS')
    for state in (f,s):
        for k in ('mS1Sq','mS2Sq'):
            if not math.isfinite(float(state[k])) or float(state[k])<=0:
                raise ValueError(f'nonpositive {k}');
    return doc


def calculate_candidate_quartics(ch, state_s: dict):
    """Uses source v4's exact CG parser and explicitly provisional v5 aliases."""
    missing_alias=sorted(set(ch.by_name)-set(CANDIDATE_ALIAS))
    missing_running=sorted({CANDIDATE_ALIAS[x] for x in ch.by_name if x in CANDIDATE_ALIAS}
                           - set(state_s))
    if missing_alias or missing_running:
        return None, (), {}, missing_alias, missing_running
    from Numerical.joint_scalar_matching.source_quartic_v4 import evaluate_channels
    mapping={k:CANDIDATE_ALIAS[k] for k in sorted(ch.by_name)}
    vals={k:state_s[v] for k,v in mapping.items()}
    arr,names=evaluate_channels(ch,vals)
    return arr,names,mapping,[],[]


def _running_at_f_snapshot(snapshot, state_f):
    """Reuse the same source-verified special-mass builder at *both* endpoints."""
    return replace(snapshot, scalar_mass1_sq=float(state_f['mS1Sq']),
                   scalar_mass2_sq=float(state_f['mS2Sq']),couplings=state_f)


def calculate_candidate_masses(model: str, boundary, snap, f_state: dict):
    """muF and muS basis matrices are constructed independently; never R_F=R_S.

    The mixed prototype is conditional because its native epsilon/phase map
    has not been independently certified against the Matchete LLSS source.
    """
    from Numerical.joint_scalar_matching.numerical_benchmark_v5 import SPECIAL
    from Numerical.joint_scalar_matching.special_mass_v6 import build_special_basis
    n=len(boundary.labels)
    if model not in SPECIAL:
        d1,_,_=MODEL_DIMS[model[0]]
        m_f=np.diag([float(f_state['mS1Sq'])]*(2*d1)+
                    [float(f_state['mS2Sq'])]*(n-2*d1))
        m_s=np.diag([snap.scalar_mass1_sq]*(2*d1)+
                    [snap.scalar_mass2_sq]*(n-2*d1))
        rf=rs=np.eye(n)
        interpretation='same original interleaved S1/S2 real basis; unbroken EW illustrative masses'
        phase_open=False
    else:
        data_f=build_special_basis(model,_running_at_f_snapshot(snap,f_state),boundary)
        data_s=build_special_basis(model,snap,boundary)
        rf=np.asarray(data_f['R_original_from_boundary'],float)
        rs=np.asarray(data_s['R_original_from_boundary'],float)
        m_f=np.asarray(data_f['mass2_boundary'],float)
        # Independent instantaneous muS rotation; express muS mass in fixed muF Wilson basis.
        m_original_s=rs@np.asarray(data_s['mass2_boundary'],float)@rs.T
        m_s=transport_mass2(m_original_s,rf)
        interpretation=data_s['mass_basis_assumption']
        phase_open=True
    mass2_matrix(m_f,n,'muF scalar Hessian')
    mass2_matrix(m_s,n,'muS scalar Hessian in muF basis')
    rel=relative_basis_rotation(rf,rs)
    vals,u=checked_eigenbasis(m_s)
    if not np.allclose(np.sort(vals),np.sort(np.linalg.eigvalsh(m_s)),rtol=1e-12):
        raise ValueError('Scalar eigenvalues changed under rotation')
    return {'R_original_from_F':rf,'R_original_from_S':rs,
            'U_F_from_S':rel,'mass2_at_F_in_F_basis':m_f,
            'mass2_at_S_in_F_basis':m_s,'diagonalizing_U_at_S_in_F_basis':u,
            'mass2_eigenvalues_at_S_gev2':vals,
            'interpretation':interpretation,'special_phase_certificate_missing':phase_open}


def check_source_versions(root: Path) -> dict:
    report={}
    for rel,sha in SOURCE_BLOBS.items():
        path=root/rel
        observed=git_blob_id(path) if path.is_file() else None
        report[rel]={'git_blob_sha1_observed':observed,'git_blob_sha1_inspected':sha,
                     'sha256_observed':digest(path) if observed else None,
                     'matches_inspected_main':observed==sha}
    return report


def model_paths(root: Path, model: str, release: Path, native_cache: Path) -> dict[str,Path]:
    from Numerical.joint_scalar_matching.fermion_boundary_bridge import slug
    from Numerical.joint_scalar_matching.source_quartic_v4 import source_seed_path
    t=slug(model)
    return {'boundary':release/t/'boundary.json',
            'quartic_seed':source_seed_path(root/'output/full16_symbolic',model),
            'quartic_cache':native_cache/f'{t}.json',
            'running':root/'output/general_observable/all16_v1'/t/'running.json'}


def inspect_model(root: Path, out: Path, model: str, release: Path, cache: Path,
                  registry: dict, *, source_versions:dict) -> dict:
    from Numerical.joint_scalar_matching.fermion_boundary_bridge import load_source
    from Numerical.joint_scalar_matching.numerical_benchmark_v5 import (
        load_running_snapshot, read_v4_cached_channels)
    from RGE.running.intermediate.IntegratedEFT1JointScalarBridge import contracts
    p=model_paths(root, model, release, cache)
    out.mkdir(parents=True, exist_ok=False)
    blockers=[]; stages={}
    hashes={name:{'path':str(path.relative_to(root)).replace('\\','/'),
                  'sha256':digest(path) if path.is_file() else None}
            for name,path in p.items()}
    if any(h['sha256'] is None for h in hashes.values()):
        raise FileNotFoundError('missing source files: '+str([k for k,v in hashes.items() if not v['sha256']]))
    source=load_source(p['boundary'],expected_model=model)
    ch=read_v4_cached_channels(p['quartic_cache'],p['quartic_seed'],source)
    stages['matchete_LLSS_tree_source']='verified_formal_polynomial_not_physical_normalized'
    stages['native_quartic_V4_CG_tensor']='verified_exact_seed_and_cache'
    if source.normalization_certified is False:
        blockers.append('GammaCC_to_master_Weyl_signed_normalization_not_certified')
    reg=registry['models'][model]
    extra=nonzero_sources(reg,set(ch.by_name))
    if extra:
        blockers.append('additional_symmetry_allowed_HHSS_tensor_channels_not_normalized_or_present_in_seed')
    if reg.get('quadratics_not_implemented'):
        blockers.append('additional_RG_closed_quadratic_couplings_not_implemented_in_source_registry')
    aliases={k:{'RGBeta_symbol_candidate':CANDIDATE_ALIAS.get(k),
                'status':'unverified_relative_tensor_normalization_and_phase',
                'matchete_source_cg':True} for k in sorted(ch.by_name)}
    if any(row['RGBeta_symbol_candidate'] is None for row in aliases.values()):
        blockers.append('native_Matchete_to_RGBeta_quartic_alias_missing')
    blockers.extend(('full_SM_Weyl_Yukawa_tensor_and_complete_gauge_generator_embedding_not_certified',
                     'muF_muS_same_EFT1_renormalisable_trajectory_is_illustrative_only',
                     'additional_one_loop_EFT1_operator_and_hard_fermion_boundary_closure_not_certified'))
    contract=contracts(read_json(p['boundary']),source.hash)
    write_json(out/'existing_solver_contract.json',contract)
    write_json(out/'symbolic_hhss.json',{
        'format':'t3-source-exact-hhss-components-v1','model':model,
        'potential_convention':'V4=lambda_abcd phi_a phi_b phi_c phi_d/4!',
        'basis':'Matchete original S1,S2 interleaved (Re,Im)',
        'source_files':hashes,'native_channels':scalar_native_catalog(ch),
        'missing_independent_HHSS_channels':extra,
        'source_operator_inventory_complete':len(extra)==0 and reg.get('quadratics_not_implemented')==0,
        'native_coupling_alias_status':aliases})
    diag={'model':model,'format':FORMAT,'models_in_EFT1':['S1','S2'],
          'thresholds':['all heavy BSM fermions together','both heavy scalar multiplets together'],
          'mu_F_gev':source.mu_f_gev,'source_files':hashes,
          'source_compatibility':{rel:x['matches_inspected_main'] for rel,x in source_versions.items()},
          'scalar_real_labels_F':list(source.labels),
          'native_quartic_channels':sorted(ch.by_name),
          'extra_HHSS_registry_channels':extra,
          'native_vs_RGBeta':aliases,
          'quartic_tensor_status':stages['native_quartic_V4_CG_tensor'],
          'physical_spinor_factor':None, 'blockers':blockers}
    numeric=None
    try:
        doc=observed_running_states(p['running'],model)
        snap=load_running_snapshot(p['running'],model)
        if abs(snap.mu_f/source.mu_f_gev-1)>1e-11:
            raise ValueError('muF running/LLSS source mismatch')
        diag['mu_S_gev']=snap.mu_s
        diag['running_status']='verified_repository_illustrative_EFT1_MSbar_renormalizable_snapshot'
        write_json(out/'running_couplings_muF_muS.json',{
            'format':'t3-source-running-EFT1-couplings-v1','model':model,
            'MSbar_mu_F_gev':snap.mu_f,'MSbar_mu_S_gev':snap.mu_s,
            'file_sha256':snap.sha256,
            'muF_renormalizable_state':doc['eft1_at_fermion_threshold'],
            'muS_renormalizable_state':doc['eft1_at_scalar_threshold'],
            'gauge_variables':'gY is SM hypercharge coupling (not GUT rescaled), g2 and g3 as stored',
            'EFT1_SM_yukawa_matrices_available':['ye','yu','yd'],
            'mass_parameters_available':['mS1Sq','mS2Sq'],
            'physical_mass_parameters_available':False,
            'heavy_fermion_Yukawas_absent_below_muF':True,
            'complete_master_yukawa_and_gauge_tensors_constructed':False})
        native,names,mapping,ma,mc=calculate_candidate_quartics(ch,snap.couplings)
        if ma: blockers.append('unmapped_native_couplings:'+','.join(ma))
        if mc: blockers.append('missing_running_couplings:'+','.join(mc))
        if native is not None:
            b=calculate_candidate_masses(model,source,snap,doc['eft1_at_fermion_threshold'])
            c=transport_hhss(native,b['R_original_from_F'])
            sumq=c.sum(axis=0)
            from Numerical.joint_scalar_matching.core import check_tensors
            check_tensors(source.raw, sumq, len(source.labels))
            for individual in c:
                check_tensors(source.raw, individual, len(source.labels))
            # Inverse rotation is an independent coordinate-covariance check.
            inv=transport_hhss(c,b['R_original_from_F'].T)
            if not np.allclose(inv,native,atol=1e-9,rtol=1e-9):
                raise ValueError('Quaric round-trip R_F covariance failure')
            if b['special_phase_certificate_missing']:
                blockers.append('special_Matchete_scalar_phase_and_muF_basis_map_unverified')
            stages['same_muS_mass_hessian_positive']='verified_for_illustrative_running_inputs'
            stages['muF_muS_relative_basis']='derived_from_separate_running_endpoint_quadratic_inputs_conditional'
            np.savez_compressed(out/'conditional_numeric_tensors.npz',
                source_to_F_real_R=b['R_original_from_F'],
                source_to_S_real_R=b['R_original_from_S'],
                relative_F_from_S=b['U_F_from_S'],
                mass_squared_gev2_muF_fixed_F=b['mass2_at_F_in_F_basis'],
                mass_squared_gev2_muS_fixed_F=b['mass2_at_S_in_F_basis'],
                diagonalizing_U_muS_in_F=b['diagonalizing_U_at_S_in_F_basis'],
                lambda_HHSS_channels_candidate_muS_fixed_F=c,
                lambda_HHSS_candidate_muS_fixed_F=sumq,
                C_LLSS_formal_muF_in_F=source.raw)
            write_json(out/'conditional_numeric_metadata.json',{
                'format':'t3-provisional-not-physical-solver-background-v1',
                'model':model,'mu_F_gev':snap.mu_f,'mu_S_gev':snap.mu_s,
                'scheme':'MSbar', 'mass_origin':'unbroken_electroweak',
                'two_grouped_thresholds_only':True,
                'arrays_npz':'conditional_numeric_tensors.npz',
                'quartic_channels':list(names),'candidate_native_to_RGBeta_symbols':mapping,
                'joint_scalar_check_tensors_interface_passed':True,
                'extra_HHSS_not_included':extra,
                'scalar_basis_original_equals_R_F_times_boundary':True,
                'F_basis_fixed_for_EFT1_and_joint_scalar_matching':True,
                'muF_muS_basis_rotation':'relative_F_from_S = R_F.T @ R_S',
                'muS_mass_eigenvalues_gev2':b['mass2_eigenvalues_at_S_gev2'].tolist(),
                'muS_mass_eigenvalues_gev':np.sqrt(b['mass2_eigenvalues_at_S_gev2']).tolist(),
                'scalar_basis_assumption':b['interpretation'],
                'scalar_basis_certified':not b['special_phase_certificate_missing'],
                'quark_and_lepton_Yukawa_running_not_embedded_in_master_Weyl_tensor':True,
                'GammaCC_to_Weyl_factor':None,
                'quartic_aliases_physically_certified':False,
                'not_safe_to_use_for_complete_seven_term_or_physical_C5':True})
            numeric={'basis':b,'channels':c,'sum':sumq,'names':names}
            diag['candidate_numerical_tensor_shapes']={'HHSS_channels':list(c.shape),
                 'HHSS_sum':list(sumq.shape),'scalar_mass2_muS':list(b['mass2_at_S_in_F_basis'].shape)}
            diag['muS_mass_eigenvalues_gev2']=b['mass2_eigenvalues_at_S_gev2'].tolist()
        else:
            diag['candidate_numerical_tensor_shapes']=None
    except (ValueError,TypeError,KeyError,ImportError,IndexError,AssertionError) as exc:
        blockers.append('numerical_EFT1_snapshot_or_mass_basis_blocked:'+str(exc))
        diag['numerical_error']=str(exc)
    # Physical certification is an AND of independent conditions, not an average.
    if any(not row['matches_inspected_main'] for row in source_versions.values()):
        blockers.append('local_source_versions_differ_from_inspected_GitHub_main_20261102')
    diag['blockers']=list(dict.fromkeys(blockers))
    diag['stages']=stages
    diag['coupling_dictionary_status']='conditional' if numeric is not None else 'blocked'
    diag['physical_status']='blocked' if diag['blockers'] else 'verified'
    diag['physically_normalized_complete_master_couplings_available']=False
    diag['joint_matching_physical_C5_available']=False
    write_json(out/'status.json',diag)
    return {'model':model,'physical_status':diag['physical_status'],
            'native_CG_status':'verified','coupling_status':diag['coupling_dictionary_status'],
            'positive_illustrative_mass':numeric is not None,
            'missing_independent_HHSS':len(extra),
            'blocking_dependencies':diag['blockers'],
            'mu_F_gev':source.mu_f_gev,'mu_S_gev':diag.get('mu_S_gev')}


def run_batch(root: Path, output: Path, models: tuple[str,...] | None = None,
              *, reject_drift: bool = True) -> dict:
    """16 independent model audits; no solver reimplementation or UV/EFT writes."""
    root=Path(root).resolve(); output=Path(output).resolve()
    owned=(root/'output/general_observable/lo_certification_couplings').resolve()
    if output != owned and owned not in output.parents:
        raise ValueError('Output must be inside owned output/general_observable/lo_certification_couplings/')
    from Numerical.joint_scalar_matching.adapter import ALL_MODELS
    from Numerical.joint_scalar_matching.fermion_boundary_bridge import load_release_manifest
    if output.exists():
        raise FileExistsError(f'Refusing to overwrite existing result directory: {output}')
    chosen=tuple(ALL_MODELS if models is None else models)
    if not chosen or len(set(chosen))!=len(chosen) or not set(chosen)<=set(ALL_MODELS):
        raise ValueError('Invalid model selection')
    release=root/'output/general_observable/lo_fermion_boundary/releases/source_021dcc1e8fe8ba0f'
    cache=root/'output/general_observable/lo_joint_scalar_matching/source_quartic_v4/native_quartic_channels'
    manifest=load_release_manifest(release)
    registry_path=root/'RGE/general/generated/general_rgbeta_quartic_registry.json'
    registry=read_json(registry_path)
    if registry.get('status')!='GeneratedCandidateBasisNotWolframVerified' or set(chosen)-set(registry.get('models',{})):
        raise ValueError('Unexpected RGBeta candidate registry; verify again')
    versions=check_source_versions(root)
    drift=[path for path,x in versions.items() if not x['matches_inspected_main']]
    if reject_drift:
        # Fails closed on a changed tensor, missing old invariant, or changes
        # outside the two explicitly re-audited files. No silent SHA bypass.
        recheck = local_source_recheck(root, chosen, versions, release, cache, registry)
    else:
        # Explicitly requested exploratory mode; never labelled reverified.
        recheck = {'status': 'explicit_unaudited_drift_override', 'drift': drift,
                   'parser_reparsed_models': [], 'registry_models_checked': [],
                   'registry_new_HHSS_candidates': {}, 'physical_source_certified': False}
    output.mkdir(parents=True,exist_ok=False)
    write_json(output/'source_versions.json',{
        'github_remote':'https://github.com/Frankyang2002/T3neutrinomassmodelpipeline',
        'github_read_only':True,'inspected_tree_sha':'2c29e5ba2c0350b7e9be292e02ceb8bc7c52c0d3',
        'files':versions,'drift':drift,'local_source_recheck':recheck,
        'registry_source_sha256':digest(registry_path)})
    rows=[]
    for model in chosen:
        from Numerical.joint_scalar_matching.fermion_boundary_bridge import slug
        target=output/slug(model)
        try:
            rows.append(inspect_model(root,target,model,release,cache,registry,source_versions=versions))
        except Exception as exc:
            # A failed model is never replaced by a zero and never suppresses other models.
            target.mkdir(parents=True,exist_ok=True)
            failure={'model':model,'physical_status':'blocked','coupling_status':'blocked',
                     'native_CG_status':'blocked','positive_illustrative_mass':False,
                     'blocking_dependencies':['source_integrity_or_interface_error: '+str(exc)]}
            write_json(target/'failure.json',failure)
            rows.append(failure)
    summary={'format':FORMAT,'model_count':len(rows),
             'all16_requested':len(chosen)==16,
             'two_matching_thresholds_only':True,
             'fermions_integrated_together_at_muF':True,
             'both_scalars_integrated_together_at_muS':True,
             'scheme':'MSbar',
             'sourced_matchete_quartic_convention':'V4=lambda_abcd phi_a phi_b phi_c phi_d/4!',
             'full_physical_eft1_couplings_certified':False,
             'physical_C5_calculated':False,
             'local_source_recheck':recheck,
             'models':rows}
    write_json(output/'summary.json',summary)
    return summary
