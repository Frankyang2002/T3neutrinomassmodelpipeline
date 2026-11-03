"""B2 G2 stage 1: exact *conditional* SU(2) pseudoreal-intertwiner audit.

Inputs are live-captured G1 RGBeta registry evidence and already-recorded
Matchete native source-CG audits. This does NOT determine RGBeta component
axes, its physical action, an absolute signed coupling dictionary, or loop
normalisation. No source, benchmarks, or existing reports are modified.
"""
from __future__ import annotations

import hashlib
import json
import re
from itertools import product
from pathlib import Path

import sympy as sp

from b2_gauge_coupling import (
    AUDITS, AuditBlocked, FLAGS, SOURCE_BLOBS, git_result, models, native_tensor,
    require, source_rotations,
)
from b2_rgbeta_ward_witness import (
    ANTI, FUND, EPS, all16_charge_warning, fit_tensor,
    formal_evidence, literal_B_y1, literal_B_y2,
    literal_rgbeta_22_mix, repaired_22_mix,
    source_provenance, ward_residual,
)
from b2_rgbeta_internal_action import (
    ALLOWED_B2_PATHS, native_snapshot_guard, read_report, validate_native_data,
)

FORMAT = 'b2-su2-pseudoreal-conditional-intertwiner-v1'
FORMAL_DEFAULT = Path(
    'output/general_observable/physical_normalisation/'
    'b2_rgbeta_native_01/formal_invariants.json')
NATIVE_DEFAULT = Path(
    'output/general_observable/physical_normalisation/'
    'b2_rgbeta_action_native_05/internal_action_registry.json')
OUTPUT_DEFAULT = Path(
    'output/general_observable/physical_normalisation/'
    'b2_su2_intertwiner_local_01')


def exact_zero(v: sp.Expr) -> bool:
    return sp.simplify(v) == 0


def transported_mix(ix: tuple[int, int, int, int], intertwiner: sp.Matrix = EPS) -> sp.Expr:
    """Transport the formal S1 *anti*-fundamental slot into a fund slot.

    M'(h,h2,s,b) = sum_k M_literal(h,h2,k,b) J[k,s].
    This is a representation isomorphism hypothesis, NOT RGBeta's native
    numerical conversion. Its overall sign is a convention.
    """
    h, h2, s, b = ix
    return sp.simplify(sum(literal_rgbeta_22_mix((h,h2,k,b)) *
                           intertwiner[k,s] for k in range(2)))


def transported_y1(ix: tuple[int, int], intertwiner: sp.Matrix = EPS) -> sp.Expr:
    """Transport the Bar[S1] formal *fund* slot into an anti-fund slot."""
    s, l = ix
    return sp.simplify(sum(intertwiner[s,k] * literal_B_y1((k,l))
                           for k in range(2)))


def core_witness(intertwiner: sp.Matrix = EPS) -> dict:
    require(intertwiner.shape == (2,2), 'intertwiner must be 2x2')
    require(intertwiner.det() != 0, 'singular intertwiner')
    # These are mathematical checks for an ASSUMED interpretation of the
    # symbolic RGBeta index slots. No native component values enter here.
    wrong_mix = ward_residual(literal_rgbeta_22_mix, (FUND,FUND,FUND,ANTI))
    assumed_mix = ward_residual(literal_rgbeta_22_mix, (FUND,FUND,ANTI,ANTI))
    fixed_mix = ward_residual(lambda ix: transported_mix(ix,intertwiner),
                              (FUND,FUND,FUND,ANTI))
    wrong_y1 = ward_residual(literal_B_y1, (ANTI,FUND))
    assumed_y1 = ward_residual(literal_B_y1, (FUND,FUND))
    fixed_y1 = ward_residual(lambda ix: transported_y1(ix,intertwiner),
                             (ANTI,FUND))
    y2 = ward_residual(literal_B_y2, (FUND,FUND))
    require(wrong_mix == (3,0,3) and wrong_y1 == (2,0,2),
            'literal negative-control Ward norms changed')
    require(all(v == 0 for v in (*assumed_mix,*fixed_mix,*assumed_y1,*fixed_y1,*y2)),
            'assumed/transported representation not a Ward singlet')
    require(all(exact_zero(transported_mix(ix,intertwiner) - repaired_22_mix(ix))
                for ix in product(range(2),repeat=4)),
            'transported B/C mixing differs from chosen delta-epsilon candidate')
    require(all(exact_zero(transported_y1(ix,intertwiner) + int(ix[0]==ix[1]))
                for ix in product(range(2),repeat=2)),
            'transported B y1 differs from negative delta')
    require(all(exact_zero(transported_mix((h,k,s,b),intertwiner) -
                           transported_mix((k,h,s,b),intertwiner))
                for h,k,s,b in product(range(2),repeat=4)),
            'transformed mixing lost Higgs exchange symmetry')
    transpose_flipped = all(exact_zero(
        transported_mix(ix,intertwiner.T) + transported_mix(ix,intertwiner))
        for ix in product(range(2),repeat=4))
    require(transpose_flipped,
            'epsilon transpose phase negative control unexpectedly changed')
    return {
        'epsilon_chosen': [[str(z) for z in row] for row in intertwiner.tolist()],
        'epsilon_is_RGBeta_native_component_witness': False,
        'alternative_formal_index_assignment': {
            'H':FUND,'S1':ANTI,'Bar[S1]':FUND,
            'S2':FUND,'Bar[S2]':ANTI,'l':FUND,
        },
        'literal_standard_BC_mix_ward_squared':[int(z) for z in wrong_mix],
        'alternative_formal_BC_mix_ward_squared':[int(z) for z in assumed_mix],
        'transported_standard_BC_mix_ward_squared':[int(z) for z in fixed_mix],
        'literal_standard_B_y1_ward_squared':[int(z) for z in wrong_y1],
        'alternative_formal_B_y1_ward_squared':[int(z) for z in assumed_y1],
        'transported_standard_B_y1_ward_squared':[int(z) for z in fixed_y1],
        'B_y2_fund_fund_ward_squared':[int(z) for z in y2],
        'transported_mix_equals_previous_delta_epsilon_candidate': True,
        'transported_B_y1_equals_negative_delta': True,
        'transpose_epsilon_changes_overall_sign_not_ward_status': True,
        'interpretation': 'conditional pseudoreality intertwiner; actual RGBeta barred-index conventions unknown',
        'physical_coupling_conversion': None,
    }


def triplet_witness(source_u: dict[int, sp.Matrix]) -> dict:
    """Source-Matchete triplet basis metric; NOT a measured RGBeta triplet."""
    u = source_u[3]
    require(u.shape == (3,3) and sp.simplify(u.H*u) == sp.eye(3),
            'Matchete triplet rotation not unitary')
    metric = sp.simplify(u.T*u)
    ts = [sp.Matrix(3,3, lambda b,c: -sp.I*sp.LeviCivita(a,b,c))
          for a in range(3)]
    tw = [sp.simplify(u.H*t*u) for t in ts]
    require(all(t.T*metric + metric*t == sp.zeros(3) for t in tw),
            'Matchete triplet metric failed Ward identity')
    require(metric.det()!=0 and metric.T==metric and metric!=sp.eye(3),
            'triplet metric not a nontrivial symmetric weight metric')
    return {
        'Matchete_triplet_weight_metric':
            [[str(z) for z in row] for row in metric.tolist()],
        'Matchete_triplet_metric_ward_generators_verified':3,
        'source_basis_rotated_using_native_Matchete_generator_audit':True,
        'metric_equals_RGBeta_delS2_axes':False,
        'RGBeta_numeric_triplet_axis_map':None,
    }


def _json_file(path: Path) -> dict:
    require(path.is_file(),f'missing required evidence: {path}')
    try:
        data = json.loads(path.read_text(encoding='utf-8-sig'))
    except (UnicodeError, json.JSONDecodeError) as exc:
        raise AuditBlocked(f'malformed JSON {path}: {exc}') from exc
    require(isinstance(data,dict),f'evidence not an object: {path}')
    return data


def g1_capture_reuse_guard(native: dict, root: Path, provenance: dict) -> dict:
    """Allow reuse of previously captured G1 on descendant commits only.

    The native G1 JSON retains its *actual* capture HEAD. When the checkout
    has moved forward for unrelated work, verify a real Git ancestry link,
    all six pinned physics-source blobs in the capture commit's tree, the
    three capture-inventory blobs, and the current canonical worktree hashes.
    This does not recapture RGBeta or relax source_provenance(root), which
    separately checks all six current HEAD/index/clean-filter source paths.
    """
    capture = native.get('git_head_at_capture')
    current = provenance.get('repository_head_measured')
    require(isinstance(capture, str) and
            re.fullmatch(r'[0-9a-f]{40}', capture) is not None,
            'G1 report missing valid native capture HEAD')
    require(isinstance(current, str) and
            re.fullmatch(r'[0-9a-f]{40}', current) is not None and
            provenance.get('all_six_source_blobs_unchanged') is True,
            'G2 current six-source provenance not verified')
    require(set(ALLOWED_B2_PATHS).issubset(SOURCE_BLOBS),
            'G1 source inventory not included in six pinned source paths')
    captured_blobs = native.get('git_source_blobs_at_capture')
    require(isinstance(captured_blobs, dict) and
            set(captured_blobs) == set(ALLOWED_B2_PATHS),
            'G1 captured source inventory missing or unexpected')

    if capture == current:
        native_snapshot_guard(native, root, provenance)
        return {'capture_head': capture, 'current_head': current,
                'capture_is_ancestor': True, 'same_head': True,
                'all_six_source_blobs_verified_in_capture_tree': True,
                'all_six_source_blobs_verified_in_current_checkout': True,
                'g1_native_reran_at_current_head': False,
                'physical_source_certified': False}

    rc, _ = git_result(root, 'merge-base', '--is-ancestor', capture, 'HEAD')
    require(rc == 0,
            'G1 native capture HEAD is not a verified ancestor of current HEAD')
    for rel, sha in SOURCE_BLOBS.items():
        rc, record = git_result(root, 'ls-tree', capture, '--', rel)
        # Git ls-tree: '<mode> blob <SHA>\t<path>'. Check the path as well.
        parts = record.split('\t', 1)
        left = parts[0].split() if parts else []
        require(rc == 0 and len(parts) == 2 and parts[1] == rel and
                len(left) == 3 and left[1] == 'blob' and left[2] == sha,
                f'G1 capture commit has different/missing pinned source: {rel}')
        rc, live = git_result(root, 'hash-object', f'--path={rel}', '--', rel)
        require(rc == 0 and live == sha,
                f'G1 source drift in current canonical checkout: {rel}')
    for rel in ALLOWED_B2_PATHS:
        require(captured_blobs[rel] == SOURCE_BLOBS[rel],
                f'G1 native recorded different source blob: {rel}')
    return {'capture_head': capture, 'current_head': current,
            'capture_is_ancestor': True, 'same_head': False,
            'all_six_source_blobs_verified_in_capture_tree': True,
            'all_six_source_blobs_verified_in_current_checkout': True,
            'g1_native_reran_at_current_head': False,
            'physical_source_certified': False}


def evidence_guard(cg: dict, covariance: dict) -> dict[tuple[int,int,int,int],dict]:
    """Require all 16 recorded Matchete source components/covariance records."""
    require(covariance.get('status')=='cg_generator_covariance_verified' and
            covariance.get('count')==16,
            'Matchete CG covariance audit not verified for 16 cases')
    require(len(cg.get('cases',[]))==16 and
            len(covariance.get('models',[]))==16,
            'Matchete source CG/covariance inventory incomplete')
    expected = {(m['dS1'],m['dS2'],m['dF'],m['alpha']) for m in models()}
    require(len(expected)==16, 'restricted model inventory changed')
    cases = {}
    for case in cg['cases']:
        m=case['model']
        key=tuple(int(m[k]) for k in ('dS1','dS2','dF','alpha'))
        require(key in expected and key not in cases, 'unknown/duplicate Matchete CG model')
        require(set(('T3Y1CG','T3Y2CG','T3MixCG')).issubset(case['vertices']),
                'native Matchete CG vertices missing')
        cases[key]=case
    covkeys=[]
    for row in covariance['models']:
        m=row['model']
        key=tuple(int(m[k]) for k in ('dS1','dS2','dF','alpha'))
        require(row.get('all_covariant') is True and
                all(row['vertices'][v]['covariant'] is True
                    for v in ('T3Y1CG','T3Y2CG','T3MixCG')),
                f'Matchete covariance not verified for {key}')
        covkeys.append(key)
    require(set(cases)==expected and len(set(covkeys))==16 and
            set(covkeys)==expected,'Matchete CG/covariance models differ')
    return cases


def compare_b_c_mixing(cases: dict, source_u: dict[int,sp.Matrix]) -> list[dict]:
    rows=[]
    for m in models():
        if m['family'] not in ('B','C'):
            continue
        key=(m['dS1'],m['dS2'],m['dF'],m['alpha'])
        cg=native_tensor(cases[key]['vertices']['T3MixCG'],
                         FLAGS['T3MixCG'],source_u)
        require(tuple(cg.shape)==(2,2,2,2),
                f'B/C mixing source CG has incorrect dimensions: {key}')
        obs=lambda ix: sp.simplify(cg[ix])
        candidate=fit_tensor(repaired_22_mix,obs,4)
        literal=fit_tensor(literal_rgbeta_22_mix,obs,4)
        require(candidate['is_proportional_exact'] is True and
                candidate['exact_residual_squared_norm']=='0' and
                literal['is_proportional_exact'] is False,
                f'B/C native Matchete CG disagrees with conditional intertwiner: {key}')
        rows.append({
            'family':m['family'],'alpha':m['alpha'],
            'Matchete_over_transported_candidate':candidate['complex_ratio_to_chosen_candidate'],
            'Matchete_candidate_exact_residual_squared_norm':'0',
            'literal_untransported_delta_delta_proportional_to_Matchete':False,
            'actual_RGBeta_to_Matchete_coupling_ratio':None,
        })
    require(len(rows)==6,'B/C native component comparison not six models')
    return rows


def run(root: Path, native_path: Path, formal_path: Path) -> dict:
    root = root.resolve()
    provenance = source_provenance(root)
    native = read_report(native_path)
    g1_reuse = g1_capture_reuse_guard(native,root,provenance)
    validated = validate_native_data(native)
    require(validated['families_checked']==5 and
            validated['native_registry_ratios_measured']==30,
            'G1 native registry gate not accepted')
    formal = _json_file(formal_path)
    validated_formal = formal_evidence(formal)
    require(validated_formal['five_classes_formal_records_validated'] and
            validated_formal['B_C_formal_delta_delta_confirmed'],
            'prior RGBeta native formal classes not verified')
    gen_path=root/AUDITS/'matchete_generator_analysis.json'
    cg_path=root/AUDITS/'component_support_audit.json'
    cov_path=root/AUDITS/'cg_covariance_audit.json'
    gen=_json_file(gen_path)
    cases=evidence_guard(_json_file(cg_path),_json_file(cov_path))
    source_u=source_rotations(gen)
    ward=core_witness()
    triplet=triplet_witness(source_u)
    bc=compare_b_c_mixing(cases,source_u)
    charge=all16_charge_warning()
    require(len(charge)==6 and all(
        all(x['electric_charge_preserved'] for x in row['explicit_Ward_invariant_candidate'])
        for row in charge), 'B/C candidate did not conserve electric charge')
    source_paths=[native_path, formal_path, gen_path, cg_path, cov_path]
    paths_sha={p.relative_to(root).as_posix() if p.is_relative_to(root) else str(p):
               hashlib.sha256(p.read_bytes()).hexdigest() for p in source_paths}
    return {
        'format':FORMAT,
        'status':'CONDITIONAL_BC_INTERTWINER_EXACT_NATIVE_RGBETA_AXES_PENDING',
        'source_provenance':provenance,
        'g1_capture_reuse_provenance':g1_reuse,
        'evidence_sha256':paths_sha,
        'g1_validated_families':5,
        'g1_validated_internal_ratios':30,
        'RGBeta_formal_families_verified':5,
        'Matchete_native_CG_cases_accounted_for':16,
        'Matchete_native_BC_mix_comparisons':bc,
        'B_C_charge_cases':len(charge),
        'ward_and_intertwiner':ward,
        'Matchete_triplet_basis':triplet,
        'conditional_orientation_tested_for':['B','C'],
        'classes_without_RGBeta_numerical_axis_comparison':['A','D','E'],
        'RGBeta_native_component_values_measured':False,
        'RGBeta_native_S1_orientation_identified':False,
        'RGBeta_native_triplet_basis_identified':False,
        'crosspackage_physical_coupling_dictionary':None,
        'physical_source_certified_models':0,
        'safe_to_integrate_into_production':False,
        'EFT_performed':False,
        'next_gate':'Independent live-native RGBeta component/field-orientation witness or action-level SU2 generators; then G3 signed physical-action dictionary',
    }


def write_report(output_dir: Path, result: dict) -> None:
    require(not output_dir.exists(), f'refusing existing output: {output_dir}')
    require(result['format']==FORMAT and
            result['physical_source_certified_models']==0 and
            result['safe_to_integrate_into_production'] is False,
            'invalid G2 conditional output/certification')
    output_dir.mkdir(parents=True,exist_ok=False)
    (output_dir/'intertwiner_witness.json').write_text(
        json.dumps(result,sort_keys=True,indent=2)+'\n',encoding='utf-8')
    (output_dir/'HANDOFF.md').write_text(
        '# B2 G2 stage 1: conditional SU(2) intertwiner witness\n\n'
        f"- Local HEAD: `{result['source_provenance']['repository_head_measured']}`.\n"
        '- Six tracked source blobs checked with HEAD/index/clean-filter guard.\n'
        f"- G1 native capture HEAD: `{result['g1_capture_reuse_provenance']['capture_head']}`; "
        f"same HEAD: `{result['g1_capture_reuse_provenance']['same_head']}`; "
        'verified ancestor and unchanged historical six-source tree.\n'
        '- G1: 5/5 native classes and 30/30 internal ratios.\n'
        '- Exact Ward: the B/C mixing literal (fund,fund,fund,anti) has (3,0,3); '
        'the alternative (fund,fund,anti,anti) has (0,0,0).\n'
        '- Under a **chosen** epsilon intertwiner on the formal S1 index, '
        'the mixing tensor equals the delta-epsilon candidate; '
        'the B y1 formal epsilon becomes negative delta.\n'
        '- Native Matchete B/C mixing comparison: 6/6 exact proportionalities '
        'to the **conditional** candidate; all 16 Matchete source CG covariance '
        'records guarded.\n'
        '- Triplet weight metric: from Matchete generator rotation only. '
        '**Not a measured RGBeta triplet metric.**\n'
        '- This does NOT establish RGBeta S1/barred-index orientation, '
        'numerical eps or delS2 components, signed physical dictionary, or '
        'B1 spinor phase. A sign reversal of epsilon is a remaining convention.\n'
        '- Physical source certifications: **0/16**. No production change, RGE or EFT.\n',
        encoding='utf-8')
