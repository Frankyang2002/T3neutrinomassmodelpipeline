"""Workstream 3 closeout: unit-reference HHSS sensitivity versus conditional C5.

Only saved Workstream 3 results are read. Reference matrices remain unit-strength
projected diagnostics in *unverified* component/normalization conventions. No
named RGBeta coupling is assigned, no physical Wilson coefficient is constructed.

For Re/Im response matrices A and B and a hypothetical complex reference
coefficient c=x+iy, |c|<=k, the tight Frobenius response bound is

    max ||x*A+y*B||_F = k*sqrt(lambda_max(Gram_R(A,B)))

where Gram_R has entries Re Tr(A^dagger B). Summing these single-reference
bounds is a conservative bound, only *within the arbitrary unit reference
convention*; it is NOT a bound on omitted physical UV terms.
"""
from __future__ import annotations

from collections import Counter
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from .inventory import ALL_MODELS

FORMAT = 't3-ws3-unit-hhss-sensitivity-closeout-v1'
SOURCE_FORMAT = 't3-ws3-conditional-portal-completion-and-reference-response-v1'
ACCEPTED_PARTS = frozenset(('CONDITIONAL_UNIT_REFERENCE_RESPONSE', 'IDENTICALLY_ZERO_REFERENCE_TENSOR'))
DEFAULT_CAPS = (1e-3, 1e-2, 1e-1)


class CloseoutError(ValueError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise CloseoutError(message)


def read_json(path: Path) -> dict:
    require(path.is_file(), f'Missing saved Workstream 3 file: {path}')
    with path.open('r', encoding='utf-8-sig') as f:
        result = json.load(f)
    require(isinstance(result, dict), f'{path}: expected JSON object')
    return result


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def slug(model: str) -> str:
    return model.replace('-', 'm').replace('+', 'p')


def complex_matrix(data: Mapping[str, Any], what: str) -> np.ndarray:
    require(isinstance(data, dict) and set(data) >= {'real', 'imag'},
            f'{what}: missing real/imag matrix')
    re = np.asarray(data['real'], dtype=float)
    im = np.asarray(data['imag'], dtype=float)
    require(re.shape == (3, 3) and im.shape == (3, 3) and
            np.isfinite(re).all() and np.isfinite(im).all(),
            f'{what}: nonfinite or non-3x3 flavour matrix')
    mat = re + 1j*im
    n = float(np.linalg.norm(mat))
    require(np.linalg.norm(mat - mat.T) <= 1e-8*max(n, 1e-300),
            f'{what}: matrix is not flavour symmetric')
    return mat


def unit_response_norm(re: np.ndarray, im: np.ndarray) -> float:
    """Real-parameter operator norm from R^2 to complex 3x3 Frobenius space."""
    a = float(np.vdot(re, re).real)
    b = float(np.vdot(im, im).real)
    cross = float(np.vdot(re, im).real)
    eigenmax = 0.5*(a+b+math.hypot(a-b, 2*cross))
    require(math.isfinite(eigenmax) and eigenmax >= 0,
            'Invalid unit complex-reference response Gram matrix')
    return math.sqrt(eigenmax)


def _one_part(part: Any, description: str) -> tuple[str, np.ndarray | None, float | None]:
    require(isinstance(part, dict), f'{description}: missing part')
    require(part.get('physical_C5_certified') is False,
            f'{description}: part must remain physically uncertified')
    status = part.get('status')
    if status == 'IDENTICALLY_ZERO_REFERENCE_TENSOR':
        return status, np.zeros((3, 3), complex), 0.0
    if status == 'CONDITIONAL_UNIT_REFERENCE_RESPONSE':
        require(part.get('physical_C5_certified') is False,
                f'{description}: unit response unexpectedly claims physical status')
        m = complex_matrix(part.get('dim5_matrix_projected_diagnostic_only'), description)
        residuals = part.get('projection_residuals')
        require(isinstance(residuals, dict) and bool(residuals),
                f'{description}: missing Ward residuals')
        vals = [float(v) for v in residuals.values()]
        require(all(math.isfinite(v) and v >= 0 for v in vals),
                f'{description}: invalid Ward residuals')
        require(max(vals) <= 1e-8,
                f'{description}: saved status passes but Ward residual fails')
        return status, m, max(vals)
    if status in ('WARD_REJECTED_DO_NOT_USE_AS_C5', 'REFERENCE_NUMERICAL_BLOCKED'):
        return status, None, None
    raise CloseoutError(f'{description}: unknown/missing status {status!r}')


def source_gate(doc: Mapping[str, Any]) -> None:
    require(doc.get('format') == SOURCE_FORMAT, 'Unexpected portal completion run format')
    require(doc.get('model_count') == 16 and doc.get('validated_input_models') == 16 and
            doc.get('blocked_source_models') == 0, 'Portal completion requires 16 validated models')
    require(doc.get('baseline_directions_verified') == 14 and
            doc.get('extra_named_symbols_verified') == 24 and
            doc.get('new_reference_unit_response_count') == 24,
            'Expected 14 baseline and 24 computed extra reference directions')
    require(doc.get('conditional_baseline_ward_passed') == 16 and
            doc.get('conditional_baseline_ward_rejected') == 0,
            'Completed baseline Ward checks are not all passing')
    require(doc.get('physical_C5_certified') is False and
            doc.get('complete_one_loop_C5_certified') is False and
            doc.get('physical_named_HHSS_coupling_maps_certified') == 0 and
            doc.get('extra_reference_coefficients_assigned') is False,
            'Cannot closeout a source that asserts unverified physical normalization')
    require(doc.get('source_GammaCC_to_Weyl_signed_factor') is None,
            'Unknown GammaCC normalization unexpectedly supplied')
    required_sequence = ('UV -> ALL BSM fermions @ ONE muF -> EFT1(H,S1,S2,SM fermions)'
                         ' -> BOTH scalars @ ONE muS -> SMEFT')
    require(doc.get('threshold_sequence') == required_sequence,
            'Threshold sequence differs from UV -> all F -> both S -> SMEFT')
    models = doc.get('models')
    require(isinstance(models, list) and len(models) == 16 and
            Counter(x.get('model') for x in models) == Counter(ALL_MODELS),
            'Incomplete or duplicated all-16 model inventory')
    require(sum(row.get('baseline_completion',{}).get('added_baseline_directions',-100)
                for row in models) == 14,
            'Per-model completed baseline directions do not sum to 14')
    for row in models:
        require(row.get('status') == 'COMPLETED_BASELINE_CONDITIONAL_C5' and
                row.get('portal_ward_projection_passed') is True,
                f'{row.get("model")}: baseline not complete or Ward rejected')
        require(row.get('physical_C5_certified') is False and
                row.get('conditional_sensitivity_only') is True,
                f'{row.get("model")}: invalid conditional source flags')


def inspect_model(model: str, meta: Mapping[str, Any], baseline: dict,
                  extras: dict, caps: tuple[float, ...]) -> dict:
    """Inspect one saved model, not a recomputation of an unknown UV coupling."""
    require(baseline.get('model') == model and
            baseline.get('status') == 'CONDITIONAL_COMPLETED_BASELINE_PORTAL' and
            baseline.get('physical_C5_certified') is False and
            baseline.get('source_normalization_certified') is False and
            baseline.get('all_independent_HHSS_included') is False,
            f'{model}: invalid baseline source/status flags')
    data = baseline.get('response', {})
    require(data.get('projection_passed') is True and
            data.get('physical_C5_certified') is False,
            f'{model}: baseline source Ward or physical status invalid')
    base = complex_matrix(data.get('matrices', {}).get('dim5_sum'), f'{model}: baseline dim5')
    base_norm = float(np.linalg.norm(base))
    require(base_norm > 0 and math.isfinite(base_norm), f'{model}: zero/invalid baseline C5 norm')
    saved_norm = float(data.get('dim5_norm', float('nan')))
    require(math.isfinite(saved_norm) and abs(base_norm-saved_norm) <=
            1e-10*base_norm, f'{model}: baseline matrix norm inconsistent with saved norm')
    require(extras.get('model') == model and
            extras.get('status') == 'REFERENCE_UNIT_SENSITIVITIES_NOT_PHYSICAL' and
            extras.get('named_Matchete_RGBeta_maps_certified') == 0 and
            extras.get('physical_coefficients_inserted') is False and
            extras.get('physical_C5_certified') is False,
            f'{model}: extra unit source falsely certified or wrong model')
    rows = extras.get('linear_response_basis')
    require(isinstance(rows, list) and
            len(rows) == meta.get('unit_reference_responses_computed'),
            f'{model}: unit reference response count differs from release summary')
    symbols = set()
    reference_keys = set()
    details = []
    valid_opnorms = []
    unavailable = 0
    ward_rejected = 0
    for ref in rows:
        refid = ref.get('reference_id')
        sector = ref.get('sector_counts')
        require(isinstance(refid, str) and refid and isinstance(sector, list) and
                len(sector) == 6 and all(type(x) is int and x >= 0 for x in sector),
                f'{model}: invalid reference identifier or six-field sector')
        key = (tuple(sector), refid)
        require(key not in reference_keys, f'{model}: repeated sector/reference identifier')
        reference_keys.add(key)
        names = ref.get('unassigned_source_symbols')
        require(isinstance(names, list) and names and
                all(isinstance(n, str) and n for n in names),
                f'{model}: reference missing named source evidence')
        symbols.update(names)
        require(ref.get('physical_coefficient_assigned') is False and
                ref.get('reference_basis_to_named_RGBeta_map_certified') is False,
                f'{model}: reference coefficient unexpectedly assigned')
        parts = ref.get('unit_parts')
        require(isinstance(parts, dict) and set(parts) == {'Re','Im'},
                f'{model}: reference must have both Re and Im unit responses')
        result = {}
        matrices = {}
        residuals = {}
        for partname in ('Re','Im'):
            status, matrix, residual = _one_part(parts[partname],
                                                 f'{model} {key} {partname}')
            result[partname] = status
            matrices[partname] = matrix
            residuals[partname] = residual
        valid = all(v in ACCEPTED_PARTS for v in result.values())
        if valid:
            m_re, m_im = matrices['Re'], matrices['Im']
            sigma = unit_response_norm(m_re, m_im)
            re_norm, im_norm = float(np.linalg.norm(m_re)), float(np.linalg.norm(m_im))
            valid_opnorms.append(sigma)
            sensitivity = sigma / base_norm
            status = ('ZERO_UNIT_RESPONSE' if sigma == 0 else
                      'VALID_CONDITIONAL_UNIT_SENSITIVITY')
        else:
            sigma = re_norm = im_norm = sensitivity = None
            status = ('WARD_REJECTED_NO_C5' if
                      'WARD_REJECTED_DO_NOT_USE_AS_C5' in result.values()
                      else 'BLOCKED_REFERENCE_RESPONSE')
            unavailable += 1
            if status == 'WARD_REJECTED_NO_C5': ward_rejected += 1
        details.append({
            'model':model,'sector_counts':sector,'reference_id':refid,
            'unassigned_source_symbols':names,'status':status,
            're_status':result['Re'],'im_status':result['Im'],
            're_dim5_frobenius_per_unit':re_norm,
            'im_dim5_frobenius_per_unit':im_norm,
            'max_dim5_frobenius_per_unit_complex_coefficient':sigma,
            'max_relative_to_baseline_per_unit_coefficient':sensitivity,
            'single_reference_coefficient_for_1pct_baseline':(0.01/sensitivity if sensitivity else None),
            'single_reference_coefficient_for_10pct_baseline':(0.1/sensitivity if sensitivity else None),
            'single_reference_coefficient_for_100pct_baseline':(1.0/sensitivity if sensitivity else None),
            'max_saved_projection_residual':max((x for x in residuals.values() if x is not None),default=None),
            'physical_named_coupling_assigned':False,
            'reference_to_named_orientation_certified':False,
            'same_unknown_GammaCC_factor_as_baseline_assumed_for_ratios':True,
        })
    expected_named = meta.get('extra_named_inventory',{}).get('extra_named_symbols')
    require(type(expected_named) is int and len(symbols) == expected_named,
            f'{model}: extra named symbol count mismatch')
    require(len(rows) == meta.get('extra_reference_basis_count'),
            f'{model}: source reference basis count differs from the saved model record')
    original_inventory=meta.get('extra_named_inventory',{}).get('source_named_symbol_evidence')
    require(isinstance(original_inventory,list) and
            {x.get('rgbeta_symbol') for x in original_inventory} == symbols,
            f'{model}: named source symbols differ from saved V4 inventory')
    upper_sum = float(sum(valid_opnorms))/base_norm
    complete = unavailable == 0
    envelope = {f'{cap:g}':(cap*upper_sum if complete else None) for cap in caps}
    return {
        'model':model,
        'status':('NO_EXTRA_INVARIANT_DIRECTIONS' if not rows else
                  'ALL_REFERENCE_SENSITIVITIES_COMPUTED' if complete else
                  'REFERENCE_SENSITIVITIES_PARTIALLY_BLOCKED'),
        'baseline_dim5_frobenius_per_unknown_N_gamma_gev_inverse':base_norm,
        'baseline_ward_passed':True,
        'reference_directions':len(rows),
        'named_source_symbols':len(symbols),
        'valid_reference_directions':len(rows)-unavailable,
        'ward_rejected_reference_directions':ward_rejected,
        'other_blocked_reference_directions':unavailable-ward_rejected,
        'max_single_direction_relative_response_per_unit':(
            max((x['max_relative_to_baseline_per_unit_coefficient'] for x in details
                 if x['max_relative_to_baseline_per_unit_coefficient'] is not None),default=0.0)),
        'sum_single_direction_relative_response_per_unit':upper_sum,
        'illustrative_max_relative_frobenius_shift_at_caps':envelope,
        'all_reference_directions_usable_as_conditional_C5':complete,
        'hypothetical_single_reference_coefficients_only_not_physical':True,
        'physical_C5_certified':False,
        'references':details,
    }


def inspect_run(directory: Path, caps: tuple[float, ...] = DEFAULT_CAPS) -> dict:
    require(len(caps) >= 1 and len(caps) <= 8 and
            all(isinstance(x, (float,int)) and math.isfinite(x) and 0 < x <= 1 for x in caps) and
            len(set(caps)) == len(caps), 'Scenario caps must be distinct positive values <=1')
    root = Path(directory)
    source_path = root / 'summary.json'
    source = read_json(source_path)
    source_gate(source)
    summaries=[]
    paths={}
    for meta in source['models']:
        model=meta['model']
        code=slug(model)
        path1 = root / 'conditional_portal_response' / (code+'.json')
        path2 = root / 'extra_reference_unit_responses' / (code+'.json')
        baseline = read_json(path1)
        extras = read_json(path2)
        result = inspect_model(model,meta,baseline,extras,caps)
        result['input_sha256']={'baseline_response':sha256(path1),'unit_reference_response':sha256(path2)}
        summaries.append(result)
        paths[model]=[str(path1),str(path2)]
    require(sum(x['reference_directions'] for x in summaries) == 24,
            'Reference direction count changed from V3/V4 24-entry release')
    require(sum(x['named_source_symbols'] for x in summaries) == 24,
            'Named source symbol inventory does not total 24')
    num_valid = sum(x['valid_reference_directions'] for x in summaries)
    num_bad = sum(x['ward_rejected_reference_directions'] for x in summaries)
    num_other = sum(x['other_blocked_reference_directions'] for x in summaries)
    overall = ('ALL_16_CONDITIONAL_BASELINES_AND_24_UNIT_RESPONSES_AUDITED'
               if num_valid == 24 else
               'ALL_16_BASELINES_AUDITED_SOME_UNIT_RESPONSES_BLOCKED')
    ranks = [dict(ref) for model in summaries for ref in model['references']
             if ref['max_relative_to_baseline_per_unit_coefficient'] is not None]
    ranks.sort(key=lambda x:x['max_relative_to_baseline_per_unit_coefficient'],reverse=True)
    return {
        'format':FORMAT,'status':overall,'model_count':16,
        'baseline_completed_models':16,'baseline_added_portal_directions':14,
        'extra_named_symbols_in_inventory':24,'extra_reference_directions':24,
        'valid_unit_reference_directions':num_valid,
        'ward_rejected_reference_directions':num_bad,
        'other_blocked_reference_directions':num_other,
        'models_with_extra_references':sum(x['reference_directions']>0 for x in summaries),
        'models_with_nonzero_valid_sensitivities':sum(any(
            r['max_relative_to_baseline_per_unit_coefficient'] not in (None,0)
            for r in x['references']) for x in summaries),
        'illustrative_unit_coefficient_caps':list(caps),
        'reference_unit_is_declared_V3_potential_coefficient_not_named_RGBeta_value':True,
        'scenario_bounds_are_conditional_not_physical_predictions':True,
        'unknown_shared_spinor_factor_cancels_only_in_consistently_defined_ratios':True,
        'source_run_summary_sha256':sha256(source_path),
        'source_run_ws2_releases':source.get('ws2_releases'),
        'physical_C5_certified':False,
        'physical_full_one_loop_certified':False,
        'exploratory_workstream3_closeout_report_produced':True,
        'physical_extra_named_coupling_maps_certified':0,
        'missing_physics_not_resolved':[
            'GammaCC-to-Weyl signed normalization',
            'Matchete/RGBeta reference coefficient normalization, phases, orientation',
            'special scalar original-to-boundary physical basis maps',
            '24 named extra HHSS coefficients and reliable full EFT1 evolution',
            'full one-loop hard finite/counterterm/operator closure',
            'earlier independent scalar projection failures in A_-2, D_+0, E_-2, E_+0',
        ],
        'eft_sequence':'UV -> integrate out ALL heavy fermions @muF -> EFT1(H,S1,S2) -> integrate out BOTH heavy scalars @muS -> SMEFT',
        'reference_sensitivity_ranking':[
            {'model':x['model'],'reference_id':x['reference_id'],
             'sector_counts':x['sector_counts'],
             'unassigned_source_symbols':x['unassigned_source_symbols'],
             'max_relative_to_baseline_per_unit_coefficient':x['max_relative_to_baseline_per_unit_coefficient'],
             'single_reference_coefficient_for_10pct_baseline':x['single_reference_coefficient_for_10pct_baseline']}
            for x in ranks],
        'models':summaries,
    }
