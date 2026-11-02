"""Batch 2: read-only consolidation of v6 native-Ward conditional T3 C5 outputs.

No new matching, running, couplings, normalization, or phenomenology is done.
Only the existing 16-model LLSS-induced leading log + one joint scalar bubble
is consolidated. Absolute GammaCC-to-Weyl normalization is UNKNOWN.

UV -> ALL heavy BSM fermions at muF -> EFT1(S1,S2)
   -> BOTH heavy BSM scalar multiplets together at muS -> SMEFT.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import sys
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SAFE = ROOT / 'output' / 'general_observable' / 'lo_joint_scalar_matching'
DEFAULT_INPUT = SAFE / 'exploratory_all16_v6_native_batch1_local'
DEFAULT_OUTPUT = SAFE / 'conditional_all16_batch2_v1'
MODELS = ('A_-4', 'A_-2', 'A_+0', 'B_-3', 'B_-1', 'B_+1',
          'C_-3', 'C_-1', 'C_+1', 'D_-2', 'D_+0', 'D_+2',
          'E_-4', 'E_-2', 'E_+0', 'E_+2')
V6_FORMAT = 't3-exploratory-all16-conditional-c5-benchmark-v6'
V6_RESULT_FORMAT = 't3-conditional-c5-per-unknown-gammacc-mixed-split-v6'
V6_INNER_FORMAT = 't3-joint-scalar-integrated-conditional-batch-v2'
SUCCESS = 'conditional_C5_per_unknown_NGamma_exploratory'
LABEL = 'EXPLORATORY / CONDITIONAL — NOT A COMPLETE PHYSICAL ONE-LOOP NEUTRINO-MASS PREDICTION'
VEV_GEV = 246.22
EV_PER_GEV = 1.0e9
VEV_TO_EV = VEV_GEV**2 * EV_PER_GEV / 2.0
C5_KEY = 'C5_LLSS_total_one_loop_per_N_gev_inverse'
LOG_KEY = 'C5_EFT1_induced_leading_log_per_N_gev_inverse'
BUBBLE_KEY = 'C5_joint_scalar_finite_bubble_per_N_gev_inverse'
CSV_COLS = (
    'model', 'status', 'scalar_case', 'blocker', 'mu_f_gev', 'mu_s_gev',
    'mS1Sq_input_gev2', 'mS2Sq_input_gev2', 'mu2_real_gev2', 'mu2_imag_gev2',
    'scalar_eigenstate_masses_gev_JSON',
    'C5_per_N_real_GeV_inv_JSON', 'C5_per_N_imag_GeV_inv_JSON',
    'mnu_per_N_real_eV_JSON', 'mnu_per_N_imag_eV_JSON',
    'mass_singular_1_eV_per_abs_N', 'mass_singular_2_eV_per_abs_N',
    'mass_singular_3_eV_per_abs_N', 'C5_frobenius_GeV_inv_per_N',
    'weinberg_projection_residual', 'matching_scale_cancellation_residual',
    'full_SU2_source_beta_parity', 'running_snapshot_sha256',
    'quartic_seed_sha256', 'fermion_boundary_sha256', 'result_sha256',
    'physical_C5_certified'
)


def require(ok: Any, message: str) -> None:
    if not ok:
        raise ValueError(message)


def read_json(path: Path) -> dict:
    obj = json.loads(path.read_text(encoding='utf-8-sig'))
    require(isinstance(obj, dict), f'expected JSON object: {path}')
    return obj


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def slug(model: str) -> str:
    require(model in MODELS, f'unknown model: {model}')
    return model.replace('-', 'm').replace('+', 'p')


def matrix_from_json(obj: Any, field: str) -> np.ndarray:
    require(isinstance(obj, dict) and 'real' in obj and 'imag' in obj,
            f'{field}: complex real/imag entries missing')
    r, i = np.asarray(obj['real'], dtype=float), np.asarray(obj['imag'], dtype=float)
    require(r.shape == (3, 3) and i.shape == (3, 3), f'{field}: expected 3x3')
    z = r + 1j*i
    require(np.isfinite(z).all(), f'{field}: nonfinite entries')
    return z


def matrix_to_json(z: np.ndarray) -> dict:
    require(z.shape == (3, 3) and np.isfinite(z).all(), 'invalid complex 3x3')
    return {'real': z.real.tolist(), 'imag': z.imag.tolist()}


def close_small(a: np.ndarray, b: np.ndarray, rtol: float = 2e-10,
                atol: float = 1e-27) -> bool:
    # Avoid numpy.allclose's default atol=1e-8: C5 is O(1e-15) GeV^-1.
    return bool(np.max(np.abs(a - b)) <=
                atol + rtol * max(float(np.max(np.abs(a))),
                                  float(np.max(np.abs(b))), 1e-30))


def close_number(a: Any, b: Any, label: str, rtol: float = 2e-9,
                 atol: float = 1e-26) -> None:
    a, b = float(a), float(b)
    require(math.isfinite(a) and math.isfinite(b), f'{label}: nonfinite metric')
    require(abs(a-b) <= atol + rtol*max(abs(a), abs(b)),
            f'{label}: stored metric differs from recomputation')


def check_v6_input(summary: dict) -> dict[str, dict]:
    require(summary.get('format') == V6_FORMAT, 'wrong v6 input summary format')
    require(summary.get('models_checked') == 16 and
            isinstance(summary.get('models'), list) and
            len(summary['models']) == len(MODELS),
            'Batch 2 requires 16 individually reported v6 models')
    models = summary['models']
    require(set(row.get('model') for row in models) == set(MODELS),
            'missing, unexpected, or duplicate model labels')
    require(summary.get('allow_exploratory_map') is True and
            summary.get('allow_exploratory_special_bases') is True,
            'v6 conditional calculation opt-ins not recorded')
    require(summary.get('physical_C5_certified') is False and
            summary.get('full_one_loop_neutrino_mass_calculation') is False and
            summary.get('absolute_GammaCC_to_Weyl_normalization') is None,
            'input makes incompatible physical certification claim')
    require(summary.get('all_scalars_integrated_at_one_muS') is True,
            'input does not preserve joint scalar threshold')
    n = sum(row.get('status') == SUCCESS for row in models)
    e = sum(str(row.get('status', '')).startswith('ERROR_') for row in models)
    require(summary.get('models_with_conditional_C5') == n and
            summary.get('errors') == e,
            'v6 input summary counts do not match model statuses')
    return {row['model']: row for row in models}


def evaluate_valid_model(model: str, row: dict, input_dir: Path) -> dict:
    """Recompute numerical diagnostics only, without recomputing physics tensors."""
    name = slug(model)
    result_file = input_dir / 'C5_per_unknown_NGamma_exploratory' / (name + '.json')
    background_file = input_dir / 'numeric_backgrounds' / (name + '.json')
    background_npz = input_dir / 'numeric_backgrounds' / (name + '.npz')
    for path in (result_file, background_file, background_npz):
        require(path.is_file(), f'{model}: required v6 input absent: {path}')
    obj = read_json(result_file)
    bg = read_json(background_file)
    inner = obj.get('response')
    require(isinstance(inner, dict), f'{model}: v6 response missing')
    require(obj.get('model') == inner.get('model') == row.get('model') == model,
            f'{model}: model provenance mismatch')
    require(obj.get('format') == V6_RESULT_FORMAT and
            inner.get('format') == V6_INNER_FORMAT,
            f'{model}: unexpected response format')
    require(obj.get('physical_C5_certified') is False and
            obj.get('complete_one_loop_C5_certified') is False and
            inner.get('physical_C5_certified') is False,
            f'{model}: unexpected physical-C5 certification')
    require(obj.get('native_quartic_alias_mapping_certified') is False and
            obj.get('muF_muS_mass_eigenfield_alignment_certified') is False,
            f'{model}: changed exploratory provenance flags')
    require(bg.get('model') == model and bg.get('no_separate_scalar_threshold') is True
            and bg.get('same_muS_quartic_and_mass') is True,
            f'{model}: background has wrong model or scalar matching sequence')
    require(bg.get('physical_C5_certified') is False,
            f'{model}: background erroneously certified')
    require(inner.get('model_source_boundary_sha256') ==
            row.get('fermion_boundary_sha256') == bg.get('source_boundary_sha256'),
            f'{model}: fermion boundary SHA mismatch')
    require(inner.get('background_json_sha256') == sha256_file(background_file) and
            inner.get('background_npz_sha256') == sha256_file(background_npz),
            f'{model}: local background files differ from matching inputs')
    require(Path(bg.get('tensor_npz', '')).name == name + '.npz',
            f'{model}: local background NPZ name inconsistent')
    mf, ms = float(inner['mu_f_gev']), float(inner['mu_s_gev'])
    require(abs(mf-100000.0) < 1e-6 and abs(ms-1000.0) < 1e-7 and
            abs(mf-float(row['mu_f_gev'])) < 1e-6 and
            abs(ms-float(row['mu_s_gev'])) < 1e-7 and
            abs(ms-float(bg['mu_s_gev'])) < 1e-7,
            f'{model}: incompatible matching scales')
    c5 = matrix_from_json(inner[C5_KEY], C5_KEY)
    log = matrix_from_json(inner[LOG_KEY], LOG_KEY)
    bubble = matrix_from_json(inner[BUBBLE_KEY], BUBBLE_KEY)
    require(close_small(c5, log + bubble), f'{model}: C5 != log + scalar bubble')
    require(close_small(c5, c5.T), f'{model}: C5 is not symmetric')
    metric = row.get('numerical_metrics', {})
    require(isinstance(metric, dict), f'{model}: v6 numerical metrics missing')
    require(obj.get('numerical_metrics') == metric,
            f'{model}: result/summary metric provenance differs')
    projection = max(float(inner['weinberg_beta_projection_residual']),
                     float(inner['weinberg_bubble_projection_residual']))
    cancellation = float(inner['relative_matching_scale_cancellation_residual'])
    require(0 <= projection <= 1e-8 and 0 <= cancellation <= 2e-9,
            f'{model}: Ward or scale cancellation validation failed')
    close_number(metric['projection_residual'], projection, 'Ward residual', atol=1e-18)
    close_number(metric['matching_scale_cancellation_residual'], cancellation,
                 'matching-scale residual', atol=1e-18)
    source_parity = row.get('full_su2_beta_source_parity', {})
    require(source_parity.get('full_SU2_unordered_pair_beta_parity') is True,
            f'{model}: full SU2 source parity missing or false')
    require(row.get('v2_pivot_verified') is True and
            row.get('v3_full_SU2_verified') is True,
            f'{model}: upstream source parity checks missing')
    require(source_parity == obj.get('full_SU2_beta_source_parity'),
            f'{model}: full SU2 source parity records disagree')
    m_eigen = np.asarray(inner['scalar_eigenstate_masses_gev'], dtype=float)
    require(m_eigen.ndim == 1 and len(m_eigen) >= 2 and
            np.isfinite(m_eigen).all() and (m_eigen > 0).all(),
            f'{model}: invalid scalar masses')
    if row.get('scalar_case') != 'ordinary':
        bgrow = row.get('background', {})
        require(bgrow.get('muS_mass_in_fixed_muF_source_basis') is True,
                f'{model}: special scalar basis not recorded as fixed native basis')
    for item in ('mS1Sq', 'mS2Sq'):
        require(float(row['scalar_mass_parameters_gev2'][item]) > 0,
                f'{model}: missing scalar mass input {item}')
    singular = np.linalg.svd(c5, compute_uv=False)
    masses_eV = singular*VEV_TO_EV
    for real, expected in zip(singular,
            metric['conditional_C5_singular_values_per_N_gev_inverse']):
        close_number(real, expected, 'C5 singular value')
    for real, expected in zip(masses_eV,
            metric['conditional_neutrino_mass_singular_values_eV_per_abs_N']):
        close_number(real, expected, 'mass singular value (eV)')
    close_number(np.linalg.norm(c5), metric['conditional_C5_frobenius_per_N_gev_inverse'],
                 'C5 Frobenius norm')
    close_number(np.linalg.norm(log), metric['log_C5_frobenius_per_N_gev_inverse'],
                 'log C5 norm')
    close_number(np.linalg.norm(bubble), metric['finite_bubble_C5_frobenius_per_N_gev_inverse'],
                 'bubble C5 norm')
    mu2 = row.get('background', {}).get('source_mu2')
    if row.get('scalar_case') == 'ordinary':
        require(mu2 is None, f'{model}: unexplained extra quadratic input')
    else:
        require(isinstance(mu2, dict) and 'real' in mu2 and 'imag' in mu2,
                f'{model}: special scalar quadratic input missing')
    c5_json = matrix_to_json(c5)
    mnu_json = matrix_to_json(c5*VEV_TO_EV)
    return {
        'model': model, 'status': SUCCESS,
        'scalar_case': row['scalar_case'],
        'conditional_C5_per_N_GeV_inv': c5_json,
        'conditional_mnu_per_N_eV': mnu_json,
        'mass_singular_values_eV_per_abs_N': masses_eV.tolist(),
        'scalar_mass_inputs_gev2': row['scalar_mass_parameters_gev2'],
        'scalar_quadratic_mu2_gev2': mu2,
        'scalar_eigenstate_masses_gev': m_eigen.tolist(),
        'mu_f_gev': mf, 'mu_s_gev': ms,
        'diagnostics': {
            'weinberg_projection_residual': projection,
            'matching_scale_cancellation_residual': cancellation,
            'C5_frobenius_GeV_inv_per_N': float(np.linalg.norm(c5)),
            'leading_log_frobenius_GeV_inv_per_N': float(np.linalg.norm(log)),
            'finite_bubble_frobenius_GeV_inv_per_N': float(np.linalg.norm(bubble)),
            'SU2_source_beta_parity': True,
            'C5_symmetry_max_abs_GeV_inv': float(np.max(abs(c5-c5.T))),
            'C5_log_plus_bubble_max_abs_GeV_inv': float(np.max(abs(c5-log-bubble))),
            'all_entries_finite': True
        },
        'provenance': {
            'original_v6_result_sha256': sha256_file(result_file),
            'source_background_json_sha256': sha256_file(background_file),
            'source_background_npz_sha256': sha256_file(background_npz),
            'running_snapshot_sha256': row['running_snapshot_sha256'],
            'source_quartic_seed_sha256': row['source_quartic_seed_sha256'],
            'fermion_boundary_sha256': row['fermion_boundary_sha256'],
            'candidate_quartic_aliases_NOT_CERTIFIED': row['candidate_coupling_map'],
            'native_special_basis_map_certified': False,
            'mass_input_source': row['mass_source'],
            'source_eft1_running_illustrative': True,
        },
        'absolute_GammaCC_to_Weyl_factor': None,
        'physical_C5_certified': False,
        'complete_one_loop_prediction': False
    }


def csv_row(result: dict) -> dict:
    row = {col: '' for col in CSV_COLS}
    row.update({'model': result['model'], 'status': result['status'],
                'scalar_case': result.get('scalar_case', ''),
                'blocker': result.get('blocker', ''),
                'physical_C5_certified': 'False'})
    if result['status'] != SUCCESS:
        return row
    mass = result['scalar_mass_inputs_gev2']
    mu = result['scalar_quadratic_mu2_gev2']
    diag = result['diagnostics']
    c5 = result['conditional_C5_per_N_GeV_inv']
    mnu = result['conditional_mnu_per_N_eV']
    prov = result['provenance']
    def as_cell(v: Any) -> str:
        return json.dumps(v, separators=(',', ':'), allow_nan=False)
    row.update({
        'mu_f_gev': result['mu_f_gev'], 'mu_s_gev': result['mu_s_gev'],
        'mS1Sq_input_gev2': mass['mS1Sq'], 'mS2Sq_input_gev2': mass['mS2Sq'],
        'mu2_real_gev2': '' if mu is None else mu['real'],
        'mu2_imag_gev2': '' if mu is None else mu['imag'],
        'scalar_eigenstate_masses_gev_JSON': as_cell(result['scalar_eigenstate_masses_gev']),
        'C5_per_N_real_GeV_inv_JSON': as_cell(c5['real']),
        'C5_per_N_imag_GeV_inv_JSON': as_cell(c5['imag']),
        'mnu_per_N_real_eV_JSON': as_cell(mnu['real']),
        'mnu_per_N_imag_eV_JSON': as_cell(mnu['imag']),
        'C5_frobenius_GeV_inv_per_N': diag['C5_frobenius_GeV_inv_per_N'],
        'weinberg_projection_residual': diag['weinberg_projection_residual'],
        'matching_scale_cancellation_residual': diag['matching_scale_cancellation_residual'],
        'full_SU2_source_beta_parity': 'True',
        'running_snapshot_sha256': prov['running_snapshot_sha256'],
        'quartic_seed_sha256': prov['source_quartic_seed_sha256'],
        'fermion_boundary_sha256': prov['fermion_boundary_sha256'],
        'result_sha256': prov['original_v6_result_sha256']
    })
    for j, value in enumerate(result['mass_singular_values_eV_per_abs_N'], 1):
        row[f'mass_singular_{j}_eV_per_abs_N'] = value
    return row


def consolidate(input_dir: Path) -> dict:
    summary_path = input_dir / 'summary.json'
    summary = read_json(summary_path)
    by_model = check_v6_input(summary)
    models = []
    for model in MODELS:
        row = by_model[model]
        if row.get('status') != SUCCESS:
            models.append({'model': model, 'status': row.get('status'),
                           'scalar_case': row.get('scalar_case', ''),
                           'blocker': row.get('error', 'v6 model produced no conditional C5'),
                           'physical_C5_certified': False})
        else:
            try:
                models.append(evaluate_valid_model(model, row, input_dir))
            except (OSError, ValueError, KeyError, TypeError, IndexError) as exc:
                models.append({'model': model, 'status': 'ERROR_batch2_input_validation',
                               'scalar_case': row.get('scalar_case', ''),
                               'blocker': f'{type(exc).__name__}: {exc}',
                               'physical_C5_certified': False})
    successes = [m for m in models if m['status'] == SUCCESS]
    by_case = {}
    for category in ('ordinary', 'mixed', 'self_s1', 'self_s2'):
        selection = [m for m in successes if m['scalar_case'] == category]
        largest = [m['mass_singular_values_eV_per_abs_N'][0] for m in selection]
        by_case[category] = {
            'valid_models': len(selection),
            'largest_mass_singular_value_range_eV_per_abs_N':
                [min(largest), max(largest)] if largest else None,
        }
    return {
        'format': 't3-conditional-all16-consolidated-batch2-v1',
        'label': LABEL, 'status': 'All16ConditionalC5Exploratory' if len(successes)==16
                              else 'PartialConditionalC5Exploratory',
        'input_source': str(summary_path),
        'input_summary_sha256': sha256_file(summary_path),
        'models_checked': 16, 'models_with_conditional_C5': len(successes),
        'models_blocked': 16-len(successes),
        'mu_f_gev': 1e5, 'mu_s_gev': 1e3,
        'matching_sequence': 'UV -> ALL BSM fermions -> EFT1(S1,S2) -> BOTH BSM scalars -> SMEFT',
        'matrix_convention': 'mnu_per_N_eV = (246.22 GeV)^2 C5_per_N_GeV_inv /2 *1e9 eV/GeV',
        'singular_value_convention': 'SVD/Takagi magnitudes of complex symmetric mnu per |N_GammaCC|',
        'physical_C5_certified': False,
        'absolute_GammaCC_to_Weyl_factor': None,
        'complete_one_loop_neutrino_mass_prediction': False,
        'physics_included': [
            'existing source-resolved LLSS conditional fermion-threshold operator',
            'existing fixed-order EFT1 quartic-induced leading log',
            'existing simultaneous finite joint scalar bubble at muS'
        ],
        'physics_omitted_or_uncertified': [
            'signed absolute GammaCC-to-Weyl source normalization',
            'Matchete invariant-to-RGBeta quartic numerical normalization/alias certification',
            'independent muF/muS scalar eigenfield identification and coupling normalization',
            'hard one-loop heavy-fermion matching and EFT1 soft/counterterm subtraction',
            'complete EFT1 operator mixing and full one-loop physical C5',
            'physical fit and external phenomenological constraints',
            'electroweak broken-phase resummed scalar mass shifts'
        ],
        'numerical_diagnostics_not_physical_uncertainty':
            'Ward residual, finite entries, matrix symmetry, matching-scale cancellation',
        'scalar_case_comparison': by_case,
        'models': models
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir', type=Path, default=DEFAULT_INPUT,
                        help='existing corrected v6 Batch 1 release (read only)')
    parser.add_argument('--output-dir', type=Path, default=DEFAULT_OUTPUT,
                        help='fresh Batch 2 release under lo_joint_scalar_matching')
    args = parser.parse_args(argv)
    input_dir, output_dir, safe = (p.resolve() for p in
                                  (args.input_dir, args.output_dir, SAFE))
    if not input_dir.is_dir():
        parser.error(f'input directory not found: {input_dir}')
    if output_dir == safe or safe not in output_dir.parents:
        parser.error(f'output must be in a new directory below {safe}')
    if output_dir == input_dir or output_dir in input_dir.parents or input_dir in output_dir.parents:
        parser.error('input and output directories cannot overlap')
    if output_dir.exists() and any(output_dir.iterdir()):
        parser.error(f'output directory is non-empty; refusing overwrite: {output_dir}')
    try:
        document = consolidate(input_dir)
    except (OSError, ValueError, KeyError, TypeError, IndexError) as exc:
        parser.error(f'cannot consolidate v6 input: {type(exc).__name__}: {exc}')
    output_dir.mkdir(parents=True, exist_ok=True)
    json_path = output_dir / 'summary.json'
    csv_path = output_dir / 'all16_conditional_mass_summary.csv'
    json_path.write_text(json.dumps(document, indent=2, allow_nan=False)+'\n',
                         encoding='utf-8')
    with csv_path.open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_COLS, extrasaction='raise')
        writer.writeheader()
        writer.writerows(csv_row(m) for m in document['models'])
    for record in document['models']:
        print(record['model'] + ': ' + record['status'])
        if record['status'] != SUCCESS:
            print('  BLOCKER: ' + record.get('blocker', 'unspecified'))
    print(LABEL)
    print(f"Batch 2: {document['models_with_conditional_C5']}/16 conditional C5; "
          f"{document['models_blocked']} blocked")
    print('JSON:', json_path)
    print('CSV: ', csv_path)
    return 0 if document['models_blocked'] == 0 else 2

if __name__ == '__main__':
    raise SystemExit(main())
