"""Source-checked Batch 2 loop-function/Laurent audit and native MF diagnostics.

No physical Matchete/RGBeta dictionary, finite MSbar Wilson matching, or
anomalous-dimension evolution is implied by this module.  The two candidate
matrix contractions below are *diagnostics*, not effective operators.
"""
from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from Numerical.parallel_t3.background.general_batch1_api import decode_native, encode_native

REPORT_FORMAT = 't3-batch2-laurent-operator-audit-v2'
PIPELINE_FORMAT = 't3-thesis-full-general-batch2-loop-audit-v1'
RAW_FORMAT = 't3-thesis-full-general-batch2-raw-fermion-v1'
RAW_SUCCESS = 'RawFullGeneralOneLoopFThresholdSourceExported_NotFiniteWilsons'
MODELS = ('A_-4', 'D_-2')


class SourceAuditError(ValueError):
    """Evidence or an explicitly required stage is absent or inconsistent."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise SourceAuditError('Batch2 loop audit: ' + message)


def sha(path: Path) -> str:
    path = Path(path)
    require(path.is_file(), 'file not found: ' + str(path))
    return sha256(path.read_bytes()).hexdigest()


def read_json(path: Path) -> dict:
    value = json.loads(Path(path).read_text(encoding='utf-8-sig'))
    require(isinstance(value, dict), 'expected JSON object: ' + str(path))
    return value


def save_json(path: Path, data: Mapping[str, Any]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + '.partial')
    tmp.write_text(json.dumps(data, sort_keys=True, indent=2, allow_nan=False) + '\n',
                   encoding='utf-8')
    tmp.replace(path)


def check_saved_raw(model: str, batch1_dir: Path, raw_dir: Path) -> dict:
    """Validate the exact source/raw ledger and hashes against local Batch 1.

    This intentionally does not attempt to validate a physical basis dictionary.
    The archived matched EFT must not be silently overwritten.
    """
    from Numerical.parallel_t3.background.general_batch1_api import load_batch1
    from Numerical.parallel_t3.wilson.general_batch2_source import verify_ledger

    require(model in MODELS, 'unknown model')
    b1 = load_batch1(Path(batch1_dir))
    require(b1.model == model, 'Batch 1 model mismatch')
    manifest_file = Path(raw_dir) / 'manifest.json'
    manifest = read_json(manifest_file)
    require(manifest.get('format') == RAW_FORMAT and manifest.get('model') == model,
            'raw matching manifest schema/model mismatch')
    require(manifest.get('status') == RAW_SUCCESS, 'raw matching has not succeeded')
    checks = {
        'batch1_manifest': 'manifest.json',
        'batch1_fermion_matching_input': 'fermion_matching_input.json',
        'batch1_scalar_matching_input': 'scalar_matching_input.json',
        'batch1_trajectory': 'trajectories.npz',
    }
    for key, filename in checks.items():
        require(manifest.get('input_sha256', {}).get(key) == sha(Path(batch1_dir) / filename),
                'Batch 1 source changed since matching: ' + filename)
    one_loop = Path(raw_dir) / 'native_matching/data/eft1_after_F_one_loop_inputform.txt'
    require(one_loop.is_file() and one_loop.stat().st_size > 0,
            'matching one-loop expression absent')
    native = Path(raw_dir) / 'native_matching/native_status.json'
    ledger_file = Path(raw_dir) / 'one_loop_operator_ledger.json'
    source_hashes = manifest.get('output_files_sha256', {})
    require(source_hashes.get('native_status') == sha(native),
            'native status changed after successful matching')
    require(source_hashes.get('operator_ledger') == sha(ledger_file),
            'ledger changed after successful matching')
    native_data = read_json(native)
    require(native_data.get('status') == 'FThresholdTreeWilsonExported'
            and native_data.get('model') == model
            and native_data.get('loop_order') == 1
            and native_data.get('integrated_fields') == ['F']
            and native_data.get('light_bsm_fields') == ['S1', 'S2'],
            'F-only matching did not retain both light scalars')
    ledger = verify_ledger(read_json(ledger_file), one_loop)
    require(ledger['term_count'] > 0, 'empty matched one-loop correction')
    require(manifest.get('native_running_values_inserted_into_Matchete_Lagrangian') is False,
            'unexpected claim of physical Matchete/RGBeta numerical mapping')
    # Verify the saved numerical handover is still exact, not only the hash of Batch 1.
    fsrc = b1.fermion_matching_input()
    saved_numeric = read_json(Path(raw_dir) / 'numerical_fermion_boundary.json')
    expected_native = fsrc['native_uv_rgbeta_couplings_at_muF']
    require(set(expected_native) == set(saved_numeric.get('native_parameter_values', {})),
            'native parameter inventory not complete at F')
    for name, value in expected_native.items():
        stored = saved_numeric['native_parameter_values'][name]
        require(np.allclose(decode_native(value), decode_native(stored),
                            rtol=1e-14, atol=1e-15),
                'saved F-boundary value changed: ' + name)
    return {
        'model': model,
        'one_loop_source': one_loop,
        'one_loop_source_sha256': sha(one_loop),
        'batch1_manifest_sha256': sha(Path(batch1_dir) / 'manifest.json'),
        'raw_manifest_sha256': sha(manifest_file),
        'raw_ledger_sha256': sha(ledger_file),
        'native_matching_status_sha256': sha(native),
        'raw_one_loop_term_count': ledger['term_count'],
        'candidate_sector_counts': ledger['raw_one_loop_operator_sector_counts'],
        'numeric': saved_numeric,
    }


def matrix_spurion_diagnostics(numeric: Mapping[str, Any], model: str) -> dict:
    """Two explicitly convention-UNRESOLVED candidate flavour contractions.

    Both coincide for diagonal real MF. The actual Matchete ΓCC/spinor
    normalization and chiral-index orientation are not inferred here.
    Not a numerical EFT1 Wilson boundary.
    """
    require(model in MODELS and numeric.get('model') == model, 'numeric model mismatch')
    values = numeric.get('native_parameter_values')
    require(isinstance(values, dict), 'missing native parameters')
    for name in ('MF', 'y1', 'y2'):
        require(name in values, 'missing numeric matrix ' + name)
    matrices = {k: np.asarray(decode_native(values[k]), dtype=np.complex128)
                for k in ('MF', 'y1', 'y2')}
    MF = matrices['MF']; y1 = matrices['y1']; y2 = matrices['y2']
    require(all(a.shape == (3, 3) and np.isfinite(a).all() for a in matrices.values()),
            'all heavy-flavour matrices must be finite 3x3')
    cond = float(np.linalg.cond(MF))
    require(np.isfinite(cond) and cond < 1e14,
            'heavy MF cannot safely be inverted in native coordinates')
    masses = np.linalg.svd(MF, compute_uv=False)
    diag = np.diag(MF)
    require(np.min(np.abs(diag)) > 0, 'mass diagonal vanishes')
    offdiag = float(np.linalg.norm(MF - np.diag(diag)) / np.linalg.norm(MF))

    def candidates(rhs: np.ndarray) -> dict:
        require(rhs.shape == (3, 3), 'right Yukawa shape mismatch')
        # Permissible coordinate algebraic forms only. The chiral map determines
        # which one corresponds to actual Matchete coefficient (or neither).
        wa = y1 @ np.linalg.solve(MF, rhs.T)
        wb = y1 @ np.linalg.solve(MF.T, rhs.T)
        wd = (y1 / diag[None, :]) @ rhs.T
        return {
            'Y1_MFinverse_Yright_T_GeVinv': encode_native(wa),
            'Y1_MFtransposeinverse_Yright_T_GeVinv': encode_native(wb),
            'diagonal_MF_reference_GeVinv': encode_native(wd),
            'orientation_max_abs_difference_GeVinv': float(np.max(np.abs(wa - wb))),
            'diagonal_limit_max_abs_difference_GeVinv':
                float(max(np.max(np.abs(wa - wd)), np.max(np.abs(wb - wd)))),
            'candidate_matrix_norms_GeVinv': {
                'inverse': float(np.linalg.norm(wa)),
                'inverse_transpose': float(np.linalg.norm(wb)),
            },
        }
    ans = {
        'model': model,
        'status': 'AlgebraicHeavyFlavourDiagnosticsOnly',
        'heavy_fermion_masses_singular_GeV': sorted(float(m) for m in masses),
        'MF_condition_number': cond,
        'MF_off_diagonal_relative_Frobenius': offdiag,
        'baseline_Y1_Y2': candidates(y2),
        'native_to_Matchete_heavy_chirality_map_verified': False,
        'physical_Wilson_coefficient': False,
        'matchete_spinor_and_CG_normalisation_applied': False,
    }
    if model == 'D_-2':
        require('y2Conj' in values and 'muT3GeneralS2Sq' in values,
                'D(-2) extra Yukawa / self-bilinear omitted')
        extra = np.asarray(decode_native(values['y2Conj']), dtype=np.complex128)
        require(extra.shape == (3, 3) and np.isfinite(extra).all(),
                'invalid D(-2) extra Yukawa')
        ans['extra_Y1_Y2Conj'] = candidates(extra)
        ans['D_self_bilinear_at_F_GeV2'] = values['muT3GeneralS2Sq']
        ans['D_five_extra_quartics_at_F'] = {
            k: v for k, v in values.items() if k.startswith('lambdaT3GeneralQ')}
        require(len(ans['D_five_extra_quartics_at_F']) == 5,
                'missing D(-2) general quartic at F')
    return ans


def check_audit_report(report: Mapping[str, Any], expected_hash: str,
                       output_folder: Path) -> dict:
    require(report.get('format') == REPORT_FORMAT, 'wrong Laurent audit report format')
    require(report.get('status') in {
        'LoopFunctionsEvaluated_ExplicitLaurentSeparated_TreeRemainderPreserved_AuditOnly',
        'LoopFunctionsEvaluated_NoExplicitRegulator_TreeRemainderPreserved_AuditOnly',
    }, 'Wolfram loop-function/Laurent audit did not finish: '
            + str(report.get('reason', report.get('status'))))
    require(report.get('source_sha256') == expected_hash,
            'Wolfram read different one-loop source bytes')
    require(report.get('full_expression_hbar_decomposition_exact') is True and
            report.get('hbar_higher_order_absent') is True and
            report.get('hbar1_nonzero') is True and
            report.get('all_operators_preserved') is True,
            'exact hbar^0 + hbar^1 decomposition or full source preservation failed')
    require(report.get('hbar0_only_recognised_quadratic_mass_terms') is True,
            'hbar^0 remainder has not been checked against the known mass inventory')
    require(report.get('unreconciled_tree_mass_coordinate_identity_verified') is False
            and report.get('tree_mass_remainder_dropped_from_source') is False
            and report.get('tree_mass_remainder_included_in_one_loop_coefficient') is False,
            'unverified tree mass coordinate relation must not be imposed')
    tree_count = report.get('unreconciled_hbar0_mass_term_count')
    loop_count = report.get('raw_hbar1_term_count')
    source_count = report.get('original_source_term_count')
    require(all(type(n) is int and n >= 0 for n in (tree_count,loop_count,source_count))
            and tree_count in (0,4,6) and loop_count > 0
            and source_count == tree_count + loop_count,
            'unexpected hbar^0 inventory or source operator accounting')
    require(report.get('unreconciled_hbar0_mass_remainder_present') is (tree_count > 0)
            and report.get('one_loop_hbar_linearity_exact') is (tree_count == 0),
            'incorrect original hbar-linearity classification')
    require(report.get('native_fermions_absent') is True,
            'heavy fermion left in matched EFT1')
    require(report.get('renormalised_MSbar_Wilson_boundary_computed') is False and
            report.get('source_scale_substitution_performed') is False,
            'raw audit incorrectly claims physical numerical Wilsons')
    files = report.get('output_sha256', {})
    required = {'raw_hbar1_coefficient.wxf',
                'unreconciled_hbar0_mass_remainder.wxf',
                'evaluated_full_hbar1_coefficient.wxf'}
    if 'ExplicitLaurentSeparated' in report['status']:
        required |= {'explicit_pole_residue_unsubtracted.wxf',
                     'epsilon_zero_regular_unsubtracted.wxf'}
        require(report.get('explicit_laurent_reconstruction_exact') is True,
                'the explicit Laurent decomposition is not exact')
    require(isinstance(files, dict) and required <= set(files),
            'audit omitted full Wilson one-loop/zero-order/Laurent WXF expressions')
    for filename, digest in files.items():
        require(isinstance(filename, str) and '/' not in filename and '\\' not in filename
                and filename.endswith('.wxf'), 'invalid Wolfram export name')
        require(sha(output_folder / filename) == digest,
                'Wolfram symbolic output checksum mismatch: ' + filename)
    return {
        'status': report['status'],
        'source_sha256': expected_hash,
        'wolfram_output_sha256': files,
        'regulator_symbol_inputform': report.get('regulator_symbol_inputform'),
        'pole_expression_formal_only': True,
        'unreconciled_hbar0_mass_term_count': tree_count,
        'raw_hbar1_term_count': loop_count,
        'original_expression_hbar_linear': tree_count == 0,
        'full_expression_hbar_decomposition_exact': True,
        'numerical_finite_Wilson_matching_complete': False,
    }
