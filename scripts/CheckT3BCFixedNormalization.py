"""Read-only B/C check of FIXED sqrt(3),2,2 EFT1 channel normalizations.

Uses the previous signed B/C sequential ablation and independent Matchete
common-log audit. No fitted parameters, RGEs, or saved reports are modified.
This tests numerical consistency only, not physical coupling-map certification.

From the repository root:
 python -X utf8 -m scripts.CheckT3BCFixedNormalization
"""
from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np

from Numerical.running.GeneralT3ConditionalBCQuarticNormalization import (
    CHANNEL_FACTORS,
)
from Numerical.running.GeneralT3MatcheteLogCoefficientAudit import (
    read_json, decode_matrix,
)

ROOT = Path(__file__).resolve().parents[1]
CASE_ROOT = ROOT / 'output/studies/general_scalar_extended/general_sequential_rge_ablation'
AUDIT_ROOT = ROOT / 'output/studies/general_scalar_extended/general_log_coefficient_audit'
MODELS = ('B_-1', 'C_-1')


def require(value: bool, message: str) -> None:
    if not value:
        raise ValueError('B/C fixed-channel log check: ' + message)


def fixed_spurion_from_signed_case(record: dict, model: str) -> np.ndarray:
    """Build common-log derivative with NO fitted or per-model normalizations."""
    require(record.get('model') == model and record.get('status') == 'Computed',
            'expected signed B/C model record')
    case = record['cases']['1111']
    provenance = case.get('quartic_provenance', {})
    require(provenance.get('bc_direct_quartic_channel_signs') == {
        'general:lambdaHHdagS2S2': -1,
        'general:lambdaHHdagS1barS1bar': -1,
    }, 'input ablation is not the signed direct-quartic calculation')
    require('bc_channel_normalization_factors' not in provenance,
            'input ablation is already channel-normalized; use the OLD signed run')
    channels = [x for x in case['per_general_interaction']
                if x.get('category') == 'quartic_decomposition']
    names = [x['name'] for x in channels]
    require(len(names) == len(set(names)) and all(names.count(x) == 1 for x in CHANNEL_FACTORS),
            'missing/duplicate expected quartic channel')
    converted = {row['name']: decode_matrix(row['c5_leading_log']) for row in channels}
    unweighted_total = sum(converted.values(), np.zeros((3, 3), complex))
    saved_total = decode_matrix(case['c5_induced_wilson_log'])
    rel = np.linalg.norm(unweighted_total-saved_total) / max(np.linalg.norm(saved_total), 1e-300)
    require(rel < 1e-8, f'{model}: stored EFT1 channel sum mismatch ({rel:.3g})')
    other_channels = {name: mat for name, mat in converted.items()
                      if name not in CHANNEL_FACTORS}
    other_norm = np.linalg.norm(sum(other_channels.values(),np.zeros((3,3),complex)))
    require(other_norm < 1e-8 * max(np.linalg.norm(saved_total), 1e-300),
            f'{model}: additional nonzero HHSS channels require individual audit')
    corrected_log = sum((factor*converted[name] for name,factor in CHANNEL_FACTORS.items()),
                        np.zeros((3, 3), complex))
    return corrected_log / (2 * math.log(1e5 / 1e3))


def compare_fixed(candidate: np.ndarray, matchete: np.ndarray) -> dict:
    a = np.asarray(candidate, dtype=complex)
    m = np.asarray(matchete, dtype=complex)
    require(a.shape == m.shape == (3,3) and np.isfinite(a).all() and
            np.isfinite(m).all() and np.linalg.norm(a) > 0 and
            np.linalg.norm(m) > 0, 'invalid nonzero finite C5 derivatives')
    best = np.vdot(a, m) / np.vdot(a, a)
    den = np.linalg.norm(m)
    return {'fixed_relative_difference': float(np.linalg.norm(m-a)/den),
            'best_fit_real_diagnostic_only': float(best.real),
            'best_fit_imag_diagnostic_only': float(best.imag),
            'best_fit_matrix_residual_diagnostic_only':
                float(np.linalg.norm(m-best*a)/den)}


def verify_model(model: str, signed_path: Path, audit_path: Path, *,
                 max_relative_difference: float) -> list[dict]:
    require(model in MODELS, 'only B_-1/C_-1 supported')
    signed, audit = read_json(signed_path), read_json(audit_path)
    require(audit['model'] == model and audit['status'] == 'Computed',
            'wrong Matchete reference log audit')
    require(audit['one_scale_reproduction_relative_residual'] < 1e-4,
            'Matchete one-scale baseline failed its parity audit')
    candidate = fixed_spurion_from_signed_case(signed, model)
    results = []
    for item in audit['hierarchy_scan']:
        h = item['heavy_mass_scale_factor']
        require(isinstance(h, int) and h >= 1,
                'invalid fermion mass scaling in Matchete audit')
        actual = decode_matrix(item['matchete_common_log_coefficient'])
        comparison = compare_fixed(candidate/h, actual)
        result = {'model':model,'heavy_mass_factor':h,**comparison}
        results.append(result)
        print(f"{model:4} MF x{h:<2} fixed_diff="
              f"{100*comparison['fixed_relative_difference']:.6f}% "
              f"best_fit={comparison['best_fit_real_diagnostic_only']:.8f} "
              f"shape={100*comparison['best_fit_matrix_residual_diagnostic_only']:.6f}%",
              flush=True)
        require(comparison['fixed_relative_difference'] < max_relative_difference,
                f'{model} MF x{h}: fixed factors do not reproduce common log '
                f'(relative difference {comparison["fixed_relative_difference"]:.6g})')
    return results


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case-run',default='general_seq_bc_signed_01')
    parser.add_argument('--reference-audit',default='log_bc_signed_01')
    parser.add_argument('--max-relative-difference',type=float,default=0.001)
    args = parser.parse_args(argv)
    require(math.isfinite(args.max_relative_difference) and
            0 < args.max_relative_difference < 1,
            'invalid comparison tolerance')
    for model in MODELS:
        tag = model.replace('-','m')
        case=CASE_ROOT/args.case_run/tag/'cases.json'
        audit=AUDIT_ROOT/args.reference_audit/tag/'log_audit.json'
        require(case.is_file() and audit.is_file(),
                f'{model}: missing old signed run or reference audit')
        verify_model(model,case,audit,
                     max_relative_difference=args.max_relative_difference)
    print('B/C FIXED factor comparison PASSED (numerical convention check only)')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
