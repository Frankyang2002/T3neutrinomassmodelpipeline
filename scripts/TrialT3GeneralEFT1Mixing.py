"""Exploratory general EFT1 C5 mixing trial; never modifies authoritative results.

Intended project path: scripts/TrialT3GeneralEFT1Mixing.py

Uses the repository's GeneralEFT1C5Mixing.solve_boundary with authority='synthetic'
for an explicit, user-provided master-Weyl boundary. This turns off *claims of*
certification, not the mathematical structural checks. It is not a validated
conversion of legacy Matchete ordered tensors, nor complete sequential matching.

Expected input: one JSON file already in the format accepted by
RGE/running/intermediate/GeneralEFT1C5Mixing.py (FORMAT), with complete
fermion/scalar field blocks, Wilson operator components, renormalisable Yukawa
and quartic components, and SU(2)L/U(1)Y generators. The trial explicitly
sets 'authority' to 'synthetic' and all certificate flags to False.

Example:
 python -X utf8 -m scripts.TrialT3GeneralEFT1Mixing \
   --boundary path/to/A_m4_master_boundary.json \
   --reference-run corrected_llss16_v4 --case 1111 \
   --trial-id general_eft1_A_m4_trial_v1
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
import math
from pathlib import Path

import numpy as np

from RGE.running.intermediate.GeneralEFT1C5Mixing import (
    FORMAT, CERTIFICATES, solve_boundary, load_boundary,
)
from Numerical.running.GeneralT3SequentialLeadingLog import decode, encode

ROOT = Path(__file__).resolve().parents[1]
INPUT_ROOT = ROOT / 'output/studies/general_scalar_extended/corrected_llss_eft1_freeze'
OUTPUT_ROOT = ROOT / 'output/studies/general_scalar_extended/general_eft1_uncertified_trials'


def read_json(p: Path) -> dict:
    data = json.loads(p.read_text(encoding='utf-8-sig'))
    if not isinstance(data, dict):
        raise ValueError(f'Expected JSON object: {p}')
    return data


def write_json(p: Path, value: dict) -> None:
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n', encoding='utf-8')


def non_scaling(full: np.ndarray, other: np.ndarray) -> tuple[float, dict]:
    denom = np.vdot(other, other)
    if abs(denom) < 1e-100:
        return float('nan'), {'real': 0., 'imag': 0.}
    factor = np.vdot(other, full) / denom
    residual = np.linalg.norm(full - factor * other) / max(np.linalg.norm(full), 1e-100)
    return float(residual), {'real': float(factor.real), 'imag': float(factor.imag)}


def run(args: argparse.Namespace) -> None:
    boundary = args.boundary.resolve()
    if not boundary.is_file():
        raise FileNotFoundError(f'Missing master-Weyl input: {boundary}')
    raw = boundary.read_bytes()
    payload = read_json(boundary)
    if payload.get('format') != FORMAT:
        raise ValueError('Input is not a master-Weyl EFT1 boundary. Do not pass a legacy ordered tensor.')
    if payload.get('model') != args.model:
        raise ValueError(f'Input model {payload.get("model")} differs from --model {args.model}')
    if not args.trial_id or not all(c.isalnum() or c in '_-' for c in args.trial_id):
        raise ValueError('--trial-id must contain alphanumerics, _ or -')
    if args.case not in ('1111', '0111', '1011', '1101'):
        raise ValueError('Select an active W case: 1111, 0111, 1011, or 1101')

    # This is an explicit, transparently uncertified trial. We NEVER set
    # certificates to True and NEVER patch the general evaluator's physics.
    payload['authority'] = 'synthetic'
    payload['certificates'] = {name: False for name in CERTIFICATES}
    # Structural input checks and actual tensor algebra remain intact.
    load_boundary(payload, allow_synthetic=True)
    general = solve_boundary(payload, allow_synthetic=True,
                             trajectory_points=args.trajectory_points)
    general_log = decode(general['total_delta_C5_log_gev_inverse'])

    tag = args.model.replace('+', 'p').replace('-', 'm')
    original_path = INPUT_ROOT / args.reference_run / tag / 'cases.json'
    if not original_path.is_file():
        raise FileNotFoundError(f'Missing corrected LLSS source: {original_path}')
    original = read_json(original_path)
    if original.get('model') != args.model or original.get('status') != 'Computed':
        raise ValueError('Reference case is not a successfully computed result for this model')
    cases = original.get('cases', {})
    full_case = cases.get(args.case)
    if not isinstance(full_case, dict):
        raise ValueError(f'Reference case {args.case} absent')
    base_w_off = args.case[:2] + '0' + args.case[3]
    frozen_case = cases.get(base_w_off)
    if not isinstance(frozen_case, dict):
        raise ValueError(f'Reference W=0 case {base_w_off} absent')
    original_log = decode(full_case['c5_induced_wilson_log'])
    bubble = decode(frozen_case['c5_joint_scalar_bubble'])
    current_c5 = decode(full_case['c5_at_scalar_threshold'])
    if not np.allclose(current_c5, bubble + original_log, atol=1e-29, rtol=1e-9):
        raise ValueError('Reference C5 is not the corrected independent LLSS sum')
    # IMPORTANT: sum only once. The general result is used INSTEAD OF the LLSS
    # log, never added on top of it. The bubble remains from the old provisional
    # convention and thus the sum is an UNCERTIFIED comparison only.
    trial_c5 = 0.5 * (bubble + general_log + (bubble + general_log).T)
    discrepancy = float(np.linalg.norm(general_log - original_log)
                        / max(np.linalg.norm(original_log), 1e-100))
    nonscaling, complex_ratio = non_scaling(trial_c5, current_c5)

    dest = OUTPUT_ROOT / args.trial_id
    if dest.exists():
        raise FileExistsError(f'Trial ID already exists: {dest}')
    result = {
        'format': 't3-general-eft1-uncertified-side-by-side-v1',
        'model': args.model,
        'status': 'EXPLORATORY_UNCERTIFIED_NOT_PHYSICAL_C5',
        'reference_run': args.reference_run,
        'reference_case': args.case,
        'reference_w_off_case': base_w_off,
        'master_boundary_path': str(boundary),
        'master_boundary_sha256': sha256(raw).hexdigest(),
        'certification_bypassed_for_trial': True,
        'actual_master_tensor_structural_checks_kept': True,
        'input_gauge_sector_names': [g.get('name') for g in payload.get('gauge_sectors', [])],
        'input_operator_source_count': len(payload.get('operators', [])),
        'input_yukawa_component_count': len(payload.get('yukawa_components', [])),
        'input_quartic_component_count': len(payload.get('quartic_components', [])),
        'seven_term_general_leading_log': encode(general_log),
        'legacy_LLSS_leading_log': encode(original_log),
        'finite_scalar_bubble_legacy': encode(bubble),
        'general_log_plus_legacy_bubble_EXPLORATORY': encode(trial_c5),
        'legacy_reference_C5': encode(current_c5),
        'general_vs_LLSS_log_relative_Frobenius': discrepancy,
        'trial_vs_legacy_C5_nonscaling': nonscaling,
        'trial_vs_legacy_C5_best_scalar': complex_ratio,
        'general_solver': general,
        'limitations': [
            'Not certified: master Wilson normalization and Matchete-to-Weyl conversion must be established.',
            'Not certified: basis alignment and all missing Wilson sources must be audited.',
            'Finite bubble originates in the previous provisional LLSS basis.',
            'This is fixed-order leading-log, not numerical evolution of all Wilson operators.',
            'Ordinary Weinberg self-running is not separately included unless the supplied operator inventory represents it.',
            'Not a physical neutrino prediction or a replacement for the validated dataset.',
        ],
    }
    write_json(dest / 'comparison.json', result)
    print(f'Model: {args.model} / case {args.case}')
    print(f'General seven-term relative to LLSS log: {discrepancy:.8g}')
    print(f'Candidate C5 non-scaling relative to old LLSS C5: {nonscaling:.8g}')
    print(f'Ward residual: {general["ward_residual_max"]:.8g}')
    print(f'Output: {dest / "comparison.json"}')
    print('UNCERTIFIED TRIAL: no physical mixing-angle claim')


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--boundary', type=Path, required=True,
                    help='Explicit full general master-Weyl boundary JSON')
    ap.add_argument('--model', required=True)
    ap.add_argument('--reference-run', default='corrected_llss16_v4')
    ap.add_argument('--case', default='1111')
    ap.add_argument('--trial-id', required=True)
    ap.add_argument('--trajectory-points', type=int, default=3)
    run(ap.parse_args())
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
