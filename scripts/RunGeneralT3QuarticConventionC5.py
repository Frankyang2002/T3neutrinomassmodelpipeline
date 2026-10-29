"""Source-derived quartic sign mapping and conditional B/C(-1) extra-C5 probe.

Reuses Stage 11, Stage 13 and Stage 14 artifacts; does NOT rerun RGBeta or
Matchete. The source convention is checked each run, and the calculation is
blocked if the expected field words, signs or tensor identities change.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from Numerical.experimental.GeneralT3C5Pilot import load_kernels
from Numerical.experimental.GeneralT3LagrangianConvention import derive_lagrangian_convention
from Numerical.experimental.GeneralT3SourceMappedBridge import source_mapped_uv_to_c5


def _load(path: Path) -> dict:
    value = json.loads(Path(path).read_text(encoding='utf-8-sig'))
    if not isinstance(value, dict):
        raise ValueError(f'{path} must be a JSON object')
    return value


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--input', type=Path, default=Path('docs/general_c5_toy_benchmark.json'))
    ap.add_argument('--matchete-source', type=Path, default=Path('Lagrangian/interactions/GeneralScalarOperators.wl'))
    ap.add_argument('--rgbeta-source', type=Path, default=Path('RGE/running/rgbeta/T3RGBetaModel.wl'))
    ap.add_argument('--trial-root', type=Path, default=Path('output/general_t3_trial'))
    ap.add_argument('--stage13-root', type=Path, default=Path('output/general_t3_stage13'))
    ap.add_argument('--stage14-root', type=Path, default=Path('output/general_t3_stage14'))
    ap.add_argument('--matchete-root', type=Path, default=Path('output/general_t3_numerical_pilot'))
    ap.add_argument('--output', type=Path, default=Path('output/general_t3_stage15'))
    ap.add_argument('--convention-only', action='store_true')
    args = ap.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    errors: list[dict] = []
    rows: list[dict] = []
    convention = None
    try:
        convention = derive_lagrangian_convention(
            args.matchete_source, args.rgbeta_source, args.stage13_root / 'doublet_tensor_map.json'
        )
        (args.output / 'source_lagrangian_convention.json').write_text(
            json.dumps(convention, indent=2, allow_nan=False) + '\n', encoding='utf-8'
        )
        print('Source Lagrangian convention:', convention['status'])
        print('  lambda_Matchete = -lambda_RGBeta (same named scalar fields; no conjugation)')
        print('Report:', args.output / 'source_lagrangian_convention.json')
    except (ValueError, KeyError, OSError, TypeError) as e:
        errors.append({'stage': 'source_convention', 'error': str(e)})
        print('BLOCKED source convention:', e, file=sys.stderr)

    if not args.convention_only and convention is not None:
        try:
            pars = _load(args.input)
            tensor = _load(args.stage13_root / 'doublet_tensor_map.json')
        except (ValueError, OSError) as e:
            errors.append({'stage': 'shared_input', 'error': str(e)})
            pars = None
            tensor = None
        if pars is not None and tensor is not None:
            for model in ('B_m1', 'C_m1'):
                try:
                    s13 = _load(args.stage13_root / f'{model}_coupled_uv.json')
                    s14 = _load(args.stage14_root / f'{model}_conditional_induced_c5.json')
                    ref = _load(args.matchete_root / f'{model}_matchete_hard_c5.json')
                    kernels = load_kernels(
                        args.trial_root / model / 'general_c5_algebraic' / 'delta_c5_atoms.json', model
                    )
                    result = source_mapped_uv_to_c5(
                        model, pars, s13, tensor, convention, kernels, ref, stage14=s14
                    )
                    file = args.output / f'{model}_source_mapped_extra_c5.json'
                    file.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n', encoding='utf-8')
                    mass = result['delta_mnu_ev']
                    rows.append({
                        'model': model,
                        'extra_only_mnu_norm_ev': mass['frobenius_norm'],
                        'smallest_singular_value_ev': mass['takagi_singular_values_ev'][0],
                        'middle_singular_value_ev': mass['takagi_singular_values_ev'][1],
                        'largest_singular_value_ev': mass['takagi_singular_values_ev'][2],
                        'stage14_negative_case_reproduced': result['stage14_negative_case_reproduced'],
                    })
                    print(f'{model} conditional source-mapped ||mnu||_F = {mass["frobenius_norm"]:.8e} eV')
                    print('Report:', file)
                except (ValueError, KeyError, OSError, RuntimeError, TypeError) as e:
                    errors.append({'stage': model, 'error': str(e)})
                    print('BLOCKED', model, e, file=sys.stderr)

    summary = {
        'status': 'Success' if not errors else 'ReviewRequired',
        'source_lagrangian_convention': convention['status'] if convention else 'Blocked',
        'matchete_rgbeta_convention': 'lambda_matchete=-lambda_rgbeta for common scalar basis',
        'models': rows,
        'errors': errors,
        'threshold_matching': False,
        'complete_two_loop_prediction': False,
        'production_code_modified': False,
    }
    (args.output / 'stage15_summary.json').write_text(
        json.dumps(summary, indent=2, allow_nan=False) + '\n', encoding='utf-8'
    )
    if rows:
        with (args.output / 'stage15_summary.csv').open('w', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        print('CSV:', args.output / 'stage15_summary.csv')
    return 0 if not errors else 2


if __name__ == '__main__':
    raise SystemExit(main())
