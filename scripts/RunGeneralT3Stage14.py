"""Run conditional UV-induced extra C5 response estimates for B/C(-1)."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from Numerical.experimental.GeneralT3C5Pilot import load_kernels
from Numerical.experimental.GeneralT3ConditionalUVBridge import conditional_uv_to_c5


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--input', type=Path, default=Path('docs/general_c5_toy_benchmark.json'))
    ap.add_argument('--trial-root', type=Path, default=Path('output/general_t3_trial'))
    ap.add_argument('--coupled-root', type=Path, default=Path('output/general_t3_stage13'))
    ap.add_argument('--matchete-root', type=Path, default=Path('output/general_t3_numerical_pilot'))
    ap.add_argument('--output', type=Path, default=Path('output/general_t3_stage14'))
    args = ap.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    errors = []
    summary_rows = []
    try:
        parameters = json.loads(args.input.read_text(encoding='utf-8-sig'))
        tensor_map = json.loads((args.coupled_root / 'doublet_tensor_map.json').read_text(encoding='utf-8-sig'))
        if not isinstance(parameters, dict):
            raise ValueError('Toy input must be a JSON object')
    except (OSError, ValueError, KeyError, TypeError) as exc:
        errors.append({'stage': 'shared_inputs', 'error': str(exc)})
        parameters = None
        tensor_map = None
    if parameters is not None and tensor_map is not None:
        for model in ('B_m1', 'C_m1'):
            try:
                stage13 = json.loads((args.coupled_root / f'{model}_coupled_uv.json').read_text(encoding='utf-8-sig'))
                reference = json.loads((args.matchete_root / f'{model}_matchete_hard_c5.json').read_text(encoding='utf-8-sig'))
                kernels = load_kernels(args.trial_root / model / 'general_c5_algebraic' / 'delta_c5_atoms.json', model)
                result = conditional_uv_to_c5(model, parameters, stage13, tensor_map, kernels, reference)
                path = args.output / f'{model}_conditional_induced_c5.json'
                path.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n', encoding='utf-8')
                limits = result['range']
                summary_rows.append({'model': model, **limits})
                print(f'{model}: hypothetical induced extra-only ||mnu||_F between '
                      f'{limits["minimum_extra_only_mnu_norm_ev"]:.6e} and '
                      f'{limits["maximum_extra_only_mnu_norm_ev"]:.6e} eV')
                print('Report:', path)
            except (ValueError, OSError, KeyError, TypeError, RuntimeError) as exc:
                errors.append({'stage': model, 'error': str(exc)})
                print('BLOCKED:', model, exc, file=sys.stderr)
    summary = {
        'status': 'Success' if not errors else 'ReviewRequired',
        'scope': 'conditional sensitivity estimate; no full RG-to-Matchete matching',
        'full_lagrangian_coupling_mapping_verified': False,
        'threshold_consistent_matching': False,
        'rows': summary_rows,
        'errors': errors,
    }
    (args.output / 'stage14_summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    if summary_rows:
        with (args.output / 'stage14_summary.csv').open('w', encoding='utf-8', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=list(summary_rows[0]))
            writer.writeheader()
            writer.writerows(summary_rows)
    return 0 if not errors else 2


if __name__ == '__main__':
    raise SystemExit(main())
