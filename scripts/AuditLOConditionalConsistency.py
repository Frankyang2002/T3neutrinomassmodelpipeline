"""Run the all-16 read-only conditional-portal consistency and external evidence audit.

No source modifications. New outputs restricted to the Workstream 3 subtree.
Default inputs are the previously generated scalar_hard_01 and published
upstream WS2/legacy v6 evidence, when available. Read-only and provisional.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from Numerical.lo_certification.completeness.consistency_audit import (
    audit_model, assemble, crossworkstream_evidence, ConsistencyError,
)
from Numerical.lo_certification.completeness.inventory import ALL_MODELS

OWNED = Path('output/general_observable/lo_certification_completeness')
DEFAULT_RUN = OWNED / 'scalar_hard_01'
DEFAULT_PORTAL = Path('output/general_observable/lo_certification_couplings/ws2_portal_normalization_v1/summary.json')
DEFAULT_GAMMACC = Path('output/audits/t3_bc_y_gammacc_convention_v1/summary.json')
DEFAULT_PREVIOUS = Path('output/general_observable/lo_joint_scalar_matching/exploratory_all16_v6/summary.json')


def resolved(p: Path) -> Path:
    return (p if p.is_absolute() else ROOT / p).resolve()


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-run-dir', type=Path, default=DEFAULT_RUN)
    parser.add_argument('--portal-summary', type=Path, default=DEFAULT_PORTAL)
    parser.add_argument('--gammacc-report', type=Path, default=DEFAULT_GAMMACC)
    parser.add_argument('--previous-scalar-summary', type=Path, default=DEFAULT_PREVIOUS)
    parser.add_argument('--output-dir', type=Path,
                        default=OWNED / 'conditional_consistency_01')
    parser.add_argument('--check-only', action='store_true',
                        help='Compute audit and print compact result without writing files')
    args = parser.parse_args(argv)
    source = resolved(args.source_run_dir)
    summary_path = source / 'summary.json'
    original = json.loads(summary_path.read_text(encoding='utf-8-sig'))
    if (original.get('format') != 't3-lo-scalar-hard-source-bridge-v1' or
            original.get('model_count') != 16 or
            original.get('conditional_hard_models') != 16 or
            original.get('physical_C5_certified') is not False or
            {r.get('model') for r in original.get('models', [])} != set(ALL_MODELS)):
        parser.error('Expected a complete 16/16 provisional scalar_hard release')
    target = resolved(args.output_dir)
    allowed = resolved(OWNED)
    if target == allowed or allowed not in target.parents:
        parser.error('Output must be a new child directory of ' + str(allowed))
    if target.exists() and not args.check_only:
        parser.error('Refusing to overwrite an existing output directory: ' + str(target))
    evidence = crossworkstream_evidence(resolved(args.portal_summary),
                                        resolved(args.gammacc_report),
                                        resolved(args.previous_scalar_summary))
    results = []
    errors = []
    for model in ALL_MODELS:
        try:
            result = audit_model(model, source)
            results.append(result)
            print(f'{model}: {result["status"]}, UV remainder={result["relative_exact_triangle_truncation_error"]:.3e}',
                  flush=True)
        except (OSError, ValueError, KeyError, TypeError, ArithmeticError) as exc:
            errors.append({'model': model, 'error': f'{type(exc).__name__}: {exc}'})
            print(f'{model}: FAILED ({exc})', flush=True)
    if errors:
        outcome = {
            'format': 't3-lo-conditional-portal-crossworkstream-audit-v1',
            'model_count': 16, 'models_audited_successfully': len(results),
            'failed_checks': errors, 'physical_C5_certified': False,
        }
        if not args.check_only:
            target.mkdir(parents=True, exist_ok=False)
            (target / 'summary.json').write_text(json.dumps(outcome, indent=2)+'\n', encoding='utf-8')
        print(json.dumps(outcome, indent=2))
        return 2
    outcome = assemble(results, evidence)
    compact = {k: outcome[k] for k in ('model_count', 'portal_numeric_consistency_passed',
        'portal_numeric_incomplete_or_failed', 'physical_C5_certified',
        'previous_scalar_projection_failures', 'models_with_missing_native_portal_directions')}
    if not args.check_only:
        target.mkdir(parents=True, exist_ok=False)
        (target / 'summary.json').write_text(json.dumps(outcome, indent=2, allow_nan=False)+'\n',
                                             encoding='utf-8')
        with (target / 'numerical_table.csv').open('w', newline='', encoding='utf-8') as f:
            fields = ['model', 'status', 'hard_frobenius_per_unknown_N_gamma_gev_inverse',
                      'dim5_frobenius_per_unknown_N_gamma_gev_inverse',
                      'hard_to_dim5_norm_ratio', 'relative_exact_triangle_truncation_error',
                      'max_projection_residual',
                      'matching_scale_muF_relative_drift_fixed_inputs',
                      'matching_scale_muS_relative_drift_fixed_inputs',
                      'scalar_mass2_power_exponent_1_to_1quarter',
                      'scalar_mass2_power_exponent_1quarter_to_1sixteenth',
                      'UV_power_scan_status', 'physical_C5_certified']
            writer = csv.DictWriter(f, fieldnames=fields)
            writer.writeheader()
            for row in outcome['models']:
                writer.writerow({k: row.get(k) for k in fields})
        (target / 'blockers.json').write_text(json.dumps({
            'physical_C5_certified': False,
            'models': [{'model': r['model'], 'physical_blockers': r['physical_blockers'],
                        'portal_span': r['portal_span'],
                        'earlier_scalar_comparison': r['earlier_scalar_comparison']}
                       for r in outcome['models']]}, indent=2, allow_nan=False)+'\n',
            encoding='utf-8')
    print(json.dumps({'output': None if args.check_only else str(target), **compact}, indent=2))
    return 0 if outcome['portal_numeric_consistency_passed'] == 16 else 2


if __name__ == '__main__':
    raise SystemExit(main())
