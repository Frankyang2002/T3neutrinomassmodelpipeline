"""Non-destructive T3 full-scalar readiness audit.

Run at project root: python scripts/AuditFullScalarReadiness.py
Structural inventory only; it does NOT certify RG coefficients or C5 matching.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REGISTRY = ROOT / 'RGE/general/generated/general_rgbeta_quartic_registry.json'
FILES = {
    'operator_basis': 'Lagrangian/interactions/GeneralScalarOperators.wl',
    'verified_quartics': 'Lagrangian/interactions/GeneralScalarVerifiedQuartics.wl',
    'ordinary_potential': 'Lagrangian/interactions/ScalarPotential.wl',
    'rgbeta_model': 'RGE/running/rgbeta/T3RGBetaModel.wl',
    'general_rgbeta_runner': 'RGE/running/rgbeta/RunGeneralT3RGBeta.wl',
    'quadratic_rgbeta_runner': 'RGE/running/rgbeta/RunGeneralT3QuadraticRGBeta.wl',
    'quadratic_mass_module': 'Numerical/experimental/GeneralT3ScalarQuadraticMasses.py',
    'mixed_mass_module': 'Numerical/experimental/GeneralT3MixedScalarMassBasis.py',
    'scalar_rg_state': 'Numerical/experimental/GeneralT3ScalarRGState.py',
    'general_campaign': 'scripts/RunT3GeneralScalar48Campaign.py',
}
EXPECTED = {
    'A_+0', 'A_-2', 'A_-4', 'B_+1', 'B_-1', 'B_-3',
    'C_+1', 'C_-1', 'C_-3', 'D_+0', 'D_+2', 'D_-2',
    'E_+0', 'E_+2', 'E_-2', 'E_-4',
}


def audit(root: Path) -> dict:
    errors = []
    present = {name: (root / rel).is_file() for name, rel in FILES.items()}
    missing = [name for name, yes in present.items() if not yes]
    if missing:
        errors.append('Missing source files: ' + ', '.join(missing))
    path = root / 'RGE/general/generated/general_rgbeta_quartic_registry.json'
    if not path.is_file():
        return {'status': 'BLOCKED', 'errors': errors + ['Missing quartic registry'],
                'source_files_present': present, 'models': {}}
    registry = json.loads(path.read_text(encoding='utf-8-sig'))
    models = registry.get('models', {})
    if set(models) != EXPECTED:
        errors.append(f'Model inventory differs: missing={sorted(EXPECTED-set(models))}; '
                      f'unexpected={sorted(set(models)-EXPECTED)}')
    report = {}
    for key, entry in sorted(models.items()):
        added = entry.get('new_rgbeta_symbols', [])
        reused = entry.get('preexisting_rgbeta_symbols', [])
        total = entry.get('total_additional_quartics')
        quadratics = entry.get('quadratics_not_implemented')
        if total != len(added) + len(reused):
            errors.append(f'{key}: quartic registry inventory inconsistent')
        names = [v.get('name') for v in added] + [v.get('existing_name') for v in reused]
        if len(names) != len(set(names)) or any(not isinstance(s, str) or not s for s in names):
            errors.append(f'{key}: duplicate or invalid symbol names')
        if quadratics is None:
            errors.append(f'{key}: missing quadratic status')
        report[key] = {
            'additional_quartics': total,
            'new_rgbeta_symbols': len(added),
            'reused_rgbeta_symbols': len(reused),
            'quadratic_structures_not_implemented': quadratics,
            'quartic_normalisations_verified': all(
                v.get('mapping_status') in ('Verified', 'Normalised') for v in added
            ) and len(added) > 0,
            'physical_c5_full_model_certified': False,
        }
    totals = {
        'models': len(report),
        'additional_quartics_across_models': sum(int(v['additional_quartics'] or 0) for v in report.values()),
        'quadratic_structures_not_implemented_across_models': sum(
            int(v['quadratic_structures_not_implemented'] or 0) for v in report.values()),
    }
    return {
        'status': 'BLOCKED' if errors else 'STRUCTURAL_AUDIT_ONLY',
        'errors': errors,
        'source_files_present': present,
        'totals': totals,
        'models': report,
        'interpretation': {
            'structural_registry_verified': not errors,
            'all_model_full_potential_rg_and_matching_verified': False,
            'reason': 'Quartic inventory does not establish CG normalisation, closed coupled RG flow, scalar mixing, threshold matching or physical C5 fitting.',
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path,
                        default=Path('output/audits/full_scalar_readiness.json'))
    args = parser.parse_args()
    result = audit(ROOT)
    output = args.output if args.output.is_absolute() else ROOT / args.output
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
    print('Status:', result['status'])
    print('Totals:', result.get('totals'))
    print('Errors:', result['errors'])
    print('Report:', output)
    return 2 if result['errors'] else 0


if __name__ == '__main__':
    raise SystemExit(main())
