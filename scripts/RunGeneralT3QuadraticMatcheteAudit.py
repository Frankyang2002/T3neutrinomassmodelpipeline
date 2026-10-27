"""Fail-closed Matchete general scalar audit with live-baseline-calibrated diagnostics.

For B_-1 and C_-1 only, the off-diagonal scalar bilinear is gauge allowed,
but Matchete's HeavyMassBasis diagnostic fails. We REQUIRE every *other*
Matchete diagnostic to pass and report the heavy-basis check as blocked.
Do NOT interpret this as validated scalar matching.

    python -m scripts.RunGeneralT3QuadraticMatcheteAudit --all
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess
from typing import Any

from scripts.AuditGeneralT3QuadraticBasis import audit, ROOT

RUNNER = ROOT / 'Lagrangian/RunGeneralQuadraticMatcheteAudit.wl'
OUT = ROOT / 'output/general_t3_quadratic_eft1/matchete'
MIXED_MODELS = frozenset({'B_-1', 'C_-1'})
HEAVY_PASS = 'MatcheteFullDetailedChecksPassed'
NON_MASS_PASS = 'MatcheteNonMassBasisDetailedChecksPassedHeavyMassBasisBlocked'


def _diag(data: Any, *, mixed: bool,
          baseline_opaque: dict[str, str] | None = None) -> dict[str, str]:
    """Validate actual live Matchete core checks and baseline-anchored metadata.

    Installed Matchete versions can expose non-Boolean diagnostic payloads
    such as Charges, GaugeAnomalies and ZeroModeMixing. A passing ordinary
    baseline anchors these values; no arbitrary True/{} coercion is permitted.
    """
    if not isinstance(data, dict):
        raise ValueError('Matchete detailed check absent')
    expectations = {
        'diagnostic_returned': True,
        'standard_check_boolean': True,
        'standard_check_passed': not mixed,
        'non_mass_basis_checks_passed': True,
        'heavy_mass_basis_reported': True,
        'heavy_mass_basis_passed': not mixed,
        'opaque_metadata_baseline_match': True,
        'failed_checks': [],
        'comparison_method': 'ExplicitCoreChecksPlusOrdinaryBaselineAnchoredOpaqueMetadata',
    }
    for key, want in expectations.items():
        got = data.get(key)
        if type(got) is not type(want) or got != want:
            raise ValueError(f'Matchete detailed check {key}: expected {want!r}, got {got!r}')
    checks = data.get('diagnostics')
    if not isinstance(checks, dict):
        raise ValueError('Matchete per-check summary absent')
    required = {
        'Open/Complex/InconsistentSpinChains': 'PassEmpty',
        'Hermiticity': 'PassTrue',
        'UncontractedIndices': 'PassEmpty',
        'CanonicallyNormalized': 'PassTrue',
        'HeavyTadpoles': 'PassEmpty',
        'AllObjectsDefined': 'PassTrue',
        'HeavyMassBasis': 'FailFalse' if mixed else 'PassTrue',
    }
    for name, expected in required.items():
        if checks.get(name) != expected:
            raise ValueError(f'Matchete {name} invalid or missing: {checks.get(name)!r}')
    if not any(checks.get(key) == 'PassTrue' for key in
               ('GaugeInvariant', 'FreeOfGaugeFields')):
        raise ValueError('Matchete gauge-invariance core check missing or failed')
    if checks.get('ChargeNeutral') != 'PassTrue' and 'Charges' not in checks:
        raise ValueError('Matchete neither verifies ChargeNeutral nor reports baseline-comparable Charges')
    if 'GaugeAnomalies' not in checks:
        raise ValueError('Matchete GaugeAnomalies diagnostic missing')

    opaque = data.get('opaque_metadata_values_inputform')
    if not isinstance(opaque, dict) or any(
        not isinstance(k, str) or not isinstance(v, str)
        for k, v in opaque.items()
    ):
        raise ValueError('Matchete opaque values absent or not preserved verbatim')
    if any(checks.get(k) != 'OpaqueMatchesValidBaseline' for k in opaque):
        raise ValueError('Matchete opaque diagnostic status differs from baseline')
    if checks.get('ChargeNeutral') != 'PassTrue' and 'Charges' not in opaque:
        raise ValueError('Matchete Charges diagnostic was not baseline-anchored')
    if 'GaugeAnomalies' not in opaque:
        raise ValueError('Matchete GaugeAnomalies was not baseline-anchored')
    if baseline_opaque is not None and opaque != baseline_opaque:
        raise ValueError('Matchete native opaque diagnostic values differ from verified baseline')
    return opaque


def _validate_report(model: str, data: Any, expected_quadratic_slots: int,
                     expected_quartics: int | None = None) -> dict:
    if not isinstance(data, dict):
        raise ValueError(f'{model}: Matchete output is not an object')
    mixed = model in MIXED_MODELS
    expected = {
        'status': NON_MASS_PASS if mixed else HEAVY_PASS,
        'quadratic_slots': expected_quadratic_slots,
        'baseline_heavy_lagrangian_valid': True,
        'combined_non_mass_basis_checks_passed': True,
        'baseline_calibrated_diagnostic_parser': True,
        'opaque_diagnostics_uninterpreted_but_baseline_identical': True,
        'combined_heavy_mass_basis_passed': not mixed,
        'checked_in_original_heavy_field_definitions': True,
        'operator_consistency_validated_except_mass_basis': mixed,
        'matchete_full_lagrangian_valid': not mixed,
        'mass_basis_diagonalized': False,
        'heavy_scalar_matching_mass_basis_usable': False,
        'mixed_heavy_mass_term_requires_diagonalization': mixed,
        'operator_normalization_certified': False,
        'rgbeta_physical_mapping_certified': False,
        'mass_eigenbasis_and_matching_done': False,
        'matching_executed': False,
    }
    for key, want in expected.items():
        got = data.get(key)
        if type(got) is not type(want) or got != want:
            raise ValueError(f'{model}: report {key}: expected {want!r}, got {got!r}')
    cls, alpha = model.split('_')
    if data.get('class') != cls or data.get('alpha') != int(alpha):
        raise ValueError(f'{model}: Matchete representation label mismatch')
    n = data.get('additional_quartics')
    if type(n) is not int or n < 0 or (expected_quartics is not None and n != expected_quartics):
        raise ValueError(f'{model}: quartic inventory inconsistent')
    baseline_opaque = _diag(data.get('baseline_diagnostic'), mixed=False)
    for kind, count, ismixed in (
        ('quartic', n, False), ('quadratic', expected_quadratic_slots, mixed)
    ):
        names = data.get(kind + '_operators')
        rows = data.get(kind + '_individual_checks')
        if not isinstance(names, list) or not isinstance(rows, list) or len(names) != count or len(rows) != count:
            raise ValueError(f'{model}: {kind} operator list incomplete')
        if len(set(names)) != count:
            raise ValueError(f'{model}: duplicate {kind} operator names')
        for name, row in zip(names, rows):
            if not isinstance(name, str) or not name or not isinstance(row, dict) or row.get('name') != name or row.get('valid') is not True:
                raise ValueError(f'{model}: {kind} operator validation failed')
            _diag(row.get('diagnostic'), mixed=ismixed, baseline_opaque=baseline_opaque)
    _diag(data.get('combined_diagnostic'), mixed=mixed, baseline_opaque=baseline_opaque)
    return data


def run(model: str, outdir: Path, executable: str) -> dict:
    basis = audit()
    if model not in basis['models']:
        raise ValueError('Unknown ordinary T3 model: ' + model)
    if not RUNNER.is_file():
        raise FileNotFoundError('Missing Matchete audit: ' + str(RUNNER))
    outdir.mkdir(parents=True, exist_ok=True)
    tag = model.replace('+', 'p').replace('-', 'm')
    path = outdir / f'{tag}_full_scalar_operator_audit.json'
    path.unlink(missing_ok=True)
    cls, alpha = model.split('_')
    token = f'm{abs(int(alpha))}' if int(alpha) < 0 else str(int(alpha))
    cmd = [executable, '-file', str(RUNNER), str(path), cls, token]
    try:
        proc = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True, check=False)
    except OSError as exc:
        raise RuntimeError(f'{model}: could not launch Matchete: {exc}') from exc
    log = outdir / f'{tag}_matchete.log'
    log.write_text(proc.stdout + '\n' + proc.stderr, encoding='utf-8')
    if proc.returncode or not path.is_file():
        detail = ''
        if path.is_file():
            try:
                failed = json.loads(path.read_text(encoding='utf-8'))
                stage = failed.get('status', 'Unknown')
                diag = (failed.get('baseline_diagnostic') or
                        failed.get('combined_diagnostic') or
                        failed.get('quadratic_individual_checks') or
                        failed.get('quartic_individual_checks'))
                detail = f'\nWolfram result status: {stage}\nDiagnostic: {diag!r}'
            except (OSError, ValueError, TypeError):
                detail = '\nWolfram wrote an unreadable JSON diagnostic.'
        raise RuntimeError(f'{model}: Matchete audit failed (exit={proc.returncode}); see {log}'
                           + detail + '\n' + proc.stdout[-3500:] + proc.stderr[-3500:])
    data = json.loads(path.read_text(encoding='utf-8'))
    record = basis['models'][model]
    slots = record['independent_complex_mass_parameters_in_declared_basis']
    # Manifest independent quartic slot count derives from the same source as
    # the independent quadratic audit; do not trust a status string alone.
    manifest_path = ROOT / 'docs/general_scalar_trial_operator_manifest.json'
    manifest = json.loads(manifest_path.read_text(encoding='utf-8'))
    quartics = manifest['models'][model]['extra_quartic_complex_couplings']
    return _validate_report(model, data, slots, quartics)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--model')
    group.add_argument('--all', action='store_true')
    parser.add_argument('--wolframscript', default='wolframscript')
    parser.add_argument('--output-dir', type=Path, default=OUT)
    args = parser.parse_args()
    basis = audit()
    names = sorted(basis['models']) if args.all else [args.model]
    rows = []
    for model in names:
        report = run(model, args.output_dir, args.wolframscript)
        rows.append(report)
        if model in MIXED_MODELS:
            print(f'{model}: all non-mass-basis Matchete checks PASSED; '
                  'HeavyMassBasis BLOCKED (scalar diagonalization required)', flush=True)
        else:
            print(f'{model}: full Matchete detailed checks PASSED; '
                  f'quartics={report["additional_quartics"]}, '
                  f'quadratics={report["quadratic_slots"]}', flush=True)
    if args.all:
        counts = {
            'models': len(rows),
            'quartics': sum(row['additional_quartics'] for row in rows),
            'quadratic_slots': sum(row['quadratic_slots'] for row in rows),
            'non_mass_basis_only_checks': sum(model in MIXED_MODELS for model in names),
            'full_matchete_checks': sum(model not in MIXED_MODELS for model in names),
        }
        required = {'models': 16, 'quartics': 50, 'quadratic_slots': 8,
                    'non_mass_basis_only_checks': 2, 'full_matchete_checks': 14}
        if counts != required:
            raise RuntimeError(f'Matchete 16-model coverage mismatch: {counts}')
        summary = {
            'status': 'OperatorAuditDetailedDiagnostics14Full2MassBasisBlocked',
            **counts,
            'physical_matchete_rgbeta_mapping_certified': False,
            'mixed_scalar_threshold_matching_done': False,
            'mass_eigenbasis_diagonalization_done': False,
            'model_statuses': {name: row['status'] for name, row in zip(names, rows)},
        }
        dest = args.output_dir / 'matchete_16_model_summary.json'
        dest.write_text(json.dumps(summary, indent=2) + '\n', encoding='utf-8')
        print('16/16 operator diagnostic audits: 14 full, 2 non-mass-basis only; '
              'mixed heavy-mass matching remains BLOCKED', flush=True)
        print('Summary:', dest)


if __name__ == '__main__':
    main()
