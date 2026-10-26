"""Stage 21: audit all 16 general-scalar UV/EFT1 RGBeta exports.

Read-only verification of the candidate quartic RGE export contract. This is
NOT a Matchete/RGBeta coefficient matching proof and does not run an ODE.

Run from the repository root: python -m scripts.RunGeneralT3Stage21
Optional diagnostic:          python -m scripts.RunGeneralT3Stage21 --numeric-probe
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import re
from typing import Any, Callable, Mapping

ROOT = Path(__file__).resolve().parents[1]
REGISTRY_PATH = ROOT / 'RGE/general/generated/general_rgbeta_quartic_registry.json'
EXPORT_DIR = ROOT / 'output/general_t3_stage20'
REPORT_DIR = ROOT / 'output/general_t3_stage21'
STAGES = ('UV', 'EFT1')
UNRESOLVED = re.compile(
    r'(?<![A-Za-z0-9_$])(BetaTerm|BetaFunction|QuarticTensors|UpsilonQuarticTensors|'
    r'CheckQuarticMixing|ProjectionToUnmixedBetas|Missing|Indeterminate|ComplexInfinity|'
    r'Infinity|Undefined|DirectedInfinity)(?=\s*\[|\b)|\$Failed|\$Aborted'
)
_CONVENTION = '16*pi^2*dX/dln(mu); gauge BetaTerm divided by 2*g'


def _require(condition: bool, text: str) -> None:
    if not condition:
        raise ValueError(text)


def _unique_strings(value: Any, context: str) -> list[str]:
    _require(isinstance(value, list), f'{context}: expected list')
    _require(all(isinstance(x, str) and x for x in value), f'{context}: invalid list element')
    _require(len(value) == len(set(value)), f'{context}: repeated coupling')
    return value


def _registry(path: Path) -> dict:
    payload = json.loads(path.read_text(encoding='utf-8-sig'))
    _require(payload.get('status') == 'GeneratedCandidateBasisNotWolframVerified',
             'Unexpected candidate registry status')
    models = payload.get('models')
    _require(isinstance(models, dict) and len(models) == 16,
             'Expected 16 distinct ordinary T3 model assignments')
    counts = {
        'matchete_extra_quartics': sum(int(v['total_additional_quartics']) for v in models.values()),
        'rgbeta_preexisting_structural': sum(len(v['preexisting_rgbeta_symbols']) for v in models.values()),
        'rgbeta_new_vertices': sum(len(v['new_rgbeta_symbols']) for v in models.values()),
        'quadratic_couplings_blocked': sum(int(v['quadratics_not_implemented']) for v in models.values()),
    }
    _require(counts == {
        'matchete_extra_quartics': 50,
        'rgbeta_preexisting_structural': 12,
        'rgbeta_new_vertices': 38,
        'quadratic_couplings_blocked': 8,
    }, f'Unexpected general scalar inventory: {counts}')
    _require(payload.get('totals') == counts, 'Registry total counters disagree with model inventory')
    return payload


def _audit_one(path: Path, model: str, stage: str, entry: dict) -> dict:
    _require(path.is_file(), f'Missing exported RGBeta JSON: {path}')
    payload = json.loads(path.read_text(encoding='utf-8-sig'))
    _require(payload.get('status') == 'Success', f'{path}: status is not Success')
    meta = payload.get('metadata')
    _require(isinstance(meta, dict), f'{path}: missing metadata')

    expected_meta = {
        'GeneralScalarModel': model,
        'GeneralScalarStage': stage,
        'GeneralScalarBasisStatus': 'CandidateNativeRGBetaUnmappedToMatchete',
        'GeneralScalarQuarticBetaBasis': 'RGBetaFullQuarticProjection',
        'ReportBetaConvention': _CONVENTION,
        'GeneralScalarQuadraticMassMixing': 'NotIncluded',
        'dS1': entry['dimensions'][0], 'dS2': entry['dimensions'][1],
        'dF': entry['dimensions'][2], 'alpha': entry['alpha'],
        'GeneralScalarAdditionalQuartics': entry['total_additional_quartics'],
        'GeneralScalarExistingStructuralCount': len(entry['preexisting_rgbeta_symbols']),
    }
    for key, expected in expected_meta.items():
        _require(meta.get(key) == expected,
                 f'{path}: {key} mismatch: expected {expected!r}; got {meta.get(key)!r}')

    added = [item['name'] for item in entry['new_rgbeta_symbols']]
    reused = [item['existing_name'] for item in entry['preexisting_rgbeta_symbols']]
    _require(len(added) == len(set(added)) and len(reused) == len(set(reused)),
             f'{model}: duplicate symbol registration')
    _require(not set(added).intersection(reused), f'{model}: new/reused coupling collision')
    _require(len(added) + len(reused) == entry['total_additional_quartics'],
             f'{model}: quartic inventory mismatch')
    _require(sorted(_unique_strings(meta.get('GeneralScalarNewBetaKeys'), f'{path}: new beta metadata')) == sorted(added),
             f'{path}: new RGBeta beta keys do not match registry')
    projected_old = _unique_strings(meta.get('GeneralScalarExistingQuarticBetasProjected'),
                                    f'{path}: projected old quartics')
    _require(set(reused).issubset(projected_old),
             f'{path}: reused quartics omitted from full projection: {sorted(set(reused)-set(projected_old))}')

    betas = payload.get('report_betas')
    raw = payload.get('betas')
    _require(isinstance(betas, dict) and isinstance(raw, dict),
             f'{path}: must export both betas and report_betas')
    _require(set(betas) == set(raw), f'{path}: raw/report beta keys differ')
    _require(all(isinstance(k, str) and k for k in betas), f'{path}: invalid beta key')
    _require(set(added + reused).issubset(betas),
             f'{path}: missing required general quartic betas: {sorted(set(added + reused)-set(betas))}')
    _require(set(projected_old).issubset(betas),
             f'{path}: projected baseline quartic missing from report betas')
    unresolved: list[str] = []
    for key, expression in betas.items():
        if not isinstance(expression, str) or not expression.strip():
            raise ValueError(f'{path}: empty or non-string beta {key}')
        if UNRESOLVED.search(expression):
            unresolved.append(key)
        if not isinstance(raw[key], str) or not raw[key].strip():
            raise ValueError(f'{path}: empty or non-string raw beta {key}')
    _require(not unresolved, f'{path}: unresolved RGBeta expressions: {unresolved}')
    _require(meta.get('NativeRGBetaCheckInvarianceUsed') is False,
             f'{path}: unsupported native gauge-check claim')

    return {
        'model': model, 'stage': stage, 'status': 'CandidateProjectedQuarticExportStructureVerified',
        'new_quartics': len(added), 'reused_quartics': len(reused),
        'projected_original_quartics': len(projected_old),
        'all_report_beta_keys': sorted(betas),
        'report_beta_count': len(betas),
        'matchete_normalisation_verified': False,
        'physical_matching_complete': False,
    }


def _smoke_env(betas: Mapping[str, str]) -> dict[str, Any]:
    """A deliberately artificial numerical environment for parser diagnostics.

    Using 3x3 test Yukawas is NOT a representation-correct flavor benchmark.
    """
    import numpy as np

    matrix = np.array([[0.1+0.002j, 0.003, 0.005],
                       [0.002, 0.11+0.003j, 0.007],
                       [0.001, 0.004, 0.09-0.006j]], dtype=complex)
    env: dict[str, Any] = {}
    text = '\n'.join(betas.values())
    for tok in set(re.findall(r'\b[A-Za-z_$][A-Za-z0-9_$]*\b', text)) | set(betas):
        env[tok] = 0.02+0.003j
    env.update({'gY':0.35,'g2':0.65,'g3':1.17,
                'yu':matrix,'yd':matrix*.13,'ye':matrix*.01,
                'y1':matrix*.16,'y2':matrix*.18, 'MF':np.eye(3,dtype=complex)*5000,
                'h':matrix*.16, 'I':1j, 'Pi':math.pi, 'E':math.e,
                'mS1Sq':1.6e6,'mS2Sq':2.0e6, 'mSSq':1.6e6,
                'lambdaH':0.13})
    return env


def _numerical_probe(export_path: Path, evaluator: Callable[..., Any]) -> dict:
    """Check InputForm evaluator compatibility; do not assert physical RG running."""
    import numpy as np
    payload = json.loads(export_path.read_text(encoding='utf-8-sig'))
    betas = payload['report_betas']
    env = _smoke_env(betas)
    errors: dict[str, str] = {}
    for key, expr in betas.items():
        try:
            val = evaluator(expr, env)
            arr = np.asarray(val, dtype=complex)
            if not np.all(np.isfinite(arr)):
                raise ValueError('Nonfinite synthetic evaluation')
        except Exception as exc:
            errors[key] = f'{type(exc).__name__}: {str(exc)[:260]}'
    return {
        'status': 'ExperimentalInputFormParserDiagnostic',
        'expression_count': len(betas),
        'parser_passed': len(betas)-len(errors),
        'parser_failed': len(errors),
        'errors': errors,
        'synthetic_inputs_only': True,
        'representation_correct_flavor_shapes': False,
    }


def audit_catalogue(*, registry_path: Path = REGISTRY_PATH,
                    export_dir: Path = EXPORT_DIR, output_dir: Path | None = REPORT_DIR,
                    numeric_probe: bool = False, evaluator: Callable[..., Any] | None = None) -> dict:
    registry = _registry(Path(registry_path))
    models = registry['models']
    export_dir = Path(export_dir)
    rows: list[dict] = []
    for model in sorted(models):
        entry = models[model]
        stage_rows = []
        for stage in STAGES:
            path = export_dir / f'{model}_{stage.lower()}_rgbeta.json'
            row = _audit_one(path, model, stage, entry)
            if numeric_probe:
                if evaluator is None:
                    try:
                        from Numerical.core.RGBetaEvaluator import evaluate_inputform_expression
                    except ImportError as exc:
                        raise RuntimeError('Production RGBeta InputForm evaluator unavailable') from exc
                    evaluator = evaluate_inputform_expression
                row['numerical_parser_probe'] = _numerical_probe(path, evaluator)
            stage_rows.append(row)
            rows.append(row)
        _require(set(stage_rows[0]['all_report_beta_keys']) - set(stage_rows[1]['all_report_beta_keys']) <= {'y1', 'y2', 'MF'},
                 f'{model}: unexpected UV-only beta keys')
        _require(set(stage_rows[1]['all_report_beta_keys']) - set(stage_rows[0]['all_report_beta_keys']) <= set(),
                 f'{model}: unexpected EFT1-only beta keys')

    totals = registry['totals']
    _require(len(rows) == 32, 'Expected 32 distinct UV/EFT1 beta exports')
    _require(sum(row['new_quartics'] for row in rows) == 76,
             'Expected 38 distinct new model-specific quartics per UV/EFT1 stage')
    _require(sum(row['reused_quartics'] for row in rows) == 24,
             'Expected 12 reused quartic slots per UV/EFT1 stage')
    parser_failures = sum(row.get('numerical_parser_probe',{}).get('parser_failed',0) for row in rows)
    report = {
        'status': 'Stage21CandidateQuarticExportStructureVerified',
        'model_count': len(models), 'export_count': len(rows),
        'inventory': totals, 'stage_export_counts': {'UV':16,'EFT1':16},
        'report_total_beta_entries': sum(row['report_beta_count'] for row in rows),
        'numeric_probe_requested': numeric_probe,
        'numeric_parser_failures': parser_failures if numeric_probe else None,
        'numeric_probe_status': ('ParserDiagnosticsRequireReview' if parser_failures else 'SyntheticParserDiagnosticPassed') if numeric_probe else 'NotRequested',
        'matchete_rgbeta_normalisation_certified': False,
        'quadratic_mass_mixing_included': False,
        'end_to_end_16_model_numerics_validated': False,
        'rows': rows,
    }
    if output_dir is not None:
        path = Path(output_dir)
        path.mkdir(parents=True, exist_ok=True)
        (path/'stage21_summary.json').write_text(json.dumps(report, indent=2)+'\n',encoding='utf-8')
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--numeric-probe', action='store_true',
                        help='Use synthetic 3x3 matrix inputs to diagnose the current InputForm parser (not physics)')
    parser.add_argument('--registry', type=Path, default=REGISTRY_PATH)
    parser.add_argument('--exports', type=Path, default=EXPORT_DIR)
    parser.add_argument('--output', type=Path, default=REPORT_DIR)
    args = parser.parse_args()
    report = audit_catalogue(registry_path=args.registry,export_dir=args.exports,
                             output_dir=args.output,numeric_probe=args.numeric_probe)
    print('Stage21:',report['status'])
    print('Models:', report['model_count'], 'UV/EFT1 beta exports:', report['export_count'])
    print('Additional quartics:',report['inventory']['matchete_extra_quartics'],
          '(38 new, 12 reused across distinct models)')
    if args.numeric_probe:
        print('Synthetic parser evaluation failures:', report['numeric_parser_failures'])
        if report['numeric_parser_failures']:
            print('Inspect output/general_t3_stage21/stage21_summary.json for expression errors.')
    print('Matchete/RGBeta normalisation: NOT CERTIFIED')
    print('Report: output/general_t3_stage21/stage21_summary.json')
    return 0 if not report['numeric_parser_failures'] else 2


if __name__ == '__main__':
    raise SystemExit(main())
