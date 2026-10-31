"""Fail-closed verifier for live Matchete B/C 110101 action-term evidence.

This verifies native term construction and CheckLagrangian in the original
scalar basis. It does not compare Matchete with RGBeta or certify matching.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

MODELS = ('B_-1', 'C_-1')
FORMAT = 't3-bc-live-matchete-110101-action-terms-v1'
SOURCE = 'Matchete live T3ExtendedBuildExtraQuarticCandidates and CheckLagrangian'
EXPECTED = (
    ('T3ExtraQ_1_1_0_1_0_1_1', 'lambdaT3ExtraQ110101I1'),
    ('T3ExtraQ_1_1_0_1_0_1_2', 'lambdaT3ExtraQ110101I2'),
)
SOURCE_MARKERS = (
    '"Class" -> If[degree === 4, "ExtraQuartic", "ExtraQuadratic"]',
    '"Expression" -> PlusHc[contribution]',
    'contribution = coupling[] product If[cg === None, 1, CG[cg, labels]]',
)


def require(ok: bool, msg: str) -> None:
    if not ok:
        raise ValueError('Matchete B/C 110101 action evidence: ' + msg)


def verify_source(path: Path) -> None:
    raw = path.read_text(encoding='utf-8-sig')
    compact = ''.join(raw.split())
    for marker in SOURCE_MARKERS:
        require(''.join(marker.split()) in compact,
                'Matchete source-plus-h.c. convention changed: ' + marker)


def verify(doc: dict[str, Any]) -> dict[str, Any]:
    require(isinstance(doc, dict) and doc.get('format') == FORMAT, 'unexpected format')
    require(doc.get('status') == 'NativeMatcheteInteractionConstructionDiagnostic',
            'unexpected stage status')
    require(doc.get('source') == SOURCE, 'unexpected evidence source')
    for flag in ('matchete_scalar_off_diagonal_mass_included',
                 'cross_tool_action_equality_verified',
                 'physical_coupling_mapping_certified', 'physical_matching_certified'):
        require(doc.get(flag) is False, 'unexpected claim: ' + flag)
    models = doc.get('models')
    require(isinstance(models, dict) and set(models) == set(MODELS),
            'expected precisely B_-1 and C_-1')
    results = []
    for name in MODELS:
        rec = models[name]
        require(isinstance(rec, dict) and rec.get('model') == name,
                name + ': missing/incorrect model')
        require(rec.get('status') == 'LiveMatchete110101ActionTermsValidated' and
                rec.get('matchete_action_validated') is True and
                rec.get('native_model_built') is True,
                name + ': live Matchete action not validated')
        require(rec.get('scalar_extra_quartics_count') == 7 and
                rec.get('selected_110101_couplings_count') == 2 and
                rec.get('combined_check_lagrangian_passed') is True and
                rec.get('combined_check_lagrangian_return_inputform') == 'True',
                name + ': incomplete quartic inventory or failed combined check')
        for flag in ('cross_tool_action_equality_verified',
                     'physical_coupling_mapping_certified', 'physical_matching_certified'):
            require(rec.get(flag) is False, name + ': unexpected physical claim')
        require(rec.get('source_operator_uses_plus_hc') ==
                'verified separately from GeneralScalarOperators.wl',
                name + ': undocumented plus-h.c. source')
        terms = rec.get('selected_terms')
        require(isinstance(terms, list) and len(terms) == 2,
                name + ': expected two native terms')
        for row, (term_name, coupling) in zip(terms, EXPECTED):
            require(isinstance(row, dict) and row.get('name') == term_name and
                    row.get('coupling') == coupling,
                    name + ': incorrect native coupling name')
            require(row.get('coupling_present_in_native_expression') is True and
                    row.get('check_lagrangian_passed') is True and
                    row.get('check_lagrangian_return_inputform') == 'True',
                    name + '/' + coupling + ': failed native single-term validation')
            expr = row.get('native_expression_inputform')
            require(isinstance(expr, str) and coupling in expr and len(expr) < 24001,
                    name + '/' + coupling + ': missing native expression')
        results.append({'model': name, 'validated_terms': 2,
                        'individual_checks_passed': 2, 'combined_check_passed': True})
    return {
        'format': 't3-bc-matchete-110101-action-term-verification-v1',
        'status': 'LiveMatchete110101ActionTermsVerifiedCrossToolEqualityPending',
        'results': results,
        'matchete_native_terms_validated': 4,
        'matchete_combined_checks_passed': 2,
        'cross_tool_action_equality_verified': False,
        'physical_coupling_mapping_certified': False,
        'physical_matching_certified': False,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--input', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--matchete-source', type=Path,
                    default=Path('Lagrangian/interactions/GeneralScalarOperators.wl'))
    args = ap.parse_args()
    require(args.input.is_file(), 'missing probe output: ' + str(args.input))
    require(not args.output.exists(), 'refusing to overwrite ' + str(args.output))
    verify_source(args.matchete_source)
    data = json.loads(args.input.read_text(encoding='utf-8-sig'))
    result = verify(data)
    result['input_file'] = args.input.as_posix()
    result['input_sha256'] = hashlib.sha256(args.input.read_bytes()).hexdigest()
    result['matchete_operator_source_sha256'] = hashlib.sha256(
        args.matchete_source.read_bytes()).hexdigest()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n',
                           encoding='utf-8')
    for row in result['results']:
        print(f"{row['model']}: live Matchete 110101 individual 2/2; combined CheckLagrangian=True")
    print('Cross-tool action/potential equality and physical matching: NOT CERTIFIED')
    print('Summary:', args.output)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
