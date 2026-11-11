"""Read-only audit of scalar quartics feeding the simultaneous T3 EFT1 Wilson RGE.

Project path: scripts/AuditT3SimultaneousQuarticReadiness.py

Runs without Wolfram; does not alter benchmark inputs or numerical source.
This checks *structural readiness*, not physical tensor normalization.
"""
from __future__ import annotations

import argparse
import ast
import json
from pathlib import Path
import re
import sys

ROOT = Path(__file__).resolve().parents[1]
REGISTRY = ROOT / 'RGE/general/generated/general_rgbeta_quartic_registry.json'
STAGE22 = ROOT / 'output/general_t3_stage22/stage22_summary.json'
ORDINARY = ROOT / 'scripts/RunT3SimultaneousEFT1All16.py'
MIXED = ROOT / 'scripts/RunT3SimultaneousEFT1MixedBC.py'
GENERAL_STATE = ROOT / 'Numerical/experimental/GeneralT3ScalarRGState.py'
REPORT = ROOT / 'output/general_t3_stage22/simultaneous_quartic_readiness.json'


def read_json(path: Path) -> dict:
    value = json.loads(path.read_text(encoding='utf-8-sig'))
    if not isinstance(value, dict):
        raise ValueError(f'{path}: JSON root must be an object')
    return value


def function_node(path: Path, name: str) -> ast.FunctionDef:
    tree = ast.parse(path.read_text(encoding='utf-8'), filename=str(path))
    funcs = [node for node in ast.walk(tree)
             if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name == name]
    if len(funcs) != 1:
        raise ValueError(f'{path}: expected exactly one function {name}, found {len(funcs)}')
    return funcs[0]


def names_in_function(path: Path, name: str) -> set[str]:
    node = function_node(path, name)
    symbols = set()
    for child in ast.walk(node):
        if isinstance(child, ast.Constant) and isinstance(child.value, str):
            if re.fullmatch(r'lambda[A-Za-z0-9_]+', child.value):
                symbols.add(child.value)
    return symbols


def source_has_call(path: Path, name: str) -> bool:
    tree = ast.parse(path.read_text(encoding='utf-8'), filename=str(path))
    return any(isinstance(x, ast.Call) and
               ((isinstance(x.func, ast.Name) and x.func.id == name) or
                (isinstance(x.func, ast.Attribute) and x.func.attr == name))
               for x in ast.walk(tree))


def audit(root: Path, registry_path: Path, stage22_path: Path) -> dict:
    ordinary = root / ORDINARY.relative_to(ROOT)
    mixed = root / MIXED.relative_to(ROOT)
    state_file = root / GENERAL_STATE.relative_to(ROOT)
    for file in (registry_path, ordinary, mixed, state_file):
        if not file.is_file():
            raise FileNotFoundError(file)
    registry = read_json(registry_path)
    models = registry.get('models', {})
    if len(models) != 16:
        raise ValueError(f'Expected 16 models, found {len(models)}')
    declared_baseline = names_in_function(ordinary, 'general_quartic')
    # Distinguish explicit state keys from the HHSS tensor, which is reconstructed separately.
    model_rows = {}
    for model, entry in sorted(models.items()):
        new = {x['name'] for x in entry['new_rgbeta_symbols']}
        reused = {x['existing_name'] for x in entry['preexisting_rgbeta_symbols']}
        if len(new) + len(reused) != entry['total_additional_quartics']:
            raise ValueError(f'{model}: inconsistent registry slot count')
        entries = [(x['name'], x['counts']) for x in entry['new_rgbeta_symbols']] + [
            (x['existing_name'], x['counts']) for x in entry['preexisting_rgbeta_symbols']]
        non_hhss = sorted(name for name, counts in entries if counts[0] + counts[1] != 2)
        hhss = sorted(name for name, counts in entries if counts[0] + counts[1] == 2)
        missing_extra = sorted((new | reused) - declared_baseline)
        model_rows[model] = {
            'additional_slots': entry['total_additional_quartics'],
            'new_symbols': sorted(new),
            'reused_symbols': sorted(reused),
            'not_explicit_in_general_quartic': missing_extra,
            'additional_hhss_symbols': hhss,
            'additional_non_hhss_symbols': non_hhss,
            'additional_non_hhss_not_explicit': sorted(set(non_hhss) - declared_baseline),
            'mixed': model in ('B_-1', 'C_-1'),
            'ready_for_complete_wilson_tensor': False,
        }

    stage22 = read_json(stage22_path) if stage22_path.is_file() else None
    normalized = bool(stage22 and stage22.get('normalisation_certified') is True and
                      stage22.get('certification_status') == 'ExactTensorAndExplicitPrefactorAlgebraPassed')
    # This audit never infers that a source-level assembly exists from a certificate alone.
    state_source = state_file.read_text(encoding='utf-8')
    ordinary_source = ordinary.read_text(encoding='utf-8')
    mixed_source = mixed.read_text(encoding='utf-8')
    state_support = ('pack_general_scalar_state' in state_source and
                     'unpack_general_scalar_state' in state_source)
    ordinary_extended = source_has_call(ordinary, 'pack_general_scalar_state')
    mixed_extended = ('unpack_general_scalar_state' in mixed_source and
                      'pack_general_scalar_state' in mixed_source)
    provisional = ('general_quartic' in ordinary_source and
                   source_has_call(ordinary, 'add_pair_contraction'))
    # A structural status cannot certify tensor correctness.  In particular, an
    # HHSS-dependent term may exist without the complete SSSS and HSSS blocks.
    totals = registry.get('totals', {})
    return {
        'format': 't3-simultaneous-quartic-structural-readiness-v1',
        'scope': 'read-only, not a physics or normalization certificate',
        'registry_totals': totals,
        'models': model_rows,
        'observed_implementation': {
            'ordinary_explicit_provisional_scalar_keys': sorted(declared_baseline),
            'ordinary_packs_extended_quartics': ordinary_extended,
            'mixed_uses_extended_scalar_state': mixed_extended,
            'extended_state_helpers_exist': state_support,
            'ordinary_calls_pair_contraction': provisional,
            'mixed_reuses_general_initial_context':
                'general_initial_context' in mixed_source,
        },
        'native_evidence': {
            'stage22_report_found': stage22 is not None,
            'stage22_status': stage22.get('certification_status') if stage22 else None,
            'stage22_normalisation_certified': normalized,
        },
        'blockers': [
            'Ordinary simultaneous IVP does not pack extended RGBeta quartic state'
            if not ordinary_extended else None,
            'Existing general_quartic uses provisional scalar-potential contraction assembly'
            if provisional else None,
            'Stage 22 exact native tensor normalization not certified'
            if not normalized else None,
            'No demonstrated map from all model-specific RGBeta quartics to the Wilson real-scalar tensor',
        ],
        'ready_for_numerical_replacement': False,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=ROOT)
    parser.add_argument('--registry', type=Path)
    parser.add_argument('--stage22', type=Path)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--require-ready', action='store_true',
                        help='Exit 2 if complete Wilson tensor readiness is not demonstrated')
    args = parser.parse_args(argv)
    root = args.root.resolve()
    output = args.output or root / REPORT.relative_to(ROOT)
    try:
        report = audit(root, args.registry or root / REGISTRY.relative_to(ROOT),
                       args.stage22 or root / STAGE22.relative_to(ROOT))
        report['blockers'] = [x for x in report['blockers'] if x]
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
        print('Structural scalar quartic audit (not certification)')
        print('Additional registry slots:', sum(x['additional_slots']
                                                for x in report['models'].values()))
        for model, row in report['models'].items():
            print(f"{model:5} extra={row['additional_slots']:2} "
                  f"non-HHSS-absent={len(row['additional_non_hhss_not_explicit']):2} "
                  f"{'MIXED' if row['mixed'] else 'ordinary'}")
        print('Stage 22:', report['native_evidence']['stage22_status'] or 'no report')
        print('Blockers:')
        for blocker in report['blockers']:
            print(' -', blocker)
        print('JSON:', output)
        return 2 if args.require_ready and not report['ready_for_numerical_replacement'] else 0
    except (OSError, ValueError, SyntaxError, KeyError, TypeError) as exc:
        print(f'ERROR: {exc}', file=sys.stderr)
        return 1


if __name__ == '__main__':
    raise SystemExit(main())
