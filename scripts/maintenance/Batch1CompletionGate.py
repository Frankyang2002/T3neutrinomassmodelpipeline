"""Read-only T3 honours numerical completion gate. Run from repository root.

This is an execution/provenance gate, not proof of physical equivalence between
RGBeta and Matchete. Does not alter model code or scientific outputs.
"""
from __future__ import annotations
import argparse
import ast
import datetime as dt
import importlib
import json
import os
from pathlib import Path
import subprocess
import sys
import traceback

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'output' / 'audit' / 'batch1'
TARGETS = {
    'general_uv': 'Numerical/experimental/T3GeneralYukawaCoupledUV.py',
    'extended_c5': 'Numerical/experimental/T3RGAssistedExtendedC5.py',
    'low_energy': 'Numerical/experimental/T3ExtendedLowEnergyStudy.py',
    'c5_observables': 'Numerical/experimental/T3ExtendedC5Observables.py',
    'production': 'Numerical/orchestration/PipelineNumericalResults.py',
    'coupled_b': 'Numerical/orchestration/T3BIntegratedCoupledPipeline.py',
    'mass_observables': 'physics/NeutrinoObservables.py',
    'charged_basis': 'physics/NeutrinoTrajectory.py',
    'eft1_wilson': 'Numerical/experimental/T3BCCoupledWilsonTensor.py',
    'general_state': 'Numerical/experimental/T3GeneralMatrixBeta.py',
}
FOCUS_TESTS = [
    'tests/python/test_general_t3_eft1_hard_pole_audit.py',
    'tests/python/test_general_t3_quartic_action_identity.py',
    'tests/python/test_general_t3_triplet_intertwiner.py',
    'tests/python/test_general_t3_triplet_native_structure.py',
    'tests/python/test_general_t3_coupled_uv_running.py',
    'tests/python/test_general_t3_eft1_scalar_pair_mixing.py',
]


def parse_source(path):
    source = path.read_text(encoding='utf-8-sig')
    tree = ast.parse(source, filename=str(path))
    imports, funcs, constants = [], [], {}
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom):
            imports.append(node.module or '')
        elif isinstance(node, ast.Import):
            imports.extend(x.name for x in node.names)
        elif isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            funcs.append(node.name)
    for node in tree.body:
        if isinstance(node, (ast.Assign, ast.AnnAssign)):
            names = [x.id for x in node.targets if isinstance(x, ast.Name)] if isinstance(node, ast.Assign) else ([node.target.id] if isinstance(node.target, ast.Name) else [])
            value = node.value
            for name in names:
                if name in {'SPECS', 'EXTRAS', 'MODELS', 'STATUS', 'SCHEMA', 'OUTPUT', 'OUT'}:
                    try:
                        v = ast.literal_eval(value)
                        if isinstance(v, dict):
                            constants[name] = {'type': 'dictionary', 'keys': list(map(str, v.keys()))}
                        elif isinstance(v, (tuple, list, str, bool, int, float)):
                            constants[name] = v
                    except (ValueError, TypeError, MemoryError):
                        pass
    return {'functions': sorted(set(funcs)), 'imports': sorted(set(imports)), 'constants': constants,
            'lines': len(source.splitlines()), 'source_has_fixed_reference': 'fitted_reference' in source,
            'source_uses_solve_ivp': 'solve_ivp(' in source}


def run_command(label, args, timeout):
    try:
        cp = subprocess.run(args, cwd=ROOT, capture_output=True, text=True,
                            encoding='utf-8', errors='replace', timeout=timeout)
        return {'label': label, 'command': args, 'returncode': cp.returncode,
                'stdout_tail': cp.stdout[-14000:], 'stderr_tail': cp.stderr[-5000:]}
    except subprocess.TimeoutExpired:
        return {'label': label, 'command': args, 'error': f'timeout ({timeout}s)'}
    except Exception as exc:
        return {'label': label, 'command': args, 'error': repr(exc)}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run-tests', action='store_true')
    p.add_argument('--run-self-tests', action='store_true', help='Call *only* self_test() in selected numeric modules, if callable without arguments')
    p.add_argument('--timeout', type=int, default=300, help='Timeout per test command')
    p.add_argument('--output-dir', type=Path, default=OUT)
    ns = p.parse_args()
    output = ns.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    report = {'timestamp_utc': dt.datetime.now(dt.timezone.utc).isoformat(),
              'repo_root': str(ROOT), 'scope': 'Executable source gates and targeted tests; no physical-basis certification',
              'components': {}, 'test_runs': [], 'self_tests': [], 'explicit_unverified': []}
    for label, rel in TARGETS.items():
        path = ROOT / rel
        entry = {'path': rel, 'exists': path.is_file()}
        if path.is_file():
            try: entry.update(parse_source(path))
            except Exception as exc: entry['parse_error'] = str(exc)
        report['components'][label] = entry
    report['explicit_unverified'] = [
        'Whether the extended-C5 bridge uses a physically derived RGBeta-to-Matchete coupling basis',
        'Whether the six-model one-scale path covers the input representation',
        'Whether one full state is passed continuously from UV to final mass and PMNS',
        'Whether all heavy thresholds are matched consistently in hierarchical cases',
        'Numerical benchmark validation against independent known limits',
    ]
    if ns.run_tests:
        for rel in FOCUS_TESTS:
            if (ROOT / rel).is_file():
                report['test_runs'].append(run_command(rel,
                    [sys.executable, '-m', 'pytest', '-q', '-rA', rel], ns.timeout))
            else: report['test_runs'].append({'label': rel, 'missing': True})
    if ns.run_self_tests:
        for label in ('general_uv', 'extended_c5', 'low_energy', 'c5_observables'):
            obj = report['components'][label]
            if not obj['exists'] or 'self_test' not in obj.get('functions', []):
                report['self_tests'].append({'module':label, 'status':'no zero-argument self_test discovered'})
                continue
            module = TARGETS[label][:-3].replace('/', '.')
            # Run the source-owned self-test as a separate process, protecting this audit.
            code = ('import importlib; m=importlib.import_module(' + repr(module) + '); m.self_test()')
            report['self_tests'].append(run_command(module, [sys.executable, '-c', code], ns.timeout))
    failures = [e['label'] for e in report['test_runs'] if e.get('returncode') not in (None, 0)]
    failures += [e['label'] for e in report['self_tests'] if e.get('returncode') not in (None, 0)]
    report['failing_commands'] = failures
    report['decision'] = 'NOT_CERTIFIED' if failures else 'SOURCE_AND_SELECTED_CHECKS_ONLY_NOT_PHYSICS_CERTIFIED'
    json_path = output / 'batch1_execution_report.json'
    json_path.write_text(json.dumps(report, indent=2, default=str) + '\n', encoding='utf-8')
    rows = ['# Batch 1 execution gate', '', 'No source files modified.', '', '| Component | Present |', '|---|---|']
    for key, item in report['components'].items():
        rows.append(f"| `{key}` | {'yes' if item['exists'] else 'NO'} |")
    rows += ['', '## Focused test result', '']
    if not ns.run_tests: rows.append('Not run; pass `--run-tests`.')
    for r in report['test_runs']:
        rows.append(f"- `{r['label']}`: return code {r.get('returncode', 'missing/error')}")
    rows += ['', '## Scientific certification', '', 'Not certified by this gate. Inspect JSON for explicit gaps.', '']
    (output / 'batch1_execution_summary.md').write_text('\n'.join(rows), encoding='utf-8')
    print('Report:', json_path)
    print('Summary:', output / 'batch1_execution_summary.md')
    print('Failing selected test commands:', len(failures))
    return 1 if failures else 0


if __name__ == '__main__':
    raise SystemExit(main())
