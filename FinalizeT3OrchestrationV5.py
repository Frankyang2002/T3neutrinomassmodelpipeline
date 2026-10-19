r"""Finish an existing T3 orchestration move using the files actually on disk.

Run at repository root:
    python .\FinalizeT3OrchestrationV5.py --plan
    python .\FinalizeT3OrchestrationV5.py --apply
    python .\FinalizeT3OrchestrationV5.py --restore

No organiser manifest/backup is required. Four recovered originals are supplied beside this script under .t3_orchestration_recovery/. Changes are transactional where possible:
all source files to be changed are copied to a local backup before writing;
pytest failure triggers restoration. No original implementation is deleted.
"""
from __future__ import annotations

import argparse
import ast
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys

BASE = Path('Numerical/orchestration')
STATE = Path('.t3_orchestration_v5_state.json')
BACKUP = Path('.t3_orchestration_v5_backup')
RECOVERY = Path('.t3_orchestration_recovery')
MARKER = '# AUTO-GENERATED T3 ORCHESTRATION FORWARDER'
GROUPS = {
    'presentation': 'BuildThreePanelInteractive CollectChapter4Interactive ExportChapter4Figures',
    'studies/common_inputs': 'OptimizeCommonYukawas RunCommonInputCouplingScan RunCommonInputEnvelopes RunCommonInputFermionScale RunCommonInputRGE RunCommonParameterComparisonV2 RunFourBenchmarkCorrelations RunScalarCouplingCorrelations',
    'studies/fermion_scales': 'RunFittedFermionScaleScan RunHeavyFermionScaleComparison',
    'studies/fitted_rge': 'RunFittedRGEDisplay RunFittedSensitivityComparison RunFittedStagedRGComparison RunFittedWeinbergRunningComparison',
    'studies/scalar_rge': 'RunCommonScalarBetaDominance RunCommonScalarComparison RunCommonScalarRGControls',
    'support': 'ResumePipelineNumericalResults',
    'archive_candidates': 'BuildComparisonReports CollectThesisFigures ModelReportLayout RunCommonParameterComparison',
}
SKIP_DIRS = {'.git', '.venv', 'venv', '__pycache__', '.pytest_cache', 'output', 'Reports', 'display', 'dist', 'build', '.mypy_cache', '.ruff_cache', '.t3_orchestration_recovery'}


def root() -> Path:
    cwd = Path.cwd().resolve()
    if not (cwd / BASE / 'PipelineNumericalResults.py').is_file():
        raise SystemExit('Run from repository root (containing Numerical/orchestration).')
    return cwd


def files(rootdir: Path):
    for p in rootdir.rglob('*.py'):
        if any(x in SKIP_DIRS or x.startswith('.orchestration_') or x.startswith('.t3_orchestration_') for x in p.relative_to(rootdir).parts[:-1]):
            continue
        if p.name in {'FinalizeT3OrchestrationV5.py', 'FinalizeT3OrchestrationV3.py', 'FinalizeT3OrchestrationV2.py', 'FinalizeT3Orchestration.py', 'OrganizeT3Orchestration.py'}:
            continue
        yield p


def plan(rootdir: Path):
    mapping = {name: group for group, names in GROUPS.items() for name in names.split()}
    valid, unmoved, broken, missing = [], [], [], []
    for name, group in mapping.items():
        old = rootdir / BASE / f'{name}.py'
        new = rootdir / BASE / group / f'{name}.py'
        old_text = old.read_text(encoding='utf-8-sig') if old.is_file() else None
        new_text = new.read_text(encoding='utf-8-sig') if new.is_file() else None
        if new_text is not None and MARKER not in new_text:
            ast.parse(new_text, filename=str(new))
            valid.append((name, old, new, old_text))
        elif old_text is not None and MARKER not in old_text:
            ast.parse(old_text, filename=str(old))
            unmoved.append((name, old, new))
        elif old_text is not None and MARKER in old_text:
            recovered = rootdir / RECOVERY / f'{name}.py'
            if recovered.is_file():
                recovered_text = recovered.read_text(encoding='utf-8-sig')
                if MARKER in recovered_text:
                    raise ValueError(f'Recovery source is only a forwarder: {recovered}')
                ast.parse(recovered_text, filename=str(recovered))
                unmoved.append((name, old, new))
            else:
                broken.append((name, old, new))
        else:
            missing.append(name)
    replacements = {f'Numerical.orchestration.{name}': 'Numerical.orchestration.' + group.replace('/', '.') + f'.{name}' for name, group in mapping.items() if name in {r[0] for r in valid + unmoved}}
    # A forwarder without the actual implementation cannot be removed or repaired
    # from a file tree alone. No attempt is made to reconstruct lost code.
    changed = {}
    for path in files(rootdir):
        current = path.read_text(encoding='utf-8-sig')
        updated = current
        for before, after in replacements.items():
            updated = re.sub(r'(?<![\w.])' + re.escape(before) + r'(?=\b)', after, updated)
        if updated != current:
            ast.parse(updated, filename=str(path))
            changed[path] = updated
    return valid, unmoved, broken, missing, changed, replacements


def show(rootdir, p):
    valid, unmoved, broken, missing, changed, _ = p
    print(f'Already relocated implementations: {len(valid)}')
    print(f'Original implementations awaiting move: {len(unmoved)}')
    for name, _, new in unmoved:
        print('  MOVE', name, '->', new.relative_to(rootdir))
    print(f'Old forwarding files eligible for removal: {sum(x[3] is not None and MARKER in x[3] for x in valid)}')
    print(f'Files with Python references to rewrite: {len(changed)}')
    print(f'Missing implementations (only broken forwarders): {len(broken)}')
    for name, _, new in broken:
        print('  WARNING: Cannot recover', name, 'because', new.relative_to(rootdir), 'does not exist; leaving forwarder untouched')
    if missing:
        print('Missing both old and new:', ', '.join(missing))


def restore(rootdir):
    statefile = rootdir / STATE
    if not statefile.is_file():
        raise SystemExit('No V5 state found to restore.')
    state = json.loads(statefile.read_text(encoding='utf-8'))
    backup = rootdir / BACKUP
    for rel in state['original_files']:
        original = backup / rel
        if not original.is_file():
            raise SystemExit('Backup missing; refusing incomplete restore: ' + rel)
    for rel in state['created_files']:
        target = rootdir / rel
        if target.is_file():
            target.unlink()
    for rel in state['original_files']:
        target = rootdir / rel
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(backup / rel, target)
    statefile.unlink()
    print('Restored all V5-edited files. Backup left at', BACKUP)


def apply(rootdir, p, pytest_args):
    valid, unmoved, broken, missing, changed, replacements = p
    if broken or missing:
        raise SystemExit('Cannot proceed: source implementations missing: ' + ', '.join(x[0] for x in broken) + ', '.join(missing))
    if (rootdir / STATE).exists():
        raise SystemExit('Existing V5 state found; restore first.')
    if (rootdir / BACKUP).exists():
        raise SystemExit('Existing V5 backup directory found; rename it before applying.')
    # Back up the original files and prepare rollback state BEFORE modifying anything.
    remove_forwarders = [old for _, old, _, text in valid if text is not None and MARKER in text]
    moved_source = [old for _, old, _ in unmoved]
    modified = sorted(set(changed) | set(remove_forwarders) | set(moved_source))
    original_rels = [p.relative_to(rootdir).as_posix() for p in modified]
    created = [new for _, _, new in unmoved]
    created_rels = [p.relative_to(rootdir).as_posix() for p in created]
    backup = rootdir / BACKUP
    backup.mkdir()
    for path, rel in zip(modified, original_rels):
        dest = backup / rel
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, dest)
    (rootdir / STATE).write_text(json.dumps({'original_files': original_rels, 'created_files': created_rels}, indent=2), encoding='utf-8')
    try:
        # Relocate untouched originals and adjust only known __file__.parents[N] expressions.
        for name, old, new in unmoved:
            source = old.read_text(encoding='utf-8-sig')
            if MARKER in source:
                recovered = rootdir / RECOVERY / f'{name}.py'
                source = recovered.read_text(encoding='utf-8-sig')
            depth = len(new.relative_to(rootdir / BASE).parts) - 1
            source = re.sub(r'(Path\(\s*__file__\s*\)\.resolve\(\)\.parents\[)(\d+)(\])', lambda m: m.group(1) + str(int(m.group(2)) + depth) + m.group(3), source)
            for before, after in replacements.items():
                source = re.sub(r'(?<![\w.])' + re.escape(before) + r'(?=\b)', after, source)
            ast.parse(source, filename=str(new))
            new.parent.mkdir(parents=True, exist_ok=True)
            (new.parent / '__init__.py').touch(exist_ok=True)
            new.write_text(source, encoding='utf-8')
            old.unlink()
        for path, text in changed.items():
            # An originally unmoved file is now at its new location.
            if path in moved_source or path in remove_forwarders:
                continue
            path.write_text(text, encoding='utf-8')
        for old in remove_forwarders:
            old.unlink()
        print('Changes applied. Running pytest:', sys.executable, '-m pytest', *pytest_args, flush=True)
        result = subprocess.run([sys.executable, '-m', 'pytest', *pytest_args], cwd=rootdir)
        if result.returncode:
            print('pytest failed (exit code', result.returncode, '); restoring source files.', flush=True)
            restore(rootdir)
            raise SystemExit(result.returncode)
        print('pytest passed. Run `python -c "import pipeline"` as an additional smoke check.')
        print('To undo: python FinalizeT3OrchestrationV5.py --restore')
    except BaseException:
        if (rootdir / STATE).exists():
            restore(rootdir)
        raise


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument('--plan', action='store_true')
    mode.add_argument('--apply', action='store_true')
    mode.add_argument('--restore', action='store_true')
    ap.add_argument('--pytest-args', nargs='*', default=['-q'])
    args = ap.parse_args()
    r = root()
    if args.restore:
        restore(r)
        return
    p = plan(r)
    show(r, p)
    if args.apply:
        apply(r, p, args.pytest_args)


if __name__ == '__main__':
    main()
