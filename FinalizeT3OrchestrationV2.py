"""Finalise a partially applied T3 orchestration reorganisation.

Run from project root:
  python .\FinalizeT3OrchestrationV2.py --plan
  python .\FinalizeT3OrchestrationV2.py --apply
  python .\FinalizeT3OrchestrationV2.py --restore

Only removes verifiable organiser-generated forwarders that have a real moved
implementation. Unmoved original files (including archive candidates) remain.
Rewrites Python module references, runs full pytest, and rolls back on failure.
No GitHub writes; no permanent source deletion. Retains backup after success.
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

ORIG = '.orchestration_organizer_manifest.json'
STATE = '.orchestration_finalizer_v2_manifest.json'
BACKUP = '.orchestration_finalizer_v2_backup'
MARKER = '# AUTO-GENERATED T3 ORCHESTRATION FORWARDER'
SKIP = {'.git', '.venv', 'venv', '__pycache__', '.pytest_cache', 'output',
        'Reports', 'display', 'dist', 'build', '.mypy_cache', '.ruff_cache',
        '.orchestration_organizer_backup', '.orchestration_finalizer_backup', BACKUP}


def root_dir():
    root = Path.cwd().resolve()
    if not (root / 'Numerical' / 'orchestration' / 'PipelineNumericalResults.py').is_file():
        raise SystemExit('Run from the project root containing Numerical/orchestration.')
    return root


def read_operations(root):
    path = root / ORIG
    if not path.is_file():
        raise SystemExit('Missing organiser manifest: ' + ORIG)
    ops = json.loads(path.read_text(encoding='utf-8'))['operations']
    relocated, untouched, missing = [], [], []
    for op in ops:
        old, new = root / op['old'], root / op['new']
        if new.is_file():
            if old.is_file():
                data = old.read_text(encoding='utf-8-sig')
                if MARKER not in data:
                    raise SystemExit('Refusing to delete non-forwarder at ' + op['old'])
            elif old.exists():
                raise SystemExit('Old path is not a regular file: ' + op['old'])
            relocated.append(op)
        elif old.is_file():
            if MARKER in old.read_text(encoding='utf-8-sig'):
                raise SystemExit('Forwarder has no moved implementation: ' + op['new'])
            untouched.append(op)
        else:
            missing.append(op)
    if missing:
        raise SystemExit('Missing original and moved implementation for: ' +
                         ', '.join(o['name'] for o in missing))
    return relocated, untouched


def py_sources(root):
    for file in root.rglob('*.py'):
        if any(part in SKIP for part in file.relative_to(root).parts):
            continue
        if file.name in {'OrganizeT3Orchestration.py', 'FinalizeT3Orchestration.py',
                         'FinalizeT3OrchestrationV2.py'}:
            continue
        yield file


def changes(root, relocated):
    substitutions = {
        '.'.join(Path(op['old']).with_suffix('').parts):
        '.'.join(Path(op['new']).with_suffix('').parts)
        for op in relocated
    }
    old_files = {root / op['old'] for op in relocated}
    edits = []
    for path in py_sources(root):
        if path in old_files:
            continue
        current = path.read_text(encoding='utf-8-sig')
        updated = current
        # Replace whole module identifiers only (works for imports and dotted strings).
        for old, new in sorted(substitutions.items(), key=lambda t: -len(t[0])):
            updated = re.sub(r'(?<![A-Za-z0-9_.])' + re.escape(old) +
                             r'(?![A-Za-z0-9_])', new, updated)
        if updated != current:
            ast.parse(updated, filename=str(path))
            edits.append((path, updated))
    return edits


def restore(root):
    statepath, backup = root / STATE, root / BACKUP
    if not statepath.is_file() or not backup.is_dir():
        raise SystemExit('No V2 finalisation backup to restore.')
    state = json.loads(statepath.read_text(encoding='utf-8'))
    for name in state['files']:
        if not (backup / name).is_file():
            raise SystemExit('Missing backup: ' + name)
    for name in state['files']:
        dst = root / name
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(backup / name, dst)
    shutil.rmtree(backup)
    statepath.unlink()
    print('Restored all original forwarders and edited source files. Moved files retained.')


def apply(root, pytest_args):
    if (root / STATE).exists() or (root / BACKUP).exists():
        raise SystemExit('V2 finalisation state exists; use --restore first.')
    relocated, untouched = read_operations(root)
    edits = changes(root, relocated)
    forwarders = [root / op['old'] for op in relocated if (root / op['old']).is_file()]
    touched = sorted(set(forwarders + [p for p, _ in edits]))
    backup = root / BACKUP
    backup.mkdir()
    try:
        for path in touched:
            dst = backup / path.relative_to(root)
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(path, dst)
        (root / STATE).write_text(json.dumps({
            'files': [str(p.relative_to(root)).replace('\\', '/') for p in touched],
            'removed_forwarders': [str(p.relative_to(root)).replace('\\', '/') for p in forwarders],
            'unmoved_preserved': [o['old'] for o in untouched],
        }, indent=2) + '\n', encoding='utf-8')
        for path, data in edits:
            path.write_text(data, encoding='utf-8')
        for path in forwarders:
            path.unlink()
        print(f'Removed {len(forwarders)} moved-file forwarders; rewritten {len(edits)} Python files.', flush=True)
        print(f'Preserved {len(untouched)} never-moved originals.', flush=True)
        # Avoid accepting stale bytecode as evidence of a working old import.
        print('Running:', sys.executable, '-m pytest', *pytest_args, flush=True)
        result = subprocess.run([sys.executable, '-m', 'pytest', *pytest_args], cwd=root)
        if result.returncode:
            print('pytest failed; restoring before finalisation.', file=sys.stderr)
            restore(root)
            return result.returncode
        print('pytest passed. Backup retained. Use --restore to undo this finalisation.')
        return 0
    except BaseException:
        if (root / STATE).is_file():
            restore(root)
        else:
            shutil.rmtree(backup, ignore_errors=True)
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument('--plan', action='store_true')
    action.add_argument('--apply', action='store_true')
    action.add_argument('--restore', action='store_true')
    parser.add_argument('--pytest-arg', action='append', default=[],
                        help='Extra pytest argument; repeat if necessary (default: -q)')
    args = parser.parse_args()
    root = root_dir()
    if args.restore:
        restore(root)
        return 0
    relocated, untouched = read_operations(root)
    edits = changes(root, relocated)
    if args.plan:
        print(f'Moved implementations found: {len(relocated)}')
        print(f'Unmoved originals preserved: {len(untouched)}')
        for op in untouched:
            print('  KEEP', op['old'], '(not relocated)')
        for op in relocated:
            if (root / op['old']).is_file():
                print('  REMOVE FORWARDER', op['old'])
        print(f'Python files to update: {len(edits)}')
        print('pytest:', sys.executable, '-m pytest', *(args.pytest_arg or ['-q']))
        print('No files changed.')
        return 0
    return apply(root, args.pytest_arg or ['-q'])


if __name__ == '__main__':
    raise SystemExit(main())
