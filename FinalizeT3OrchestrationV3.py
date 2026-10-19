r"""Finish orchestration relocation, recovering missing archive implementations.

From repository root:
  python .\FinalizeT3OrchestrationV3.py --plan
  python .\FinalizeT3OrchestrationV3.py --apply
  python .\FinalizeT3OrchestrationV3.py --restore

All modified and removed source files are backed up. Any missing moved implementation
is recovered only from the original organiser's verified backup. Full pytest is
run automatically, with rollback on failure.
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
STATE = '.orchestration_finalizer_v3_manifest.json'
BACKUP = '.orchestration_finalizer_v3_backup'
MARKER = '# AUTO-GENERATED T3 ORCHESTRATION FORWARDER'
SKIP = {'.git', '.venv', 'venv', '__pycache__', '.pytest_cache', 'output',
        'Reports', 'display', 'dist', 'build', '.mypy_cache', '.ruff_cache',
        '.orchestration_organizer_backup', '.orchestration_finalizer_backup', BACKUP}


def root_dir():
    root = Path.cwd().resolve()
    if not (root / 'Numerical' / 'orchestration' / 'PipelineNumericalResults.py').is_file():
        raise SystemExit('Run from the project root containing Numerical/orchestration.')
    return root


def original_from_backup(root, op):
    """Verified original source, recovered only from the initial organiser backup."""
    original = root / '.orchestration_organizer_backup' / op['name']
    if not original.is_file():
        raise SystemExit('Missing original organiser backup for: ' + op['name'])
    import hashlib
    data = original.read_bytes()
    if hashlib.sha256(data).hexdigest() != op['original_sha256']:
        raise SystemExit('Organiser backup checksum differs: ' + op['name'])
    return data.decode('utf-8-sig')


def relocate_from_backup(root, op):
    """Replicate the original organiser's path-depth remapping."""
    source = original_from_backup(root, op)
    depth = len(Path(op['group']).parts)
    return re.sub(
        r'(Path\(\s*__file__\s*\)\.resolve\(\)\.parents\[)(\d+)(\])',
        lambda match: match.group(1) + str(int(match.group(2)) + depth) + match.group(3),
        source,
    )


def read_operations(root):
    path = root / ORIG
    if not path.is_file():
        raise SystemExit('Missing organiser manifest: ' + ORIG)
    ops = json.loads(path.read_text(encoding='utf-8'))['operations']
    relocated, recover = [], []
    for op in ops:
        old, new = root / op['old'], root / op['new']
        if new.is_file():
            if not old.is_file() or MARKER not in old.read_text(encoding='utf-8-sig'):
                raise SystemExit('Refusing to delete modified or missing forwarder: ' + op['old'])
            relocated.append(op)
        elif old.is_file() and MARKER in old.read_text(encoding='utf-8-sig'):
            recovered = relocate_from_backup(root, op)
            ast.parse(recovered, filename=str(new))
            recover.append((op, recovered))
            relocated.append(op)
        else:
            raise SystemExit('Cannot safely locate the source and forwarder: ' + op['name'])
    return relocated, recover


def py_sources(root):
    for file in root.rglob('*.py'):
        if any(part in SKIP for part in file.relative_to(root).parts):
            continue
        if file.name in {'OrganizeT3Orchestration.py', 'FinalizeT3Orchestration.py',
                         'FinalizeT3OrchestrationV2.py', 'FinalizeT3OrchestrationV3.py'}:
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
        raise SystemExit('No V3 finalisation backup to restore.')
    state = json.loads(statepath.read_text(encoding='utf-8'))
    for name in state['files']:
        if not (backup / name).is_file():
            raise SystemExit('Missing backup: ' + name)
    for name in state.get('created', []):
        created = root / name
        if created.is_file():
            created.unlink()
    for name in state['files']:
        dst = root / name
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(backup / name, dst)
    shutil.rmtree(backup)
    statepath.unlink()
    print('Restored all original forwarders and edited source files. Moved files retained.')


def apply(root, pytest_args):
    if (root / STATE).exists() or (root / BACKUP).exists():
        raise SystemExit('V3 finalisation state exists; use --restore first.')
    relocated, recover = read_operations(root)
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
            'created': [o['new'] for o, _ in recover],
        }, indent=2) + '\n', encoding='utf-8')
        for op, data in recover:
            moved = root / op['new']
            moved.parent.mkdir(parents=True, exist_ok=True)
            moved.write_text(data, encoding='utf-8')
        for path, data in edits:
            path.write_text(data, encoding='utf-8')
        for path in forwarders:
            path.unlink()
        print(f'Removed {len(forwarders)} moved-file forwarders; rewritten {len(edits)} Python files.', flush=True)
        print(f'Recovered {len(recover)} missing moved files from verified organiser backups.', flush=True)
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
    relocated, recover = read_operations(root)
    edits = changes(root, relocated)
    if args.plan:
        print(f'Moved implementations found: {len(relocated)}')
        print(f'Missing implementations to recover: {len(recover)}')
        for op, _ in recover:
            print('  RECOVER FROM VERIFIED BACKUP', op['new'])
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
