#!/usr/bin/env python3
"""One-time LOCAL cleanup of non-thesis numerical research implementations.

Preview: python -X utf8 tools/PruneLegacyNumerical.py
Apply:   python -X utf8 tools/PruneLegacyNumerical.py --apply

Never reads GitHub refs, HEAD, commits, or remote state. Never edits production
files; only archives then deletes paths from a GitHub-reviewed candidate scope.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import shutil
import sys
from datetime import datetime
from pathlib import Path
from uuid import uuid4

# The four directory inventories were inspected from the repository's default
# branch via GitHub Contents API. FNV-1a checksums cover newline-joined,
# sorted immediate filenames. This prevents deletion of unexpected local files.
REVIEWED_DIRECTORIES = {
    'Numerical/experimental': (64, '4f4325ea'),
    'Numerical/direct_one_loop': (8, '0cbc0129'),
    'Numerical/lo_certification/direct_uv': (24, 'd7e50a5d'),
    'Numerical/lo_patterns_16': (32, '5d1bfdd6'),
}

# Scripts identified alongside those independent / provisional projects.
# Missing paths are harmless: earlier cleanup rounds may already have removed
# them. Present paths are treated as *candidates*, not unconditional deletions.
REVIEWED_SCRIPTS = (
    'scripts/AnalyseLOPatternStudy16.py',
    'scripts/AuditT3DirectUVPortalBasis.py',
    'scripts/BuildLOFinalSynthesis16.py',
    'scripts/BuildT3DirectOneLoopInputs.py',
    'scripts/ExplainLOPatterns16.py',
    'scripts/ResolveT3DirectUVPortalIdentity.py',
    'scripts/RunLODirectUVAmplitudeBridge.py',
    'scripts/RunLODirectUVCertification.py',
    'scripts/RunLODirectUVFermionSpinorBridge.py',
    'scripts/RunLODirectUVMixingDictionary.py',
    'scripts/RunLODirectUVNativeAmplitudeContract.py',
    'scripts/RunLODirectUVNativeFermionAction.wl',
    'scripts/RunLODirectUVNativeFermionRegistry.wl',
    'scripts/RunLODirectUVNativeLoopReference.py',
    'scripts/RunLODirectUVNativeMixingTensorExport.wl',
    'scripts/RunLODirectUVNativeOperatorExport.wl',
    'scripts/RunLODirectUVNativePortalTensorExport.wl',
    'scripts/RunLODirectUVPortalDictionary.py',
    'scripts/RunLODirectUVSpinorMassDictionary.py',
    'scripts/RunLOPatternStudy16.py',
    'scripts/RunT3DirectOneLoop.py',
    'scripts/RunT3DirectOneLoopPortalSignSensitivity.py',
    'scripts/SynthesiseLOPatterns16.py',
    'Numerical/diagnostics/worker5_scalar_matching.py',
)
SOURCE_ROOTS = (
    'common', 'model', 'Lagrangian', 'RGE', 'Numerical', 'physics',
    'studies', 'validation', 'tests', 'configs', 'input', 'scripts',
    'tools', 'display', 'Reports',
)
SOURCE_SUFFIXES = {
    '.py', '.pyw', '.wl', '.wls', '.m', '.ps1', '.psm1', '.sh',
    '.json', '.toml', '.yaml', '.yml', '.ini', '.cfg', '.ipynb',
}
EXCLUDED_PARTS = {
    '.git', '.venv', '__pycache__', 'output', '.cleanup_quarantine',
    'node_modules', '.pytest_cache', '.mypy_cache', '.ruff_cache',
}


def fnv1a_inventory(names: list[str]) -> str:
    value = 2166136261
    for byte in '\n'.join(sorted(names)).encode('utf-8'):
        value = ((value ^ byte) * 16777619) & 0xffffffff
    return f'{value:08x}'


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as fp:
        for chunk in iter(lambda: fp.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def as_rel(root: Path, p: Path) -> str:
    return p.relative_to(root).as_posix()


def add_reviewed_candidates(root: Path) -> tuple[set[str], dict[str, str]]:
    candidates: set[str] = set()
    groups: dict[str, str] = {}
    for directory, (count, expected_hash) in REVIEWED_DIRECTORIES.items():
        target = root / directory
        if not target.exists():
            print(f'ALREADY ABSENT: {directory}/')
            continue
        if target.is_symlink() or not target.is_dir():
            raise RuntimeError(f'Unexpected directory type: {directory}')
        entries = list(target.iterdir())
        # Untracked Python caches and local subdirectories may exist after a
        # scientific run. They are outside the reviewed top-level file list:
        # preserve them, and never descend into them to select deletions.
        # We still reject symlinks and non-regular entries explicitly.
        links = [p.name for p in entries if p.is_symlink()]
        if links:
            raise RuntimeError(
                f'Symlink(s) in {directory} are not safe to inventory: '
                + ', '.join(sorted(links))
            )
        folders = sorted(p for p in entries if p.is_dir())
        unknown = [p.name for p in entries if not p.is_file() and not p.is_dir()]
        if unknown:
            raise RuntimeError(
                f'Non-file entry/entries in {directory}: '
                + ', '.join(sorted(unknown))
            )
        for folder in folders:
            print(f'PRESERVED unreviewed subdirectory: {as_rel(root, folder)}/')
        # File counts and hashes still must agree with the exact reviewed
        # inventory. Unexpected top-level files remain a fail-closed condition.
        regular_files = [p for p in entries if p.is_file()]
        names = [p.name for p in regular_files]
        actual = fnv1a_inventory(names)
        if len(names) != count or actual != expected_hash:
            raise RuntimeError(
                f'{directory}/ differs from reviewed inventory: '
                f'count={len(names)}, FNV1a={actual}; '
                f'expected {count}, {expected_hash}. No deletion attempted.'
            )
        for p in regular_files:
            rel = as_rel(root, p)
            candidates.add(rel)
            groups[rel] = directory
    for rel in REVIEWED_SCRIPTS:
        p = root / rel
        if p.is_symlink():
            raise RuntimeError(f'Refusing symlink: {rel}')
        if p.is_file():
            candidates.add(rel)
            groups[rel] = 'script:' + rel
    return candidates, groups


def read_text(path: Path) -> str:
    # All scanned files are plausible executable/configuration sources.
    # Invalid text encodings are decoded conservatively rather than ignored.
    return path.read_text(encoding='utf-8-sig', errors='replace')


def source_files(root: Path, candidates: set[str], self_path: str):
    for folder in SOURCE_ROOTS:
        base = root / folder
        if not base.is_dir():
            continue
        for p in base.rglob('*'):
            if not p.is_file() or p.is_symlink():
                continue
            if any(x in EXCLUDED_PARTS for x in p.relative_to(root).parts):
                continue
            rel = as_rel(root, p)
            if rel in candidates or rel == self_path:
                continue
            if p.suffix.lower() not in SOURCE_SUFFIXES:
                continue
            if p.stat().st_size > 5_000_000:
                raise RuntimeError(f'Executable/config file too large to scan safely: {rel}')
            yield rel, read_text(p)
    for name in ('pipeline.py', 'full16_symbolic.py'):
        if (root / name).is_file():
            yield name, read_text(root / name)


def token_patterns(rel: str, group: str) -> list[re.Pattern]:
    p = Path(rel)
    if group.startswith('script:'):
        tokens = [re.escape(p.name), re.escape(p.stem)]
        return [re.compile(r'(?<!\w)' + token + r'(?!\w)') for token in tokens]
    if group == 'Numerical/experimental':
        if p.name == '__init__.py':
            return [re.compile(r'(?m)^\s*(?:import\s+Numerical\.experimental\s*$|from\s+Numerical\s+import\s+experimental\b)')]
        stem = re.escape(p.stem)
        return [
            re.compile(r'(?<!\w)' + stem + r'(?!\w)'),
            re.compile(re.escape(rel)),
        ]
    # The direct UV and LO-pattern directories are mutually dependent packages.
    # A reference to any module or package path protects the full directory.
    return [re.compile(re.escape(group.replace('/', '.'))),
            re.compile(re.escape(group))]


def scan_references(root: Path, candidates: set[str], groups: dict[str, str], self_path: str):
    patterns = {rel: token_patterns(rel, groups[rel]) for rel in sorted(candidates)}
    anchors = list(source_files(root, candidates, self_path))
    protected: dict[str, str] = {}
    for source, text in anchors:
        for rel, regexes in patterns.items():
            if rel in protected:
                continue
            if any(pattern.search(text) for pattern in regexes):
                protected[rel] = source

    # For packaged experimental modules, dependencies within candidate files
    # must also be kept whenever any of those files is independently required.
    # Whole direct-UV / LO-pattern package groups are kept together.
    changed = True
    while changed:
        changed = False
        for rel in tuple(protected):
            group = groups[rel]
            if group in REVIEWED_DIRECTORIES and group != 'Numerical/experimental':
                for member, member_group in groups.items():
                    if member_group == group and member not in protected:
                        protected[member] = rel + ' (package consistency)'
                        changed = True
            if group == 'Numerical/experimental' and Path(rel).name != '__init__.py':
                init = 'Numerical/experimental/__init__.py'
                if init in candidates and init not in protected:
                    protected[init] = rel + ' (package __init__)'
                    changed = True
            if Path(rel).suffix.lower() not in SOURCE_SUFFIXES:
                continue
            content = read_text(root / rel)
            for target, regexes in patterns.items():
                if target not in protected and any(pattern.search(content) for pattern in regexes):
                    protected[target] = rel
                    changed = True
    return anchors, protected


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--apply', action='store_true', help='Backup and delete unreferenced candidates')
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    myself = as_rel(root, Path(__file__).resolve())
    print(f'Repository: {root}')
    print('Read-only source inventory; NO GitHub/HEAD/commit access or git commands.')
    candidates, groups = add_reviewed_candidates(root)
    anchors, protected = scan_references(root, candidates, groups, myself)
    remove = sorted(candidates.difference(protected))
    total = sum((root / rel).stat().st_size for rel in remove)
    print(f'Reviewed candidates found locally: {len(candidates)}')
    print(f'Retained-source files scanned: {len(anchors)}')
    print(f'Protected due to source references: {len(protected)}')
    print(f'Proposed removals: {len(remove)} files ({total / 1048576:.2f} MiB)')
    print('\nPROTECTED:')
    for rel in sorted(protected):
        print(f'  KEEP   {rel} <- {protected[rel]}')
    print('\nCANDIDATES TO REMOVE:')
    for rel in remove:
        print('  REMOVE ' + rel)
    if not args.apply:
        print('\nPREVIEW ONLY; use --apply after reviewing the list.')
        return 0
    if not remove:
        print('Nothing eligible for removal.')
        return 0

    # All local files, including uncommitted edits, are archived before deletion.
    # No original is touched until every copy is SHA-256 verified.
    archive = root.parent / f'T3_removed_numerical_round4_{datetime.now():%Y%m%d_%H%M%S}_{uuid4().hex[:8]}'
    if archive.is_relative_to(root):
        raise RuntimeError('Backup location is inside the repository')
    archive.mkdir(parents=False, exist_ok=False)
    records = []
    for rel in remove:
        src, dst = root / rel, archive / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)
        a, b = digest(src), digest(dst)
        if a != b:
            raise RuntimeError(f'SHA256 mismatch in external archive for {rel}. No deletion started.')
        records.append({'path': rel, 'sha256': a, 'bytes': src.stat().st_size})
    with (archive / 'removed_files.csv').open('w', newline='', encoding='utf-8') as fp:
        writer = csv.DictWriter(fp, fieldnames=('path', 'sha256', 'bytes'))
        writer.writeheader()
        writer.writerows(records)
    (archive / 'protected_files.json').write_text(json.dumps(protected, indent=2) + '\n', encoding='utf-8')
    (archive / 'summary.json').write_text(json.dumps({
        'repository': str(root), 'candidate_count': len(candidates),
        'protected_count': len(protected), 'removed_count': len(remove),
        'created': datetime.now().isoformat(),
        'notice': 'No GitHub or git metadata was queried.',
    }, indent=2) + '\n', encoding='utf-8')
    # Verify originals again immediately before unlink, avoiding modified files
    # between the initial scan and deletion.
    for record in records:
        rel = record['path']
        if digest(root / rel) != record['sha256']:
            raise RuntimeError(f'Local content changed before deletion: {rel}; originals untouched.')
    for record in records:
        (root / record['path']).unlink()
    for dirname in sorted(REVIEWED_DIRECTORIES, key=len, reverse=True):
        folder = root / dirname
        if folder.is_dir() and not any(folder.iterdir()):
            folder.rmdir()
    print(f'\nREMOVED {len(remove)} candidate files; external backup: {archive}')
    print('Run the core regression tests and full pytest suite before further cleanup.')
    print('If checks fail, restore the listed files from the external backup.')
    return 0


if __name__ == '__main__':
    try:
        sys.exit(main())
    except (OSError, ValueError, RuntimeError) as exc:
        print(f'ABORTED: {exc}', file=sys.stderr)
        sys.exit(2)
