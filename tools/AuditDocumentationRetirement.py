"""Read-only reference audit and opt-in safe external archival of old T3 notes.

Run from the repository root:
    python -X utf8 -m tools.AuditDocumentationRetirement

This utility only considers the explicitly enumerated historical Markdown files.
It never edits production code, surviving documentation, tests, or numerical data.
The optional --apply mode archives only documents with no retained live references,
first copying and SHA-256 verifying every byte in an external location.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
AUDIT_DIR = Path('output/cleanup_audit/documentation_v1')
ARCHIVE_DIR = Path('T3_cleanup_archive/documentation_v1')
SELF = 'tools/AuditDocumentationRetirement.py'
AUDIT_GUIDE = 'docs/REPOSITORY_CLEANUP_AUDIT.md'

# Initial retirement group only. Do not infer source-code obsolescence from names.
# The target is an authoritative doc in the checkout, not a claim that the old
# note has ceased being historically useful (all bytes are externally archived).
CANDIDATES = {
    'docs/T3_SCALAR_POWER_COUNTING.md': 'docs/SCALAR_SECTOR.md',
    'docs/T3_SCALAR_VERTEX_BATCH2.md': 'docs/SCALAR_SECTOR.md',
    'docs/T3_SCALAR_QUADRATIC_BATCH3.md': 'docs/SCALAR_SECTOR.md',
    'docs/T3_SCALAR_RG_CLOSURE_BATCH4.md': 'docs/SCALAR_SECTOR.md',
    'docs/T3_SCALAR_REDUCED_RG_BATCH5.md': 'docs/SCALAR_SECTOR.md',
    'docs/T3_LEADING_SCALAR_POTENTIAL.md': 'docs/SCALAR_SECTOR.md',
    'docs/T3_REDUCED_SCALAR_POTENTIAL_CONCLUSION.md': 'docs/SCALAR_SECTOR.md',
    'docs/T3_ZERO_EXTRA_BOUNDARY_BENCHMARK_BATCH7.md': 'docs/SCALAR_SECTOR.md',
    'docs/T3_ZERO_BOUNDARY_LOW_ENERGY_BATCH8.md': 'docs/SCALAR_SECTOR.md',
    'docs/T3_GENERAL_SCALAR_STUDY_SEPARATION.md': 'docs/SCALAR_SECTOR.md',
    'docs/T3_REDUCED_SCALAR_SOBOL_COMPARISON.md': 'docs/NUMERICAL.md',
    'docs/T3_REDUCED_SCALAR_FINAL_RESULTS.md': 'docs/RESULTS.md',
    'docs/T3_REDUCED_SCALAR_EXISTING_INTERACTIVE_REPORTS.md': 'docs/NUMERICAL.md',
    'docs/T3_REDUCED_SCALAR_LAMBDA_STYLE_REPORTS.md': 'docs/NUMERICAL.md',
}
AUTHORITATIVE_DOCS = frozenset({'docs/SCALAR_SECTOR.md', 'docs/NUMERICAL.md', 'docs/RESULTS.md'})
# SHA-256 of the delivered consolidated documents after normalising CRLF to LF.
# Prevent retirement if the old, unconsolidated NUMERICAL/RESULTS documents
# are still installed, or if the masters were subsequently rewritten.
MASTER_LF_SHA256 = {
    'docs/SCALAR_SECTOR.md': 'af587e0e1612d2cd4f7ac8c7e370332b5b1149c7cf3f9c724bf5758a70dd8253',
    'docs/NUMERICAL.md': '480e1aa8d7282720a92ec60a82cb3b60a63c20d894a8b80f752bc6ef7828f7f2',
    'docs/RESULTS.md': '3c9546179d3959d729022b0df24eebce1d8da3e5bab22c9ba5dfc586825b56e7',
}
PROVENANCE_NOTES = AUTHORITATIVE_DOCS
IGNORE_DIR_NAMES = frozenset({'.git', '.venv', 'venv', 'env', '__pycache__', '.pytest_cache',
                              '.mypy_cache', '.ruff_cache', '.tox', 'node_modules', '.idea',
                              '.vscode', '.ipynb_checkpoints'})
TEXT_EXTENSIONS = frozenset({'.md', '.txt', '.py', '.wl', '.m', '.ps1', '.psm1', '.psd1',
                             '.json', '.jsonl', '.rst', '.tex', '.toml', '.yml', '.yaml',
                             '.ini', '.cfg', '.sh', '.bat', '.cmd', '.html', '.htm', '.xml',
                             '.csv', '.tsv', '.ipynb', '.js', '.ts', '.css', '.bib', '.wls'})
TEXT_NAMES = frozenset({'README', 'Makefile', 'Dockerfile', '.gitignore', '.gitattributes'})
SKIP_REFERENCE_SOURCES = frozenset({SELF, AUDIT_GUIDE})
# Use the original path/stem, rather than guessing whether a renamed script exists.
NAVIGATION_RE = re.compile(r'\b(?:see|refer(?:ring)? to|read|described in|documented in|follow|details in|consult|as in)\b', re.I)
MARKDOWN_LINK_RE = re.compile(r'\]\([^\n)]*\)|^\s*\[[^\]]+\]:\s*\S+', re.I)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def normalised_text_sha256(path: Path) -> str:
    data = path.read_bytes().replace(b'\r\n', b'\n')
    return hashlib.sha256(data).hexdigest()


def master_status(root: Path) -> dict[str, bool]:
    return {name: (root / name).is_file() and
            normalised_text_sha256(root / name) == expected
            for name, expected in MASTER_LF_SHA256.items()}


def relative_posix(path: Path, root: Path) -> str:
    return path.relative_to(root).as_posix()


def current_head(root: Path) -> str | None:
    try:
        cmd = subprocess.run(['git', '-C', str(root), 'rev-parse', 'HEAD'],
                             text=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
                             timeout=10, check=True)
        return cmd.stdout.strip()
    except (OSError, subprocess.CalledProcessError, subprocess.TimeoutExpired):
        return None


def list_text_files(root: Path, report_path: Path, archive_path: Path):
    """Scan actual local files (tracked and untracked), including saved reports.

    Avoid reading virtualenvs, Git internals, compiled/binary assets and the
    external archive. Symlinked source trees are not followed.
    """
    for directory, folders, files in os.walk(root, followlinks=False):
        base = Path(directory)
        folders[:] = [name for name in folders
                      if name not in IGNORE_DIR_NAMES
                      and not (base / name).is_symlink()
                      and not (base / name).resolve().is_relative_to(report_path)
                      and not (base / name).resolve().is_relative_to(archive_path)]
        for name in files:
            file = base / name
            if file.is_symlink() or not file.is_file():
                continue
            if file.suffix.lower() not in TEXT_EXTENSIONS and name not in TEXT_NAMES and not name.startswith('README_'):
                continue
            relative = relative_posix(file, root)
            if relative in SKIP_REFERENCE_SOURCES:
                continue
            if file.resolve().is_relative_to(report_path) or file.resolve().is_relative_to(archive_path):
                continue
            yield file, relative


def classify_reference(source: str, line: str, target: str) -> str:
    """A master document's *plain provenance list* is not a broken link.

    All other literal references block automatic archiving. This deliberately
    treats many non-link prose mentions as blockers rather than guessing intent.
    """
    if source in PROVENANCE_NOTES:
        if MARKDOWN_LINK_RE.search(line) or NAVIGATION_RE.search(line):
            return 'blocking_authoritative_navigation'
        return 'historical_provenance_mention'
    return 'blocking_reference'


def scan(root: Path, report_path: Path, archive_path: Path):
    names = {path: (Path(path).name.casefold(), Path(path).stem.casefold())
             for path in CANDIDATES}
    results = {path: [] for path in CANDIDATES}
    scanned = 0
    errors = []
    for source_path, relative in list_text_files(root, report_path, archive_path):
        scanned += 1
        try:
            # Text streams are line-buffered, avoiding whole-file memory loads.
            with source_path.open('r', encoding='utf-8-sig', errors='replace', newline='') as stream:
                for index, line in enumerate(stream, 1):
                    folded = line.casefold()
                    hits = [target for target, (name, stem) in names.items()
                            if name in folded or stem in folded]
                    if not hits:
                        continue
                    for target in hits:
                        results[target].append({
                            'source': relative,
                            'line': index,
                            'kind': classify_reference(relative, line, target),
                            'excerpt': line.strip()[:320],
                        })
        except (OSError, UnicodeError) as error:
            errors.append({'source': relative, 'error': f'{type(error).__name__}: {error}'})
    return results, scanned, errors


def safe_archive_set(present: set[str], refs: dict):
    """Greatest safe set: no retained file may refer to an archived file.

    Internal references between files archived *together* are fine; if any
    candidate is blocked, that retained candidate may block additional ones.
    """
    allowed = set(present)
    blockers = {}
    while True:
        remove = {}
        for target in sorted(allowed):
            active = [entry for entry in refs[target]
                      if entry['source'] not in allowed and entry['kind'] != 'historical_provenance_mention']
            if active:
                remove[target] = active
        if not remove:
            break
        allowed.difference_update(remove)
        blockers.update(remove)
    for target in present - allowed:
        blockers[target] = [entry for entry in refs[target]
                            if entry['source'] not in allowed and entry['kind'] != 'historical_provenance_mention']
    return allowed, blockers


def require_destination_outside_repo(root: Path, dest: Path):
    if dest == root or dest.is_relative_to(root):
        raise ValueError('Archive destination MUST be outside the repository')
    if root.is_relative_to(dest):
        raise ValueError('Archive destination cannot contain the entire repository')


def write_report(root: Path, report_dir: Path, archive_dir: Path,
                 records: dict, refs: dict, safe: set, blockers: dict,
                 scanned: int, errors: list, head: str | None, masters: dict[str, bool]):
    report_dir.mkdir(parents=True, exist_ok=True)
    now = datetime.now(timezone.utc).isoformat()
    record = {
        'format': 't3-documentation-retirement-audit-v1',
        'generated_at_utc': now,
        'checkout_head': head,
        'repository_root': str(root),
        'archive_destination': str(archive_dir),
        'source_files_scanned': scanned,
        'source_read_errors': errors,
        'required_authoritative_docs': sorted(AUTHORITATIVE_DOCS),
        'consolidated_master_files_verified': masters,
        'selection_rule': 'Safe = candidate with no references outside the whole proposed archive group; '
                          'provenance-only mentions in authoritative documents are recorded but nonblocking',
        'candidates': {
            path: {
                **entry,
                'safe_for_external_archive': path in safe and not errors,
                'reference_count': len(refs[path]),
                'references': refs[path],
                'blocking_references': blockers.get(path, []),
            } for path, entry in records.items()
        },
        'safe_candidate_count': len(safe) if not errors else 0,
        'blocked_candidate_count': len([p for p, r in records.items() if r['exists'] and p not in safe]) + (len(safe) if errors else 0),
        'missing_candidate_count': sum(not r['exists'] for r in records.values()),
        'no_files_moved': True,
    }
    (report_dir / 'audit.json').write_text(json.dumps(record, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')
    with (report_dir / 'references.csv').open('w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=['target', 'source', 'line', 'kind', 'excerpt'])
        writer.writeheader()
        for target in CANDIDATES:
            for reference in refs[target]:
                writer.writerow({'target': target, **reference})
    lines = [
        '# Historical documentation reference audit', '',
        f'- Checkout HEAD: `{head or "not available"}`',
        f'- Text files scanned: **{scanned}**',
        f'- Candidates present: **{sum(x["exists"] for x in records.values())} / {len(CANDIDATES)}**',
        f'- Safe for external archival: **{len(safe) if not errors else 0}**',
        f'- Blocked by retained references: **{len([p for p, r in records.items() if r["exists"] and p not in safe])}**',
        f'- Source read errors: **{len(errors)}**',
        f'- Authoritative document fingerprint match: **{sum(masters.values())}/{len(masters)}**',
        '',
        'This is a dry-run. Files are **not** moved unless `--apply --confirm ARCHIVE_SAFE_ONLY` is supplied.',
        'Archival additionally requires all three installed master documents to match the delivered consolidated versions.',
        'Historical source mentions in authoritative documentation are retained as provenance, not live links.',
        '',
        '| Candidate historical document | Classification | Live blockers | Replacement |',
        '|---|---|---:|---|',
    ]
    for target, entry in records.items():
        status = 'MISSING' if not entry['exists'] else 'SAFE' if target in safe and not errors else 'BLOCKED'
        lines.append(f'| `{target}` | {status} | {len(blockers.get(target, []))} | `{entry["consolidated_in"]}` |')
    lines.extend(['', '## What prevents archiving', ''])
    for target in CANDIDATES:
        blocked = blockers.get(target, [])
        if not blocked:
            continue
        lines.append(f'### `{target}`')
        lines.append('')
        for ref in blocked[:30]:
            lines.append(f'- `{ref["source"]}:{ref["line"]}` — {ref["kind"]}: `{ref["excerpt"].replace("`", "")[:145]}`')
        if len(blocked) > 30:
            lines.append(f'- ... plus {len(blocked)-30} additional references in `references.csv`.')
        lines.append('')
    lines.extend([
        '## Safe procedure', '',
        '1. Review `audit.json`, `references.csv` and this report, including blocked files.',
        '2. Do not archive any blocked file until references are updated in full replacement files.',
        '3. To archive only the safe group, explicitly run `--apply --confirm ARCHIVE_SAFE_ONLY`.',
        '4. This copies every selected original to the external archive, verifies SHA-256, then removes its active copy.',
        '5. Rerun `python -X utf8 -m pytest -q --tb=line`; do not remove fitting data or source code.',
        '',
    ])
    if errors:
        lines.append('## Audit read failures (archiving disabled)')
        for error in errors:
            lines.append(f'- `{error["source"]}`: {error["error"]}')
    (report_dir / 'summary.md').write_text('\n'.join(lines) + '\n', encoding='utf-8')
    return record


def archive_safe(root: Path, archive_dir: Path, records: dict, allowed: set[str],
                 head: str | None):
    if not allowed:
        print('No safe candidates to archive; no files changed.')
        return
    if archive_dir.exists():
        # Never merge into an old archive: its prior manifest may describe a
        # different source tree, even if candidate names match.
        raise ValueError(f'Archive destination already exists: {archive_dir}')
    # Preflight ALL candidate inputs/destinations before making changes.
    for name in sorted(allowed):
        source = root / name
        if source.is_symlink() or not source.is_file() or sha256_file(source) != records[name]['sha256']:
            raise ValueError(f'Source changed after audit; refusing archive: {name}')
        target = archive_dir / name
        if target.exists():
            raise ValueError(f'Archive target already exists: {target}')
    archive_dir.mkdir(parents=True, exist_ok=False)
    manifest = {
        'format': 't3-external-documentation-archive-v1',
        'archived_at_utc': datetime.now(timezone.utc).isoformat(),
        'source_head': head,
        'source_checkout': str(root),
        'archived_files': [
            {'original_project_path': name,
             'archive_relative_path': name,
             'sha256': records[name]['sha256'],
             'size_bytes': records[name]['size_bytes'],
             'consolidated_in': records[name]['consolidated_in'],
             'reason': 'Historical batch-specific documentation consolidated into authoritative reference'}
            for name in sorted(allowed)
        ],
    }
    # Never unlink until ALL copies exist and their hashes are verified.
    for name in sorted(allowed):
        dest = archive_dir / name
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(root / name, dest)
        if sha256_file(dest) != records[name]['sha256']:
            raise ValueError('Copied archive did not preserve source bytes: ' + name)
    (archive_dir / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n', encoding='utf-8')
    for name in sorted(allowed):
        if sha256_file(root / name) != records[name]['sha256']:
            raise ValueError('Source changed during copy; originals retained where possible: ' + name)
    for name in sorted(allowed):
        (root / name).unlink()
        print('ARCHIVED', name)
    print('External archive:', archive_dir)
    print('Manifest:', archive_dir / 'manifest.json')


def run(root: Path, report_dir: Path, archive_dir: Path, apply: bool, confirm: str | None):
    root = root.resolve()
    report_dir = report_dir.resolve()
    archive_dir = archive_dir.resolve()
    require_destination_outside_repo(root, archive_dir)
    if not (root / 'docs').is_dir():
        raise ValueError('Not a T3 repository root (docs/ absent)')
    if apply and confirm != 'ARCHIVE_SAFE_ONLY':
        raise ValueError('--apply requires --confirm ARCHIVE_SAFE_ONLY')
    absent = sorted(doc for doc in AUTHORITATIVE_DOCS if not (root / doc).is_file())
    if absent:
        raise ValueError('Missing consolidated docs; refusing any retirement: ' + ', '.join(absent))
    records = {}
    present = set()
    for name, replacement in CANDIDATES.items():
        file = root / name
        exists = file.is_file() and not file.is_symlink()
        if exists:
            present.add(name)
        records[name] = {'exists': exists, 'sha256': sha256_file(file) if exists else None,
                         'size_bytes': file.stat().st_size if exists else None,
                         'consolidated_in': replacement}
    masters = master_status(root)
    refs, scanned, errors = scan(root, report_dir, archive_dir)
    safe, blockers = safe_archive_set(present, refs)
    write_report(root, report_dir, archive_dir, records, refs, safe, blockers, scanned, errors,
                 current_head(root), masters)
    print('Checked', len(present), 'historical documents; scanned', scanned, 'text files.')
    print('Master doc fingerprints verified:', sum(masters.values()), '/', len(masters))
    print('Safe:', len(safe) if not errors else 0,
          'Blocked:', len(present - safe), 'Read errors:', len(errors))
    print('Review:', report_dir / 'summary.md')
    print('Detailed references:', report_dir / 'references.csv')
    if apply:
        if not all(masters.values()):
            raise RuntimeError('Installed consolidated master docs differ from delivered versions; no files archived: ' + ', '.join(name for name, ok in masters.items() if not ok))
        if errors:
            raise RuntimeError('Audit has read errors; no files archived')
        archive_safe(root, archive_dir, records, safe, current_head(root))
        if safe:
            data_path = report_dir / 'audit.json'
            audited = json.loads(data_path.read_text(encoding='utf-8'))
            audited['no_files_moved'] = False
            audited['archived_files'] = sorted(safe)
            data_path.write_text(json.dumps(audited, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')
            summary_path = report_dir / 'summary.md'
            summary = summary_path.read_text(encoding='utf-8')
            summary = summary.replace('This is a dry-run. Files are **not** moved unless `--apply --confirm ARCHIVE_SAFE_ONLY` is supplied.',
                                      'APPLIED: the safe set was copied to the external archive with SHA-256 verification; original active copies were then removed.')
            summary_path.write_text(summary, encoding='utf-8')
    else:
        print('DRY RUN ONLY: no files moved.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=ROOT,
                        help='Local T3 project root (default: parent of this script)')
    parser.add_argument('--report-dir', type=Path,
                        help='Audit output (default: <root>/output/cleanup_audit/documentation_v1)')
    parser.add_argument('--archive-dir', type=Path,
                        help='External archive (default: <root parent>/T3_cleanup_archive/documentation_v1)')
    parser.add_argument('--apply', action='store_true',
                        help='Copy SHA-verified, unreferenced documents outside repo and remove originals')
    parser.add_argument('--confirm', help='Required exact value ARCHIVE_SAFE_ONLY when --apply is used')
    args = parser.parse_args()
    root = args.root.resolve()
    report_dir = args.report_dir if args.report_dir else root / AUDIT_DIR
    archive_dir = args.archive_dir if args.archive_dir else root.parent / ARCHIVE_DIR
    try:
        run(root, report_dir, archive_dir, args.apply, args.confirm)
    except (OSError, ValueError, RuntimeError) as error:
        print(f'ERROR: {error}', file=sys.stderr)
        return 2
    return 0


if __name__ == '__main__':
    sys.exit(main())
