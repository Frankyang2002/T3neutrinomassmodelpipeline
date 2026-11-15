#!/usr/bin/env python3
"""Locally consolidate historical T3 technical notes into an external archive.

This program is conservative: it only archives clean, tracked Markdown files
listed in audits/technical_note_inventory.csv, after checking references from
retained source files. It never edits scientific code, numerical data, tests,
Git history, or GitHub. It does not decide whether an experimental physical
claim is valid.

From the project root:
    python -X utf8 scripts/maintenance/ConsolidateTechnicalNotes.py plan
    python -X utf8 scripts/maintenance/ConsolidateTechnicalNotes.py apply --confirm ARCHIVE_TECHNICAL_NOTES
    python -X utf8 scripts/maintenance/ConsolidateTechnicalNotes.py restore --manifest 'C:\\...\\manifest.json'

The full original text of each archived file remains in a checksum-verified
ZIP in a sibling directory outside the repository, recoverable by 'restore'.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import posixpath
import re
import subprocess
import sys
import zipfile

ROOT = Path(__file__).resolve().parents[2]
SOURCE_HEAD = '8719090f19e9f6a224b0e4eb1de1c642056277f0'
INVENTORY = ROOT / 'audits' / 'technical_note_inventory.csv'
PREFLIGHT = ROOT / 'audits' / 'technical_note_preflight.csv'
HISTORY = ROOT.parent / 'T3_technical_notes_external_archive'
INDEX_REL = 'docs/research/INDEX.md'
REQUIRED_DOCS = (
    'docs/UV_MODEL_AND_RUNNING.md',
    'docs/FERMION_MATCHING.md',
    'docs/SCALAR_EFT_AND_WEINBERG.md',
    'docs/NUMERICAL_VALIDATION.md',
)
# Conservative text formats for literal-reference checks. Binary outputs are
# not considered textual references; no physics files are ever selected here.
TEXT_EXTENSIONS = {
    '.py', '.wl', '.wls', '.m', '.md', '.rst', '.txt', '.ps1', '.sh',
    '.bat', '.cmd', '.json', '.jsonl', '.yaml', '.yml', '.toml',
    '.ini', '.cfg', '.csv', '.tsv', '.tex', '.bib', '.html', '.xml',
}
MAX_SCAN_BYTES = 8 * 1024 * 1024
CONFIRM = 'ARCHIVE_TECHNICAL_NOTES'


def git(*args: str) -> bytes:
    result = subprocess.run(['git', '-C', str(ROOT), *args],
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            check=False)
    if result.returncode:
        raise RuntimeError('git ' + ' '.join(args) + ': ' +
                           result.stderr.decode('utf-8', 'replace')[:350])
    return result.stdout


def git_head() -> str:
    return git('rev-parse', 'HEAD').decode('ascii').strip()


def tracked() -> list[str]:
    return sorted(os.fsdecode(x).replace('\\', '/') for x in
                  git('ls-files', '-z').split(b'\0') if x)


def validate_rel(rel: str) -> str:
    p = PurePosixPath(rel)
    if (not rel or p.is_absolute() or p.as_posix() != rel or
            any(s in {'', '.', '..'} for s in p.parts) or
            '\\' in rel or ':' in rel or not rel.lower().endswith('.md')):
        raise ValueError('Unexpected candidate path: ' + repr(rel))
    return rel


def local_file(rel: str) -> Path:
    validate_rel(rel)
    path = ROOT.joinpath(*PurePosixPath(rel).parts)
    if path.is_symlink() or not path.resolve().is_relative_to(ROOT.resolve()):
        raise RuntimeError('Symlink or escaping path: ' + rel)
    return path


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    return sha256_bytes(path.read_bytes())


def clean_tracked_file(rel: str, known: set[str]) -> str | None:
    if rel not in known:
        return 'NOT_TRACKED'
    path = local_file(rel)
    if not path.is_file():
        return 'MISSING'
    if git('status', '--porcelain=v1', '--', rel).strip():
        return 'LOCAL_MODIFICATION'
    expected = git('rev-parse', 'HEAD:' + rel).decode('ascii').strip()
    actual = git('hash-object', '--path=' + rel, '--', rel).decode('ascii').strip()
    if expected != actual:
        return 'BLOB_CHANGED'
    return None


def inventory_rows() -> list[dict[str, str]]:
    with INVENTORY.open(newline='', encoding='utf-8-sig') as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != 77:
        raise RuntimeError('Expected 77 documented note candidates, got ' + str(len(rows)))
    seen = set()
    for row in rows:
        rel = validate_rel(row['path'])
        if rel in seen:
            raise RuntimeError('Duplicate target: ' + rel)
        if row['topic'] not in {'uv', 'matching', 'scalar', 'numerical'}:
            raise RuntimeError('Unknown technical topic for ' + rel)
        if row['proposed_treatment'] != 'archive_if_dependency_safe':
            raise RuntimeError('Unexpected proposed treatment for ' + rel)
        seen.add(rel)
    return rows


def candidate_refs(rows: list[dict[str, str]], paths: list[str]) -> dict[str, set[str]]:
    """Search exact paths and sufficiently distinctive basenames in tracked text.

    A text reference from another file earmarks that source for preservation.
    The root docs/research index is a special case: eligible simple Markdown
    link entries can be removed from that index in the same local operation.
    """
    names = [row['path'] for row in rows]
    patterns: dict[str, set[str]] = {}
    for rel in names:
        patterns.setdefault(rel.casefold(), set()).add(rel)
        base = posixpath.basename(rel).casefold()
        if len(base) >= 12:
            patterns.setdefault(base, set()).add(rel)
    regex = re.compile('|'.join(re.escape(k) for k in sorted(patterns, key=len, reverse=True)), re.I)
    found: dict[str, set[str]] = {rel: set() for rel in names}
    for source in paths:
        if source in {'audits/technical_note_inventory.csv', 'audits/technical_note_preflight.csv',
                      'audits/front_facing_review.csv', 'audits/front_facing_reference_report.csv'}:
            continue  # planning inventories are not runtime dependencies
        path = ROOT / source
        if not path.is_file() or path.is_symlink() or path.suffix.lower() not in TEXT_EXTENSIONS:
            continue
        if path.stat().st_size > MAX_SCAN_BYTES:
            continue
        # Git-tracked numerical manifests may be large. They are read (up to
        # this byte limit) to avoid assuming a saved source is unreferenced.
        try:
            text = path.read_text(encoding='utf-8-sig', errors='replace').casefold().replace('\\', '/')
        except OSError as exc:
            raise RuntimeError('Cannot scan retained text ' + source) from exc
        for match in regex.finditer(text):
            for target in patterns[match.group(0)]:
                if target != source:
                    found[target].add(source)
    return found


def normalize_index_target(text: str) -> str:
    # A docs/research/INDEX.md Markdown link is relative to docs/research/.
    return posixpath.normpath(posixpath.join('docs/research', text.split('#')[0]))


def rewrite_index(content: str, archive: set[str]) -> tuple[str, set[str]]:
    removed: set[str] = set()
    output: list[str] = []
    link = re.compile(r'^\s*-\s+\[[^\]]+\]\(([^)]+)\)\s*$')
    for line in content.splitlines(keepends=True):
        m = link.match(line.strip('\r\n'))
        target = normalize_index_target(m.group(1)) if m else None
        if target in archive:
            removed.add(target)
            continue
        output.append(line)
    if removed:
        newline = '\r\n' if '\r\n' in content else '\n'
        # This source index otherwise points to absent old paths. Its other
        # entries remain byte-for-byte unchanged.
        prefix = ('> Consolidated technical explanations: '
                  '[UV theory](../UV_MODEL_AND_RUNNING.md), '
                  '[fermion matching](../FERMION_MATCHING.md), '
                  '[scalar EFT](../SCALAR_EFT_AND_WEINBERG.md), and '
                  '[numerical checks](../NUMERICAL_VALIDATION.md).\n\n')
        prefix = prefix.replace('\n', newline)
        if output and output[0].lstrip('\ufeff').startswith('#'):
            output.insert(1, newline + prefix)
        else:
            output.insert(0, prefix)
    return ''.join(output), removed


def plan() -> tuple[list[dict], str | None, dict]:
    rows = inventory_rows()
    paths = tracked()
    known = set(paths)
    head = git_head()
    all_candidates = {row['path'] for row in rows}
    blockers: dict[str, set[str]] = {r['path']: set() for r in rows}
    reasons: dict[str, list[str]] = {r['path']: [] for r in rows}
    if head != SOURCE_HEAD:
        for path in all_candidates:
            reasons[path].append('DIFFERENT_HEAD')
    for required in REQUIRED_DOCS:
        if not (ROOT / required).is_file():
            for path in all_candidates:
                reasons[path].append('MISSING_CONSOLIDATED_DOCS')
            break
    for rel in all_candidates:
        issue = clean_tracked_file(rel, known)
        if issue:
            reasons[rel].append(issue)

    inbound = candidate_refs(rows, paths)
    indexed: set[str] = set()
    index_content = None
    index_old = None
    if INDEX_REL in known:
        index_issue = clean_tracked_file(INDEX_REL, known)
        if index_issue:
            # Index may need rewrite; avoid overwriting local edits.
            for rel in all_candidates:
                if INDEX_REL in inbound[rel]:
                    reasons[rel].append('INDEX_' + index_issue)
        else:
            index_content = (ROOT / INDEX_REL).read_bytes().decode('utf-8')
            index_old = index_content
            _, indexed = rewrite_index(index_content, all_candidates)
    for rel in all_candidates:
        for src in inbound[rel]:
            # Index references are ignored *only if* a matching complete link
            # can be removed from the index in the same transaction.
            if src == INDEX_REL and rel in indexed:
                continue
            blockers[rel].add(src)

    ready = {r for r in all_candidates if not reasons[r]}
    changed = True
    while changed:
        changed = False
        for rel in tuple(ready):
            external = {s for s in blockers[rel] if s not in ready}
            if external:
                ready.remove(rel)
                changed = True
    result = []
    for row in rows:
        rel = row['path']
        outside = sorted(s for s in blockers[rel] if s not in ready)
        state = 'READY' if rel in ready else 'HOLD'
        issue = ';'.join(sorted(set(reasons[rel])))
        if outside:
            issue = (issue + ';' if issue else '') + 'REFERENCED_BY_RETAINED_FILE'
        result.append({**row, 'status': state, 'reason_codes': issue,
                       'referrers': ';'.join(outside[:12]),
                       'referrer_count': len(outside),
                       'sha256': sha256_file(ROOT / rel) if (ROOT / rel).is_file() else ''})
    summary = {'source_head': head, 'baseline': SOURCE_HEAD,
               'ready': sum(r['status'] == 'READY' for r in result),
               'hold': sum(r['status'] == 'HOLD' for r in result),
               'candidates': len(result),
               'index_rewrite_possible': bool(index_content is not None)}
    return result, index_old, summary


def write_preflight(rows: list[dict]) -> None:
    PREFLIGHT.parent.mkdir(parents=True, exist_ok=True)
    keys = ['path', 'topic', 'status', 'reason_codes', 'referrer_count',
            'referrers', 'sha256']
    with PREFLIGHT.open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=keys, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def archive_apply(rows: list[dict], old_index: str | None) -> None:
    selected = sorted((r for r in rows if r['status'] == 'READY'), key=lambda r: r['path'])
    if not selected:
        print('No note files pass the audit. Nothing archived.')
        return
    ready_paths = {r['path'] for r in selected}
    index_new, index_removed = (rewrite_index(old_index, ready_paths)
                                if old_index is not None else ('', set()))
    index_changed = bool(index_removed)
    if any(INDEX_REL in r['referrers'].split(';') for r in selected):
        raise RuntimeError('Unrewritten index dependency in selected paths')
    for rel in ready_paths:
        if rel.startswith('docs/research/') and old_index and rel not in index_removed:
            raise RuntimeError('Index link could not be removed for ' + rel)
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    HISTORY.mkdir(parents=True, exist_ok=True)
    archive = HISTORY / ('technical_notes_' + stamp + '.zip')
    manifest_path = HISTORY / ('technical_notes_' + stamp + '.manifest.json')
    files = []
    for item in selected:
        rel = item['path']
        data = local_file(rel).read_bytes()
        files.append({'path': rel, 'bytes': len(data), 'sha256': sha256_bytes(data),
                      'topic': item['topic']})
    meta = {'schema': 't3-consolidated-technical-notes-v1',
            'head': git_head(), 'baseline': SOURCE_HEAD,
            'archive_zip': str(archive),
            'files': files,
            'index_path': INDEX_REL if index_changed else None,
            'index_sha256_before': sha256_bytes(old_index.encode('utf-8')) if index_changed else None,
            'index_sha256_after': sha256_bytes(index_new.encode('utf-8')) if index_changed else None,
            'scientific_claim': 'Original research records, not independently physics-certified'}
    # Archive and verify all bytes before modifying the working tree.
    with zipfile.ZipFile(archive, 'x', compression=zipfile.ZIP_DEFLATED) as zout:
        for item in files:
            zout.write(ROOT / item['path'], item['path'])
        if index_changed:
            zout.writestr('__history__/research_index_before.md', old_index.encode('utf-8'))
        zout.writestr('__history__/manifest.json', json.dumps(meta, indent=2).encode('utf-8'))
    with zipfile.ZipFile(archive) as zin:
        for item in files:
            if sha256_bytes(zin.read(item['path'])) != item['sha256']:
                raise RuntimeError('Archive hash verification failed: ' + item['path'])
        if index_changed and sha256_bytes(zin.read('__history__/research_index_before.md')) != meta['index_sha256_before']:
            raise RuntimeError('Index backup verification failed')
    # Do not accept files changed between planning and archive creation.
    for item in files:
        if clean_tracked_file(item['path'], set(tracked())):
            raise RuntimeError('File changed after archive creation: ' + item['path'])
    if index_changed and (ROOT / INDEX_REL).read_bytes() != old_index.encode('utf-8'):
        raise RuntimeError('Index changed after preflight')
    done: list[str] = []
    try:
        if index_changed:
            (ROOT / INDEX_REL).write_bytes(index_new.encode('utf-8'))
        for item in files:
            (ROOT / item['path']).unlink()
            done.append(item['path'])
        manifest_path.write_text(json.dumps(meta, indent=2), encoding='utf-8')
    except Exception:
        with zipfile.ZipFile(archive) as zin:
            for rel in reversed(done):
                dst = ROOT / rel
                dst.parent.mkdir(parents=True, exist_ok=True)
                if not dst.exists():
                    dst.write_bytes(zin.read(rel))
        if index_changed:
            (ROOT / INDEX_REL).write_bytes(old_index.encode('utf-8'))
        raise
    print(f'Archived {len(files)} original Markdown files to: {archive}')
    print(f'Manifest for restoration: {manifest_path}')
    print('GitHub unchanged. Inspect local git status/diff and run regression checks.')


def restore(manifest_path: Path) -> None:
    record = json.loads(manifest_path.read_text(encoding='utf-8'))
    if record.get('schema') != 't3-consolidated-technical-notes-v1':
        raise RuntimeError('Invalid archive manifest')
    archive = Path(record['archive_zip'])
    entries = record['files']
    for item in entries:
        rel = validate_rel(item['path'])
        dst = local_file(rel)
        if dst.exists():
            raise RuntimeError('Refusing to overwrite restored original: ' + rel)
    index_rel = record.get('index_path')
    if index_rel:
        path = ROOT / index_rel
        if sha256_file(path) != record['index_sha256_after']:
            raise RuntimeError('Modified research index differs; refusing restoration')
    with zipfile.ZipFile(archive) as zin:
        for item in entries:
            data = zin.read(item['path'])
            if sha256_bytes(data) != item['sha256']:
                raise RuntimeError('Archive checksum mismatch: ' + item['path'])
        if index_rel:
            data = zin.read('__history__/research_index_before.md')
            if sha256_bytes(data) != record['index_sha256_before']:
                raise RuntimeError('Index checksum mismatch')
        for item in entries:
            dst = local_file(item['path'])
            dst.parent.mkdir(parents=True, exist_ok=True)
            with dst.open('xb') as handle:
                handle.write(zin.read(item['path']))
        if index_rel:
            (ROOT / index_rel).write_bytes(zin.read('__history__/research_index_before.md'))
    print(f'Restored {len(entries)} original Markdown files; GitHub unchanged.')


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='action', required=True)
    sub.add_parser('plan', help='Write local audit only; no original files changed')
    a = sub.add_parser('apply', help='Archive verified notes and remove local copies')
    a.add_argument('--confirm', required=True)
    b = sub.add_parser('restore', help='Restore originals from saved archive')
    b.add_argument('--manifest', type=Path, required=True)
    args = parser.parse_args()
    if args.action == 'restore':
        restore(args.manifest)
        return 0
    rows, index, summary = plan()
    write_preflight(rows)
    print('Commit:', summary['source_head'])
    print('Candidates:', summary['candidates'], 'READY:', summary['ready'], 'HOLD:', summary['hold'])
    print('Preflight CSV:', PREFLIGHT)
    if args.action == 'apply':
        if args.confirm != CONFIRM:
            parser.error('Explicit confirmation must equal ' + CONFIRM)
        if summary['source_head'] != SOURCE_HEAD:
            parser.error('Repository HEAD differs from audited baseline; stop and re-review')
        archive_apply(rows, index)
    else:
        print('AUDIT ONLY. No source or output files moved.')
    return 0


if __name__ == '__main__':
    try:
        raise SystemExit(main())
    except (OSError, RuntimeError, ValueError, zipfile.BadZipFile) as exc:
        print('ERROR:', exc, file=sys.stderr)
        raise SystemExit(2)
