#!/usr/bin/env python3
"""Thesis-scoped repository audit and reversible external archival (local only).

This utility does NOT decide scientific relevance from an unreferenced file name.
It is conservative about source code, configurations, generated matching inputs,
benchmark data, and figure provenance. `plan` performs read-only inventory;
`apply` archives only safe, approved historical documentation/build byproducts.

PowerShell from repository root:
  python -X utf8 scripts/maintenance/CurateThesisRepository.py plan
  python -X utf8 scripts/maintenance/CurateThesisRepository.py apply --plan <plan.json> --confirm ARCHIVE_THESIS_NONESSENTIAL
  # Manually reviewed REVIEW items: use --review-file reviewed.csv and a stronger confirmation
  python -X utf8 scripts/maintenance/CurateThesisRepository.py restore --manifest <manifest.json>

No GitHub operations, code execution, scientific computation or git commit.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict, deque
import csv
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import shutil
import subprocess
import sys

FORMAT = 't3-thesis-curation-v1'
CONFIRM = 'ARCHIVE_THESIS_NONESSENTIAL'
REVIEW_CONFIRM = 'ARCHIVE_MANUALLY_REVIEWED_THESIS_FILES'
SCRIPT = 'scripts/maintenance/CurateThesisRepository.py'
BASELINE = '91f8d7b0f81417df751632bacb436156a04872de'  # inspected public GitHub main, 2026-11-15

# These areas are needed to keep code+data provenance for all thesis chapters.
# A specialised audit may later retire truly unused scientific code, but not by
# inference from names like "batch", "experimental", "stage" or "worker".
PRESERVE_PREFIX = (
    'common/', 'model/', 'Lagrangian/', 'RGE/', 'Numerical/',
    'physics/', 'studies/', 'validation/', 'tests/', 'configs/',
    'input/', 'benchmarks/', 'Result Graphs/', 'display/',
    'Reports/figure/',
    'Reports/output/thesis_general_5_1/',
    'Reports/output/thesis_general_5_2/',
    'Reports/output/thesis_general_5_3/',
    'Reports/output/thesis_qualitative/',
    'output/studies/thesis_qualitative/',
    'output/studies/general_scalar_extended/zee_style_general_sensitivity/zee_general_all16_merged_v1/',
    'output/studies/fermion_start_lambda_comparison/',
    'output/parallel_t3/thesis_',
    'output/parallel_t3/eft1_llss_',
    'output/parallel_t3/background/',
    'output/full16_symbolic/',
)
CORE_ROOT = frozenset({
    'README.md', '.gitignore', 'pipeline.py', 'full16_symbolic.py',
    'requirements.txt',
})
CORE_DOCS = frozenset({
    'README.md', 'INSTALLATION.md', 'IMPLEMENTATION_STATUS.md',
    'REPRODUCIBILITY.md', 'PROJECT_GUIDE.md', 'REPOSITORY_STRUCTURE.md',
    'PHYSICS.md', 'MODEL_AND_PIPELINE.md', 'SCALAR_SECTOR.md',
    'NUMERICAL.md', 'RESULTS.md', 'TEST_STRATEGY.md',
    'REGENERATE_CURRENT_REPORTS.md', 'GENERAL_SYMBOLIC_PHYSICS_AUDIT.md',
    'THESIS_SCOPE.md', 'THESIS_ARTIFACT_MAP.md',
})
ROOT_NOTE = re.compile(
    r'^(?:README[_-]|BATCH\d+[_-]|EFT1[_-]|WORKER[_-]|STREAM\d+[_-]).*\.(?:md|txt)$', re.I)
HISTORY_DOC = re.compile(
    r'^(?:RESEARCH[_-]|GENERAL_T3[_-]|WORKER[_-]|WORKSTREAM[_-]|'
    r'T3_MATCHETE[_-]|T3_MASS_EIGEN[_-]|T3_SCALAR_.*BATCH|'
    r'T3_GENERAL_.*BATCH|REPOSITORY_MAINTENANCE_HISTORY|'
    r'HISTORICAL_DOCS_CONSOLIDATION_GUIDE).*\.md$', re.I)
BUILD_SUFFIXES = ('.compile.log', '.aux', '.toc', '.out', '.log')
TEXT_SUFFIXES = frozenset({
    '.py', '.wl', '.wls', '.m', '.md', '.rst', '.txt', '.ps1',
    '.psm1', '.json', '.jsonl', '.yaml', '.yml', '.toml', '.ini',
    '.cfg', '.csv', '.tsv', '.tex', '.sh', '.bat', '.cmd', '.html', '.htm',
    '.xml', '.ipynb', '.bib',
})
SCAN_MAX_BYTES = 8 * 1024 * 1024


def git(root: Path, *args: str) -> bytes:
    p = subprocess.run(['git', '-C', str(root), *args], capture_output=True, check=False)
    if p.returncode:
        raise RuntimeError(f"git {' '.join(args)}: {p.stderr.decode('utf8','replace')}")
    return p.stdout


def head(root: Path) -> str:
    return git(root, 'rev-parse', 'HEAD').decode('ascii').strip()


def tracked(root: Path) -> list[str]:
    return sorted(PurePosixPath(os.fsdecode(x)).as_posix() for x in
                  git(root, 'ls-files', '-z').split(b'\0') if x)


def clean(root: Path) -> None:
    if git(root, 'status', '--porcelain=v1', '--untracked-files=no').strip():
        raise RuntimeError('Tracked checkout has local changes; commit or back up before applying.')


def safe(root: Path, rel: str) -> Path:
    p = PurePosixPath(rel)
    if not rel or p.is_absolute() or p.as_posix() != rel or any(
            x in ('', '.', '..') for x in p.parts) or ':' in rel or '\\' in rel:
        raise RuntimeError(f'Unsafe relative path: {rel!r}')
    target = root.joinpath(*p.parts)
    if target.is_symlink() or not target.resolve().is_relative_to(root.resolve()):
        raise RuntimeError(f'Symlink/path escape: {rel!r}')
    return target


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as f:
        for part in iter(lambda: f.read(1024 * 1024), b''):
            h.update(part)
    return h.hexdigest()


def file_kind(rel: str, allpaths: set[str]) -> tuple[str, str]:
    """Return (KEEP | REVIEW | CANDIDATE, explanation)."""
    if rel in CORE_ROOT:
        return 'KEEP', 'Principal application entry point/configuration'
    if rel.startswith(PRESERVE_PREFIX):
        return 'KEEP', 'Thesis physics, implementation, input, result or provenance family'
    if rel.startswith('docs/'):
        name = rel.split('/')[-1]
        if rel.count('/') == 1 and name in CORE_DOCS:
            return 'KEEP', 'Current scientific documentation'
        if rel.startswith('docs/research/') or HISTORY_DOC.fullmatch(name):
            return 'CANDIDATE', 'Historical research narrative (not runtime physics)'
        return 'REVIEW', 'Scientific documentation: review thesis relevance manually'
    if '/' not in rel and (ROOT_NOTE.fullmatch(rel) or rel == 'T3_COMMAND_REFERENCE.md'):
        return 'CANDIDATE', 'Development/packaging note at repository root'
    if rel.startswith('report_cleanup_archive/'):
        return 'CANDIDATE', 'Stored intermediate LaTeX build files'
    if rel.startswith('Reports/output/'):
        if any(rel.endswith(x) for x in BUILD_SUFFIXES):
            stem = next(rel[:-len(x)] for x in BUILD_SUFFIXES if rel.endswith(x))
            if stem + '.tex' in allpaths:
                return 'CANDIDATE', 'Regenerable LaTeX build auxiliary with tracked TeX source'
        return 'REVIEW', 'Saved report/data: retain until its thesis use and provenance are reviewed'
    if rel.startswith('output/'):
        return 'REVIEW', 'Saved computation/export: may be a required source input or thesis figure provenance'
    if rel.startswith(('scripts/', 'tools/', 'Reports/')) or rel.endswith(('.py', '.wl', '.wls')):
        return 'REVIEW', 'Scientific/utility source: lack of imports cannot prove irrelevance'
    return 'REVIEW', 'Not conclusively established as unrelated to thesis'


def references(root: Path, paths: list[str], candidates: set[str]) -> dict[str, list[str]]:
    """Find literal mentions of *candidate* filenames in tracked text.

    Be deliberately conservative: an incidental mention can block archiving,
    but must never allow us to retire a required file. Once a retained file
    references a candidate, that candidate and its own dependencies are kept.
    References between files archived *together* do not block archival.
    """
    aliases: dict[str, set[str]] = defaultdict(set)
    for rel in candidates:
        aliases[rel.casefold()].add(rel)
        name = PurePosixPath(rel).name.casefold()
        # Generic names such as INDEX.md and README.md cause false positives.
        # Their references are still found through the full path string.
        if len(name) >= 13 and name not in {'readme.md'}:
            aliases[name].add(rel)
    if not aliases:
        return {}
    pattern = re.compile('|'.join(re.escape(k) for k in sorted(aliases, key=len, reverse=True)), re.I)
    edges: dict[str, set[str]] = defaultdict(set)
    for source in paths:
        p = root / source
        if p.suffix.lower() not in TEXT_SUFFIXES or not p.is_file() or p.stat().st_size > SCAN_MAX_BYTES:
            continue
        try:
            content = p.read_text(encoding='utf-8-sig', errors='replace').casefold()
        except OSError:
            raise RuntimeError(f'Cannot inspect reference source: {source}')
        for m in pattern.finditer(content):
            for dest in aliases[m.group(0)]:
                if source != dest:
                    edges[source].add(dest)
    # Traverse backwards from all files that would be retained.
    incoming: dict[str, list[str]] = defaultdict(list)
    initial = set(paths) - candidates
    kept = set(initial)
    q = deque(initial)
    while q:
        src = q.popleft()
        for target in edges.get(src, ()):
            if target not in kept:
                kept.add(target)
                q.append(target)
            incoming[target].append(src)
    return {rel: sorted(set(incoming.get(rel, []))) for rel in candidates if rel in kept}


def default_archive(root: Path) -> Path:
    return root.parent / 'T3_thesis_external_archive'


def external(root: Path, path: Path) -> Path:
    where = path.resolve()
    if where == root.resolve() or where.is_relative_to(root.resolve()):
        raise RuntimeError('Archive must be outside the Git repository.')
    return where


def make_plan(root: Path, archive: Path) -> Path:
    root = root.resolve()
    archive = external(root, archive)
    allfiles = tracked(root)
    allset = set(allfiles)
    kinds = {p: file_kind(p, allset) for p in allfiles}
    candidates = {p for p, (status, _) in kinds.items() if status == 'CANDIDATE'}
    blocked = references(root, allfiles, candidates)
    rows = []
    for rel in allfiles:
        path = safe(root, rel)
        if not path.is_file():
            raise RuntimeError(f'Tracked file missing in local working tree: {rel}')
        kind, reason = kinds[rel]
        state = 'BLOCKED_REFERENCE' if rel in blocked else 'READY_ARCHIVE' if kind == 'CANDIDATE' else kind
        rows.append({'path': rel, 'decision': state, 'category': kind,
                     'reason': reason, 'bytes': path.stat().st_size,
                     'references': blocked.get(rel, [])})
    out = archive / 'plans' / head(root)[:12]
    out.mkdir(parents=True, exist_ok=True)
    plan = {'schema': FORMAT + '-plan', 'source_head': head(root),
            'audited_baseline': BASELINE, 'repository': str(root),
            'archive_root': str(archive), 'entries': rows,
            'summary': dict(Counter(r['decision'] for r in rows)),
            'notice': 'READY_ARCHIVE only covers reference-free documentation/build intermediates. REVIEW never auto-deleted.'}
    json_path = out / 'plan.json'
    json_path.write_text(json.dumps(plan, indent=2, ensure_ascii=False) + '\n', encoding='utf8')
    with (out / 'inventory.csv').open('w', newline='', encoding='utf-8-sig') as f:
        w = csv.DictWriter(f, fieldnames=['path','decision','category','bytes','reason','references'])
        w.writeheader()
        for row in rows:
            w.writerow({**row, 'references': '; '.join(row['references'][:8])})
    print(f'HEAD: {plan["source_head"]}')
    print(f'Tracked: {len(rows):,}; classifications: {plan["summary"]}')
    print('READY_ARCHIVE means only: research note or LaTeX build intermediate, no retained literal references.')
    print('No scientific code, matching export, fitting configuration or saved output is automatically removed.')
    print('PLAN:', json_path)
    return json_path


def load_review_csv(path: Path) -> set[str]:
    """Only explicit path,decision,justification rows authorise archival.

    A lack of a literal reference does not prove unused science; these rows
    require a human check of run manifests, figure generation and provenance.
    """
    with path.open(newline='', encoding='utf-8-sig') as f:
        reader = csv.DictReader(f)
        if not {'path', 'decision', 'justification'} <= set(reader.fieldnames or []):
            raise RuntimeError('Review CSV must have path,decision,justification columns.')
        selected = set()
        for row in reader:
            if row['decision'].strip().upper() != 'ARCHIVE':
                continue
            justification = row['justification'].strip()
            if len(justification) < 15:
                raise RuntimeError('Provide a specific human-reviewed justification for ' + row['path'])
            path_name = row['path'].strip()
            if path_name in selected:
                raise RuntimeError('Duplicate review path: ' + path_name)
            selected.add(path_name)
        return selected


def apply(root: Path, plan_file: Path, confirm: str, review_file: Path | None = None) -> Path:
    expected = REVIEW_CONFIRM if review_file is not None else CONFIRM
    if confirm != expected:
        raise RuntimeError('Missing explicit confirmation token.')
    root = root.resolve()
    clean(root)
    plan = json.loads(plan_file.read_text(encoding='utf8'))
    if plan.get('schema') != FORMAT + '-plan' or plan.get('repository') != str(root):
        raise RuntimeError('Wrong plan / checkout.')
    if plan['source_head'] != head(root):
        raise RuntimeError('Git HEAD differs from plan; regenerate audit.')
    archive = external(root, Path(plan['archive_root']))
    # Require a fresh local reference check. This prevents stale plans from
    # removing files required by new tracked source and data manifests.
    files = tracked(root)
    if {e['path'] for e in plan['entries']} != set(files):
        raise RuntimeError('Tracked file list changed; regenerate audit.')
    kinds = {p: file_kind(p, set(files)) for p in files}
    candidates = {p for p, (k, _) in kinds.items() if k == 'CANDIDATE'}
    review = load_review_csv(review_file) if review_file is not None else set()
    states = {e['path']: e for e in plan['entries']}
    for rel in review:
        if rel not in states or states[rel]['decision'] != 'REVIEW':
            raise RuntimeError('Only audited REVIEW paths are eligible for reviewed archival: ' + rel)
        if kinds[rel][0] != 'REVIEW':
            raise RuntimeError('Protected file may not be manually selected: ' + rel)
    selected = [e for e in plan['entries']
                if e['decision'] == 'READY_ARCHIVE' or e['path'] in review]
    selected_names = {e['path'] for e in selected}
    blocking = references(root, files, selected_names)
    for e in selected:
        p = e['path']
        if p in blocking:
            raise RuntimeError('Selected file has retained tracked references: ' + p +
                               ' <- ' + ', '.join(blocking[p][:6]))
        if p not in candidates and p not in review:
            raise RuntimeError('Previously safe candidate is now protected: ' + p)
        source = safe(root, p)
        if not source.is_file() or e['bytes'] != source.stat().st_size:
            raise RuntimeError('Candidate changed: ' + p)
    if not selected:
        print('Nothing automatically eligible for archival.')
        return plan_file
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    dest = archive / 'archived' / (stamp + '-' + head(root)[:10])
    if dest.exists():
        raise RuntimeError('Archive destination exists unexpectedly.')
    entries = []
    # Copy + verify everything BEFORE deleting anything. Missing source or
    # failed digest leaves working tree untouched. Never follow symlinks.
    for item in selected:
        rel = item['path']
        src = safe(root, rel)
        dst = dest / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst, follow_symlinks=False)
        digest = sha(src)
        if sha(dst) != digest:
            raise RuntimeError('Checksum mismatch for ' + rel)
        entries.append({'path': rel, 'sha256': digest, 'bytes': item['bytes']})
    manifest = dest / 'manifest.json'
    document = {'schema': FORMAT + '-archive', 'source_head': head(root),
                'repository': str(root), 'entries': entries}
    manifest.write_text(json.dumps(document, indent=2) + '\n', encoding='utf8')
    # Final preflight complete. Do not change Git history; Git records deletions
    # in the working tree for the user to review/commit.
    for e in entries:
        safe(root, e['path']).unlink()
    print(f'ARCHIVED: {len(entries)} files outside Git checkout; manifest: {manifest}')
    print('Review `git status --short` and run the core tests before committing.')
    return manifest


def restore(root: Path, manifest_file: Path) -> None:
    root = root.resolve()
    doc = json.loads(manifest_file.read_text(encoding='utf8'))
    if doc.get('schema') != FORMAT + '-archive' or doc.get('repository') != str(root):
        raise RuntimeError('Archive manifest does not match checkout.')
    base = external(root, manifest_file.parent)
    for e in doc['entries']:
        source = safe(base, e['path'])
        target = safe(root, e['path'])
        if target.exists() or target.is_symlink():
            raise RuntimeError('Restore refuses to overwrite: ' + e['path'])
        if not source.is_file() or sha(source) != e['sha256']:
            raise RuntimeError('Corrupt archived file: ' + e['path'])
    for e in doc['entries']:
        source = base / e['path']
        target = root / e['path']
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
        if sha(target) != e['sha256']:
            raise RuntimeError('Restore checksum mismatch: ' + e['path'])
    print(f'RESTORED: {len(doc["entries"])} exact files.')


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=('plan','apply','restore'))
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument('--archive', type=Path, default=None)
    parser.add_argument('--plan', type=Path)
    parser.add_argument('--manifest', type=Path)
    parser.add_argument('--confirm', type=str)
    parser.add_argument('--review-file', type=Path, help='CSV listing individually reviewed REVIEW files to archive')
    args = parser.parse_args(argv)
    root = args.root.resolve()
    try:
        if args.command == 'plan':
            make_plan(root, args.archive or default_archive(root))
        elif args.command == 'apply':
            if not args.plan:
                parser.error('--apply requires --plan')
            apply(root, args.plan.resolve(), args.confirm or '',
                  args.review_file.resolve() if args.review_file else None)
        else:
            if not args.manifest:
                parser.error('restore requires --manifest')
            restore(root, args.manifest.resolve())
    except (OSError, RuntimeError, ValueError, KeyError) as exc:
        print('ERROR:', exc, file=sys.stderr)
        return 2
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
