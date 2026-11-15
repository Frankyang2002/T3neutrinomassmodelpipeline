#!/usr/bin/env python3
"""Thesis-focused, external-archive-first review of a CurateThesisRepository plan.

Does not execute scientific code, change Git history, or contact GitHub.

From the repository root (Windows PowerShell):
  python -X utf8 scripts/maintenance/ReviewThesisCleanup.py audit --plan PATH
  python -X utf8 scripts/maintenance/ReviewThesisCleanup.py apply --plan PATH --approved PATH --confirm ARCHIVE_REVIEWED_HISTORY
  python -X utf8 scripts/maintenance/ReviewThesisCleanup.py restore --manifest PATH

Audit output is outside the repo, next to the supplied plan. The default
review CSV approves NOTHING. Only explicitly reviewed historical documents
and LaTeX auxiliaries may be archived. Physics source, tests, configurations,
matching exports, numerical results and thesis figures are NEVER eligible.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
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
from typing import Any

SCHEMA = 't3-thesis-curation-v1-plan'
REVIEW_SCHEMA = 't3-thesis-history-review-v2'
CONFIRM = 'ARCHIVE_REVIEWED_HISTORY'
DOC_EXT = {'.md', '.txt', '.rst'}
AUX_EXT = ('.aux', '.log', '.toc', '.out', '.compile.log')

# These are *recognised as historical reference sources*, not evidence that
# the referenced document is unimportant to a thesis. Each archival decision
# still requires a human approval and a checksum-backed external copy.
HISTORICAL_REFERENCE_EXACT = {
    'docs/RESEARCH_HISTORY_INDEX.md',
    'docs/REPOSITORY_MAINTENANCE_HISTORY.md',
    'docs/HISTORICAL_DOCS_CONSOLIDATION_GUIDE.md',
    'tests/python/test_public_repository_cleanup.py',  # synthetic fixtures only
}
HISTORICAL_REFERENCE_PREFIX = (
    'docs/research/',
    'output/cleanup_audit/',
    'output/maintenance/',
    'output/audits/',
    'report_cleanup_archive/',
)
# Old generated provenance can mention a filename even when nothing loads it.
# This is an archival *caveat*, not a reason to pretend the link is unneeded.
SAVED_RECORD_PREFIX = ('output/', 'Reports/output/')

# No broad classification of e.g. tools/* or tests/* as harmless; those may
# read inputs or assert file existence, and are blockers until inspected.
CORE_PUBLIC_DOCUMENTS = {
    'docs/README.md', 'docs/PHYSICS.md', 'docs/SCALAR_SECTOR.md',
    'docs/NUMERICAL.md', 'docs/RESULTS.md', 'docs/MODEL_AND_PIPELINE.md',
    'docs/PROJECT_GUIDE.md', 'docs/REPOSITORY_STRUCTURE.md',
    'docs/REPRODUCIBILITY.md', 'docs/IMPLEMENTATION_STATUS.md',
    'docs/THESIS_SCOPE.md', 'docs/THESIS_ARTIFACT_MAP.md',
    'docs/THESIS_CLEANUP.md', 'docs/TEST_STRATEGY.md',
    'docs/INSTALLATION.md', 'docs/REGENERATE_CURRENT_REPORTS.md',
}
IMPORTANT_FIGURE_PREFIXES = (
    'Reports/output/thesis_general_5_1/',
    'Reports/output/thesis_general_5_2/',
    'Reports/output/thesis_general_5_3/',
    'Reports/output/thesis_qualitative/',
    'output/studies/thesis_qualitative/',
    'output/studies/fermion_start_lambda_comparison/',
    'output/full16_symbolic/',
    'output/parallel_t3/background/',
    'output/parallel_t3/eft1_llss_',
    'output/parallel_t3/thesis_',
    'output/studies/general_scalar_extended/zee_style_general_sensitivity/zee_general_all16_merged_v1/',
)
DOC_DEV_RE = re.compile(
    r'(?i)(?:^README[_-].*(?:batch|worker|stage|fix|trial|install|pilot|release)|'
    r'^BATCH\d+|^WORKER[_-]|^WORKSTREAM[_-]|^ROADMAP|^HANDOFF|'
    r'(?:batch\d+|_batch\d+|_roadmap|_handover|_dev_notes|_next_gate))'
)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for data in iter(lambda: handle.read(1 << 20), b''):
            h.update(data)
    return h.hexdigest()


def git(root: Path, *args: str) -> bytes:
    p = subprocess.run(['git', '-C', str(root), *args], capture_output=True, check=False)
    if p.returncode:
        raise RuntimeError(f'git {" ".join(args)}: {p.stderr.decode("utf8", "replace")}')
    return p.stdout


def git_head(root: Path) -> str:
    return git(root, 'rev-parse', 'HEAD').decode('ascii').strip()


def tracked(root: Path) -> set[str]:
    return {os.fsdecode(x).replace('\\', '/') for x in git(root, 'ls-files', '-z').split(b'\0') if x}


def safe_path(root: Path, rel: str) -> Path:
    posix = PurePosixPath(rel)
    if not rel or rel != posix.as_posix() or posix.is_absolute() or '\\' in rel \
            or ':' in rel or any(x in ('', '.', '..') for x in posix.parts):
        raise RuntimeError(f'Invalid project-relative path: {rel!r}')
    result = root.joinpath(*posix.parts)
    if result.is_symlink() or not result.resolve().is_relative_to(root.resolve()):
        raise RuntimeError(f'Unsafe symlink or path escape: {rel}')
    return result


def outside_repo(root: Path, directory: Path) -> Path:
    dest = directory.resolve()
    if dest == root.resolve() or dest.is_relative_to(root.resolve()):
        raise RuntimeError('The archive must be outside the repository')
    return dest


def load_plan(path: Path) -> dict[str, Any]:
    data = json.loads(path.read_text(encoding='utf-8-sig'))
    if data.get('schema') != SCHEMA or not isinstance(data.get('entries'), list):
        raise RuntimeError('Not a CurateThesisRepository v1 plan')
    paths = [e['path'] for e in data['entries']]
    if len(paths) != len(set(paths)):
        raise RuntimeError('Duplicate file paths in source plan')
    if not data.get('source_head') or not data.get('archive_root') or not data.get('repository'):
        raise RuntimeError('Missing plan provenance')
    return data


def ref_type(ref: str) -> str:
    # Do not orphan a retained research index: its Markdown links are real.
    if ref.endswith('/INDEX.md') or ref == 'docs/research/INDEX.md':
        return 'DOCUMENT_OR_TEST'
    if ref in HISTORICAL_REFERENCE_EXACT or ref.startswith(HISTORICAL_REFERENCE_PREFIX):
        return 'HISTORY'
    if ref.startswith(SAVED_RECORD_PREFIX):
        return 'SAVED_PROVENANCE'
    # Source-code references in tools and one-off scripts are NOT assumed inert.
    if ref.endswith(('.py', '.ps1', '.wl', '.wls', '.m', '.sh')):
        return 'CODE_OR_TEST'
    if ref.startswith(('docs/', 'Numerical/', 'RGE/', 'Lagrangian/', 'tests/')):
        return 'DOCUMENT_OR_TEST'
    return 'OTHER'


def eligible_kind(row: dict[str, Any]) -> bool:
    """Eligibility for REVIEW of historical files, not permission to delete."""
    rel = row['path']
    if rel in CORE_PUBLIC_DOCUMENTS or rel.startswith(IMPORTANT_FIGURE_PREFIXES):
        return False
    # A scientific source/config/result must never enter the archival approval CSV.
    if row.get('category') != 'CANDIDATE':
        return False
    if rel.startswith('report_cleanup_archive/'):
        return rel.lower().endswith(AUX_EXT)
    if rel.startswith('docs/') or '/' not in rel:
        return PurePosixPath(rel).suffix.lower() in DOC_EXT
    if rel.startswith('Reports/output/'):
        return rel.lower().endswith(AUX_EXT)
    return False


def historical_review(row: dict[str, Any]) -> dict[str, Any] | None:
    if not eligible_kind(row):
        return None
    kinds = Counter(ref_type(r) for r in row.get('references', []))
    blockers = [r for r in row.get('references', []) if ref_type(r) in {'CODE_OR_TEST', 'DOCUMENT_OR_TEST', 'OTHER'}]
    rel = row['path']
    doc = rel.startswith('docs/')
    if blockers:
        suggestion = 'BLOCKED_LIVE_REFERENCE'
    elif doc:
        suggestion = 'REVIEW_HISTORICAL_DOCUMENT'
    elif rel.startswith('report_cleanup_archive/') or rel.startswith('Reports/output/'):
        suggestion = 'REVIEW_BUILD_AUXILIARY'
    elif '/' not in rel:
        suggestion = 'REVIEW_ROOT_DEVELOPMENT_NOTE'
    else:
        suggestion = 'REVIEW_HISTORY'
    return {
        'path': rel, 'suggestion': suggestion, 'original_decision': row.get('decision', ''),
        'bytes': row.get('bytes', 0),
        'history_references': kinds['HISTORY'],
        'saved_provenance_references': kinds['SAVED_PROVENANCE'],
        'strong_reference_count': len(blockers),
        'strong_reference_examples': ' | '.join(blockers[:5]),
        'reference_examples': ' | '.join(row.get('references', [])[:5]),
        'approval': '', 'justification': '',
    }


def additional_development_docs(entries: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Protected-by-directory docs can still clutter the reader view.

    These remain inspection-only and cannot be archived by this tool.
    """
    result = []
    for row in entries:
        rel = row['path']
        name = PurePosixPath(rel).name
        if row['decision'] != 'KEEP' or rel in CORE_PUBLIC_DOCUMENTS:
            continue
        if PurePosixPath(rel).suffix.lower() in DOC_EXT and DOC_DEV_RE.search(name):
            result.append({'path': rel, 'bytes': row.get('bytes', 0),
                           'classification': 'MANUAL_SOURCE_TREE_DOC_REVIEW',
                           'note': 'Protected parent directory; inspect actual content/imports and thesis provenance'})
    return result


def group_outputs(entries: list[dict[str, Any]]) -> list[dict[str, Any]]:
    groups: dict[str, dict[str, Any]] = {}
    for e in entries:
        rel = e['path']
        if not (rel.startswith('output/') or rel.startswith('Reports/output/')):
            continue
        parts = rel.split('/')
        depth = 3 if parts[0] == 'output' else 4
        group = '/'.join(parts[:min(depth, len(parts)-1)])
        obj = groups.setdefault(group, {'group': group, 'files': 0, 'bytes': 0,
                                        'protected_files': 0, 'review_files': 0})
        obj['files'] += 1
        obj['bytes'] += int(e.get('bytes', 0))
        if e['decision'] == 'KEEP':
            obj['protected_files'] += 1
        else:
            obj['review_files'] += 1
    for value in groups.values():
        value['recommendation'] = ('MIXED_PROTECT_AND_REVIEW' if value['protected_files'] and value['review_files']
                                   else 'PRESERVE_THESIS_SOURCE' if value['protected_files']
                                   else 'REVIEW_UPSTREAM_DEPENDENCIES')
    return sorted(groups.values(), key=lambda r: (-r['bytes'], r['group']))


def write_csv(path: Path, data: list[dict[str, Any]], fields: list[str]) -> None:
    with path.open('w', newline='', encoding='utf-8-sig') as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(data)


def analyze(plan: dict[str, Any], directory: Path) -> tuple[Path, dict[str, Any]]:
    directory.mkdir(parents=True, exist_ok=True)
    items = [r for x in plan['entries'] if (r := historical_review(x)) is not None]
    extra = additional_development_docs(plan['entries'])
    groups = group_outputs(plan['entries'])
    path = directory / 'historical_document_review.csv'
    write_csv(path, items, ['path', 'suggestion', 'original_decision', 'bytes',
                            'history_references', 'saved_provenance_references',
                            'strong_reference_count', 'strong_reference_examples',
                            'reference_examples', 'approval', 'justification'])
    write_csv(directory / 'embedded_development_docs.csv', extra,
              ['path', 'bytes', 'classification', 'note'])
    write_csv(directory / 'saved_output_groups.csv', groups,
              ['group', 'files', 'bytes', 'protected_files', 'review_files', 'recommendation'])
    summary = {
        'source_head': plan['source_head'], 'source_rows': len(plan['entries']),
        'original_classification': dict(Counter(r['decision'] for r in plan['entries'])),
        'original_sizes': dict(Counter()),
        'document_review_count': len(items),
        'document_suggestions': dict(Counter(x['suggestion'] for x in items)),
        'historical_documents_with_no_strong_references': sum(x['strong_reference_count'] == 0 for x in items),
        'embedded_development_document_count': len(extra),
        'saved_output_group_count': len(groups),
        'saved_output_total_bytes': sum(x['bytes'] for x in groups),
        'archive_default_approvals': 0,
    }
    by_decision: dict[str, int] = defaultdict(int)
    for e in plan['entries']:
        by_decision[e['decision']] += int(e.get('bytes', 0))
    summary['original_sizes'] = dict(by_decision)
    (directory / 'review_summary.json').write_text(json.dumps(summary, indent=2) + '\n', encoding='utf-8')
    md = [
        '# Thesis repository cleanup: review of existing plan', '',
        f'Git HEAD: `{plan["source_head"]}`. This is a **metadata/references** review, not a live reproduction test.',
        f'Tracked entries: {len(plan["entries"]):,}. Existing categories: {summary["original_classification"]}.',
        f'Historical-document and auxiliary review entries: **{len(items):,}**.',
        f'Among them, **{summary["historical_documents_with_no_strong_references"]:,}** have no direct code/test/active-document references according to the ORIGINAL PLAN. They are candidates for review, not proof of irrelevance.',
        f'Additional development-style documents hidden inside protected source directories: **{len(extra):,}** (not auto-archivable).',
        f'Saved output groups: **{len(groups):,}**, aggregating **{summary["saved_output_total_bytes"] / (1024**2):.1f} MiB**. Every group requires thesis provenance review before removal.', '',
        '## How to use', '',
        'Review `historical_document_review.csv`. Select only rows with no strong references,',
        'confirm that none of their information is uniquely needed for the thesis or reproducibility,',
        'and set `approval=ARCHIVE` with a meaningful `justification`.',
        'Do not approve rows with `BLOCKED_LIVE_REFERENCE`.',
        'Use `apply --approved historical_document_review.csv --confirm ARCHIVE_REVIEWED_HISTORY`.',
        'All originals are SHA-256 verified in an external archive before removing the tracked files.',
        'Review `embedded_development_docs.csv` separately: it is not authorised for archive by this tool.',
        'Review `saved_output_groups.csv` for future dependency closure; this tool will not delete data.', '',
        '## Science boundaries', '',
        'Retain thesis inputs, accepted fits, 16-model data, figures and their source snapshots.',
        'The required EFT ordering is UV -> integrate out all fermions -> integrate out all scalars -> SMEFT.',
        'Existing one-scale conditional numerical results must not be relabelled as a verified general sequential matching calculation.',
    ]
    (directory / 'review_summary.md').write_text('\n'.join(md) + '\n', encoding='utf-8')
    return path, summary


def checked_archive_root(root: Path, plan: dict[str, Any]) -> Path:
    return outside_repo(root, Path(plan['archive_root']))


def preflight(root: Path, plan: dict[str, Any]) -> dict[str, dict[str, Any]]:
    if str(root.resolve()).casefold() != str(Path(plan['repository']).resolve()).casefold():
        raise RuntimeError('Repository path differs from the plan')
    if git_head(root) != plan['source_head']:
        raise RuntimeError('Git HEAD differs from the plan; regenerate the original curation plan')
    if git(root, 'status', '--porcelain=v1', '--untracked-files=no').strip():
        raise RuntimeError('Tracked working tree is dirty; commit or back up changes first')
    actual = tracked(root)
    entries = {e['path']: e for e in plan['entries']}
    if actual != set(entries):
        raise RuntimeError('Tracked paths changed since plan; regenerate it')
    return entries


def approval_rows(path: Path) -> dict[str, dict[str, str]]:
    with path.open(newline='', encoding='utf-8-sig') as f:
        reader = csv.DictReader(f)
        if not {'path', 'approval', 'justification'} <= set(reader.fieldnames or []):
            raise RuntimeError('Approval CSV must have path, approval, justification')
        result = {}
        for row in reader:
            rel = row['path'].strip()
            if rel in result:
                raise RuntimeError('Duplicate approval path: ' + rel)
            if row['approval'].strip().upper() not in ('', 'ARCHIVE', 'KEEP'):
                raise RuntimeError(f'Invalid approval value for {rel}')
            if row['approval'].strip().upper() == 'ARCHIVE':
                if len(row['justification'].strip()) < 25:
                    raise RuntimeError('Justification must be >=25 characters for ' + rel)
                result[rel] = row
        return result


def apply(root: Path, plan: dict[str, Any], approved_path: Path, confirm: str) -> Path | None:
    if confirm != CONFIRM:
        raise RuntimeError('Incorrect explicit confirmation token')
    entries = preflight(root, plan)
    approved = approval_rows(approved_path)
    if not approved:
        print('No files approved: zero removals.')
        return None
    selected: list[dict[str, Any]] = []
    for rel in sorted(approved):
        if rel not in entries:
            raise RuntimeError('Approval path missing in plan: ' + rel)
        row = entries[rel]
        review = historical_review(row)
        if review is None or review['strong_reference_count']:
            raise RuntimeError(f'Not eligible for reviewed archival: {rel}')
        src = safe_path(root, rel)
        if not src.is_file() or src.stat().st_size != row['bytes']:
            raise RuntimeError('Source file differs from the plan: ' + rel)
        selected.append({'path': rel, 'bytes': int(row['bytes']),
                         'sha256': sha256(src),
                         'reason': approved[rel]['justification'].strip()})
    # Also verify an audit operation cannot approve a core public file despite
    # its path being spoofed in a modified CSV.
    archive = checked_archive_root(root, plan)
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    dest = archive / 'phase2_archives' / (stamp + '-' + plan['source_head'][:12])
    if dest.exists():
        raise RuntimeError('Archive target already exists')
    for item in selected:
        src = safe_path(root, item['path'])
        dst = dest / item['path']
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst, follow_symlinks=False)
        if sha256(dst) != item['sha256']:
            raise RuntimeError('Archive verification failed: ' + item['path'])
    manifest = dest / 'manifest.json'
    record = {'schema': REVIEW_SCHEMA + '-archive', 'source_head': plan['source_head'],
              'repository': str(root.resolve()), 'entries': selected}
    manifest.write_text(json.dumps(record, indent=2) + '\n', encoding='utf-8')
    # Final check before mutation. Disallow concurrent tracked changes.
    preflight(root, plan)
    for item in selected:
        src = safe_path(root, item['path'])
        if sha256(src) != item['sha256']:
            raise RuntimeError('Source changed during archival: ' + item['path'])
    for item in selected:
        safe_path(root, item['path']).unlink()
    print(f'ARCHIVED {len(selected)} historical documents/auxiliaries; manifest: {manifest}')
    return manifest


def restore(root: Path, manifest: Path) -> None:
    record = json.loads(manifest.read_text(encoding='utf-8'))
    if record.get('schema') != REVIEW_SCHEMA + '-archive' or \
            Path(record.get('repository', '')).resolve() != root.resolve():
        raise RuntimeError('Manifest does not belong to this repository')
    base = outside_repo(root, manifest.parent)
    for item in record['entries']:
        src, dst = safe_path(base, item['path']), safe_path(root, item['path'])
        if not src.is_file() or sha256(src) != item['sha256'] or src.stat().st_size != item['bytes']:
            raise RuntimeError('Missing or corrupt archived source: ' + item['path'])
        if dst.exists() or dst.is_symlink():
            raise RuntimeError('Refusing to overwrite existing file: ' + item['path'])
    for item in record['entries']:
        src, dst = safe_path(base, item['path']), safe_path(root, item['path'])
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)
        if sha256(dst) != item['sha256']:
            raise RuntimeError('Restore verification failed: ' + item['path'])
    print(f'RESTORED {len(record["entries"])} files')


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('command', choices=['audit', 'apply', 'restore'])
    ap.add_argument('--plan', type=Path)
    ap.add_argument('--approved', type=Path)
    ap.add_argument('--manifest', type=Path)
    ap.add_argument('--out', type=Path)
    ap.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[2])
    ap.add_argument('--confirm')
    args = ap.parse_args(argv)
    try:
        if args.command in ('audit', 'apply'):
            if args.plan is None:
                raise RuntimeError('--plan is required')
            plan = load_plan(args.plan)
            if args.command == 'audit':
                path, summary = analyze(plan, args.out or args.plan.parent / 'phase2_review')
                print('SOURCE HEAD:', plan['source_head'])
                print('SUGGESTIONS:', summary['document_suggestions'])
                print('NO FILES REMOVED. Review CSV:', path)
            else:
                if args.approved is None:
                    raise RuntimeError('--approved is required')
                apply(args.root.resolve(), plan, args.approved, args.confirm or '')
        elif args.command == 'restore':
            if args.manifest is None:
                raise RuntimeError('--manifest is required')
            restore(args.root.resolve(), args.manifest)
    except (ValueError, KeyError, OSError, RuntimeError) as exc:
        print('ERROR:', exc, file=sys.stderr)
        return 2
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
