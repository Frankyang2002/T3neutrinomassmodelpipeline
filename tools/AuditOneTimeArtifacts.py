"""Conservative, local-only audit/archival of verified root one-time artifacts.

No scientific output/, Reports/, configs/, physics code, or historical research
scripts are deletion targets. Only five exact top-level filenames can be archived.

    python -X utf8 -m tools.AuditOneTimeArtifacts
    python -X utf8 -m tools.AuditOneTimeArtifacts --apply --confirm ARCHIVE_ONE_TIME_ROOT_FILES

An archive copy and SHA-256 manifest are created OUTSIDE the repository before
any original is unlinked. Plain text references, dirty or untracked files, and
incomplete scan coverage block automatic archival.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
REPORT = Path('output/cleanup_audit/one_time_files_v1')
CONFIRM = 'ARCHIVE_ONE_TIME_ROOT_FILES'
# Explicit allowlist. No glob-based deletion.
TARGETS = {
    'overnight_refinement.log': 'historical reduced-scalar optimisation console log',
    'reduced_interactive_reports.log': 'historical report-generation console log',
    'reduced_interactive_reports_retry.log': 'historical retry console log',
    't -q --tb=line': 'accidental pager/terminal-help capture',
    't -q --tb=lineqqq': 'accidental pager/terminal-help capture',
}
PAGER_TARGETS = frozenset(k for k in TARGETS if k.startswith('t -q'))
READABLE_SUFFIXES = {
    '.py', '.ps1', '.wl', '.md', '.txt', '.json', '.jsonl', '.csv', '.tsv',
    '.yaml', '.yml', '.toml', '.xml', '.html', '.htm', '.tex', '.log',
    '.ini', '.cfg', '.sh', '.bat', '.rst', '.bib', '.js', '.css', '.ipynb',
}
IGNORED_DIRS = frozenset({'.git', '.venv', 'venv', '__pycache__', '.pytest_cache',
                           'node_modules', '.tox', '.mypy_cache', '.ruff_cache'})
SELF_PATHS = frozenset({
    'tools/AuditOneTimeArtifacts.py',
    'tests/python/test_one_time_artifact_cleanup.py',
    'docs/ONE_TIME_ARTIFACT_CLEANUP.md',
})
REVIEW_NOTES = ('README_*.txt',)
REVIEW_ROOT_ZIPS = ('T3_*.zip',)
REVIEW_LAUNCHERS = (
    'scripts/ArchiveSupersededFits.ps1',
    'scripts/AuditAndPruneReportsFull.ps1',
    'scripts/RetireOldStudyDisplays.ps1',
    'scripts/FullRunWithWithinModelReports.ps1',
)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def git(root: Path, *args: str) -> subprocess.CompletedProcess[bytes]:
    return subprocess.run(['git', '-C', str(root), *args], stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE, check=False)


def file_git_state(root: Path, name: str) -> str:
    # `git status` includes both staged and unstaged changes, and untracked.
    q = git(root, 'ls-files', '--error-unmatch', '--', name)
    if q.returncode:
        return 'UNTRACKED_OR_NOT_IN_INDEX'
    st = git(root, 'status', '--porcelain', '--untracked-files=all', '--', name)
    if st.returncode:
        return 'GIT_STATUS_ERROR'
    return 'CLEAN' if not st.stdout else 'MODIFIED_IN_GIT'


def _has_any(haystack: bytes, needles: dict[str, bytes]) -> set[str]:
    lower = haystack.lower()
    return {name for name, needle in needles.items() if needle in lower}


def _search_file(path: Path, needles: dict[str, bytes]) -> set[str]:
    # Streaming search, includes tokens that cross a buffer boundary.
    found: set[str] = set()
    overlap = max((len(x) for x in needles.values()), default=0) - 1
    previous = b''
    with path.open('rb') as f:
        while True:
            block = f.read(128 * 1024)
            if not block:
                break
            found.update(_has_any(previous + block, needles))
            if len(found) == len(needles):
                return found
            previous = (previous + block)[-overlap:] if overlap > 0 else b''
    return found


def text_references(root: Path, names: list[str]) -> tuple[dict[str, list[str]], int, list[str]]:
    references: dict[str, list[str]] = {name: [] for name in names}
    needles = {name: name.encode('ascii').lower() for name in names}
    scanned = 0
    errors: list[str] = []
    for parent, subdirs, filenames in os.walk(root):
        basedir = Path(parent)
        subdirs[:] = sorted(d for d in subdirs if d not in IGNORED_DIRS and
                            (basedir / d).relative_to(root).as_posix() != 'output/cleanup_audit')
        for filename in sorted(filenames):
            path = basedir / filename
            rel = path.relative_to(root).as_posix()
            if (rel in names or rel in SELF_PATHS or rel.startswith('output/cleanup_audit/')
                    or path.suffix.lower() not in READABLE_SUFFIXES or not path.is_file()):
                continue
            try:
                found = _search_file(path, needles)
                scanned += 1
                for name in found:
                    references[name].append(rel)
            except (OSError, ValueError) as exc:
                errors.append(f'{rel}: {exc}')
    return references, scanned, errors


def make_plan(root: Path) -> dict:
    if git(root, 'rev-parse', '--is-inside-work-tree').returncode:
        raise RuntimeError('Project root is not a Git worktree')
    present = [name for name in TARGETS if (root / name).is_file()]
    refs, scanned, errors = text_references(root, present)
    rows: list[dict] = []
    for name in TARGETS:
        src = root / name
        if not src.is_file():
            continue
        state = file_git_state(root, name)
        content_ok = True
        if name in PAGER_TARGETS:
            try:
                # Do not remove a file with an odd name if it no longer
                # contains the observed terminal-help capture.
                content_ok = b'Commands marked with *' in src.read_bytes()
            except OSError:
                content_ok = False
        elif src.stat().st_size > 5_000_000:
            content_ok = False  # Not a small one-time console log anymore.
        why = []
        if state != 'CLEAN':
            why.append(state)
        if refs.get(name):
            why.append('REFERENCED_BY_OTHER_FILES')
        if not content_ok:
            why.append('UNEXPECTED_CONTENT_OR_SIZE')
        if errors:
            why.append('TEXT_SCAN_ERRORS')
        rows.append({
            'path': name,
            'classification': 'ARCHIVE_CANDIDATE' if not why else 'REVIEW_ONLY',
            'reason': ';'.join(why) if why else TARGETS[name],
            'bytes': src.stat().st_size,
            'sha256': sha256(src),
            'git_state': state,
            'references': sorted(refs[name]),
        })
    review = []
    for pattern in REVIEW_NOTES + REVIEW_ROOT_ZIPS:
        for p in sorted(root.glob(pattern)):
            if p.is_file():
                review.append({'path': p.relative_to(root).as_posix(),
                               'reason': 'historical note or scientific output bundle; review manually'})
    for item in REVIEW_LAUNCHERS:
        if (root / item).is_file():
            review.append({'path': item, 'reason': 'legacy utility; check manual commands before retiring'})
    for p in sorted((root / 'report_cleanup_archive').glob('*/cleanup_manifest.json')):
        review.append({'path': p.parent.relative_to(root).as_posix(),
                       'reason': 'previous cleanup archive; preserve until archival provenance reviewed'})
    return {
        'format': 't3-one-time-artifact-audit-v1',
        'root': str(root),
        'head': git(root, 'rev-parse', 'HEAD').stdout.decode(errors='replace').strip(),
        'targets': rows,
        'review_only': review,
        'text_files_scanned': scanned,
        'read_errors': errors,
        'note': 'Only exact allowlisted root files may be externally archived; no scientific artifacts touched.',
    }


def write_report(root: Path, plan: dict, report: Path) -> None:
    folder = report if report.is_absolute() else root / report
    folder.mkdir(parents=True, exist_ok=True)
    (folder / 'audit.json').write_text(json.dumps(plan, indent=2) + '\n', encoding='utf-8')
    with (folder / 'inventory.csv').open('w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=['path', 'classification', 'reason', 'bytes', 'sha256', 'git_state', 'references'])
        writer.writeheader()
        for row in plan['targets']:
            writer.writerow({**row, 'references': '; '.join(row['references'])})
    summary = ['# One-time artifact cleanup audit', '',
               f"- HEAD: `{plan['head']}`",
               f"- Text files scanned: **{plan['text_files_scanned']}**",
               f"- Read errors: **{len(plan['read_errors'])}**",
               '- Source, tests, docs and numerical output are never auto-archived.', '',
               '## Exact root-file candidates', '',
               '| Filename | Classification | Reason | References |',
               '|---|---|---|---:|']
    for row in plan['targets']:
        summary.append(f"| `{row['path']}` | {row['classification']} | {row['reason']} | {len(row['references'])} |")
    summary += ['', '## Additional manual-review candidates', '',
                '| Filename | Reason |', '|---|---|']
    for row in plan['review_only']:
        summary.append(f"| `{row['path']}` | {row['reason']} |")
    if not plan['review_only']:
        summary.append('| None | — |')
    (folder / 'summary.md').write_text('\n'.join(summary) + '\n', encoding='utf-8')


def archive(root: Path, plan: dict, destination: Path) -> list[str]:
    candidates = [row for row in plan['targets'] if row['classification'] == 'ARCHIVE_CANDIDATE']
    if not candidates:
        return []
    dest = destination.resolve()
    if dest == root or root in dest.parents:
        raise RuntimeError('Archive must be outside the repository')
    if dest.exists():
        raise RuntimeError(f'External archive path already exists: {dest}')
    # Recheck all source bytes before changing anything.
    for row in candidates:
        path = root / row['path']
        if not path.is_file() or sha256(path) != row['sha256'] or file_git_state(root, row['path']) != 'CLEAN':
            raise RuntimeError(f'Candidate changed after audit: {row["path"]}')
    dest.mkdir(parents=True, exist_ok=False)
    copied = []
    try:
        for row in candidates:
            src = root / row['path']
            target = dest / row['path']
            shutil.copy2(src, target)
            if sha256(target) != row['sha256']:
                raise RuntimeError(f'Backup hash mismatch for {row["path"]}')
            copied.append(row)
        manifest = {
            'format': 't3-one-time-file-archive-v1',
            'archive_generated_utc': datetime.now(timezone.utc).isoformat(),
            'project_root': str(root), 'head': plan['head'],
            'archived_files': [{k: row[k] for k in ('path', 'bytes', 'sha256')} for row in copied],
        }
        (dest / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n', encoding='utf-8')
        for row in copied:
            if sha256(root / row['path']) != row['sha256']:
                raise RuntimeError('Source changed before deletion: ' + row['path'])
        # Already copied and verified all files; rollback sources if any unlink fails.
        removed: list[dict] = []
        try:
            for row in copied:
                (root / row['path']).unlink()
                removed.append(row)
        except OSError:
            for row in removed:
                src = root / row['path']
                if not src.exists():
                    shutil.copy2(dest / row['path'], src)
            raise
        return [row['path'] for row in copied]
    except Exception:
        # Keep the external backup for inspection; never delete a valid archive.
        raise


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=ROOT)
    parser.add_argument('--report', type=Path, default=REPORT)
    parser.add_argument('--apply', action='store_true')
    parser.add_argument('--confirm', default='')
    parser.add_argument('--archive', type=Path,
                        help='External, nonexistent archive folder; defaults to dated sibling directory')
    ns = parser.parse_args(argv)
    root = ns.root.resolve()
    if ns.apply and ns.confirm != CONFIRM:
        parser.error(f'--apply requires --confirm {CONFIRM}')
    plan = make_plan(root)
    write_report(root, plan, ns.report)
    print(f"ONE-TIME ROOT CANDIDATES: {sum(r['classification'] == 'ARCHIVE_CANDIDATE' for r in plan['targets'])}")
    print(f"BLOCKED ROOT CANDIDATES: {sum(r['classification'] != 'ARCHIVE_CANDIDATE' for r in plan['targets'])}")
    print(f"OTHER MANUAL REVIEW FILES: {len(plan['review_only'])}")
    print(f"TEXT FILES SCANNED: {plan['text_files_scanned']}; READ ERRORS: {len(plan['read_errors'])}")
    report = ns.report if ns.report.is_absolute() else root / ns.report
    print('REPORT:', report / 'summary.md')
    if not ns.apply:
        print('DRY RUN. No files changed.')
        return 0
    date_suffix = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    destination = ns.archive or (root.parent / ('T3_one_time_archive_' + date_suffix))
    moved = archive(root, plan, destination)
    print(f'EXTERNALLY ARCHIVED: {len(moved)}')
    if moved:
        print('ARCHIVE:', destination.resolve())
        print('MANIFEST:', destination.resolve() / 'manifest.json')
    print('No numerical outputs, fitted configurations, scripts, or tests changed.')
    return 0


if __name__ == '__main__':
    try:
        raise SystemExit(main())
    except (OSError, RuntimeError) as exc:
        print(f'FAIL CLOSED: {exc}', file=sys.stderr)
        raise SystemExit(2) from exc
