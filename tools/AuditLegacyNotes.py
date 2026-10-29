"""Conservative local audit/optional external archive of one-time root README notes.

Runs from T3 project root: python -X utf8 -m tools.AuditLegacyNotes
Never deletes scripts, numerical outputs, report archives, scientific ZIPs, or thesis notes.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
REPORT = Path('output/cleanup_audit/legacy_notes_v2')
SKIP_DIRS = {'.git', '.venv', 'venv', '__pycache__', '.mypy_cache',
             '.pytest_cache', '.ruff_cache', 'node_modules', 'output',
             'Reports', 'report_cleanup_archive', '.cleanup_quarantine'}
TEXT_SUFFIXES = {'.py', '.ps1', '.wl', '.md', '.txt', '.tex', '.json', '.csv',
                 '.yml', '.yaml', '.toml', '.ini', '.cfg', '.sh', '.bat', '.xml'}
PROTECTED_NOTE_NAMES = {'README_THESIS_COMPARISONS.txt'}
SCIENTIFIC_BUNDLES = (
    'T3_EFT1_OneLoop_Raw_Outputs.zip',
    'T3_Scalar_Threshold_C5_Outputs.zip',
    'T3_symbolic_C5_outputs.zip',
)
RETAIN_SCRIPTS = (
    'scripts/ArchiveSupersededFits.ps1',
    'scripts/AuditAndPruneReportsFull.ps1',
    'scripts/RetireOldStudyDisplays.ps1',
    'scripts/FullRunWithWithinModelReports.ps1',
)
REVIEW_ARCHIVES = (
    'report_cleanup_archive/20261020T030047Z-a225d356',
    'report_cleanup_archive/20261023T231628Z-f2dcc347',
)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as source:
        for block in iter(lambda: source.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def git_output(root: Path, *args: str) -> str:
    cp = subprocess.run(['git', *args], cwd=root, capture_output=True, text=True)
    if cp.returncode:
        raise RuntimeError('git ' + ' '.join(args) + ': ' + cp.stderr.strip())
    return cp.stdout


def tracked_clean(root: Path, relative: str) -> tuple[bool, str]:
    cp = subprocess.run(['git', 'ls-files', '--error-unmatch', '--', relative],
                        cwd=root, capture_output=True, text=True)
    if cp.returncode:
        return False, 'NOT_GIT_TRACKED'
    status = git_output(root, 'status', '--porcelain=v1', '--untracked-files=all', '--', relative)
    if status.strip():
        return False, 'MODIFIED_OR_STAGED'
    return True, 'CLEAN_TRACKED'


def discover_notes(root: Path) -> list[str]:
    return sorted(p.name for p in root.glob('README_*.txt') if p.is_file())


def read_sources(root: Path, note_names: list[str]):
    """Scan project source text for literal note filenames, excluding generated data.

    All references are conservative: textual mentions may be false positives, but
    no referenced note is approved for automatic archival.
    """
    matches: dict[str, list[dict]] = {name: [] for name in note_names}
    errors: list[str] = []
    count = 0
    for path in root.rglob('*'):
        if any(part in SKIP_DIRS for part in path.relative_to(root).parts[:-1]):
            continue
        if not path.is_file() or path.suffix.lower() not in TEXT_SUFFIXES:
            continue
        if path.parent == root and path.name in note_names:
            continue
        if path.stat().st_size > 8 * 1024 * 1024:
            continue
        relative = path.relative_to(root).as_posix()
        try:
            content = path.read_text(encoding='utf-8-sig', errors='strict')
        except (UnicodeError, OSError) as exc:
            errors.append(f'{relative}: {exc}')
            continue
        count += 1
        for line_no, line in enumerate(content.splitlines(), 1):
            for name in note_names:
                if name in line:
                    matches[name].append({'file': relative, 'line': line_no,
                                          'excerpt': line.strip()[:180]})
    return matches, count, errors


def plan(root: Path) -> dict:
    notes = discover_notes(root)
    refs, scanned, errors = read_sources(root, notes)
    rows = []
    for name in notes:
        clean, state = tracked_clean(root, name)
        protected = name in PROTECTED_NOTE_NAMES
        ref_count = len(refs[name])
        if protected:
            classification = 'RETAIN_THESIS_PROVENANCE'
        elif ref_count:
            classification = 'REVIEW_REFERENCED'
        elif not clean:
            classification = 'REVIEW_LOCAL_CHANGES'
        elif errors:
            classification = 'REVIEW_INCOMPLETE_SCAN'
        else:
            classification = 'ARCHIVE_CANDIDATE'
        file = root / name
        rows.append({'path': name, 'classification': classification,
                     'references': ref_count, 'git_state': state,
                     'bytes': file.stat().st_size, 'sha256': sha256(file)})
    return {'root': str(root), 'head': git_output(root, 'rev-parse', 'HEAD').strip(),
            'scanned_text_files': scanned, 'scan_errors': errors,
            'notes': rows, 'references': refs,
            'protected_scientific_bundles': list(SCIENTIFIC_BUNDLES),
            'protected_one_time_scripts': list(RETAIN_SCRIPTS),
            'previous_archives_not_moved': list(REVIEW_ARCHIVES),
            'excluded': ['README.txt', 'README.md', 'SHA256_MANIFEST.json']}


def write_report(root: Path, report: dict, destination: Path) -> None:
    folder = destination if destination.is_absolute() else root / destination
    folder.mkdir(parents=True, exist_ok=True)
    (folder / 'audit.json').write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    with (folder / 'inventory.csv').open('w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=['path', 'classification', 'references',
                                                'git_state', 'bytes', 'sha256'])
        writer.writeheader()
        writer.writerows(report['notes'])
    lines = [
        '# Legacy root README notes: local read-only audit', '',
        f"- HEAD: `{report['head']}`",
        f"- Text source files scanned: **{report['scanned_text_files']}**",
        f"- Read errors: **{len(report['scan_errors'])}**", '',
        '| File | Classification | References |', '|---|---|---:|',
    ]
    for row in report['notes']:
        lines.append(f"| `{row['path']}` | `{row['classification']}` | {row['references']} |")
    lines += ['', '## Reference details', '']
    for row in report['notes']:
        name = row['path']
        if not report['references'][name]:
            continue
        lines.append(f'### `{name}`')
        for item in report['references'][name][:100]:
            lines.append(f"- `{item['file']}:{item['line']}` — `{item['excerpt'].replace('`', chr(39))}`")
    lines += ['', '## Explicitly retained / manual review', '',
              '- Root `README.txt`, main `README.md`, `SHA256_MANIFEST.json`.',
              '- `README_THESIS_COMPARISONS.txt` (thesis asset provenance).',
              '- All three scientific-output ZIP archives.',
              '- Four legacy PS1 utilities: they are scripts, not one-time console logs.',
              '- Both existing `report_cleanup_archive` directories and their manifests.',
              '- Historical `output/` and `Reports/` numerical/research artifacts.', '',
              'No files are changed in dry-run mode. Apply only archives unreferenced, clean, tracked notes.',
              'The external archive records each original path and SHA-256.', '']
    (folder / 'summary.md').write_text('\n'.join(lines), encoding='utf-8')


def archive_notes(root: Path, report: dict, archive_root: Path) -> int:
    eligible = [x for x in report['notes'] if x['classification'] == 'ARCHIVE_CANDIDATE']
    if not eligible:
        print('NO ELIGIBLE FILES: nothing archived.')
        return 0
    if archive_root.exists():
        raise RuntimeError(f'Archive already exists: {archive_root}')
    if archive_root.is_relative_to(root):
        raise RuntimeError('Archive must be outside project root')
    archive_root.mkdir(parents=True, exist_ok=False)
    copied: list[dict] = []
    removed: list[Path] = []
    try:
        for entry in eligible:
            source = root / entry['path']
            clean, state = tracked_clean(root, entry['path'])
            if not clean or sha256(source) != entry['sha256']:
                raise RuntimeError(f'Source changed since audit: {entry["path"]} ({state})')
            dest = archive_root / entry['path']
            shutil.copy2(source, dest)
            if sha256(dest) != entry['sha256']:
                raise RuntimeError(f'Archive verification failed: {entry["path"]}')
            copied.append({**entry, 'archive_path': str(dest)})
        manifest = {
            'format': 't3-legacy-notes-external-archive-v1',
            'head': report['head'],
            'project_root': str(root),
            'archived_files': copied,
            'retained_unchanged': ['README_THESIS_COMPARISONS.txt',
                                   *SCIENTIFIC_BUNDLES, *RETAIN_SCRIPTS, *REVIEW_ARCHIVES],
        }
        (archive_root / 'manifest.json').write_text(
            json.dumps(manifest, indent=2) + '\n', encoding='utf-8')
        for entry in eligible:
            source = root / entry['path']
            if sha256(source) != entry['sha256']:
                raise RuntimeError(f'Source changed before retirement: {entry["path"]}')
            source.unlink()
            removed.append(source)
    except Exception:
        for source in removed:
            backup = archive_root / source.name
            if backup.is_file():
                shutil.copy2(backup, source)
        raise
    print(f'EXTERNALLY ARCHIVED: {len(copied)}')
    print(f'ARCHIVE: {archive_root}')
    print(f'MANIFEST: {archive_root / "manifest.json"}')
    return len(copied)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--root', type=Path, default=ROOT)
    ap.add_argument('--report', type=Path, default=REPORT)
    ap.add_argument('--apply', action='store_true')
    ap.add_argument('--confirm', default='')
    ap.add_argument('--archive', type=Path)
    args = ap.parse_args(argv)
    root = args.root.resolve()
    if not (root / '.git').exists() or not (root / 'pipeline.py').is_file():
        ap.error('Not a T3 project checkout')
    if args.apply and args.confirm != 'ARCHIVE_UNREFERENCED_LEGACY_NOTES':
        ap.error('Apply requires --confirm ARCHIVE_UNREFERENCED_LEGACY_NOTES')
    audit = plan(root)
    write_report(root, audit, args.report)
    candidates = sum(x['classification'] == 'ARCHIVE_CANDIDATE' for x in audit['notes'])
    print(f'ROOT README NOTES: {len(audit["notes"])}')
    print(f'ELIGIBLE UNREFERENCED NOTES: {candidates}')
    print(f'BLOCKED/RETAINED NOTES: {len(audit["notes"]) - candidates}')
    print(f'TEXT SOURCES SCANNED: {audit["scanned_text_files"]}; READ ERRORS: {len(audit["scan_errors"])}')
    path = args.report if args.report.is_absolute() else root / args.report
    print(f'REPORT: {path / "summary.md"}')
    if args.apply:
        date = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
        archive = (args.archive or root.parent / f'T3_legacy_notes_archive_{date}').resolve()
        archive_notes(root, audit, archive)
    else:
        print('DRY RUN. No files moved.')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
