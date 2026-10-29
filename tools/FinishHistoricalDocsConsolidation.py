"""Finish the T3 historical-document consolidation without touching physics code.

This is a LOCAL-ONLY, archive-first migration. It embeds 27 documents retained
solely by historical one-time tool literals in the seven existing thematic
volumes, updates Markdown navigation links and the historical index, and
preserves exact originals + pre-migration volume/index bytes externally.

Use ``python -X utf8 -m tools.FinishHistoricalDocsConsolidation`` for a dry run.
Apply explicitly with ``--apply --confirm CONSOLIDATE_HISTORICAL_REFERENCES``.
"""
from __future__ import annotations

import argparse
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

# Read-only inspected project source: Git commit 7400a4c64f2731a4a2b3351ba76ae3961e1adf27.
# This migration reads ALL source contents from the user's LOCAL checkout.
BASELINE_HEAD = '7400a4c64f2731a4a2b3351ba76ae3961e1adf27'
CONFIRM = 'CONSOLIDATE_HISTORICAL_REFERENCES'
REPORT = Path('output/cleanup_audit/docs_consolidation_final/summary.json')

# Twenty-three documents mentioned by the already-completed filename migration.
RENAME_ONLY = (
    'GENERAL_T3_16_MODEL_UV_RUNNING.md',
    'GENERAL_T3_C5_FINITE_MASS_LOGARITHMS.md',
    'GENERAL_T3_C5_LINEAR_RESPONSE_CORRECTION.md',
    'GENERAL_T3_CONDITIONAL_UV_TO_C5.md',
    'GENERAL_T3_COUPLED_UV_RUNNING.md',
    'GENERAL_T3_EFT1_SCALAR_PAIR_MIXING.md',
    'GENERAL_T3_EFT1_WILSON_SOURCE_TRANSPORT.md',
    'GENERAL_T3_EXTRA_C5_PHASE_RESPONSE.md',
    'GENERAL_T3_FERMION_FIRST_THRESHOLDS.md',
    'GENERAL_T3_FINAL_C5_AND_WEINBERG_RUNNING.md',
    'GENERAL_T3_GENERAL_QUARTIC_RGBETA_PROGRAMME.md',
    'GENERAL_T3_MATCHETE_TENSOR_INITIALIZATION.md',
    'GENERAL_T3_MATCHING_RG_SOURCE_PILOTS.md',
    'GENERAL_T3_QUARTIC_CONVENTION_C5.md',
    'GENERAL_T3_RGBETA_EXPORT_COVERAGE_AUDIT.md',
    'GENERAL_T3_SCALAR_QUARTIC_VALIDATION.md',
    'GENERAL_T3_SCALAR_TENSOR_CERTIFICATION.md',
    'GENERAL_T3_TENSOR_110101_LOCAL_CORRECTION.md',
    'GENERAL_T3_TENSOR_110101_SU2_COVARIANCE.md',
    'GENERAL_T3_TENSOR_ALIGNED_GEOMETRY.md',
    'GENERAL_T3_TENSOR_EXPLICIT_GEOMETRY.md',
    'GENERAL_T3_WEINBERG_SOURCE_TRANSPORT.md',
    'GENERAL_T3_WOLFRAMSCRIPT_SELFTEST_FIX.md',
)

# Four notes referenced only by historical one-time audit tool definitions.
AUDIT_ONLY = (
    'ONE_TIME_ARTIFACT_CLEANUP.md',
    'REPOSITORY_CLEANUP_AUDIT.md',
    'SCRIPT_DEPENDENCY_REFINEMENT.md',
    'T3_ZERO_BOUNDARY_LOW_ENERGY_BATCH8.md',
)
CANDIDATES = RENAME_ONLY + AUDIT_ONLY

# These have live code/metadata references and are NOT included:
# HIGH_SCALAR_HALF_FERMION_START.md, T3_REDUCED_SCALAR_POTENTIAL_CONCLUSION.md,
# TEST_BASELINE_FIX.md. Other protected core documents also remain untouched.

VOLUMES = {
    'tensors': 'RESEARCH_GENERAL_T3_TENSORS_AND_QUARTICS.md',
    'running': 'RESEARCH_GENERAL_T3_RUNNING_AND_MATCHING.md',
    'matchete': 'RESEARCH_MATCHETE_AND_MASS_EIGEN.md',
    'eft1': 'RESEARCH_EFT1_THRESHOLDS_AND_WEINBERG.md',
    'fits': 'RESEARCH_NUMERICAL_FITS_AND_BENCHMARKS.md',
    'general': 'RESEARCH_GENERAL_SCALAR_AND_PILOTS.md',
    'maintenance': 'REPOSITORY_MAINTENANCE_HISTORY.md',
}
INDEX = 'RESEARCH_HISTORY_INDEX.md'

# These are semantic, HISTORICAL literals (not files loaded by these tools).
# Each occurrence is separately inspected by exact line structure below.
STATIC_REFERENCES = {
    'tools/RenameNumberedResearchFiles.py': set(RENAME_ONLY),
    'tests/python/test_descriptive_research_file_names.py': {
        'GENERAL_T3_16_MODEL_UV_RUNNING.md',
        'GENERAL_T3_EFT1_SCALAR_PAIR_MIXING.md',
        'GENERAL_T3_TENSOR_110101_SU2_COVARIANCE.md',
    },
    'tools/AuditOneTimeArtifacts.py': {'ONE_TIME_ARTIFACT_CLEANUP.md'},
    'tools/AuditDocumentationRetirement.py': {
        'REPOSITORY_CLEANUP_AUDIT.md', 'T3_ZERO_BOUNDARY_LOW_ENERGY_BATCH8.md',
    },
    'tools/RefineScriptDependencyAudit.py': {'SCRIPT_DEPENDENCY_REFINEMENT.md'},
}
# Exact inert historical literals in the *previous* consolidation tool's
# synthetic test fixture. The final audit previously counted these seven
# fixture/assertion lines as unknown runtime references, even though they never
# read project documents. Do NOT exempt the whole test module: other lines,
# including real file reads or unreviewed new references, remain blockers.
PREVIOUS_CONSOLIDATION_FIXTURE_PATH = 'tests/python/test_historical_docs_consolidation.py'
PREVIOUS_CONSOLIDATION_FIXTURE_LINES = frozenset({
    "(root/'docs/GENERAL_T3_TENSOR_ALIGNED_GEOMETRY.md').write_text(",
    "'# Action identity\\n\\n[Geometry](GENERAL_T3_TENSOR_ALIGNED_GEOMETRY.md)\\n\\n'",
    "'```python\\n# [keep](GENERAL_T3_TENSOR_ALIGNED_GEOMETRY.md)\\n```\\n',encoding='utf8')",
    "assert group_for('GENERAL_T3_TENSOR_ALIGNED_GEOMETRY.md')=='tensors'",
    "assert group_for('GENERAL_T3_FINAL_C5_AND_WEINBERG_RUNNING.md')=='running'",
    "assert group_for('ONE_TIME_ARTIFACT_CLEANUP.md')=='maintenance'",
    "assert rows_by_name(plan)['GENERAL_T3_TENSOR_ALIGNED_GEOMETRY.md']['classification']=='RETAIN_REFERENCED'",
})

SELF = {'tools/FinishHistoricalDocsConsolidation.py',
        'tests/python/test_final_historical_docs_consolidation.py'}
TEXT_EXT = {'.py', '.wl', '.ps1', '.md', '.rst', '.txt', '.tex', '.json',
            '.csv', '.toml', '.ini', '.cfg', '.yml', '.yaml', '.sh', '.bat'}
EXCLUDE_TOP = {'.git', '.venv', 'output', '.cleanup_quarantine',
               'report_cleanup_archive', '__pycache__'}


def sha_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha(path: Path) -> str:
    return sha_bytes(path.read_bytes())


def git(root: Path, *args: str) -> bytes:
    result = subprocess.run(['git', '-C', str(root), *args], capture_output=True)
    if result.returncode:
        raise RuntimeError('git ' + ' '.join(args) + ': ' +
                           result.stderr.decode('utf-8', 'replace').strip())
    return result.stdout


def tracked(root: Path) -> set[str]:
    return {PurePosixPath(os.fsdecode(x)).as_posix()
            for x in git(root, 'ls-files', '-z').split(b'\0') if x}


def dirty(root: Path) -> set[str]:
    values: set[str] = set()
    chunks = git(root, 'status', '--porcelain=v1', '-z', '--untracked-files=all').split(b'\0')
    i = 0
    while i < len(chunks) and chunks[i]:
        line = os.fsdecode(chunks[i])
        if len(line) < 4:
            raise RuntimeError('Unexpected git status record')
        values.add(PurePosixPath(line[3:]).as_posix())
        if 'R' in line[:2] or 'C' in line[:2]:
            i += 1
            if i < len(chunks) and chunks[i]:
                values.add(PurePosixPath(os.fsdecode(chunks[i])).as_posix())
        i += 1
    return values


def anchor(name: str) -> str:
    return 'source-' + re.sub(r'[^a-z0-9]+', '-', name.lower()).strip('-')


def volume_for(name: str) -> str:
    if name in AUDIT_ONLY[:3]:
        return VOLUMES['maintenance']
    if name == 'T3_ZERO_BOUNDARY_LOW_ENERGY_BATCH8.md':
        return VOLUMES['fits']
    if 'EFT1_' in name or 'FERMION_FIRST_' in name or 'WEINBERG_' in name:
        return VOLUMES['eft1']
    if any(s in name for s in ('QUARTIC', 'TENSOR', 'GEOMETRY', 'COVARIANCE',
                              'RGBETA', 'MATCHETE_TENSOR')):
        return VOLUMES['tensors']
    return VOLUMES['running']


def target_for(name: str) -> str:
    return volume_for(name) + '#' + anchor(name)


def static_reference_is_historical(path: str, name: str, line: str) -> bool:
    """Verify allowlisted occurrences are inert mapping/test/ignore literals."""
    if path == PREVIOUS_CONSOLIDATION_FIXTURE_PATH:
        return name in line and line.strip() in PREVIOUS_CONSOLIDATION_FIXTURE_LINES
    if name not in STATIC_REFERENCES.get(path, set()):
        return False
    s = line.strip()
    if path == 'tools/RenameNumberedResearchFiles.py':
        return bool(re.fullmatch(r"['\"]GENERAL_T3_STAGE[^'\"]+\.md['\"]\s*:\s*['\"]" +
                                 re.escape(name) + r"['\"],?", s))
    if path == 'tests/python/test_descriptive_research_file_names.py':
        return name in s and ('assert ' in s or "'" in s or '"' in s)
    if path == 'tools/AuditOneTimeArtifacts.py':
        return s == "'docs/ONE_TIME_ARTIFACT_CLEANUP.md',"
    if path == 'tools/AuditDocumentationRetirement.py':
        return s == "AUDIT_GUIDE = 'docs/REPOSITORY_CLEANUP_AUDIT.md'" or bool(
            re.fullmatch(r"'docs/T3_ZERO_BOUNDARY_LOW_ENERGY_BATCH8.md':\s*'docs/SCALAR_SECTOR.md',", s))
    if path == 'tools/RefineScriptDependencyAudit.py':
        return s == '"docs/SCRIPT_DEPENDENCY_REFINEMENT.md",'
    return False


def iter_text(root: Path):
    for current, dirnames, filenames in os.walk(root, followlinks=False):
        relative = Path(current).relative_to(root)
        dirnames[:] = [d for d in dirnames if d not in EXCLUDE_TOP and
                       not (Path(current) / d).is_symlink()]
        if relative.as_posix() == 'Reports':
            dirnames[:] = [d for d in dirnames if d != 'output']
        for filename in filenames:
            p = Path(current) / filename
            if p.is_symlink() or p.suffix.lower() not in TEXT_EXT:
                continue
            rel = p.relative_to(root).as_posix()
            if rel in SELF or rel.startswith('Reports/output/'):
                continue
            yield rel, p


def inspect_references(root: Path, names: set[str]):
    code_issues: list[dict[str, Any]] = []
    doc_mentions: list[dict[str, Any]] = []
    permitted: list[dict[str, Any]] = []
    errors: list[str] = []
    scanned = 0
    for rel, path in iter_text(root):
        try:
            lines = path.read_text(encoding='utf-8-sig').splitlines()
        except (OSError, UnicodeError) as e:
            errors.append(f'{rel}: {e}')
            continue
        scanned += 1
        for n, line in enumerate(lines, 1):
            if '.md' not in line:
                continue
            for name in names:
                if name not in line:
                    continue
                row = {'source': rel, 'line': n, 'name': name,
                       'text': line.strip()[:260]}
                if rel == 'docs/' + name:
                    continue
                if rel.startswith('docs/') and rel.endswith('.md'):
                    doc_mentions.append(row)
                elif static_reference_is_historical(rel, name, line):
                    permitted.append(row)
                else:
                    code_issues.append(row)
    return {'scanned': scanned, 'read_errors': errors,
            'blocking_external_references': code_issues,
            'markdown_mentions': doc_mentions,
            'historical_literal_references': permitted}


def index_update(text: str, named: dict[str, dict[str, str]]) -> str:
    """Promote the requested original files from 'retained' to the lookup table."""
    old_header = '## Retained individual documents'
    lookup = '## Original-file lookup'
    if text.count(old_header) != 1 or text.count(lookup) != 1:
        raise RuntimeError('Unexpected research index structure; refusing rewrite')
    prefix, retained = text.split(old_header, 1)
    if not prefix.endswith('\n'):
        raise RuntimeError('Unexpected research index lookup boundary')
    additions = []
    for name in sorted(named):
        old = f'| [{name}]({name}) | `RETAIN_REFERENCED` |'
        if retained.count(old) != 1:
            raise RuntimeError('Retained index entry not unique: ' + name)
        retained = retained.replace(old + '\n', '')
        if old in retained:
            raise RuntimeError('Could not remove retained entry: ' + name)
        vol = named[name]['volume']
        digest = named[name]['sha256']
        additions.append(f'| `docs/{name}` | [{vol}](<{target_for(name)}>) | `{digest}` |')
    # The original-file lookup ends immediately before retained-header.
    prefix = prefix.rstrip('\n') + '\n' + '\n'.join(additions) + '\n\n'
    result = prefix + old_header + retained
    for name in named:
        if f'| [{name}]({name}) |' in result:
            raise RuntimeError('Old standalone link remains in index: ' + name)
    return result


def rewrite_links(text: str, names: set[str]) -> str:
    """Redirect *Markdown links*, not narrative names or fenced historical code."""
    def new_url(url: str) -> str:
        if url.startswith(('#', '/', 'http:', 'https:', 'mailto:')) or '://' in url:
            return url
        leaf = url.split('#', 1)[0].split('?', 1)[0].replace('\\', '/').rsplit('/', 1)[-1]
        if leaf in names:
            return target_for(leaf)
        return url

    inline = re.compile(r'(?P<a>\]\(\s*<?)(?P<u>[^)\s>]+)(?P<b>[^)]*\))')
    reflink = re.compile(r'^(?P<a>\s*\[[^]]+\]:\s*<?)(?P<u>\S+?)(?P<b>>?\s*)$')
    output = []
    fence = None
    for line in text.splitlines(keepends=True):
        match = re.match(r'^\s*(`{3,}|~{3,})', line)
        if match:
            token = match.group(1)
            if fence is None:
                fence = token[0]
            elif token[0] == fence:
                fence = None
            output.append(line)
            continue
        if fence is None:
            line = inline.sub(lambda m: m['a'] + new_url(m['u']) + m['b'], line)
            ref = reflink.match(line.rstrip('\r\n'))
            if ref:
                line = ref['a'] + new_url(ref['u']) + ref['b'] + ('\n' if line.endswith('\n') else '')
        output.append(line)
    return ''.join(output)


def append_records(volume_text: str, source_docs: dict[str, bytes]) -> str:
    if not volume_text.startswith('# ') or '## Contents' not in volume_text:
        raise RuntimeError('Unexpected historical volume structure')
    if '\n---\n' not in volume_text:
        raise RuntimeError('Missing volume table-of-contents delimiter')
    toc = []
    blocks = []
    for name, original in sorted(source_docs.items()):
        raw = original.decode('utf-8-sig')
        title = next((x[2:].strip()[:180] for x in raw.splitlines()[:40]
                      if x.startswith('# ')), name.replace('_', ' '))
        toc.append(f'- [{title}](#{anchor(name)}) — `docs/{name}`')
        body = rewrite_links(raw, set(CANDIDATES)).rstrip('\r\n')
        blocks.append('\n'.join([
            f'<a id="{anchor(name)}"></a>', '',
            f'## Original document: `docs/{name}`', '',
            f'**Original SHA-256:** `{sha_bytes(original)}`', '',
            '**Status:** Historical record; verify current source and core documentation',
            'before using its old commands or provisional conclusions.', '',
            '<details>', f'<summary>Expand original: {title}</summary>', '',
            body, '', '</details>', '', '---', '',
        ]))
    # Amend existing contents rather than producing a new topic document.
    before, after = volume_text.split('\n---\n', 1)
    before = before.rstrip('\n') + '\n' + '\n'.join(toc) + '\n'
    amended = before + '\n---\n' + after
    amended = rewrite_links(amended, set(CANDIDATES))
    return amended.rstrip('\r\n') + '\n\n' + '\n'.join(blocks).rstrip('\r\n') + '\n'


def make_plan(root: Path, *, expected_head: str | None = BASELINE_HEAD):
    root = root.resolve()
    head = git(root, 'rev-parse', 'HEAD').decode().strip()
    if expected_head and head != expected_head:
        raise RuntimeError(f'Unexpected HEAD {head}; audited baseline is {expected_head}')
    tracked_paths = tracked(root)
    dirty_paths = dirty(root)
    targets = {f'docs/{n}' for n in CANDIDATES}
    targets.update(f'docs/{v}' for v in VOLUMES.values())
    targets.add(f'docs/{INDEX}')
    errors = []
    if len(CANDIDATES) != 27 or len(set(CANDIDATES)) != len(CANDIDATES):
        errors.append('Incorrect or duplicate candidate inventory')
    for rel in sorted(targets):
        p = root / rel
        if rel not in tracked_paths or not p.is_file() or p.is_symlink():
            errors.append('Missing/untracked/symlinked required path: ' + rel)
        if rel in dirty_paths:
            errors.append('Dirty required path: ' + rel)
    if errors:
        raise RuntimeError('\n'.join(errors))
    contents = {name: (root / 'docs' / name).read_bytes() for name in CANDIDATES}
    for n, payload in contents.items():
        try:
            payload.decode('utf-8-sig')
        except UnicodeError as e:
            raise RuntimeError('Non-UTF8 source document ' + n) from e
    named = {name: {'volume': volume_for(name), 'sha256': sha_bytes(contents[name]),
                    'bytes': len(contents[name])} for name in CANDIDATES}
    idx = (root / 'docs' / INDEX).read_text(encoding='utf-8-sig')
    index_update(idx, named)  # validate original-file index before any mutation
    evidence = inspect_references(root, set(CANDIDATES))
    return {'head': head, 'candidate_count': len(CANDIDATES),
            'root_markdown_before': len(list((root / 'docs').glob('*.md'))),
            'root_markdown_after': len(list((root / 'docs').glob('*.md'))) - len(CANDIDATES),
            'document_records': named,
            'volumes_changed': sorted({volume_for(n) for n in CANDIDATES}),
            'reference_audit': evidence,
            'ready': not evidence['read_errors'] and not evidence['blocking_external_references']}


def write_report(root: Path, plan: dict[str, Any], report: Path):
    dest = report if report.is_absolute() else root / report
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(json.dumps(plan, indent=2) + '\n', encoding='utf-8')
    summary = dest.with_suffix('.md')
    audit = plan['reference_audit']
    lines = [
        '# Final documentation consolidation: dry-run audit', '',
        f'- HEAD: `{plan["head"]}`',
        f'- Historical documents eligible for final archival: **{plan["candidate_count"]}**',
        f'- Current root Markdown: **{plan["root_markdown_before"]}**',
        f'- Projected root Markdown: **{plan["root_markdown_after"]}**',
        f'- Existing volumes to extend: **{len(plan["volumes_changed"])}**',
        f'- Text sources scanned: **{audit["scanned"]}**',
        f'- Historical code/test/maintenance literals: **{len(audit["historical_literal_references"])}**',
        f'- Markdown mentions/links (redirected if active): **{len(audit["markdown_mentions"])}**',
        f'- Unknown external references: **{len(audit["blocking_external_references"])}**',
        f'- Read errors: **{len(audit["read_errors"])}**',
        f'- Safe to apply: **{plan["ready"]}**', '',
        '## Protected live references', '',
        '- `HIGH_SCALAR_HALF_FERMION_START.md`: active orchestration entry point',
        '- `T3_REDUCED_SCALAR_POTENTIAL_CONCLUSION.md`: source-path metadata in a numerical runner',
        '- `TEST_BASELINE_FIX.md`: referenced by historical SHA256 manifest', '',
        '## Blockers', '',
    ]
    for row in audit['blocking_external_references']:
        lines.append(f'- `{row["source"]}:{row["line"]}` -> `{row["name"]}`')
    for s in audit['read_errors']:
        lines.append(f'- Read error: `{s}`')
    lines += ['', 'This tool never modifies physics code, tests, configs, outputs or GitHub.',
              'Archive-first apply preserves exact original docs plus pre-change index and volumes.']
    summary.write_text('\n'.join(lines) + '\n', encoding='utf-8')
    return summary


def apply(root: Path, plan: dict[str, Any], archive: Path) -> dict[str, Any]:
    root = root.resolve()
    archive = archive.resolve()
    if not plan['ready']:
        raise RuntimeError('Unknown references or read errors block the migration')
    if archive == root or root in archive.parents or archive.exists():
        raise RuntimeError('Archive must be a new directory OUTSIDE the checkout')
    current = make_plan(root, expected_head=plan['head'])
    if current != plan:
        raise RuntimeError('Local files or dependency references changed; rerun dry run')
    target_names = set(CANDIDATES)
    source_contents = {n: (root / 'docs' / n).read_bytes() for n in CANDIDATES}
    named = plan['document_records']
    volumes: dict[str, dict[str, bytes]] = {}
    for name in CANDIDATES:
        volumes.setdefault(volume_for(name), {})[name] = source_contents[name]
    new_texts: dict[str, bytes] = {}
    for vol, members in volumes.items():
        rel = f'docs/{vol}'
        old = (root / rel).read_text(encoding='utf-8-sig')
        new_texts[rel] = append_records(old, members).encode('utf-8')
    index_rel = 'docs/' + INDEX
    index_old = (root / index_rel).read_text(encoding='utf-8-sig')
    new_texts[index_rel] = rewrite_links(index_update(index_old, named), target_names).encode('utf-8')
    # Adjust links in other retained docs only, preserving code/outputs untouched.
    for p in sorted((root / 'docs').glob('*.md')):
        rel = p.relative_to(root).as_posix()
        if rel in new_texts or p.name in target_names:
            continue
        old = p.read_text(encoding='utf-8-sig')
        updated = rewrite_links(old, target_names)
        if updated != old:
            new_texts[rel] = updated.encode('utf-8')
    existing_paths = sorted(new_texts)
    archives: dict[str, bytes] = {}
    for name in CANDIDATES:
        archives['originals/docs/' + name] = source_contents[name]
    for rel in existing_paths:
        archives['before/' + rel] = (root / rel).read_bytes()
    # Ensure no untracked/modified files would be overwritten.
    if set(existing_paths) - tracked(root):
        raise RuntimeError('A target Markdown file is not tracked')
    if set(existing_paths).intersection(dirty(root)):
        raise RuntimeError('A target Markdown file is dirty; refusing overwrite')
    manifest = {
        'format': 't3-final-historical-doc-consolidation-v1',
        'head': plan['head'],
        'source_root': str(root),
        'archived_original_documents': [
            {'path': 'docs/' + n, 'sha256': named[n]['sha256'],
             'new_location': target_for(n)} for n in CANDIDATES],
        'prechange_snapshots': [
            {'path': rel, 'sha256': sha_bytes(archives['before/' + rel])}
            for rel in existing_paths],
        'replacement_hashes': [
            {'path': rel, 'sha256': sha_bytes(body)}
            for rel, body in sorted(new_texts.items())],
    }
    archive.mkdir(parents=True, exist_ok=False)
    changed: list[str] = []
    removed: list[str] = []
    try:
        for rel, payload in archives.items():
            dest = archive / rel
            dest.parent.mkdir(parents=True, exist_ok=True)
            dest.write_bytes(payload)
            if sha(dest) != sha_bytes(payload):
                raise RuntimeError('Archive SHA-256 verification failed: ' + rel)
        (archive / 'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n', encoding='utf-8')
        for rel, body in new_texts.items():
            if (root / rel).read_bytes() != archives['before/' + rel]:
                raise RuntimeError('Local file changed immediately before write: ' + rel)
            (root / rel).write_bytes(body)
            changed.append(rel)
            if sha(root / rel) != sha_bytes(body):
                raise RuntimeError('New documentation hash verification failed: ' + rel)
        # Original docs are only removed after all new volumes/index are verified.
        for n in CANDIDATES:
            src = root / 'docs' / n
            if sha(src) != named[n]['sha256']:
                raise RuntimeError('Original changed immediately before removal: ' + n)
            src.unlink()
            removed.append('docs/' + n)
        return {'originals_archived': len(CANDIDATES), 'updated_markdown_files':len(new_texts),
                'archive':str(archive),'manifest':str(archive / 'manifest.json'),
                'root_markdown_after': len(list((root / 'docs').glob('*.md')))}
    except BaseException:
        for rel in removed:
            src = archive / 'originals' / rel
            if src.is_file() and not (root / rel).exists():
                shutil.copy2(src, root / rel)
        for rel in changed:
            backup = archive / 'before' / rel
            if backup.is_file():
                shutil.copy2(backup, root / rel)
        raise


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--report', type=Path, default=REPORT)
    parser.add_argument('--apply', action='store_true')
    parser.add_argument('--confirm', default='')
    parser.add_argument('--archive', type=Path)
    args = parser.parse_args(argv)
    if args.apply and args.confirm != CONFIRM:
        parser.error('Apply requires --confirm ' + CONFIRM)
    try:
        root = args.root.resolve()
        plan = make_plan(root)
        summary = write_report(root, plan, args.report)
        print('CURRENT ROOT MARKDOWN:', plan['root_markdown_before'])
        print('HISTORICAL NOTES TO CONSOLIDATE:', plan['candidate_count'])
        print('PROJECTED ROOT MARKDOWN:', plan['root_markdown_after'])
        print('UNKNOWN REFERENCE BLOCKERS:', len(plan['reference_audit']['blocking_external_references']))
        print('READ ERRORS:', len(plan['reference_audit']['read_errors']))
        print('REPORT:', summary)
        if not args.apply:
            print('DRY RUN: no original documents changed')
            return 0 if plan['ready'] else 2
        if not plan['ready']:
            raise RuntimeError('Audit has blockers; do not apply')
        archive = args.archive or root.parent / ('T3_final_historical_docs_archive_' +
                    datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ'))
        result = apply(root, plan, archive)
        print('ORIGINAL DOCUMENTS ARCHIVED:', result['originals_archived'])
        print('DOCS MARKDOWN FILES UPDATED:', result['updated_markdown_files'])
        print('REMAINING ROOT MARKDOWN:', result['root_markdown_after'])
        print('ARCHIVE:', result['archive'])
        print('MANIFEST:', result['manifest'])
        return 0
    except (RuntimeError, OSError, UnicodeError, ValueError) as error:
        print('ERROR:', error, file=sys.stderr)
        return 2


if __name__ == '__main__':
    raise SystemExit(main())
