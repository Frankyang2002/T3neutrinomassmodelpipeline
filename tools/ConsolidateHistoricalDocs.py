"""Conservative local-only T3 historical-document consolidation.

No executable modules, fitted data, or saved result files are edited.
Only clean, Git-tracked docs/*.md with no active retained-file references
can be bundled. The complete original bytes are kept in an external archive.
"""
from __future__ import annotations

import argparse
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
import tempfile
from typing import Any

PROTECTED = frozenset({
    'PHYSICS.md', 'SCALAR_SECTOR.md', 'NUMERICAL.md', 'RESULTS.md',
    'T3_THESIS_PLAN.md', 'REGENERATE_CURRENT_REPORTS.md',
    'THESIS_FIGURE_PLACEMENT.md', 'THESIS_FIGURE_CAPTIONS.md',
    'THESIS_EXISTING_COMPARISONS_PLACEMENT.md',
    'GENERAL_SCALAR_48_FINAL_REPORTING.md',
    'T3_GENERAL_SCALAR_STUDY_SEPARATION.md',
    'FERMION_START_STUDY.md', 'T3B_YUKAWA_BASIS.md',
})

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
REPORT = Path('output/cleanup_audit/docs_consolidation_v1')
ENCODINGS = ('utf-8-sig',)
TEXT_SUFFIXES = frozenset({'.md', '.rst', '.txt', '.py', '.ps1', '.wl', '.json',
    '.yaml', '.yml', '.toml', '.ini', '.cfg', '.csv', '.tex', '.sh', '.bat'})
SCAN_SKIP_TOP = frozenset({'.git', '.venv', 'output', 'report_cleanup_archive',
    '__pycache__', '.cleanup_quarantine'})
SKIP_TOOL = 'tools/ConsolidateHistoricalDocs.py'
SKIP_TEST = 'tests/python/test_historical_docs_consolidation.py'
SKIP_GUIDE = 'docs/HISTORICAL_DOCS_CONSOLIDATION_GUIDE.md'
CONFIRM = 'CONSOLIDATE_HISTORICAL_DOCS'


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def git(root: Path, *args: str) -> bytes:
    p = subprocess.run(['git', '-C', str(root), *args], capture_output=True)
    if p.returncode:
        raise RuntimeError('git ' + ' '.join(args) + ': ' +
                           p.stderr.decode('utf-8', 'replace').strip())
    return p.stdout


def git_tracked(root: Path) -> set[str]:
    return {PurePosixPath(os.fsdecode(p)).as_posix()
            for p in git(root, 'ls-files', '-z').split(b'\0') if p}


def git_dirty(root: Path) -> set[str]:
    # --porcelain=v1 -z: rename entries have two NUL-delimited paths.
    data = git(root, 'status', '--porcelain=v1', '-z', '--untracked-files=all')
    chunks = data.split(b'\0')
    dirty: set[str] = set()
    i = 0
    while i < len(chunks) and chunks[i]:
        line = os.fsdecode(chunks[i])
        if len(line) < 4:
            raise RuntimeError('Unexpected git status -z entry')
        dirty.add(PurePosixPath(line[3:]).as_posix())
        if line[:2][0] in 'RC' or line[:2][1] in 'RC':
            i += 1
            if i < len(chunks) and chunks[i]:
                dirty.add(PurePosixPath(os.fsdecode(chunks[i])).as_posix())
        i += 1
    return dirty


def group_for(name: str) -> str:
    upper = name.upper()
    if upper.startswith('GENERAL_T3_'):
        if any(t in upper for t in ('TENSOR', 'QUARTIC', 'GEOMETRY', 'COVARIANCE',
                'INVARIANT', 'RGBETA', 'TRIPLET', 'ACTION', 'PREFACTOR',
                'SCALAR_QUARTIC', 'BC_LOCAL', 'BC_REUSED')):
            return 'tensors'
        return 'running'
    if upper.startswith(('T3_MATCHETE_', 'T3_MATCHE', 'T3_MASS_EIGEN_',
                         'T3_MIXED_SCALAR_MASS_')):
        return 'matchete'
    if upper.startswith(('T3_GENERAL_EFT1_', 'T3_GENERAL_FERMION_',
                         'T3_GENERAL_T3_', 'T3_GENERAL_TWO_THRESHOLD_')):
        return 'eft1'
    if upper.startswith(('T3_BC_', 'T3B_', 'T3_B_PLUS1_', 'FERMION_START_',
                         'HIGH_SCALAR_', 'CONSOLIDATE_FITS_',
                         'T3_REDUCED_', 'T3_ZERO_')):
        return 'fits'
    if any(t in upper for t in ('CLEANUP', 'SCRIPT_DEPENDENCY', 'DESCRIPTIVE_RESEARCH',
                                'RENAMED_QUARTIC_PROVENANCE', 'TEST_BASELINE',
                                'RETENTION_AUDIT', 'FIT_CANDIDATE_RETENTION',
                                'REPOSITORY_', 'ONE_TIME_', 'LEGACY_NOTES')):
        return 'maintenance'
    return 'general'


def anchor_for(name: str) -> str:
    return 'source-' + re.sub(r'[^a-z0-9]+', '-', name.lower()).strip('-')


def heading(path: Path) -> str:
    # Only a label/preview. Full original file is reproduced in the volume.
    try:
        content = path.read_text(encoding='utf-8-sig')
    except UnicodeError as e:
        raise RuntimeError('Could not decode Markdown ' + str(path)) from e
    for line in content.splitlines()[:40]:
        if line.startswith('# '):
            return line[2:].strip()[:180]
    return path.stem.replace('_', ' ').title()


def text_sources(root: Path):
    for current, dirs, names in os.walk(root, followlinks=False):
        rel_dir = Path(current).relative_to(root)
        dirs[:] = [d for d in dirs if not (Path(current)/d).is_symlink() and
                   (rel_dir.parts or d not in SCAN_SKIP_TOP) and
                   d not in {'.git', '.venv', '__pycache__', '.pytest_cache', '.mypy_cache'}]
        if rel_dir.as_posix() == 'Reports':
            dirs[:] = [d for d in dirs if d != 'output']
        # Never include any raw scientific dataset or generated figures in this textual audit.
        if rel_dir.parts and rel_dir.parts[0] in SCAN_SKIP_TOP:
            dirs.clear()
            continue
        for name in names:
            path = Path(current)/name
            if path.is_symlink() or path.suffix.lower() not in TEXT_SUFFIXES:
                continue
            rel = path.relative_to(root).as_posix()
            if rel in {SKIP_TOOL, SKIP_TEST, SKIP_GUIDE}:
                continue
            yield rel, path


def is_reference(line: str, name: str, suffix: str) -> tuple[bool, str] | None:
    """Classify an exact filename occurrence; avoid stem-only weak matches.

    - Executable/config/test references are blocking.
    - Active Markdown links are blocking.
    - Plain prose mentions in Markdown are nonblocking provenance.
    """
    if name not in line:
        return None
    # A Markdown link with target ending in the original filename is actionable.
    # Also cover reference-style URLs and angle autolinks, including ../docs/.
    for m in re.finditer(r'\]\(([^)]+)\)|^\s*\[[^]]+\]:\s*(\S+)|<([^<>]+)>', line):
        target = next((x for x in m.groups() if x is not None), '')
        target = target.split('#', 1)[0].split('?', 1)[0].strip().strip('<>')
        if target.replace('\\', '/').endswith('/' + name) or target == name:
            return True, 'MARKDOWN_LINK'
    if suffix in {'.md', '.rst', '.txt', '.tex'}:
        return False, 'TEXT_MENTION'
    return True, 'EXECUTABLE_OR_CONFIG_REFERENCE'


def collect_refs(root: Path, docs: dict[str, str]) -> tuple[list[dict[str, Any]], int, list[str]]:
    # Scan source tree once, and classify candidate names on exact line substrings.
    basenames = sorted(docs)
    references: list[dict[str, Any]] = []
    scanned = 0
    errors: list[str] = []
    for rel, path in text_sources(root):
        try:
            lines = path.read_text(encoding='utf-8-sig').splitlines()
        except (OSError, UnicodeError) as e:
            errors.append(f'{rel}: {e}')
            continue
        scanned += 1
        suffix = path.suffix.lower()
        for index, line in enumerate(lines, 1):
            if not line or '.md' not in line:
                continue
            for name in basenames:
                result = is_reference(line, name, suffix)
                if result is None:
                    continue
                hard, kind = result
                references.append({'source': rel, 'target': 'docs/' + name,
                    'line': index, 'kind': kind, 'blocking': hard,
                    'excerpt': line.strip()[:320]})
    return references, scanned, errors


def make_plan(root: Path) -> dict[str, Any]:
    root = root.resolve()
    if not (root/'docs').is_dir() or (root/'docs').is_symlink():
        raise RuntimeError('Missing docs/ directory or docs/ is a symlink')
    head = git(root, 'rev-parse', 'HEAD').decode().strip()
    tracked = git_tracked(root)
    dirty = git_dirty(root)
    md_files = sorted((root/'docs').glob('*.md'))
    if not md_files:
        raise RuntimeError('No root docs/*.md files found')
    candidates: dict[str, str] = {}
    rows = []
    for path in md_files:
        rel = path.relative_to(root).as_posix()
        name = path.name
        if path.is_symlink():
            status, group = 'RETAIN_SYMLINK', None
        elif name in PROTECTED or name == INDEX or name in VOLUMES.values() or rel == SKIP_GUIDE:
            status, group = 'RETAIN_CANONICAL', None
        elif rel not in tracked:
            status, group = 'RETAIN_UNTRACKED', None
        elif rel in dirty:
            status, group = 'RETAIN_LOCALLY_MODIFIED', None
        else:
            status, group = 'CANDIDATE', group_for(name)
            candidates[rel] = group
        rows.append({'path': rel, 'title': heading(path), 'size': path.stat().st_size,
                     'sha256': sha256(path), 'group': group, 'classification': status,
                     'references': 0, 'blocking_sources': []})
    refs, scan_count, errors = collect_refs(root, {Path(x).name: x for x in candidates})
    row_by_name = {x['path']: x for x in rows}
    # Fixed point: a blocked doc becomes retained, which can block docs it references.
    active = set(candidates)
    if errors:
        active.clear()  # Fail closed; never archive with incomplete reference coverage.
    else:
        while True:
            blocked = {r['target'] for r in refs if r['blocking'] and
                       r['target'] in active and r['source'] not in active}
            if not blocked:
                break
            active.difference_update(blocked)
    for ref in refs:
        if ref['target'] in row_by_name:
            row_by_name[ref['target']]['references'] += 1
            if ref['blocking'] and ref['source'] not in active:
                row_by_name[ref['target']]['blocking_sources'].append(ref['source'])
    for row in rows:
        if row['classification'] != 'CANDIDATE':
            continue
        if errors:
            row['classification'] = 'BLOCK_SCAN_ERRORS'
        elif row['path'] not in active:
            row['classification'] = 'RETAIN_REFERENCED'
        else:
            row['classification'] = 'ARCHIVE_AND_CONSOLIDATE'
    volume_members = {v: sorted([r['path'] for r in rows
         if r['classification'] == 'ARCHIVE_AND_CONSOLIDATE' and VOLUMES[r['group']] == v])
         for v in VOLUMES.values()}
    volume_members = {v: members for v, members in volume_members.items() if members}
    destinations = [f'docs/{v}' for v in volume_members] + [f'docs/{INDEX}']
    for rel in destinations:
        if (root/rel).exists() or rel in tracked:
            raise RuntimeError('Destination already exists; refuse overwrite: '+rel)
    return {'format': 't3-docs-consolidation-plan-v1', 'head': head,
            'scan_count': scan_count, 'read_errors': errors,
            'rows': rows, 'references': refs, 'volumes': volume_members,
            'new_files': destinations,
            'totals': {
                'root_markdown': len(md_files),
                'eligible': len(active),
                'protected_or_blocked': len(md_files)-len(active),
                'new_volumes': len(volume_members),
                'new_index': 1,
                'projected_root_markdown': len(md_files)-len(active)+len(volume_members)+1,
            }}


def rewrite_archived_doc_links(original: str, doc_to_location: dict[str, str]) -> str:
    """Retarget working Markdown hyperlinks, preserving code and historical prose.

    The archive preserves identical original bytes; consolidated copies change
    only hyperlinks that target a moved Markdown file. Other links are untouched.
    """
    def replacement(url: str) -> str:
        if url.startswith(('#', '/', 'http:', 'https:', 'mailto:')) or '://' in url:
            return url
        leaf = url.split('#', 1)[0].split('?', 1)[0].replace('\\', '/').rsplit('/', 1)[-1]
        return doc_to_location.get(leaf, url)

    link_re = re.compile(r'(?P<start>\]\(\s*<?)(?P<url>[^)\s>]+)(?P<end>[^)]*\))')
    def link_sub(m: re.Match) -> str:
        return m['start'] + replacement(m['url']) + m['end']

    ref_re = re.compile(r'^(?P<start>\s*\[[^]]+\]:\s*<?)(?P<url>\S+?)(?P<end>>?\s*)$')
    output: list[str] = []
    fence = None
    for line in original.splitlines(keepends=True):
        match = re.match(r'^\s*(`{3,}|~{3,})', line)
        if match:
            token = match.group(1)
            if fence is None:
                fence = token[0]
            elif token[0] == fence:
                fence = None
            output.append(line)
            continue
        if fence is not None:
            output.append(line)
            continue
        line = link_re.sub(link_sub, line)
        reference = ref_re.match(line.rstrip('\r\n'))
        if reference:
            line = reference['start'] + replacement(reference['url']) + reference['end'] + \
                    ('\n' if line.endswith('\n') else '')
        output.append(line)
    return ''.join(output)


def archived_link_targets(plan: dict[str, Any]) -> dict[str,str]:
    return {Path(path).name: volume + '#' + anchor_for(Path(path).name)
            for volume, paths in plan['volumes'].items() for path in paths}


def render_volume(root: Path, members: list[str], volume_name: str,
                  rows: dict[str, dict], link_targets: dict[str, str]) -> str:
    title = volume_name.removesuffix('.md').replace('_', ' ').title()
    lines = [f'# {title}', '',
       '> Historical research record. This document combines **original source',
       '> Markdown** under original-path headings. Local Markdown links to moved',
       '> files are updated; exact original bytes are in the external archive.',
       '> This is not a new certification',
       '> of physics, a rewritten conclusion, or a new computational workflow.',
       '> Original bytes and SHA-256 values are stored in the external archive.', '',
       '## Contents', '']
    for rel in members:
        b=Path(rel).name
        lines.append(f'- [{rows[rel]["title"]}](#{anchor_for(b)}) — `{rel}`')
    lines.extend(['', '---', ''])
    for rel in members:
        entry = rows[rel]
        name = Path(rel).name
        lines.extend([f'<a id="{anchor_for(name)}"></a>', '',
           f'## Original document: `{rel}`', '',
           f'**Original SHA-256:** `{entry["sha256"]}`', '',
           '**Status:** Historical record; check the current source and main',
           'documentation before executing commands reproduced below.', '',
           '<details>', f'<summary>Expand original: {entry["title"]}</summary>', ''])
        raw = (root/rel).read_text(encoding='utf-8-sig')
        lines.append(rewrite_archived_doc_links(raw, link_targets).rstrip('\r\n'))
        lines.extend(['', '</details>', '', '---', ''])
    return '\n'.join(lines).rstrip() + '\n'


def render_index(plan: dict[str, Any]) -> str:
    rows = {r['path']:r for r in plan['rows']}
    a = ['# T3 Historical Research Documents', '',
         'This index maps historical `docs/*.md` documents to consolidated thematic',
         'research records. The active scientific interpretation remains in',
         '`PHYSICS.md`, `SCALAR_SECTOR.md`, `NUMERICAL.md` and `RESULTS.md`.', '',
         '**Archiving does not validate an experimental claim.** Original files are',
         'recoverable from the external manifest and Git history.', '',
         f'**Source commit:** `{plan["head"]}`', '',
         '## Consolidated volumes', '']
    for volume, members in plan['volumes'].items():
        a.append(f'- [{volume}]({volume}) — {len(members)} original documents')
    a.extend(['', '## Original-file lookup', '',
       '| Original path | New location | SHA-256 |', '|---|---|---|'])
    for volume, members in plan['volumes'].items():
        for rel in members:
            e = rows[rel]
            a.append(f'| `{rel}` | [{volume}](<{volume}#{anchor_for(Path(rel).name)}>) | `{e["sha256"]}` |')
    retained=[r for r in plan['rows'] if r['classification'] != 'ARCHIVE_AND_CONSOLIDATE']
    a.extend(['', '## Retained individual documents', '',
       'Documents listed here were not moved; they may be active, referenced,',
       'locally modified, or otherwise protected.', '',
       '| File | Reason |', '|---|---|'])
    for r in retained:
        a.append(f'| [{Path(r["path"]).name}]({Path(r["path"]).name}) | `{r["classification"]}` |')
    a.extend(['', '## Boundaries', '',
        '- Scientific JSON/CSV/TEX files, fitted inputs, code, tests, and all `output/` are unchanged.',
        '- Archived Markdown may contain obsolete commands and provisional claims.',
        '- No original document text is deleted without a verified external copy.',
        '- The exact original bytes can be recovered from the external archive manifest.', ''])
    return '\n'.join(a)


def write_report(root: Path, plan: dict[str, Any], report: Path) -> None:
    out = report if report.is_absolute() else root/report
    out.mkdir(parents=True, exist_ok=True)
    (out/'plan.json').write_text(json.dumps(plan, indent=2, ensure_ascii=False)+'\n',encoding='utf-8')
    with (out/'inventory.csv').open('w',encoding='utf-8',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=['path','classification','group','size','sha256','references','blocking_sources'])
        writer.writeheader()
        for row in plan['rows']:
            writer.writerow({k: ('; '.join(row[k]) if k=='blocking_sources' else row[k])
                             for k in writer.fieldnames})
    t=plan['totals'];lines=['# Documentation consolidation plan (read-only)', '',
       f'- HEAD: `{plan["head"]}`',
       f'- Root docs Markdown: **{t["root_markdown"]}**',
       f'- Can consolidate and archive: **{t["eligible"]}**',
       f'- Retained or blocked: **{t["protected_or_blocked"]}**',
       f'- New topic volumes: **{t["new_volumes"]}**',
       f'- Expected docs Markdown count after apply: **{t["projected_root_markdown"]}**',
       f'- Text files scanned: **{plan["scan_count"]}**; read errors: **{len(plan["read_errors"])}**',
       '', '## Volume allocations', '']
    for volume,members in plan['volumes'].items():
        lines.append(f'- `{volume}`: {len(members)} source documents')
    lines.extend(['', '## Retained and blocked', ''])
    for row in plan['rows']:
        if row['classification']!='ARCHIVE_AND_CONSOLIDATE':
            reasons=', '.join(sorted(set(row['blocking_sources']))[:8])
            lines.append(f'- `{row["path"]}` — {row["classification"]}' +
                         (f'; source: {reasons}' if reasons else ''))
    if plan['read_errors']:
        lines.extend(['', '## Read errors (all moves blocked)', '']+plan['read_errors'])
    lines.extend(['', '## Safety', '',
       'Archives are created OUTSIDE the repository with original bytes and SHA-256 hashes.',
       'The planner checks current tracked status and references. No physics source or result is edited.',
       'Do not use --apply before reviewing inventory.csv and blocking reasons.', ''])
    (out/'summary.md').write_text('\n'.join(lines),encoding='utf-8')


def archive_and_apply(root: Path, plan: dict[str,Any], archive_dir: Path) -> dict[str,Any]:
    if not plan['volumes'] or not plan['totals']['eligible']:
        raise RuntimeError('No eligible historical documentation; nothing to archive')
    root=root.resolve(); archive_dir=archive_dir.resolve()
    if archive_dir==root or root in archive_dir.parents:
        raise RuntimeError('Archive must be outside the repository')
    if archive_dir.exists():
        raise RuntimeError('Archive path already exists: '+str(archive_dir))
    for rel in plan['new_files']:
        if (root/rel).exists():
            raise RuntimeError('Refusing overwrite of '+rel)
    # Immediately rescan to reject changes between report and action.
    fresh=make_plan(root)
    if json.dumps(fresh,sort_keys=True)!=json.dumps(plan,sort_keys=True):
        raise RuntimeError('Local source or reference state changed during planning; re-audit')
    rows={r['path']:r for r in plan['rows']}
    members=sorted([rel for v in plan['volumes'].values() for rel in v])
    link_targets=archived_link_targets(plan)
    prepared={f'docs/{vol}':render_volume(root,mem,vol,rows,link_targets)
              for vol,mem in plan['volumes'].items()}
    prepared[f'docs/{INDEX}']=render_index(plan)
    copied=[]; created=[]; removed=[]
    archive_dir.mkdir(parents=True,exist_ok=False)
    try:
        for rel in members:
            src=root/rel;dest=archive_dir/rel
            if sha256(src)!=rows[rel]['sha256']:
                raise RuntimeError('Source modified: '+rel)
            dest.parent.mkdir(parents=True,exist_ok=True)
            shutil.copy2(src,dest)
            if sha256(dest)!=rows[rel]['sha256']:
                raise RuntimeError('Archive verification failed: '+rel)
            copied.append(rel)
        manifest={'format':'t3-docs-consolidation-manifest-v1',
                  'head':plan['head'], 'source_root':str(root),
                  'files':[{'original_path':p,'sha256':rows[p]['sha256'],
                            'volume':next(v for v,paths in plan['volumes'].items() if p in paths)}
                           for p in members],
                  'new_files':[{ 'path':p,'sha256':hashlib.sha256(content.encode('utf-8')).hexdigest()}
                               for p,content in sorted(prepared.items())]}
        (archive_dir/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
        for rel,content in prepared.items():
            dest=root/rel
            if dest.exists():
                raise RuntimeError('Destination appeared during apply: '+rel)
            dest.write_text(content,encoding='utf-8',newline='\n')
            created.append(rel)
        # Verify complete consolidated text before retiring any original.
        for v,m in plan['volumes'].items():
            body=(root/'docs'/v).read_text(encoding='utf-8')
            for p in m:
                original=rewrite_archived_doc_links(
                    (root/p).read_text(encoding='utf-8-sig'),link_targets).rstrip('\r\n')
                if original not in body:
                    raise RuntimeError('Consolidated document omitted source text: '+p)
        for rel in members:
            if sha256(archive_dir/rel)!=rows[rel]['sha256']:
                raise RuntimeError('Archive mutated: '+rel)
            (root/rel).unlink()
            removed.append(rel)
        return {'archived':len(members),'volumes':len(plan['volumes']),
                'archive':str(archive_dir), 'manifest':str(archive_dir/'manifest.json')}
    except BaseException:
        # Restore original files, then remove newly generated docs, never touching
        # unrelated files or any saved numerical artifacts.
        for rel in removed:
            src=archive_dir/rel
            if src.exists() and not (root/rel).exists():
                shutil.copy2(src,root/rel)
        for rel in created:
            (root/rel).unlink(missing_ok=True)
        raise


def main(argv=None) -> int:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1])
    parser.add_argument('--report',type=Path,default=REPORT)
    parser.add_argument('--apply',action='store_true')
    parser.add_argument('--confirm',default='')
    parser.add_argument('--archive',type=Path)
    args=parser.parse_args(argv)
    root=args.root.resolve()
    if args.apply and args.confirm != CONFIRM:
        parser.error('Apply requires --confirm '+CONFIRM)
    try:
        plan=make_plan(root)
        write_report(root,plan,args.report)
        t=plan['totals']
        print('ROOT MARKDOWN DOCS:',t['root_markdown'])
        print('SAFE TO CONSOLIDATE:',t['eligible'])
        print('RETAINED/BLOCKED:',t['protected_or_blocked'])
        print('THEMATIC VOLUMES:',t['new_volumes'])
        print('PROJECTED DOCS:',t['projected_root_markdown'])
        print('REPORT:',(args.report if args.report.is_absolute() else root/args.report)/'summary.md')
        if not args.apply:
            print('DRY RUN: no existing documents moved or edited')
            return 0
        when=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
        archive=args.archive or root.parent/('T3_historical_docs_archive_'+when)
        outcome=archive_and_apply(root,plan,archive)
        print('CONSOLIDATED:',outcome['archived'],'historical documents')
        print('CREATED:',outcome['volumes'],'topic volumes and 1 index')
        print('ARCHIVE:',outcome['archive'])
        print('MANIFEST:',outcome['manifest'])
        return 0
    except (RuntimeError,OSError,ValueError) as exc:
        print('ERROR:',exc,file=sys.stderr)
        return 2


if __name__=='__main__':
    raise SystemExit(main())
