"""Local-only T3 repository housekeeping. Python 3.10+. Never touches GitHub.

Default: inventory/plan only. --apply copies reports/outputs to a sibling
archive; moves redundant stage docs, backup copies, and unreferenced probes
into that archive. Primary code and runtime input files remain in place.
Restore using --restore ARCHIVE_DIR. All moves have SHA-256 manifest entries.
"""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import shutil
import sys

CODE_DIRS = ('Lagrangian', 'Numerical', 'RGE', 'common', 'model', 'physics',
             'scripts', 'tests', 'tools', 'validation', 'configs', 'studies')
ROOT_DOC_NAMES = {'BATCH1_README.md', 'BATCH2_README.md', 'SMOKE_TEST_README.md',
                  'README.txt', 'README_OnePromptT3Integration.md',
                  'README_T3CompletionBatch.md', 'README_T3EFT1WilsonTensor.md',
                  'README_T3ExtendedC5ExpressionReduction.md',
                  'README_T3ExtendedFlavorIndexFix.md',
                  'README_T3ExtendedLowEnergyBatch.md',
                  'README_T3ExtendedMatcheteLoadFix.md',
                  'README_T3ExtendedPhysicalMatching.md',
                  'README_T3HierarchicalExtraYukawa.md',
                  'README_T3RGAssistedExtendedC5.md',
                  'README_T3Stage2IndexDiagnostic.md',
                  'README_THESIS_COMPARISONS.txt'}
GENERATED_ROOT_PATTERNS = ('*.log', '*_diagnostics.txt', '*probe_details.txt')
BACKUP_RE = re.compile(r'(?i)(?:\.backup[-_]|\.bak$|\.orig$|\.tmp$|~$)')
PROBE_RE = re.compile(r'(?i)(?:^Probe|Probe\.|_probe\.|^Run.*Probe\.)')
SCAN_EXTS = {'.py','.wl','.wls','.ps1','.json','.toml','.yaml','.yml','.md','.txt'}


def checksum(path:Path)->str:
    h=hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(1048576),b''):h.update(chunk)
    return h.hexdigest()


def files_at(p:Path):
    if not p.exists():return
    if p.is_file():
        yield p
    else:
        for x in p.rglob('*'):
            if x.is_file():yield x


def reference_index(root:Path):
    # Deliberately excludes generated output to avoid matching a probe's own log.
    refs=[]
    for dirname in CODE_DIRS:
        d=root/dirname
        if not d.exists():continue
        for p in d.rglob('*'):
            if p.is_file() and p.suffix.lower() in SCAN_EXTS and p.stat().st_size < 2500000:
                if BACKUP_RE.search(p.name) or PROBE_RE.search(p.name):continue
                try:refs.append((p, p.read_text(encoding='utf-8-sig',errors='replace')))
                except OSError:pass
    return refs


def build_plan(root:Path):
    refs=reference_index(root)
    actions=[]; held=[]
    def add(path,category,action,reason):
        if not path.is_file():return
        actions.append(dict(path=path.relative_to(root).as_posix(),category=category,
                            action=action,reason=reason,size=path.stat().st_size))
    for dirname in ('output','Reports/output'):
        for f in files_at(root/dirname):
            add(f,'generated-results','copy', 'Archive without removing runtime/generated inputs')
    for f in root.iterdir():
        if not f.is_file():continue
        if f.name in ROOT_DOC_NAMES:
            add(f,'historical-notes','move','Consolidated in docs/PROJECT_GUIDE.md; original archived')
        elif BACKUP_RE.search(f.name) or any(f.match(pattern) for pattern in GENERATED_ROOT_PATTERNS):
            add(f,'root-temp','move','Backup or temporary root-level artifact')
    for dname in CODE_DIRS:
        for f in files_at(root/dname):
            if f.suffix.lower() not in SCAN_EXTS and not BACKUP_RE.search(f.name):continue
            if BACKUP_RE.search(f.name):
                add(f,'backup','move','Versioned backup; active source untouched')
            elif PROBE_RE.search(f.name) and f.suffix.lower() in {'.py','.wl','.wls'}:
                # Check references by name/stem; unreferenced only is eligible.
                nm=f.name; stem=f.stem
                hits=[str(p.relative_to(root)) for p,s in refs
                      if p!=f and (nm in s or (len(stem)>9 and stem in s))]
                if hits:
                    held.append({'path':f.relative_to(root).as_posix(),
                                 'reason':'Referenced by code/config/docs', 'examples':hits[:4]})
                else:
                    add(f,'unreferenced-probe','move','No literal references found in source/config inventory')
    # Never move files if a path is included in a production config reference.
    actions.sort(key=lambda x:x['path'])
    return actions,held


def guide(root, actions, held):
    docs=[x for x in actions if x['category']=='historical-notes']
    lines=['# T3 neutrino mass project — architecture and research history','',
      'This guide consolidates the historical stage notes into a navigable record.',
      'The original documents are preserved in the external cleanup archive.',
      'This is an inventory and workflow guide, **not** a new physics certification.','',
      '## Verified production workflow (one-point smoke test, 2026-11-01)','',
      'The ordinary T3-B, alpha=-1 benchmark completed the current production numerical path:',
      '', '1. Supply a representation, SM state, and heavy-sector UV benchmark.',
      '2. Load the exported UV and intermediate RGBeta beta functions.',
      '3. Integrate the UV RGEs to the heavy-fermion threshold.',
      '4. Run the scalar-only intermediate EFT1 renormalisable state.',
      '5. Evaluate the supplied Matchete-derived final Weinberg coefficient using `FinalC5TrajectoryAdapter`.',
      '6. Run SM + Weinberg evolution to the low scale.',
      '7. Construct the neutrino mass matrix, Takagi diagonalise, extract oscillation observables.',
      '', 'The prior smoke result reported `PASS_END_TO_END_OPERATIONAL` with finite oscillation outputs.',
      'This demonstrates software operation for that benchmark, not physics validation for every model.',
      '', '## Where to look first','',
      '- `pipeline.py`: overall project entry point (consult its CLI before invoking).',
      '- `Numerical/orchestration/PipelineNumericalResults.py`: integrated numerical diagnostics.',
      '- `Numerical/fitting/ScanCLI.py`: supported fixed/scan configuration entry point.',
      '- `Numerical/fitting/ParameterScan.py`: per-point numerical workflow.',
      '- `Numerical/running/`: numerical RGE and EFT evolution.',
      '- `Numerical/core/`: state representations and beta evaluation.',
      '- `physics/NeutrinoObservables.py`: Takagi and mixing observables.',
      '- `Lagrangian/` and `RGE/`: matching, model construction and beta generation.',
      '- `configs/`: benchmark and study configuration files.',
      '- `tests/`: regression and physics checks.',
      '', '## Research progression (conceptual chronology)','',
      '1. **Model setup:** T3 representations, group invariants, field content and Lagrangians.',
      '2. **RG and matching exports:** Matchete and RGBeta source and coefficient artifacts.',
      '3. **Numerical model:** UV parameters, state vectors, coupled running.',
      '4. **EFT evolution:** heavy fermion threshold, scalar EFT, Weinberg coefficient.',
      '5. **Phenomenology:** neutrino masses, mass splittings, mixing and scans.',
      '6. **General extensions:** additional quartics/Yukawas and conditional one-step C5 studies.',
      '7. **Operational integration:** single-point production smoke test passed.',
      '', '## Important scope distinction','',
      '- **Production path:** ordinary split-scalar numerical scan, with a working one-point end-to-end run.',
      '- **Experimental extensions:** additional couplings and hypothetical RGBeta-to-Matchete maps; not automatically production-certified.',
      '- **Known restrictions:** a fully certified general coupling-basis map and sequential thresholds are not established for every extension.',
      '', '## Historical source notes preserved by this cleanup','',
      'The archive manifest gives the original and archived location for every item.',
      'The headings below are extracted from local notes when this cleanup runs.','']
    for row in docs:
        p=root/row['path']; txt=p.read_text(encoding='utf-8-sig',errors='replace')
        headings=[l.strip() for l in txt.splitlines() if l.lstrip().startswith('#')][:12]
        lines += ['### '+p.name,'','Original: `'+row['path']+'`','']
        if headings:lines+=['Sections:']+['- '+h.lstrip('# ').strip() for h in headings]+['']
        else:lines+=['Historical stage note (no Markdown headings detected).','']
    lines+=['## Cleanup conventions','',
      '- `output/` and `Reports/output/` are **copied** to external archive, not removed: they may contain required runtime artifacts.',
      '- Old stage README files are moved after being indexed above.',
      '- Unreferenced probe source files and backup copies are moved to an external archive, not irreversibly deleted.',
      '- Referenced probes stay in the project until their dependencies can be evaluated.',
      '- No primary physics source, tests or configs are changed by this cleanup.',
      '', '## Restore','',
      'Run `python scripts/maintenance/CleanRepository.py --restore "PATH_TO_ARCHIVE"` from the repository root.',
      'Restoration refuses to overwrite modified or pre-existing paths.','']
    return '\n'.join(lines)+'\n'


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[2])
    ap.add_argument('--apply',action='store_true',help='Perform archive and local cleanup after review')
    ap.add_argument('--restore',type=Path,help='Restore files moved by a prior cleanup')
    ap.add_argument('--archive-root',type=Path,help='External archive base (default sibling of repo)')
    args=ap.parse_args(); root=args.root.resolve()
    if not (root/'Numerical').is_dir() or not (root/'Lagrangian').is_dir():
        ap.error('Not a T3 project repository root')
    if args.restore:
        archive=args.restore.resolve()
        manifest=json.loads((archive/'manifest.json').read_text(encoding='utf-8'))
        copied=0
        for item in manifest['completed']:
            if item['action']!='move':continue
            src=archive/'files'/item['path']; dst=root/item['path']
            if dst.exists():raise FileExistsError(f'Cannot overwrite existing {dst}')
            if not src.is_file() or checksum(src)!=item['sha256']:
                raise RuntimeError(f'Archived file changed or missing: {src}')
            dst.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(src,dst);copied+=1
        print('Restored files:',copied);return
    actions,held=build_plan(root)
    now=datetime.now(timezone.utc).strftime('%Y%m%d_%H%M%S_UTC')
    base=(args.archive_root or (root.parent/'T3_repository_archives')).resolve()
    if base==root or root in base.parents:ap.error('Archive MUST be outside repository')
    archive=base/('cleanup_'+now)
    plan={'format':'t3-cleanup-plan-v1','root':str(root),'archive':str(archive),
          'summary':{'copy':sum(x['action']=='copy' for x in actions),
                     'move':sum(x['action']=='move' for x in actions),
                     'held_probes':len(held)},'actions':actions,'held_probes':held}
    report=root/'cleanup_preview.json'
    if not args.apply:
        report.write_text(json.dumps(plan,indent=2),encoding='utf-8')
        print('PREVIEW ONLY — no source or results moved')
        print('Archive copies:',plan['summary']['copy'],'moves:',plan['summary']['move'],
              'held/referenced probes:',len(held))
        print('Review:',report);return
    archive.mkdir(parents=True,exist_ok=False)
    guide_content=guide(root,actions,held)  # Capture note contents before moving originals.
    completed=[]
    manifest={'format':'t3-cleanup-manifest-v1','created_utc':now,
              'source_root':str(root),'completed':completed,'held_probes':held}
    out=archive/'manifest.json'
    # Manifest is rewritten after every file, making partial transfers recoverable.
    for item in actions:
        src=root/item['path']; dst=archive/'files'/item['path']
        dst.parent.mkdir(parents=True,exist_ok=True)
        digest=checksum(src)
        shutil.copy2(src,dst)
        if checksum(dst)!=digest:raise IOError('Archive integrity failure: '+str(src))
        record={**item,'sha256':digest}
        completed.append(record)
        out.write_text(json.dumps(manifest,indent=2),encoding='utf-8')
        if item['action']=='move':src.unlink()
    guidepath=root/'docs/PROJECT_GUIDE.md'
    guidepath.parent.mkdir(parents=True,exist_ok=True)
    if guidepath.exists():
        existing=archive/'preexisting_PROJECT_GUIDE.md';shutil.copy2(guidepath,existing)
    guidepath.write_text(guide_content,encoding='utf-8')
    (archive/'cleanup_plan.json').write_text(json.dumps(plan,indent=2),encoding='utf-8')
    if report.exists():report.unlink()
    print('COMPLETE — external archive:',archive)
    print('Archive integrity recorded in:',out)
    print('Guide:',guidepath)
    print('NOTE: Output artifacts remain active to preserve runtime dependencies.')

if __name__=='__main__':
    try:main()
    except Exception as exc:
        print('CLEANUP ERROR:',exc,file=sys.stderr);raise
